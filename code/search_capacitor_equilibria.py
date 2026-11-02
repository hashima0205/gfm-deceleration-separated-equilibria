"""Exact steady network/P-Q droop equations, checked against full-model roots."""
from pathlib import Path
import ast,json,re
import numpy as np
import pandas as pd
from scipy.optimize import root

R=Path(__file__).resolve().parents[1]
D=R/'basin_recalculation/capacitor_equilibrium_search';D.mkdir(exist_ok=True,parents=True)
S=pd.read_csv(R/'results/conditional_basins/data/states.csv')
settings=pd.read_csv(R/'data/references/fixed_settings.csv')
source=(R/'provenance/source_snapshot/studies/04_clearing/07_grid_strength/code/run.jl').read_text()
branches=ast.literal_eval(re.search(r'^const branches=(\[.*\])$',source,re.M).group(1))
meta=json.loads((R/'results/conditional_basins/data/definition.json').read_text())
G=2.2682301152333055;TARGET_B=.5561806200995012
buses=[1,2,3,7];nodes=[1,2,3,4,5,7,8];indices=[nodes.index(b) for b in buses]
A=1+.05*settings.Qset.to_numpy();PS=np.array([.5,.4,.4,.4]);ZG=.01+.0209j

def vc_from_full(x):
    def val(b,k):return x[int(S.loc[(S.bus==b)&(S.symbol==k),'index'].iloc[0])-1]
    return np.array([complex(val(b,f'inv{b}₊vsrc₊filter₊V_C_r'),val(b,f'inv{b}₊vsrc₊filter₊V_C_i')) for b in buses])
def packed(v):return np.r_[v.real,v.imag]
def unpacked(x):return x[:4]+1j*x[4:]
def system(B,post=True,k=None):
    k=meta["k"] if k is None else k
    Y=np.zeros((8,8),complex)
    for j,(a,b,m) in enumerate(branches,1):
        if post and j==7:continue
        y=m/(k if j in [2,3] else 1)*np.exp(-1j*np.deg2rad(75));a-=1;b-=1
        Y[a,a]+=y;Y[b,b]+=y;Y[a,b]-=y;Y[b,a]-=y
    ni=np.array(nodes)-1;M=Y[np.ix_(ni,ni)].copy();M[nodes.index(5),nodes.index(5)]+=G+1j*B
    for i in indices:M[i,i]+=1/ZG
    inverse=np.linalg.inv(M);base=-Y[ni,5]
    def network(x):
        vc=unpacked(x);rhs=base.copy();rhs[indices]+=vc/ZG
        U=inverse@rhs;ig=(vc-U[indices])/ZG
        return vc,U,ig
    def fun(x):
        vc,U,ig=network(x);power=vc*np.conj(ig)
        return np.r_[power.real-PS,np.abs(vc)+.05*power.imag-A]
    return fun,network

def main():
    known=[('pre',R/'data/seeds/nf02_1_pre0.csv',0.,False),('stable',R/'data/seeds/nf02_1_post0.csv',0.,True),('saddle',R/'data/seeds/nf02_1_post1.csv',0.,True)]
    for t in ['093','095','098']:
        known.append(('pre'+t,R/f'basin_recalculation/capacitor_probe_predictor/v{t}_tau010_pre_root.csv',{'093':.31740252775324773,'095':TARGET_B,'098':.9150354507507668}[t],False))
    validation=[]
    for name,path,B,post in known:
        x=np.loadtxt(path,delimiter=',');f,_=system(B,post)
        er=float(np.max(np.abs(f(packed(vc_from_full(x))))));assert er<1e-8,(name,er)
        validation.append(dict(root=name,max_steady_residual=er))
    pd.DataFrame(validation).to_csv(D/'reduced_validation.csv',index=False)
    stable=packed(vc_from_full(np.loadtxt(known[1][1],delimiter=',')))
    saddle=packed(vc_from_full(np.loadtxt(known[2][1],delimiter=',')))
    rng=np.random.default_rng(20261030)
    starts=[stable,saddle]
    # Coherent remote phases over a complete period, plus independent phases.
    for angle in np.linspace(-np.pi,np.pi,49)[:-1]:
        for mag in [.65,.85,1.,1.15,1.4]:
            vc=unpacked(stable).copy();vc[1:]=mag*np.exp(1j*(angle+np.array([-.012,0,.017])))
            starts.append(packed(vc))
    for i in range(400):
        v=rng.uniform(.4,1.5,4)*np.exp(1j*rng.uniform(-np.pi,np.pi,4));starts.append(packed(v))
    found=[];attempts=[];f,network=system(TARGET_B)
    for i,x in enumerate(starts):
        result=root(f,x,method='hybr',options={'xtol':1e-10,'maxfev':1200})
        er=float(max(abs(f(result.x))));accepted=er<1e-8
        if accepted and not any(max(abs(result.x-y))<1e-5 for y in found):
            found.append(result.x);np.savetxt(D/f'root_{len(found):02d}_vc.csv',result.x,delimiter=',')
        attempts.append(dict(seed=i,residual=er,accepted=accepted,solver_success=bool(result.success)))
    pd.DataFrame(attempts).to_csv(D/'multistart.csv',index=False)
    # Follow both known equilibrium branches and retain roots as extra seeds.
    continuation=[]
    for name,start in [('stable',stable),('saddle',saddle)]:
        x=start.copy()
        for B in np.linspace(0,.00428,108):
            fun,net=system(B);result=root(fun,x,tol=1e-11);er=float(max(abs(fun(result.x))))
            continuation.append(dict(branch=name,B_pu=B,residual=er,mean_remote_phase=float(np.angle(unpacked(result.x)[1:]).mean())))
            if er>1e-8:break
            x=result.x
        np.savetxt(D/f'last_{name}_vc.csv',x,delimiter=',')
    pd.DataFrame(continuation).to_csv(D/'continuation.csv',index=False)
    summary=dict(B_pu=TARGET_B,seed_count=len(starts),distinct_reduced_roots=len(found),minimum_multistart_residual=min(a['residual'] for a in attempts),
                 maximum_validation_residual=max(v['max_steady_residual'] for v in validation),nonexistence_proven=False)
    (D/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
