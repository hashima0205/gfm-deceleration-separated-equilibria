"""Exact steady-state reduction (not a dynamic approximation) of the ZIP=1 case."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import root, minimize_scalar, minimize

HERE=Path(__file__).resolve().parents[2]
OUT=HERE/'output/equilibrium'
BUSES=[1,2,3,7]
LINES=[(1,4,2.14),(4,6,3.43),(5,6,.68),(3,5,2.21),(2,3,1.91),(2,5,2.14),(4,5,1.8),(2,8,2.17),(2,7,6.39),(7,8,1.69),(5,8,8.96)]
def network(d,scale=1,opened=True):
    y=np.zeros((12,12),complex)
    for k,(a,b,m) in enumerate(LINES,1):
        if k==7 and opened:continue
        a-=1;b-=1;v=m*np.exp(-1j*np.pi*75/180)
        y[a,a]+=v;y[b,b]+=v;y[a,b]-=v;y[b,a]-=v
    y[4,4]+=d.load_conductance.iloc[0]*scale
    for k,b in enumerate(BUSES):
        a=b-1;j=8+k;v=1/(d.Rg.iloc[k]+1j*d.Lg.iloc[k])
        y[a,a]+=v;y[j,j]+=v;y[a,j]-=v;y[j,a]-=v
    free=[0,1,2,3,4,6,7];fixed=[5,8,9,10,11]
    h=-np.linalg.solve(y[np.ix_(free,free)],y[np.ix_(free,fixed)])
    red=y[np.ix_(fixed,fixed)]+y[np.ix_(fixed,free)]@h
    return red[1:,1:],red[1:,0],h

def equations(x,d,A,c):
    n=len(d);e=x[n:]*np.exp(1j*x[:n]);s=e*np.conj(A@e+c)
    return np.r_[s.real-d.Pset.to_numpy(),x[n:]-d.Vset.to_numpy()+d.Kq.to_numpy()*(s.imag-d.Qset.to_numpy())]

def certificate(d,A,c):
    # VC integrator equilibrium implies Vc_q=0 and Vc_d=Vout.
    # CC integrator equilibrium and the limiter imply |If| <= Imax.
    # Q=Im(E conj(If))+C|E|^2. This bounds all positive-voltage equilibria.
    a=d.Vset.to_numpy()+d.Kq.to_numpy()*d.Qset.to_numpy();k=d.Kq.to_numpy();cap=d.C.to_numpy();imax=1.2
    lower=2*a/(1+k*imax+np.sqrt((1+k*imax)**2+4*k*cap*a))
    upper=2*a/(1-k*imax+np.sqrt((1-k*imax)**2+4*k*cap*a))
    G=(A+A.conj().T)/2
    mineig=np.linalg.eigvalsh(G)[0]
    def bound(lam):
        H=G-lam*np.eye(len(d))
        return lam*np.sum(lower**2)-.25*np.vdot(c,np.linalg.solve(H,c)).real
    opt=minimize_scalar(lambda q:-bound(q),bounds=(0,mineig*(1-1e-9)),method='bounded',options={'xatol':1e-14})
    lam=opt.x;H=G-lam*np.eye(len(d));v=np.linalg.solve(H,c)
    # Separately exclude negative Vout: from If=(A+jC)E+c and |If|<=Imax.
    Z=np.linalg.inv(A+1j*np.diag(cap));abs_upper=abs(Z@c)+imax*np.sum(abs(Z),axis=1)
    required_negative=(1-k*imax+np.sqrt((1-k*imax)**2+4*k*cap*a))/(2*k*cap)
    return dict(lower_V=lower.tolist(),upper_V=upper.tolist(),lambda_value=float(lam),
        lambda_min_G=float(mineig),lambda_min_H=float(np.linalg.eigvalsh(H)[0]),
        solve_residual=float(np.max(abs(H@v-c))),power_lower_bound=float(bound(lam)),
        total_Pset=float(d.Pset.sum()),gap=float(bound(lam)-d.Pset.sum()),
        arbitrary_sign_voltage_upper=abs_upper.tolist(),negative_Vout_requires=required_negative.tolist(),
        negative_Vout_excluded=bool(np.all(abs_upper<required_negative)))

def main():
    rng=np.random.default_rng(48273);rows=[];certs={}
    for load in [.8,.9,1.]:
        d=pd.read_csv(OUT/f'initial_L{load}.csv');A,c,_=network(d)
        r=d.iloc[1:].reset_index(drop=True);Ar=A[1:,1:];cr=c[1:]
        certs[str(load)]=certificate(r,Ar,cr)
        solutions=[];best=1e9
        for trial in range(400):
            guess=np.r_[rng.uniform(-np.pi,np.pi,3),rng.uniform(.85,1.1,3)]
            res=root(equations,guess,args=(r,Ar,cr),tol=1e-10)
            err=max(abs(equations(res.x,r,Ar,cr)));best=min(best,err)
            if err>1e-8 or min(res.x[3:])<=0:continue
            x=res.x.copy();x[:3]=np.angle(np.exp(1j*x[:3]))
            if any(np.linalg.norm(np.exp(1j*x[:3])-np.exp(1j*q[:3]))<1e-5 for q in solutions):continue
            solutions.append(x)
        for i,x in enumerate(solutions):
            e=x[3:]*np.exp(1j*x[:3]);ig=Ar@e+cr;curr=abs(ig+1j*r.C.to_numpy()*e)
            rows.append(dict(load=load,branch=i,residual=max(abs(equations(x,r,Ar,cr))),max_If=max(curr),**{f'angle{b}':x[j] for j,b in enumerate([2,3,7])},**{f'Vc{b}':x[3+j] for j,b in enumerate([2,3,7])}))
        print(load,'roots',len(solutions),'best',best,'certificate',certs[str(load)],flush=True)
    pd.DataFrame(rows).to_csv(OUT/'reduced_roots.csv',index=False)
    (OUT/'nonexistence_bounds.json').write_text(json.dumps(certs,indent=2))
    # Follow the normal equilibrium with the MAIN CASE controller settings fixed.
    d=pd.read_csv(OUT/'initial_L1.0.csv');r=d.iloc[1:].reset_index(drop=True)
    x=np.r_[np.zeros(3),np.ones(3)];rows=[]
    for scale in np.arange(.5,1.0001,.001):
        A,c,_=network(d,scale);res=root(equations,x,args=(r,A[1:,1:],c[1:]),tol=1e-11)
        err=max(abs(equations(res.x,r,A[1:,1:],c[1:])));rows.append(dict(load_conductance_scale=scale,residual=err,success=err<1e-8,angle2=res.x[0],Vc2=res.x[3]))
        if err<1e-8:x=res.x
    pd.DataFrame(rows).to_csv(OUT/'fixed_settings_continuation.csv',index=False)

if __name__=='__main__':main()
