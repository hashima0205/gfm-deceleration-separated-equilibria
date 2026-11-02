"""Fold check and figure for the requested compensated equilibrium basin."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import root
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from search_capacitor_equilibria import R,D,system,TARGET_B

def jacobian(fun,x):
    J=np.empty((len(x),len(x)))
    for j in range(len(x)):
        h=1e-6*max(1,abs(x[j]));a=x.copy();b=x.copy();a[j]+=h;b[j]-=h
        J[:,j]=(fun(a)-fun(b))/(2*h)
    return J

def main():
    F=R/'figures/02_voltage_adjusted';F.mkdir(exist_ok=True,parents=True)
    c=pd.read_csv(D/'continuation.csv');attempts=pd.read_csv(D/'multistart.csv')
    xs=np.loadtxt(D/'last_stable_vc.csv',delimiter=',');xu=np.loadtxt(D/'last_saddle_vc.csv',delimiter=',')
    x0=(xs+xu)/2;B0=c.B_pu.max();f,_=system(B0)
    _,_,vh=np.linalg.svd(jacobian(f,x0));v0=vh[-1]
    def fold_eq(z):
        x=z[:8];B=z[8];v=z[9:];fun,_=system(B)
        return np.r_[fun(x),jacobian(fun,x)@v,np.dot(v,v)-1]
    result=root(fold_eq,np.r_[x0,B0,v0],tol=1e-10,options={'maxfev':3000})
    residual=float(np.max(abs(fold_eq(result.x))))
    assert residual<1e-7,(result.message,residual)
    foldB=float(result.x[8]);assert .004<foldB<.005
    J=jacobian(system(foldB)[0],result.x[:8])
    uu,ss,vv=np.linalg.svd(J);left=uu[:,-1];right=vv[-1];h=1e-5
    fp=system(foldB+h)[0](result.x[:8]);fm=system(foldB-h)[0](result.x[:8])
    fun=system(foldB)[0];xx=result.x[:8]
    transverse=float(left@((fp-fm)/(2*h)))
    curvature=float(left@((fun(xx+h*right)-2*fun(xx)+fun(xx-h*right))/h**2))
    assert ss[-1]<1e-7 and ss[-2]>.1 and abs(transverse)>1e-4 and abs(curvature)>1e-4
    summary=json.loads((D/'summary.json').read_text())
    summary.update(fold_smallest_singular_value=float(ss[-1]),fold_next_singular_value=float(ss[-2]),
                   fold_transversality=transverse,fold_curvature=curvature)

    summary.update(local_fold_B_pu=foldB,fold_augmented_residual=residual,basin_computed=False,
                   basin_status='No verified post-fault stable equilibrium at requested B; cannot label a recovery basin.')
    (D/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    np.savetxt(D/'fold_augmented_solution.csv',result.x,delimiter=',')
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for label,title in [('stable','Branch from original stable root'),('saddle','Branch from original saddle')]:
        d=c[c.branch==label];axes[0].plot(d.B_pu,d.mean_remote_phase,label=title)
    delta=np.angle(result.x[:4]+1j*result.x[4:8])[1:].mean()
    axes[0].plot(foldB,delta,'ko',ms=4)
    axes[0].set(xlabel='Added shunt susceptance B [pu]',ylabel='Mean internal phase / Bus6 [rad]',title=f'Local steady-state fold: B={foldB:.7f}')
    axes[0].legend(fontsize=7)
    axes[1].hist(np.log10(attempts.residual.clip(lower=1e-14)),bins=35,color='#bdbdbd',edgecolor='white')
    axes[1].axvline(-8,color='#0072b2',ls='--',label='Root acceptance: residual < 1e-8')
    axes[1].set(xlabel='log10(maximum steady-equation residual)',ylabel='Number of initial guesses',title=f'Requested B={TARGET_B:.5f}: 0 roots / {len(attempts)} guesses')
    axes[1].legend(fontsize=7)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle('Voltage-adjusted setting: equilibrium search preceding the basin calculation')
    for ext in ['png','pdf']:fig.savefig(F/f'03_postfault_equilibrium_search.{ext}',dpi=180,bbox_inches='tight')
    (F/'equilibrium_search_status.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
