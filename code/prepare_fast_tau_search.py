"""Candidate equilibria: approach the k fold at fixed B, then vary compensation."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import root
from search_capacitor_equilibria import R,system,TARGET_B,packed,vc_from_full
from plot_capacitor_equilibrium_search import jacobian
D=R/'basin_recalculation/fast_tau_search';D.mkdir(exist_ok=True,parents=True)
start=np.loadtxt(R/'basin_recalculation/capacitor_k_scan/critical_k_solution.csv',delimiter=',')
preseed=packed(vc_from_full(np.loadtxt(R/'basin_recalculation/capacitor_k_strength/k080_pre.csv',delimiter=',')))
rows=[];folds=[]
for B in [TARGET_B,.8,1.,1.2]:
    def aug(z):
        x=z[:8];k=z[8];v=z[9:];fun=system(B,k=k)[0]
        return np.r_[fun(x),jacobian(fun,x)@v,v@v-1]
    r=root(aug,start,tol=1e-10,options={'maxfev':4000});er=max(abs(aug(r.x)))
    assert er<1e-7,(B,r.message,er)
    start=r.x;critical=r.x[8];folds.append(dict(B_pu=B,critical_k=critical,residual=er))
    kvals=[.804,.806,.807,.80725,.80730,.80731] if B==TARGET_B else [critical*.99,critical*.999,critical*.9999]
    for k in kvals:
        f,net=system(B,post=False,k=k);p=root(f,preseed,tol=1e-11);assert max(abs(f(p.x)))<1e-8
        V=abs(net(p.x)[1][4]);eligible=.95<=V<=1.05
        rows.append(dict(case=f'fast_{len(rows)+1:02d}',B_pu=B,k=k,critical_k=critical,k_margin_fraction=(critical-k)/critical,pre_V5_pu=V,eligible_voltage=eligible))
        print(rows[-1],flush=True)
for row in rows:
    row['fault_G_pu']=20.;row['fault_duration_s']=.15
originals={r['case']:r.copy() for r in rows}
for name in ['fast_06','fast_15']:
    for fg in [2.,5.,10.,40.,80.]:
        row=originals[name].copy();row['case']=f'fast_{len(rows)+1:02d}';row['fault_G_pu']=fg;rows.append(row)
pd.DataFrame(rows).to_csv(D/'candidates.csv',index=False);pd.DataFrame(folds).to_csv(D/'folds.csv',index=False)
print('prepared',len(rows),'candidates')
