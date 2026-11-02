"""Steady-state k scan at the unchanged capacitor/load settings of case 02."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import root
from search_capacitor_equilibria import R,system,TARGET_B,packed,unpacked,vc_from_full
D=R/'basin_recalculation/capacitor_k_scan';D.mkdir(exist_ok=True,parents=True)
base=vc_from_full(np.loadtxt(R/'data/seeds/nf02_1_post0.csv',delimiter=','))
rows=[];found_by_k={}
for k in [1.0299,1.,.95,.9,.85,.8,.75,.7,.6,.5]:
    fun,net=system(TARGET_B,k=k);starts=[packed(base)]
    for phase in np.linspace(-np.pi,np.pi,37)[:-1]:
        v=base.copy();v[1:]=np.exp(1j*(phase+np.array([-.012,0,.017])));starts.append(packed(v))
    roots=[]
    for seed in starts:
        sol=root(fun,seed,tol=1e-11);er=max(abs(fun(sol.x)))
        if er<1e-8 and not any(max(abs(sol.x-x))<1e-5 for x in roots):roots.append(sol.x)
    roots.sort(key=lambda x:np.angle(unpacked(x)[1:]).mean(),reverse=True)
    found_by_k[str(k)]=len(roots)
    for i,x in enumerate(roots):
        vc,U,ig=net(x);phase=np.angle(vc[1:]).mean()
        np.savetxt(D/f'k{int(round(k*10000)):05d}_root{i+1}_vc.csv',x,delimiter=',')
        rows.append(dict(k=k,root=i+1,mean_phase=phase,V5=abs(U[4]),residual=max(abs(fun(x)))))
    print('k',k,'roots',len(roots),[(round(np.angle(unpacked(x)[1:]).mean(),4),round(abs(net(x)[1][4]),4)) for x in roots],flush=True)
pd.DataFrame(rows).to_csv(D/'roots.csv',index=False)
(D/'summary.json').write_text(json.dumps(found_by_k,indent=2)+'\n')
