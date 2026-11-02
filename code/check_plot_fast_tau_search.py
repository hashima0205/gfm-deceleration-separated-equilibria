"""Verify the structured tau=0.1 search and put one summary in case 02."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from search_capacitor_equilibria import system,vc_from_full,packed,R,G,buses,branches

P=R/'basin_recalculation/fast_tau_search';D=P/'runs';F=R/'figures/02_voltage_adjusted'
candidates=pd.read_csv(P/'candidates.csv');results=pd.read_csv(D/'results.csv');eq=pd.read_csv(D/'equilibria.csv')
assert len(results)==len(candidates) and results['case'].is_unique and len(eq)==2*len(results)
assert eq.residual.max()<1e-8 and eq.max_real_eigen.lt(0).all() and eq.unstable_modes.eq(0).all()
ref=R/'basin_recalculation/capacitor_k_strength';old=R/'basin_recalculation/capacitor_probe_predictor'
refs={stage:np.loadtxt(ref/f'k080_tau010_{stage}_parameters.csv',delimiter=',') for stage in ['pre','fault','post']}
orig=np.loadtxt(old/'v095_tau010_pre_parameters.csv',delimiter=',')
zero=np.loadtxt(R/'checkpoints/tau010_T015_h150/pre_parameters.csv',delimiter=',')
lineidx=np.flatnonzero(abs(orig-refs['pre'])>1e-14);qidx=np.flatnonzero(abs(orig-zero)>1e-14)
assert len(lineidx)==4 and len(qidx)==1
faultidx=np.flatnonzero(refs['fault']==-20.);assert len(faultidx)==1
independent=[];nodal=[]
S=pd.read_csv(R/'results/conditional_basins/data/states.csv')
def val(x,b,s):return x[int(S.loc[(S.bus==b)&(S.symbol==s),'index'].iloc[0])-1]
def check_nodal(x,k,B,post):
 U={b:1+0j if b==6 else complex(val(x,b,'busbar₊u_r'),val(x,b,'busbar₊u_i')) for b in range(1,9)}
 I={b:0j for b in U}
 for j,(a,b,y) in enumerate(branches,1):
  z=(k if j in [2,3] else 1)/y*np.exp(1j*np.deg2rad(75));ii=0j if post and j==7 else (U[a]-U[b])/z
  I[a]+=ii;I[b]-=ii
 target={b:0j for b in U};target[5]=-(G+1j*B)*U[5]
 for b in buses:target[b]=complex(val(x,b,f'inv{b}₊vsrc₊filter₊i_g_r'),val(x,b,f'inv{b}₊vsrc₊filter₊i_g_i'))
 error=max(abs(I[b]-target[b]) for b in U if b!=6);assert error<1e-7
 return float(error)
late=[]
for row in results.itertuples(index=False):
 candidate=candidates[candidates['case']==row.case].iloc[0]
 assert abs(row.pre_V5_pu-candidate.pre_V5_pu)<1e-6 and .95<=row.pre_V5_pu<=1.05
 assert row.horizon_s in [30,60,150,300]
 for stage in ['pre','fault','post']:
  actual=np.loadtxt(D/f'{row.case}_{stage}_parameters.csv',delimiter=',');expected=refs[stage].copy()
  expected[lineidx]*=row.k/.8;expected[qidx]*=row.B_pu/.5561806200995012
  if stage=='fault':expected[faultidx]=-row.fault_G_pu
  assert np.allclose(actual,expected,rtol=0,atol=1e-12)
 for stage in ['pre','post']:
  xx=np.loadtxt(D/f'{row.case}_{stage}_root.csv',delimiter=',');fun,_=system(row.B_pu,post=stage=='post',k=row.k)
  er=float(max(abs(fun(packed(vc_from_full(xx))))));assert er<1e-8;independent.append(er)
  nodal.append(check_nodal(xx,row.k,row.B_pu,stage=='post'))
 for stage in ['clear','end']:
  xx=np.loadtxt(D/f'{row.case}_{stage}.csv',delimiter=',');nodal.append(check_nodal(xx,row.k,row.B_pu,True))
 tr=pd.read_csv(D/f'{row.case}.csv');f=pd.read_csv(D/f'{row.case}_fault.csv')
 assert tr.t_post_s.iloc[0]==0 and tr.t_post_s.iloc[-1]==row.horizon_s and np.isfinite(tr.to_numpy()).all()
 assert np.isclose(f.t_absolute_s.iloc[0],1.) and np.isclose(f.t_absolute_s.iloc[-1],1+row.fault_duration_s)
 assert np.allclose(f[['theta_infinite_bus_rad','df_mean_Hz']].iloc[-1],tr[['theta_infinite_bus_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-8)
 assert abs((tr.theta_infinite_bus_rad.iloc[-1]-tr.theta_infinite_bus_rad.iloc[0])/(2*np.pi)-row.remote_turns)<1e-10
 tail=tr[tr.t_post_s>=.7*row.horizon_s]
 late.append(float(np.polyfit(tail.t_post_s,tail.theta_infinite_bus_rad,1)[0]/(2*np.pi)))
results['late_mean_df_Hz']=late
results.to_csv(F/'fast_tau_search_results.csv',index=False)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
folds=pd.read_csv(P/'folds.csv').sort_values('B_pu')
axes[0].plot(folds.B_pu,folds.critical_k,'k--o',ms=4,label='Tracked equilibrium fold')
fixed=results[np.isclose(results.fault_G_pu,20.) & np.isclose(results.fault_duration_s,.15)]
axes[0].scatter(fixed.B_pu,fixed.k,s=22,facecolors='none',edgecolors='#0072b2',label='Tested stable cases')
axes[0].set(xlabel='Compensation susceptance B [pu]',ylabel='External-grid impedance multiplier k',title='Equilibrium constraint and sampled settings')
nearest=fixed.sort_values('k_margin_fraction').groupby('B_pu',as_index=False).first().sort_values('B_pu')
axes[1].plot(nearest.B_pu,nearest.pre_V5_pu,'o-',label='Pre-fault equilibrium')
axes[1].plot(nearest.B_pu,nearest.post_eq_V5_pu,'s--',label='Post-fault stable equilibrium')
axes[1].set(xlabel='Compensation susceptance B [pu]',ylabel='Bus5 voltage [pu]',title='Nearest tested setting to each fold')
for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle(r'$\tau_p=0.1$ s: structured k / compensation / fault-strength search')
fig.text(.5,-.035,f'{len(results)} cases; outcome counts: '+', '.join(f'{key}={value}' for key,value in results.classification.value_counts().items())+'; 30-300 s finite observation.',ha='center',fontsize=8)
for ext in ['png','pdf']:fig.savefig(F/f'06_fast_tau_search.{ext}',dpi=180,bbox_inches='tight')
negative=results[results.classification=='negative_running'].sort_values('k_margin_fraction',ascending=False)
if len(negative):
 winner=negative.iloc[0]
 fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
 tr=pd.read_csv(D/f'{winner["case"]}.csv');f=pd.read_csv(D/f'{winner["case"]}_fault.csv')
 root_eq=eq[(eq['case']==winner['case'])&(eq.stage=='post')].iloc[0]
 axes[0].plot(tr.theta_infinite_bus_rad,tr.df_mean_Hz,color='black');axes[0].plot(f.theta_infinite_bus_rad,f.df_mean_Hz,color='#762a83')
 axes[0].plot(root_eq.mean_phase_rad,0,'*',ms=10,mfc='white',mec='black')
 axes[0].set(xlabel='Mean internal phase / Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title='Phase-plane projection')
 axes[1].plot(tr.t_post_s,tr.theta_infinite_bus_rad);axes[1].axhline(root_eq.mean_phase_rad,color='black',ls=':',lw=.8)
 axes[1].set(xlabel='Time after clearing [s]',ylabel='Mean internal phase / Bus6 [rad]',title='Internal phase')
 axes[2].plot(tr.t_post_s,tr.df_mean_Hz);axes[2].set(xlabel='Time after clearing [s]',ylabel='Mean frequency deviation [Hz]',title='Frequency')
 for ax in axes:ax.grid(alpha=.2)
 fig.suptitle(f'taup=0.1 s; k={winner.k:.8f}, B={winner.B_pu:.5f}, pre V5={winner.pre_V5_pu:.4f}; stable post-fault equilibrium exists')
 for ext in ['png','pdf']:fig.savefig(F/f'07_fast_tau_example.{ext}',dpi=180,bbox_inches='tight')
report=dict(case_count=len(results),classifications=results.classification.value_counts().to_dict(),all_pre_post_roots_stable=True,
            max_root_residual=float(eq.residual.max()),max_independent_steady_residual=max(independent),max_nodal_residual=max(nodal),
            negative_example=negative.iloc[0]['case'] if len(negative) else None)
(F/'fast_tau_search_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(results.to_string(index=False));print(json.dumps(report,indent=2));print('PASS: full roots, physical parameters, nodal balances, event continuity, finite observation, and voltage targets')
