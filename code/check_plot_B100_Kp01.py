"""Audit mp=1% and compare gain effects without confusing frequency/state coordinates."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

R=Path(__file__).resolve().parents[1]
BASE=R/'basin_recalculation/B100_tau_comparison';D=R/'basin_recalculation/B100_Kp01_tau_comparison';F=R/'figures/02_voltage_adjusted'
def plot_events():
    events=pd.read_csv(D/'events.csv');eq=pd.read_csv(D/'equilibria.csv')
    baseline=pd.read_csv(D/'baseline_validation.csv')
    assert len(events)==4 and len(eq)==12 and eq.residual.max()<1e-8
    assert eq.loc[eq.root!='saddle','max_real_eigen'].lt(0).all()
    assert len(baseline)==2 and baseline.post150_max_full_state_difference.max()<2e-7
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(10,6.5),constrained_layout=True)
    for col,tau in enumerate([.5,1.]):
        for factor,color in [(1.,'#d55e00'),(.5,'#0072b2')]:
            event=events[np.isclose(events.taup,tau)&np.isclose(events.gain_factor,factor)].iloc[0]
            tr=pd.read_csv(D/f'{event["case"]}.csv');assert tr.t_post_s.iloc[-1]==300
            label=f'mp={2*factor:g}%';axes[0,col].plot(tr.t_post_s,tr.theta_infinite_bus_rad,label=label,color=color)
            axes[1,col].plot(tr.t_post_s,tr.df_mean_Hz,label=label,color=color)
        axes[0,col].set(title=rf'$\tau_p={tau:g}$ s',ylabel='Mean internal phase / Bus6 [rad]',xlabel='Time after clearing [s]')
        axes[1,col].set(ylabel='Mean frequency deviation [Hz]',xlabel='Time after clearing [s]')
        for ax in axes[:,col]:ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('B=1, same pre-fault equilibrium and 0.15-s fault: mp=1% versus 2%')
    for ext in ['png','pdf']:fig.savefig(F/f'12_trajectories_B100_Kp_comparison.{ext}',dpi=190,bbox_inches='tight')
    plt.close(fig)
    events.to_csv(F/'B100_Kp_comparison_events.csv',index=False)
    return events

if '--events-only' in __import__('sys').argv:
    print(plot_events().to_string(index=False));print('PASS: event equilibria, common horizon, and replay validation')
    raise SystemExit(0)

# Reuse the original full audit and grid plots with explicit gain-dependent edits.
source=(R/'code/check_plot_B100_tau_comparison.py').read_text()
source=source.replace("D=R/'basin_recalculation/B100_tau_comparison'","D=R/'basin_recalculation/B100_Kp01_tau_comparison'")
source=source.replace('kp=np.full(3,7.5398223686155035)','kp=np.full(3,3.7699111843077517)')
source=source.replace("eq=pd.read_csv(D/'equilibria.csv');events=pd.read_csv(D/'events.csv');grid=pd.read_csv(D/'grid.csv')", "eq_all=pd.read_csv(D/'equilibria.csv');events_all=pd.read_csv(D/'events.csv');grid=pd.read_csv(D/'grid.csv')\neq=eq_all[np.isclose(eq_all.gain_factor,.5)].copy();events=events_all[np.isclose(events_all.gain_factor,.5)].copy()")
source=source.replace("tag=f'B100_tau{suffix}'","tag=f'B100_Kp01_tau{suffix}'")
old="""  old=np.loadtxt(R/f'basin_recalculation/fast_tau_search/runs/fast_12_{stage}_parameters.csv',delimiter=',')
  changed=np.flatnonzero(abs(actual-old)>1e-13)
  assert len(changed)==4 and np.allclose(actual[changed]/old[changed],tau/.1,rtol=0,atol=1e-12)"""
new="""  old=np.loadtxt(R/f'basin_recalculation/B100_tau_comparison/B100_tau{suffix}_{stage}_parameters.csv',delimiter=',')
  changed=np.flatnonzero(abs(actual-old)>1e-13)
  assert len(changed)==4 and np.allclose(actual[changed]/old[changed],.5,rtol=0,atol=1e-12)
  assert np.allclose(old[changed],7.5398223686155035,rtol=0,atol=1e-12)"""
assert old in source;source=source.replace(old,new)
source=source.replace('track.t_post_s.iloc[-1]==150','track.t_post_s.iloc[-1]==300')
source=source.replace('xlim=(-.15,150)','xlim=(-.15,300)')
for a,b in [('07_basins_B100_tau050_tau100','09_basins_B100_Kp01_tau050_tau100'),('08_trajectories_B100_tau050_tau100','10_trajectories_B100_Kp01_tau050_tau100'),('B100_tau_comparison_validation.json','B100_Kp01_tau_comparison_validation.json'),('B100_tau_comparison_events.csv','B100_Kp01_tau_comparison_events.csv'),('B100_tau_comparison_equilibria.csv','B100_Kp01_tau_comparison_equilibria.csv'),('B100_basin_tau','B100_Kp01_basin_tau')]:source=source.replace(a,b)
source=source.replace("fig.suptitle(f'B=1.0 pu, k=", "fig.suptitle(f'mp=1% on common base; B=1.0 pu, k=")
source=source.replace("fig.suptitle(f'B=1.0 pu, k=", "fig.suptitle(f'mp=1% on common base; B=1.0 pu, k=")
ns={'__file__':str(R/'code/check_plot_B100_tau_comparison.py'),'__name__':'gain01_audit'}
exec(compile(source,str(R/'code/check_plot_B100_tau_comparison.py'),'exec'),ns)
S=ns['S'];dyn=ns['dyn'];pf=ns['pf'];cols=ns['cols'];xs=ns['xs'];xu=ns['xu'];phase=ns['phase'];gap=ns['gap']
newsets=ns['sets'];newinitials=ns['initials'];report=ns['report'];oldsets={};checks={};pair_rows=[]
baseline=pd.read_csv(D/'baseline_validation.csv');assert len(baseline)==2 and baseline.post150_max_full_state_difference.max()<2e-7
all_eq=pd.read_csv(D/'equilibria.csv');assert len(all_eq)==12 and all_eq.residual.max()<1e-8
assert all_eq.loc[all_eq.root!='saddle','max_real_eigen'].lt(0).all()
for name in ['pre','stable','saddle']:assert np.array_equal(np.loadtxt(D/f'{name}.csv',delimiter=','),np.loadtxt(BASE/f'{name}.csv',delimiter=','))
def category(row):
 if row.classification.startswith('recovery'):
  return 'same' if all(getattr(row,f'branch_b{b}')==0 for b in [1,2,3,7]) else 'shift'
 return {'negative_running':'negative','positive_running':'positive','saddle_initial':'saddle'}.get(row.classification,'unknown')
for tau in [.5,1.]:
 suffix=f'{int(round(100*tau)):03d}';old=pd.read_csv(BASE/f'basin_tau{suffix}.csv');old['category']=[category(r) for r in old.itertuples(index=False)];oldsets[tau]=old
 oi=pd.read_csv(BASE/f'basin_tau{suffix}_initial.csv',float_precision='round_trip').set_index('case');ni=newinitials[tau]
 new=newsets[tau];assert len(new)==len(old)==83
 assert np.allclose(new[['theta_infinite_bus_rad','df_Hz']],old[['theta_infinite_bus_rad','df_Hz']],rtol=0,atol=1e-10)
 other=np.setdiff1d(dyn,pf[1:]);assert np.max(abs(ni[cols].to_numpy()[:,other]-oi[cols].to_numpy()[:,other]))<1e-12
 matched=[]
 for j,r in enumerate(new.itertuples(index=False)):
  x=ni.loc[r.case,cols].to_numpy(float)
  chord=xs[pf[1:]]+r.alpha*(xu[pf[1:]]-xs[pf[1:]])
  ox=oi.iloc[j][cols].to_numpy(float)
  assert np.allclose(x[pf[1:]]-chord,2*(ox[pf[1:]]-chord),rtol=0,atol=1e-12)
  candidates=old[(abs(old.alpha-r.alpha)<1e-10)&(abs(old.df_offset_Hz-2*r.df_offset_Hz)<1e-12)]
  if len(candidates):
   assert len(candidates)==1;rr=candidates.iloc[0];err=float(max(abs(x-oi.loc[rr['case'],cols].to_numpy(float))));assert err<1e-10
   matched.append(err);pair_rows.append(dict(taup=tau,mp01_case=r.case,mp02_case=rr['case'],max_full_state_difference=err,mp01_category=r.category,mp02_category=rr.category))
 checks[str(tau)]=dict(mp01_counts=new.category.value_counts().to_dict(),mp02_counts=old.category.value_counts().to_dict(),matched_full_state_points=len(matched),max_matched_full_state_difference=max(matched))
pd.DataFrame(pair_rows).to_csv(F/'B100_Kp01_matched_full_state_comparison.csv',index=False)
# Gain comparison in an invariant state coordinate: Pset-Pf, not gain-dependent frequency.
fig,axes=plt.subplots(2,2,figsize=(11,8))
colors=ns['colors'];KREF=7.5398223686155035
for row,(mp,sets) in enumerate([(.02,oldsets),(.01,newsets)]):
 gain=mp*(2*np.pi*60)
 for col,tau in enumerate([.5,1.]):
  ax=axes[row,col];samples=sets[tau];xy=samples[['theta_infinite_bus_rad','df_Hz']].to_numpy().copy();xy[:,1]*=2*np.pi/gain
  tri=Delaunay(np.c_[xy[:,0]/2.,xy[:,1]/.2]).simplices;labels=samples.category.to_numpy()
  for face in tri:
   cats=set(labels[face]);same=len(cats)==1 and 'saddle' not in cats
   ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(cats))] if same else '#eeeeee',edgecolors='none' if same else '#aaaaaa',hatch=None if same else '////',lw=.2,alpha=.7))
  ax.scatter(xy[:,0],xy[:,1],c=[colors[x] for x in labels],s=12,edgecolors='#444444',lw=.25)
  ax.plot(phase(xs),0,'*',mfc='white',mec='black',ms=11);ax.plot(phase(xu),0,'D',mfc='none',mec='black',ms=5)
  ax.set(xlim=(phase(xu)-.55*abs(gap),phase(xs)+.55*abs(gap)),ylim=(-.01*2*np.pi/KREF,.01*2*np.pi/KREF),
         xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Common Pset - Pf offset [pu]',title=rf'$m_p={100*mp:g}\%$, $\tau_p={tau:g}$ s')
  ax.grid(alpha=.2)
handles=[Patch(facecolor=colors[k],label=v) for k,v in [('same','Recovery: same branch'),('shift','Recovery: changed branch'),('negative','Negative running'),('unknown','Unresolved')]]
handles.append(Patch(facecolor='#eeeeee',edgecolor='#aaaaaa',hatch='////',label='Mixed sampled outcomes'))
fig.legend(handles=handles,loc='lower center',ncol=3,fontsize=8,bbox_to_anchor=(.5,.02))
fig.suptitle('Gain comparison in P-filter state coordinates; identical other differential states on the chord')
fig.text(.5,.005,'Matched state pairs are audited separately. Sampling grids differ in this coordinate; counts are not basin volumes.',ha='center',fontsize=7)
fig.tight_layout(rect=(0,.12,1,.95));ns['save'](fig,'11_basins_B100_Kp_state_comparison')
events=plot_events()
report.update(mp_common=.01,reference_mp_common=.02,baseline_post150_max_full_state_difference=float(baseline.post150_max_full_state_difference.max()),
              gain_comparison=checks,coordinate_convention='Same phase/frequency grid across gains requires different Pf offsets; invariant-state comparison uses Pset-Pf.',
              all_event_outcomes=[dict(case=r.case,mp_common=r.mp_common,taup=r.taup,classification=r.classification,horizon_s=r.horizon_s) for r in events.itertuples(index=False)])
(F/'B100_Kp01_tau_comparison_validation.json').write_text(json.dumps(report,indent=2)+'\n')
events.to_csv(F/'B100_Kp_comparison_events.csv',index=False)
print(json.dumps(report,indent=2));print('PASS: only Kp changes, common equilibria, replayed baseline, phase/frequency grids, and matched full-state comparisons')
