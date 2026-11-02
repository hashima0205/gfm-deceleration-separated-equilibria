"""Validate and plot B=1 basins for taup=0.5 and 1.0 on identical initial states."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from search_capacitor_equilibria import system,vc_from_full,packed

R=Path(__file__).resolve().parents[1];D=R/'basin_recalculation/B100_tau_comparison';F=R/'figures/02_voltage_adjusted';F.mkdir(exist_ok=True,parents=True)
K=.6812730259441179;B=1.;TAUS=[.5,1.]
S=pd.read_csv(R/'results/conditional_basins/data/states.csv');dyn=S.loc[S.mass!=0,'index'].to_numpy(int)-1
ix=lambda b,s:int(S.loc[(S.bus==b)&(S.symbol==f'inv{b}₊droop₊{s}'),'index'].iloc[0])-1
delta=np.array([ix(b,'δ') for b in [1,2,3,7]]);pf=np.array([ix(b,'Pfilt') for b in [1,2,3,7]])
kp=np.full(3,7.5398223686155035);cols=[f'x{i}' for i in range(1,67)]
xs=np.loadtxt(D/'stable.csv',delimiter=',');xu=np.loadtxt(D/'saddle.csv',delimiter=',');xp=np.loadtxt(D/'pre.csv',delimiter=',')
phase=lambda x:float(x[delta[1:]].mean());gap=phase(xu)-phase(xs)
eq=pd.read_csv(D/'equilibria.csv');events=pd.read_csv(D/'events.csv');grid=pd.read_csv(D/'grid.csv')
assert len(eq)==6 and len(events)==2 and len(grid)==83
assert eq.residual.max()<1e-8 and eq.loc[eq.root!='saddle','max_real_eigen'].lt(0).all()
assert eq.loc[eq.root=='saddle','unstable_modes'].eq(1).all()
root_checks=[]
for name,x,post in [('pre',xp,False),('stable',xs,True),('saddle',xu,True)]:
 f,_=system(B,post=post,k=K);er=float(max(abs(f(packed(vc_from_full(x))))));assert er<1e-8
 root_checks.append(dict(root=name,independent_residual=er))
sets={};initials={};tracks={};faults={};reports={};section_checks=[]
for tau in TAUS:
 suffix=f'{int(round(100*tau)):03d}';tag=f'B100_tau{suffix}'
 pinit=pd.read_csv(D/f'basin_tau{suffix}_initial.csv',float_precision='round_trip').set_index('case')
 samples=pd.read_csv(D/f'basin_tau{suffix}.csv');history=pd.read_csv(D/f'basin_tau{suffix}_history.csv')
 audit=pd.read_csv(D/f'basin_tau{suffix}_audit.csv');assert len(samples)==83 and samples['case'].is_unique
 assert samples.max_algebraic_residual.max()<1e-7 and audit.algebraic_residual.max()<1e-8 and audit.differential_change.eq(0).all()
 assert np.allclose(samples.k,K,rtol=0,atol=1e-14) and samples.B_pu.eq(B).all() and samples.taup.eq(tau).all()
 rowevent=events[np.isclose(events.taup,tau)].iloc[0]
 for stage in ['pre','fault','post']:
  actual=np.loadtxt(D/f'{tag}_{stage}_parameters.csv',delimiter=',')
  old=np.loadtxt(R/f'basin_recalculation/fast_tau_search/runs/fast_12_{stage}_parameters.csv',delimiter=',')
  changed=np.flatnonzero(abs(actual-old)>1e-13)
  assert len(changed)==4 and np.allclose(actual[changed]/old[changed],tau/.1,rtol=0,atol=1e-12)
 categories=[]
 for row in samples.itertuples(index=False):
  x=pinit.loc[row.case,cols].to_numpy(float);expected=xs[dyn]+row.alpha*(xu[dyn]-xs[dyn])
  for j,pp in enumerate(pf[1:]):expected[np.flatnonzero(dyn==pp)[0]]-=row.df_offset_Hz*2*np.pi/kp[j]
  assert np.allclose(x[dyn],expected,rtol=0,atol=1e-12)
  assert abs(phase(x)-row.theta_infinite_bus_rad)<1e-10 and abs(phase(x)-phase(xs)-row.alpha*gap)<1e-10
  h=history[history['case']==row.case];assert h.horizon_s.tolist()==[30,60,150,300][:len(h)]
  assert h.classification.iloc[:-1].isin(['unresolved','other']).all() and h.horizon_s.iloc[-1]==row.horizon_s
  if row.alpha in [0.,1.] and row.df_offset_Hz==0:
   assert np.max(abs(x-(xs if row.alpha==0 else xu)))<1e-10
  end=np.loadtxt(D/f'{row.case}_end.csv',delimiter=',')
  if row.classification.startswith('recovery'):
   assert max(abs(np.angle(np.exp(1j*(end[delta]-xs[delta])))))<.002
   category='same' if all(getattr(row,f'branch_b{b}')==0 for b in [1,2,3,7]) else 'shift'
  else:category={'negative_running':'negative','positive_running':'positive','saddle_initial':'saddle'}.get(row.classification,'unknown')
  categories.append(category)
 samples['category']=categories;sets[tau]=samples;initials[tau]=pinit
 track=pd.read_csv(D/f'{tag}.csv');fault=pd.read_csv(D/f'{tag}_fault.csv')
 assert track.t_post_s.iloc[0]==0 and track.t_post_s.iloc[-1]==150 and np.isfinite(track.to_numpy()).all()
 assert np.allclose(fault[['theta_infinite_bus_rad','df_mean_Hz']].iloc[-1],track[['theta_infinite_bus_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-8)
 assert np.isclose(fault.t_absolute_s.iloc[0],1.) and np.isclose(fault.t_absolute_s.iloc[-1],1.15)
 assert abs((track.theta_infinite_bus_rad.iloc[-1]-track.theta_infinite_bus_rad.iloc[0])/(2*np.pi)-rowevent.remote_turns)<1e-10
 tracks[tau]=track;faults[tau]=fault
 xc=np.loadtxt(D/f'{tag}_clear.csv',delimiter=',');alpha=(phase(xc)-phase(xs))/gap
 section=xs[dyn]+alpha*(xu[dyn]-xs[dyn]);df=track.df_mean_Hz.iloc[0]
 for j,pp in enumerate(pf[1:]):section[np.flatnonzero(dyn==pp)[0]]-=df*2*np.pi/kp[j]
 section_checks.append(dict(taup=tau,clearing_alpha=alpha,clearing_df_Hz=df,unequal_differential_states=int(np.sum(abs(xc[dyn]-section)>1e-8))))
 reports[str(tau)]=dict(points=83,counts=samples.category.value_counts().to_dict(),horizons={str(k):int(v) for k,v in samples.horizon_s.value_counts().sort_index().items()},
                       max_algebraic_residual=float(samples.max_algebraic_residual.max()))
 samples.to_csv(F/f'B100_basin_tau{suffix}_points.csv',index=False)
common_error=float(np.max(abs(initials[.5][cols].to_numpy()-initials[1.][cols].to_numpy())))
assert common_error<1e-12
xy=sets[.5][['theta_infinite_bus_rad','df_Hz']].to_numpy()
assert np.allclose(xy,sets[1.][['theta_infinite_bus_rad','df_Hz']],rtol=0,atol=1e-12)
tri=Delaunay(np.c_[xy[:,0]/2.,xy[:,1]/.2]).simplices
colors={'same':'#9ecae1','shift':'#fdbb84','negative':'#d73027','positive':'#9970ab','unknown':'#d9d9d9','saddle':'#bdbdbd'}
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(F/f'{name}.{ext}',dpi=190,bbox_inches='tight')
 plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8.5))
xglobal=(phase(xs)-.42,phase(xs)+1.62);xlocal=(min(phase(xs),phase(xu))-.55*abs(gap),max(phase(xs),phase(xu))+.55*abs(gap))
ymin=-.16
for tau in TAUS:
 t=tracks[tau];visible=t[t.theta_infinite_bus_rad.between(*xglobal)]
 if len(visible):ymin=min(ymin,float(visible.df_mean_Hz.min())-.005)
for col,tau in enumerate(TAUS):
 s=sets[tau];labels=s.category.to_numpy();track=tracks[tau];fault=faults[tau]
 for row,local in enumerate([False,True]):
  ax=axes[row,col]
  for face in tri:
   cc=set(labels[face]);same=len(cc)==1 and 'saddle' not in cc
   ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(cc))] if same else '#eeeeee',edgecolors='none' if same else '#aaaaaa',hatch=None if same else '////',lw=.2,alpha=.7))
  ax.scatter(xy[:,0],xy[:,1],c=[colors[c] for c in labels],s=11,edgecolors='#444444',lw=.25,zorder=3)
  ax.plot(track.theta_infinite_bus_rad,track.df_mean_Hz,color='black',lw=1.1,zorder=4)
  ax.plot(fault.theta_infinite_bus_rad,fault.df_mean_Hz,color='#762a83',lw=1.7,zorder=5)
  ax.plot(track.theta_infinite_bus_rad.iloc[0],track.df_mean_Hz.iloc[0],'o',mfc='white',mec='black',ms=4,zorder=6)
  ax.plot(phase(xs),0,'*',mfc='white',mec='black',ms=11,zorder=7);ax.plot(phase(xu),0,'D',mfc='none',mec='black',ms=5,zorder=7)
  ax.set(xlim=xlocal if local else xglobal,ylim=(-.011,.011) if local else (ymin,.06),xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',
         title=rf'$\tau_p={tau:g}$ s: '+('equilibrium / saddle detail' if local else 'wide view and event projection'))
  ax.grid(alpha=.15)
legends=[Patch(facecolor=colors[k],label=v) for k,v in [('same','Recovery: same branch'),('shift','Recovery: changed branch'),('negative','Negative running'),('unknown','Unresolved')]]
legends += [Patch(facecolor='#eeeeee',edgecolor='#aaaaaa',hatch='////',label='Mixed sampled outcomes'),Line2D([],[],color='black',label='Post-fault projection'),Line2D([],[],color='#762a83',label='Fault-on projection'),Line2D([],[],marker='*',mfc='white',mec='black',ls='',ms=9,label='Stable root'),Line2D([],[],marker='D',mfc='none',mec='black',ls='',label='Saddle on section')]
fig.legend(handles=legends,loc='lower center',ncol=5,fontsize=8,bbox_to_anchor=(.5,.025))
fig.suptitle(f'B=1.0 pu, k={K:.8f}, pre-fault V5={events.pre_V5_pu.iloc[0]:.4f} pu; identical 83 initial states for both taus')
fig.text(.5,.005,'Chord section contains both equilibria. Actual trajectories are projections. Classification at 30/60/150/300 s; hatching is not a verified boundary.',ha='center',fontsize=7)
fig.tight_layout(rect=(0,.12,1,.95));save(fig,'07_basins_B100_tau050_tau100')
fig,axes=plt.subplots(1,3,figsize=(12,3.7),constrained_layout=True)
for tau in TAUS:
 t=tracks[tau];f=faults[tau];label=rf'$\tau_p={tau:g}$ s';time=np.r_[f.t_absolute_s.to_numpy()-1.15,t.t_post_s.to_numpy()]
 color={.5:'#0072b2',1.:'#d55e00'}[tau]
 axes[0].plot(t.theta_infinite_bus_rad,t.df_mean_Hz,label=label,color=color);axes[0].plot(f.theta_infinite_bus_rad,f.df_mean_Hz,ls='--',color=color)
 axes[1].plot(time,np.r_[f.theta_infinite_bus_rad,t.theta_infinite_bus_rad],label=label,color=color)
 axes[2].plot(time,np.r_[f.df_mean_Hz,t.df_mean_Hz],label=label,color=color)
axes[0].plot(phase(xs),0,'*',ms=10,mfc='white',mec='black');axes[0].set(xlabel='Mean internal phase / Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title='Phase-plane projection')
axes[1].set(xlabel='Time after clearing [s]',ylabel='Mean internal phase / Bus6 [rad]',title='Internal phase',xlim=(-.15,150));axes[1].axhline(phase(xs),color='black',ls=':',lw=.7)
axes[2].set(xlabel='Time after clearing [s]',ylabel='Mean frequency deviation [Hz]',title='Frequency',xlim=(-.15,150))
for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle(f'B=1.0 pu, k={K:.8f}: same 0.15-s fault and line opening');save(fig,'08_trajectories_B100_tau050_tau100')
report=dict(B_pu=B,k=K,pre_V5_pu=float(events.pre_V5_pu.iloc[0]),post_eq_V5_pu=float(events.post_eq_V5_pu.iloc[0]),
            root_checks=root_checks,max_full_root_residual=float(eq.residual.max()),between_tau_initial_state_difference=common_error,
            basins=reports,clearing_section_checks=section_checks,event_classifications={str(r.taup):r.classification for r in events.itertuples(index=False)})
(F/'B100_tau_comparison_validation.json').write_text(json.dumps(report,indent=2)+'\n');events.to_csv(F/'B100_tau_comparison_events.csv',index=False);eq.to_csv(F/'B100_tau_comparison_equilibria.csv',index=False)
print(events.to_string(index=False));print(json.dumps(report,indent=2));print('PASS: roots, parameters, identical full-state grids, staged histories, endpoints, and event continuity')
