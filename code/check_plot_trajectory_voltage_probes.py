"""Audit and plot trajectory-centered sections and load-conductance probes."""
from pathlib import Path
import argparse
import json
import ast
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from plot_saved_basins import save

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--neighborhood-only',action='store_true')
parser.add_argument('--voltage-directory',type=Path)
args=parser.parse_args()
P=ROOT/'basin_recalculation/trajectory_voltage_probes'
VP=args.voltage_directory or P
F=ROOT/'figures/trajectory_voltage_probes'
F.mkdir(parents=True,exist_ok=True)
samples=pd.read_csv(P/'neighborhood_samples.csv')
anchors=pd.read_csv(P/'anchors.csv',float_precision='round_trip').set_index('anchor_id')
initial=pd.read_csv(P/'neighborhood_initial_states.csv',float_precision='round_trip').set_index('case')
history=pd.read_csv(P/'neighborhood_history.csv')
centers=pd.read_csv(P/'center_validation.csv')
states=pd.read_csv(ROOT/'results/conditional_basins/data/states.csv')
meta=json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text())
held=np.array(meta['held_differential_indices'])-1
idx=lambda b,k:int(states.loc[(states.bus==b)&(states.symbol==f'inv{b}₊droop₊{k}'),'index'].iloc[0])-1
delta=np.array([idx(b,'δ') for b in [2,3,7]])
pf=np.array([idx(b,'Pfilt') for b in [2,3,7]])
kp=pd.read_csv(ROOT/'basin_recalculation/tau010_chord/parameters.csv').set_index('bus').loc[[2,3,7],'Kp'].to_numpy()
cols=[f'x{i}' for i in range(1,67)]
assert len(samples)==100 and samples['case'].is_unique and len(anchors)==4
assert set(anchors.anchor_post_s)=={0,1,3,10}
assert len(centers)==4 and centers.max_full_state_difference.max()<2e-7
assert pd.read_csv(P/'reference_validation.csv').max_full_state_difference.max()<2e-7
for row in samples.itertuples(index=False):
    x=initial.loc[row.case,cols].to_numpy(dtype=float)
    base=anchors.loc[row.anchor_id,cols].to_numpy(dtype=float)
    assert np.array_equal(x[held],base[held])
    assert np.allclose(x[delta]-base[delta],row.angle_offset_rad,rtol=0,atol=1e-12)
    assert np.allclose(x[pf]-base[pf],-row.df_offset_Hz*2*np.pi/kp,rtol=0,atol=1e-12)
    h=history[history['case']==row.case]
    assert h.horizon_s.tolist()==[30,60,150][:len(h)]
    assert h.classification.iloc[:-1].isin(['unresolved','other']).all()
    assert h.horizon_s.iloc[-1]==row.horizon_s and h.classification.iloc[-1]==row.classification
assert samples.max_algebraic_residual.max()<1e-7
assert pd.read_csv(P/'neighborhood_initialization.csv').algebraic_residual.max()<1e-8

def category(row):
    if row.classification.startswith('recovery'):
        return 'same' if all(getattr(row,f'branch_b{b}')==0 for b in [1,2,3,7]) else 'shift'
    return {'negative_running':'negative','positive_running':'positive'}.get(row.classification,'unknown')

samples['category']=[category(r) for r in samples.itertuples(index=False)]
colors={'same':'#9ecae1','shift':'#fdbb84','negative':'#d73027','positive':'#9970ab','unknown':'#d9d9d9'}
reference=pd.read_csv(P/'reference_post160.csv')
faces=[]
for j in range(4):
    for i in range(4):
        a=j*5+i;faces += [[a,a+1,a+6],[a,a+6,a+5]]
faces=np.array(faces)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(9,7))
for ax,(anchor_id,anchor) in zip(axes.flat,anchors.iterrows()):
    d=samples[samples.anchor_id==anchor_id].sort_values(['df_offset_Hz','angle_offset_rad'])
    assert len(d)==25
    xy=d[['theta_infinite_bus_rad','df_Hz']].to_numpy();labels=d.category.to_numpy()
    for face in faces:
        label=set(labels[face])
        ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(label))] if len(label)==1 else '#eeeeee',
                                        edgecolors='none' if len(label)==1 else '#bbbbbb',
                                        linewidths=.3,hatch=None if len(label)==1 else '////',alpha=.7))
    ax.plot(reference.theta_infinite_bus_rad,reference.df_mean_Hz,color='#555555',lw=1,zorder=2)
    ax.scatter(xy[:,0],xy[:,1],c=[colors[c] for c in labels],s=18,edgecolors='#444444',linewidths=.4,zorder=3)
    ax.plot(anchor.theta_infinite_bus_rad,anchor.df_Hz,'o',mfc='black',mec='white',ms=7,zorder=5)
    ax.set(xlim=(xy[:,0].min()-.008,xy[:,0].max()+.008),ylim=(xy[:,1].min()-.003,xy[:,1].max()+.003),
           xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',
           title=f'Anchor: {anchor.anchor_post_s:g} s after clearing')
    ax.grid(alpha=.15)
legend=[Patch(facecolor=colors[k],label=v) for k,v in [('same','Recovery: reference branch'),('shift','Recovery: changed branch'),('negative','Negative running'),('unknown','Unresolved / other')]]
legend += [Line2D([],[],color='#555555',label='Reference trajectory projection'),Line2D([],[],marker='o',color='black',ls='',label='Full-state anchor')]
fig.legend(handles=legend,loc='lower center',ncol=3,fontsize=8,bbox_to_anchor=(.5,.025))
fig.suptitle(r'Trajectory-centered conditional recovery sections, $\tau_p=0.1$ s')
fig.text(.5,.008,'Each panel fixes the other 46 initial differential states at its own anchor; 30/60/150 s finite observation.',ha='center',fontsize=7)
fig.tight_layout(rect=(0,.13,1,.96))
save(fig,F,'tau010_trajectory_neighborhood')

if args.neighborhood_only:
    report=dict(counts={f't={t:g}s/{c}':int(n) for (t,c),n in samples.groupby(['anchor_post_s','category']).size().items()},
                horizons={str(k):int(v) for k,v in samples.horizon_s.value_counts().sort_index().items()},
                max_center_full_state_difference=float(centers.max_full_state_difference.max()))
    (P/'neighborhood_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));print('PASS: trajectory-centered initial states, histories, and centers')
    raise SystemExit(0)

voltage=pd.read_csv(VP/'voltage_cases.csv')
roots=pd.read_csv(VP/'voltage_equilibria.csv')
params=pd.read_csv(VP/'voltage_parameters.csv')
assert len(voltage)==8 and set(voltage.taup)=={.1,1.2}
assert roots.residual.max()<1e-8 and roots.max_real_eigen.max()<0 and roots.unstable_modes.eq(0).all()
assert len(params)==8*3*4 and params.tauq.eq(.02).all()
assert np.allclose(params.Pset,params.bus.map({1:.5,2:.4,3:.4,7:.4}),rtol=0,atol=0)
baseline_errors={}
for name,reference_case in [('load100_tau010','tau010_T015_h150'),('load100_tau120','smoke_tau120')]:
    actual=np.loadtxt(VP/f'{name}_clear.csv',delimiter=',')
    expected=np.loadtxt(ROOT/'checkpoints'/reference_case/'post030_initial.csv',delimiter=',')
    baseline_errors[name+'_clearing']=float(np.max(abs(actual-expected)))
    assert baseline_errors[name+'_clearing']<2e-7
actual=np.loadtxt(VP/'load100_tau010_post_end.csv',delimiter=',')
expected=np.loadtxt(ROOT/'checkpoints/tau010_T015_h150/post150_end.csv',delimiter=',')
baseline_errors['load100_tau010_post150']=float(np.max(abs(actual-expected)))
assert baseline_errors['load100_tau010_post150']<2e-7
reference120=pd.read_csv(ROOT/'results/conditional_basins/data/trajectories/fault_tau1.20.csv.gz')
reference120=reference120.iloc[(reference120.t_post_s-150).abs().argmin()]
actual120=pd.read_csv(VP/'load100_tau120.csv').iloc[-1]
assert np.allclose([actual120.theta_infinite_bus_rad,actual120.df_mean_Hz],
                   [reference120.theta_relative_rad+meta['origin_rad'],reference120.df_mean_Hz],rtol=0,atol=1e-7)
source=(ROOT/'provenance/source_snapshot/studies/04_clearing/07_grid_strength/code/run.jl').read_text()
branches=ast.literal_eval(re.search(r'^const branches=(\[.*\])$',source,re.M).group(1))
def network(x,post):
    def v(b,k):return x[int(states.loc[(states.bus==b)&(states.symbol==k),'index'].iloc[0])-1]
    U={b:1+0j if b==6 else complex(v(b,'busbar₊u_r'),v(b,'busbar₊u_i')) for b in range(1,9)}
    I={b:0j for b in U};loss=0.
    for j,(a,b,y) in enumerate(branches,1):
        z=(meta['k'] if j in [2,3] else 1)/y*np.exp(1j*np.deg2rad(75))
        current=0j if post and j==7 else (U[a]-U[b])/z
        I[a]+=current;I[b]-=current;loss+=abs(current)**2*z.real
    S={b:U[b]*np.conj(I[b]) for b in U}
    assert abs(sum(s.real for s in S.values())-loss)<1e-8
    assert abs(S[4])<1e-8 and abs(S[8])<1e-8 and abs(S[5].imag)<1e-8
    for b,Pset in zip([1,2,3,7],meta['Pset']):
        ig=complex(v(b,f'inv{b}₊vsrc₊filter₊i_g_r'),v(b,f'inv{b}₊vsrc₊filter₊i_g_i'))
        vc=complex(v(b,f'inv{b}₊vsrc₊filter₊V_C_r'),v(b,f'inv{b}₊vsrc₊filter₊V_C_i'))
        assert abs(S[b]-U[b]*np.conj(ig))<1e-8 and abs((vc*np.conj(ig)).real-Pset)<1e-8
    return dict(load_P_pu=-S[5].real,line_loss_P_pu=loss,slack_net_P_pu=S[6].real,load_G_pu=-S[5].real/abs(U[5])**2)
flow_rows=[]
for row in voltage.itertuples(index=False):
    record=row._asdict()
    for stage in ['pre','post']:
        root=np.loadtxt(VP/f'{row.case}_{stage}_root.csv',delimiter=',')
        record.update({stage+'_'+k:v for k,v in network(root,stage=='post').items()})
    assert abs(record['pre_load_G_pu']-record['post_load_G_pu'])<1e-8
    assert abs(record['pre_load_G_pu']/2.2682301152333055-row.load_factor)<1e-8
    flow_rows.append(record)
voltage=pd.DataFrame(flow_rows)
voltage.to_csv(VP/'voltage_cases_with_power_balance.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(9,3.7),constrained_layout=True)
for tau,style,marker in [(.1,'--','s'),(1.2,'-','o')]:
    d=voltage[np.isclose(voltage.taup,tau)].sort_values('pre_V5_pu')
    axes[0].plot(d.pre_V5_pu,d.remote_internal_turns,style,marker=marker,label=rf'$\tau_p={tau:g}$ s')
axes[0].set(xlabel='Pre-fault Bus5 voltage [pu]',ylabel='Remote post-fault angle change [turns]',title='Same dispatch and 0.15-s fault; 150-s observation')
axes[0].legend(fontsize=8)
d=voltage[np.isclose(voltage.taup,.1)].sort_values('load_factor')
axes[1].plot(d.load_factor,d.pre_V5_pu,'o-',label='Pre-fault equilibrium')
axes[1].plot(d.load_factor,d.post_eq_V5_pu,'s--',label='Post-fault stable equilibrium')
axes[1].set(xlabel='Load conductance / original conductance',ylabel='Bus5 voltage [pu]',title='Operating points change with conductance')
axes[1].legend(fontsize=8)
for ax in axes:ax.grid(alpha=.2)
save(fig,F,'voltage_load_probe')
report=dict(neighborhood_counts=samples.groupby(['anchor_post_s','category']).size().to_dict(),
            neighborhood_horizons=samples.horizon_s.value_counts().sort_index().to_dict(),
            max_center_full_state_difference=float(centers.max_full_state_difference.max()),
            max_neighborhood_algebraic_residual=float(samples.max_algebraic_residual.max()),
            voltage_case_count=len(voltage),all_voltage_equilibria_locally_stable=True,
            baseline_full_state_errors=baseline_errors)
# Tuple keys are rendered explicitly for portable JSON.
report['neighborhood_counts']={f't={t:g}s/{c}':int(n) for (t,c),n in report['neighborhood_counts'].items()}
report['neighborhood_horizons']={str(k):int(v) for k,v in report['neighborhood_horizons'].items()}
(P/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
print(voltage[['case','load_factor','taup','pre_V5_pu','post_eq_V5_pu','classification','remote_internal_turns','pre_load_P_pu','post_load_P_pu']].to_string(index=False))
print('PASS: initial-state sections, staged histories, trajectory centers, equilibria, parameters, and power balances')
