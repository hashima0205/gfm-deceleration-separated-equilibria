"""Audit constant-B shunt cases and plot only this capacitor study."""
import ast
import json
import re
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parents[1]
D=R/'basin_recalculation/capacitor_probe_predictor'
F=R/'figures/capacitor_probe';F.mkdir(exist_ok=True,parents=True)
T=R/'data/capacitor_probe';T.mkdir(exist_ok=True,parents=True)
cases=pd.read_csv(D/'cases.csv')
eq=pd.read_csv(D/'equilibria.csv')
params=pd.read_csv(D/'parameters.csv')
machines=pd.read_csv(D/'machines.csv')
meta=json.loads((R/'results/conditional_basins/data/definition.json').read_text())
states=pd.read_csv(R/'results/conditional_basins/data/states.csv')
source=(R/'provenance/source_snapshot/studies/04_clearing/07_grid_strength/code/run.jl').read_text()
branches=ast.literal_eval(re.search(r'^const branches=(\[.*\])$',source,re.M).group(1))
assert len(cases)==8 and cases['case'].is_unique
assert np.allclose(sorted(cases.target_pre_V5_pu.unique()),[.9036434059107407,.93,.95,.98],rtol=0,atol=1e-10)
assert set(cases.taup)=={.1,1.2}
assert len(params)==96 and params.tauq.eq(.02).all()
assert np.allclose(params.Pset,params.bus.map({1:.5,2:.4,3:.4,7:.4}),rtol=0,atol=0)
assert np.allclose(params.load_G_pu,2.2682301152333055,rtol=0,atol=1e-12)
assert np.max(abs(cases.pre_V5_pu-cases.target_pre_V5_pu))<1e-6
assert eq.residual.max()<1e-8 and eq.max_real_eigen.max()<0 and eq.unstable_modes.eq(0).all()
assert len(machines)==32
baseline_errors={}
for tau,old in [(.1,'tau010_T015_h150'),(1.2,'smoke_tau120')]:
    row=cases[(cases.B_pu==0)&(cases.taup==tau)].iloc[0]
    expected=np.loadtxt(R/'checkpoints'/old/'post030_initial.csv',delimiter=',')
    actual=np.loadtxt(D/f'{row["case"]}_clear.csv',delimiter=',')
    baseline_errors[f'tau{tau:g}_clearing']=float(max(abs(actual-expected)))
    assert baseline_errors[f'tau{tau:g}_clearing']<2e-7
    if tau==.1:
        expected=np.loadtxt(R/'checkpoints'/old/'post150_end.csv',delimiter=',')
        actual=np.loadtxt(D/f'{row["case"]}_end.csv',delimiter=',')
        baseline_errors['tau0.1_end150']=float(max(abs(actual-expected)))
        assert baseline_errors['tau0.1_end150']<2e-7
    else:
        saved=pd.read_csv(R/'results/conditional_basins/data/trajectories/fault_tau1.20.csv.gz')
        saved=saved.iloc[(saved.t_post_s-150).abs().argmin()]
        actual=pd.read_csv(D/f'{row["case"]}.csv').iloc[-1]
        assert np.allclose([actual.theta_infinite_bus_rad,actual.df_mean_Hz],
                           [saved.theta_relative_rad+meta['origin_rad'],saved.df_mean_Hz],rtol=0,atol=1e-7)

def network(x,post,B):
    def v(b,k):return x[int(states.loc[(states.bus==b)&(states.symbol==k),'index'].iloc[0])-1]
    U={b:1+0j if b==6 else complex(v(b,'busbar₊u_r'),v(b,'busbar₊u_i')) for b in range(1,9)}
    I={b:0j for b in U};loss=0.
    for j,(a,b,y) in enumerate(branches,1):
        z=(meta['k'] if j in [2,3] else 1)/y*np.exp(1j*np.deg2rad(75))
        current=0j if post and j==7 else (U[a]-U[b])/z
        I[a]+=current;I[b]-=current;loss+=abs(current)**2*z.real
    S={b:U[b]*np.conj(I[b]) for b in U}
    assert abs(sum(s.real for s in S.values())-loss)<1e-8
    assert abs(S[4])<1e-8 and abs(S[8])<1e-8
    assert abs(S[5]-complex(-2.2682301152333055,B)*abs(U[5])**2)<1e-8
    for b,ps in zip([1,2,3,7],meta['Pset']):
        ig=complex(v(b,f'inv{b}₊vsrc₊filter₊i_g_r'),v(b,f'inv{b}₊vsrc₊filter₊i_g_i'))
        vc=complex(v(b,f'inv{b}₊vsrc₊filter₊V_C_r'),v(b,f'inv{b}₊vsrc₊filter₊V_C_i'))
        assert abs(S[b]-U[b]*np.conj(ig))<1e-8
        assert abs((vc*np.conj(ig)).real-ps)<1e-8
    return dict(load_P_pu=-S[5].real,cap_Q_pu=S[5].imag,line_loss_P_pu=loss,grid_net_P_pu=S[6].real)

records=[]
for row in cases.itertuples(index=False):
    rec=row._asdict()
    audit=params[params['case']==row.case]
    assert np.allclose(audit.cap_B_pu,row.B_pu,rtol=0,atol=1e-12)
    assert np.allclose(audit.taup,row.taup,rtol=0,atol=0)
    assert len(eq[eq['case']==row.case])==(2 if row.post_root_found else 1)
    for stage in ['pre','post']:
        if stage=='post' and not row.post_root_found:continue
        x=np.loadtxt(D/f'{row.case}_{stage}_root.csv',delimiter=',')
        result=network(x,stage=='post',row.B_pu)
        for key in ['load_P_pu','cap_Q_pu']:
            assert abs(rec[stage+'_'+key]-result[key])<1e-8
        rec.update({stage+'_'+k:v for k,v in result.items()})
    track=pd.read_csv(D/f'{row.case}.csv')
    assert track.t_post_s.iloc[0]==0 and track.t_post_s.iloc[-1]==150
    assert np.isfinite(track.to_numpy()).all()
    assert abs((track.theta_infinite_bus_rad.iloc[-1]-track.theta_infinite_bus_rad.iloc[0])/(2*np.pi)-row.remote_turns)<1e-10
    records.append(rec)
summary=pd.DataFrame(records)
summary.to_csv(T/'cases_with_power_balance.csv',index=False)
for name in ['equilibria.csv','parameters.csv','machines.csv']:
    shutil.copyfile(D/name,T/name)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(10,6.5),constrained_layout=True)
for col,tau in enumerate([.1,1.2]):
    for row in cases[cases.taup==tau].sort_values('B_pu').itertuples(index=False):
        d=pd.read_csv(D/f'{row.case}.csv');f=pd.read_csv(D/f'{row.case}_fault.csv')
        label=f'Pre V5={row.pre_V5_pu:.3f}, B={row.B_pu:.3f}'
        time=np.r_[f.t_absolute_s.to_numpy()-1.15,d.t_post_s.to_numpy()]
        angle=np.r_[f.theta_infinite_bus_rad.to_numpy(),d.theta_infinite_bus_rad.to_numpy()]
        freq=np.r_[f.df_mean_Hz.to_numpy(),d.df_mean_Hz.to_numpy()]
        axes[0,col].plot(time,angle,lw=1.2,label=label)
        axes[1,col].plot(time,freq,lw=1.2)
    axes[0,col].set(title=rf'$\tau_p={tau:g}$ s',ylabel='Mean internal phase / Bus6 [rad]',xlabel='Time after clearing [s]',xlim=(-.15,150))
    axes[1,col].set(ylabel='Mean frequency deviation [Hz]',xlabel='Time after clearing [s]',xlim=(-.15,15))
    axes[0,col].legend(fontsize=7)
    for ax in axes[:,col]:
        ax.axvline(0,color='black',lw=.6,ls=':');ax.grid(alpha=.2)
fig.suptitle('Bus5 constant shunt capacitor; unchanged load conductance and GFM P setpoints')
for ext in ['png','pdf']:
    fig.savefig(F/f'capacitor_responses.{ext}',dpi=180)
shutil.copyfile(F/'capacitor_responses.png',R/'docs/figures/capacitor_responses.png')
report=dict(case_count=len(cases),root_count=len(eq),all_found_roots_locally_stable=True,
            classifications=cases.classification.value_counts().to_dict(),max_root_residual=float(eq.residual.max()),
            max_voltage_target_error=float(max(abs(cases.pre_V5_pu-cases.target_pre_V5_pu))),
            baseline_full_state_errors=baseline_errors)
(T/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(summary[['case','B_pu','taup','pre_V5_pu','post_eq_V5_pu','pre_load_P_pu','pre_cap_Q_pu','post_root_found','classification','remote_turns']].to_string(index=False))
print(json.dumps(report,indent=2))
print('PASS: capacitor sign and V-squared law, unchanged G and Pset, steady power balances, baseline states, and trajectories')
