"""Audit strengthened-grid equilibria/events/chord basins, keep figures in case 02."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from search_capacitor_equilibria import system,vc_from_full,packed,TARGET_B

R=Path(__file__).resolve().parents[1];D=R/'basin_recalculation/capacitor_k_strength';F=R/'figures/02_voltage_adjusted';F.mkdir(parents=True,exist_ok=True)
a=argparse.ArgumentParser();a.add_argument('--events-only',action='store_true');args=a.parse_args()
states=pd.read_csv(R/'results/conditional_basins/data/states.csv');dyn=states.loc[states.mass!=0,'index'].to_numpy(int)-1
ix=lambda b,s:int(states.loc[(states.bus==b)&(states.symbol==f'inv{b}₊droop₊{s}'),'index'].iloc[0])-1
delta=np.array([ix(b,'δ') for b in [1,2,3,7]]);pf=np.array([ix(b,'Pfilt') for b in [1,2,3,7]])
kp=np.full(3,7.5398223686155035)
events=pd.read_csv(D/'events.csv');roots=pd.read_csv(D/'equilibria.csv')
assert len(events)==6 and len(roots)==18
assert roots.residual.max()<1e-8
assert roots.loc[roots.root!='saddle','max_real_eigen'].lt(0).all()
assert roots.loc[roots.root=='saddle','unstable_modes'].ge(1).all()
checks=[]
for row in events.itertuples(index=False):
    assert row.B_pu==TARGET_B
    suffix=f'{int(round(row.taup*100)):03d}'
    for stage in ['pre','fault','post']:
        p=np.loadtxt(D/f'{row.case}_{stage}_parameters.csv',delimiter=',')
        old=np.loadtxt(R/f'basin_recalculation/capacitor_probe_predictor/v095_tau{suffix}_{stage}_parameters.csv',delimiter=',')
        changed=np.flatnonzero(abs(p-old)>1e-14)
        assert len(changed)==4 and np.allclose(p[changed]/old[changed],row.k/1.0299,rtol=0,atol=1e-12)
    for name,post in [('pre',False),('stable',True),('saddle',True)]:
        x=np.loadtxt(D/f'k{int(round(100*row.k)):03d}_{name}.csv',delimiter=',')
        fun,_=system(TARGET_B,post=post,k=row.k);er=float(max(abs(fun(packed(vc_from_full(x))))))
        assert er<1e-8
        checks.append(dict(case=row.case,root=name,reduced_residual=er))
    tr=pd.read_csv(D/f'{row.case}.csv');fault=pd.read_csv(D/f'{row.case}_fault.csv')
    assert tr.t_post_s.iloc[0]==0 and tr.t_post_s.iloc[-1]==150
    assert np.isfinite(tr.to_numpy()).all()
    assert np.allclose(fault[['theta_infinite_bus_rad','df_mean_Hz']].iloc[-1],tr[['theta_infinite_bus_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-9)
    assert abs((tr.theta_infinite_bus_rad.iloc[-1]-tr.theta_infinite_bus_rad.iloc[0])/(2*np.pi)-row.remote_turns)<1e-10
pd.DataFrame(checks).to_csv(D/'independent_root_validation.csv',index=False)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(F/f'{name}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(10,6.5),constrained_layout=True)
for col,tau in enumerate([.1,1.2]):
    for row in events[np.isclose(events.taup,tau)].itertuples(index=False):
        d=pd.read_csv(D/f'{row.case}.csv');f=pd.read_csv(D/f'{row.case}_fault.csv')
        tt=np.r_[f.t_absolute_s.to_numpy()-1.15,d.t_post_s.to_numpy()]
        axes[0,col].plot(tt,np.r_[f.theta_infinite_bus_rad,d.theta_infinite_bus_rad],label=f'k={row.k:g}')
        axes[1,col].plot(tt,np.r_[f.df_mean_Hz,d.df_mean_Hz],label=f'k={row.k:g}')
    axes[0,col].set(title=rf'$\tau_p={tau:g}$ s',ylabel='Mean internal phase / Bus6 [rad]',xlabel='Time after clearing [s]',xlim=(-.15,150))
    axes[1,col].set(ylabel='Mean frequency deviation [Hz]',xlabel='Time after clearing [s]',xlim=(-.15,150))
    for ax in axes[:,col]:ax.legend(fontsize=8);ax.grid(alpha=.2);ax.axvline(0,color='black',ls=':',lw=.5)
fig.suptitle('Fixed capacitor B=0.55618 pu and load conductance; stronger external grid')
save(fig,'04_stronger_grid_trajectories')
events.to_csv(F/'stronger_grid_events.csv',index=False);roots.to_csv(F/'stronger_grid_equilibria.csv',index=False)
report=dict(event_count=6,verified_full_equilibria=18,maximum_root_residual=float(roots.residual.max()),
            independent_steady_residual=float(max(x['reduced_residual'] for x in checks)),
            changed_parameters='Only four R/X entries of lines 4-6 and 5-6; ratios k/1.0299',basins={})
if not args.events_only:
    xs=np.loadtxt(D/'k080_stable.csv',delimiter=',');xu=np.loadtxt(D/'k080_saddle.csv',delimiter=',')
    phases=lambda x:float(x[delta[1:]].mean())
    gap=phases(xu)-phases(xs)
    first=pd.read_csv(D/'basin_tau010_initial.csv',float_precision='round_trip')
    second=pd.read_csv(D/'basin_tau120_initial.csv',float_precision='round_trip')
    assert len(first)==len(second)==55
    statecols=[f'x{i}' for i in range(1,67)]
    common_error=float(np.max(abs(first[statecols].to_numpy()-second[statecols].to_numpy())))
    assert common_error<1e-12
    report['between_tau_initial_state_difference']=common_error
    section_checks=[]
    for tau in [.1,1.2]:
        suffix=f'{int(round(100*tau)):03d}'
        xc=np.loadtxt(D/f'k080_tau{suffix}_clear.csv',delimiter=',')
        dd=pd.read_csv(D/f'k080_tau{suffix}.csv').iloc[0]
        alpha=(phases(xc)-phases(xs))/gap
        section=xs[dyn]+alpha*(xu[dyn]-xs[dyn])
        for b,pp in enumerate(pf[1:]):section[np.flatnonzero(dyn==pp)[0]]-=dd.df_mean_Hz*2*np.pi/kp[b]
        diff=xc[dyn]-section
        section_checks.append(dict(taup=tau,clearing_alpha=alpha,clearing_df_Hz=dd.df_mean_Hz,
                                   unequal_differential_states=int(np.sum(abs(diff)>1e-8)),
                                   max_raw_mixed_unit_difference=float(max(abs(diff)))))
    pd.DataFrame(section_checks).to_csv(F/'stronger_grid_clearing_section_checks.csv',index=False)
    report['clearing_section_checks']=section_checks

    alpha_grid=np.array([-5.,-3.,-1.,-.5,0.,.5,.9,1.,1.1,1.5,3.]);df_grid=np.array([-.15,-.1,-.05,0.,.05])
    faces=[]
    for j in range(4):
        for i in range(10):
            n=j*11+i;faces.extend([[n,n+1,n+12],[n,n+12,n+11]])
    colors={'same':'#9ecae1','shift':'#fdbb84','negative':'#d73027','positive':'#9970ab','unknown':'#d9d9d9','saddle':'#bdbdbd'}
    for tau in [.1,1.2]:
        suffix=f'{int(round(tau*100)):03d}'
        samples=pd.read_csv(D/f'basin_tau{suffix}.csv')
        initial=pd.read_csv(D/f'basin_tau{suffix}_initial.csv',float_precision='round_trip').set_index('case')
        audit=pd.read_csv(D/f'basin_tau{suffix}_audit.csv');history=pd.read_csv(D/f'basin_tau{suffix}_history.csv')
        assert len(samples)==55 and samples['case'].is_unique
        assert set(samples.alpha)==set(alpha_grid) and set(samples.df_offset_Hz)==set(df_grid)
        assert audit.algebraic_residual.max()<1e-8 and audit.differential_change.eq(0).all()
        assert samples.max_algebraic_residual.max()<1e-7
        categories=[]
        for row in samples.itertuples(index=False):
            x=initial.loc[row.case,[f'x{i}' for i in range(1,67)]].to_numpy(float)
            expect=xs[dyn]+row.alpha*(xu[dyn]-xs[dyn])
            for b,pp in enumerate(pf[1:]):expect[np.flatnonzero(dyn==pp)[0]]-=row.df_offset_Hz*2*np.pi/kp[b]
            assert np.allclose(x[dyn],expect,rtol=0,atol=1e-12)
            assert abs(phases(x)-row.theta_infinite_bus_rad)<1e-10
            assert abs(phases(x)-phases(xs)-row.alpha*gap)<1e-10
            h=history[history['case']==row.case];assert h.horizon_s.tolist()==[30,60,150][:len(h)]
            assert h.classification.iloc[:-1].isin(['unresolved','other']).all()
            if row.alpha in [0.,1.] and row.df_offset_Hz==0:
                assert np.max(abs(x-(xs if row.alpha==0. else xu)))<1e-10
            if row.classification.startswith('recovery'):
                end=np.loadtxt(D/f'{row.case}_end.csv',delimiter=',')
                assert max(abs(np.angle(np.exp(1j*(end[delta]-xs[delta])))))<.002
                cat='same' if all(getattr(row,f'branch_b{b}')==0 for b in [1,2,3,7]) else 'shift'
            else:cat={'negative_running':'negative','positive_running':'positive','saddle_initial':'saddle'}.get(row.classification,'unknown')
            categories.append(cat)
        samples['category']=categories
        tr=pd.read_csv(D/f'k080_tau{suffix}.csv');fault=pd.read_csv(D/f'k080_tau{suffix}_fault.csv')
        xy=samples[['theta_infinite_bus_rad','df_Hz']].to_numpy()
        fig,axes=plt.subplots(1,2,figsize=(11,4.7))
        for ax,local in zip(axes,[False,True]):
            for face in faces:
                cc={categories[i] for i in face};same=len(cc)==1 and 'saddle' not in cc
                ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(cc))] if same else '#eeeeee',edgecolors='none' if same else '#aaaaaa',hatch=None if same else '////',lw=.25,alpha=.7))
            ax.scatter(xy[:,0],xy[:,1],c=[colors[c] for c in categories],s=13,edgecolors='#444444',lw=.3,zorder=3)
            ax.plot(tr.theta_infinite_bus_rad,tr.df_mean_Hz,color='black',lw=1.2,zorder=4)
            ax.plot(fault.theta_infinite_bus_rad,fault.df_mean_Hz,color='#762a83',lw=1.8,zorder=5)
            ax.plot(tr.theta_infinite_bus_rad.iloc[0],tr.df_mean_Hz.iloc[0],'o',mfc='white',mec='black',ms=5,zorder=6)
            ax.plot(phases(xs),0,'*',mfc='white',mec='black',ms=12,zorder=7)
            ax.plot(phases(xu),0,'D',mfc='none',mec='black',ms=6,zorder=7)
            ax.set(xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title='Equilibrium / saddle neighborhood' if local else 'Section and fault trajectory projection')
            ax.set_xlim((phases(xu)-.12,phases(xs)+.12) if local else (xy[:,0].min()-.03,xy[:,0].max()+.03))
            inview=tr[tr.theta_infinite_bus_rad.between(xy[:,0].min()-.03,xy[:,0].max()+.03)]
            low=min(-.16,float(inview.df_mean_Hz.min())-.005)
            ax.set_ylim((-.055,.055) if local else (low,.06));ax.grid(alpha=.2)
        legends=[Patch(facecolor=colors[k],label=lab) for k,lab in [('same','Recovery: same branch'),('shift','Recovery: changed branch'),('negative','Negative running'),('positive','Positive running'),('unknown','Unresolved')] if k in ['same','shift','negative'] or k in categories]
        legends.append(Patch(facecolor='#eeeeee',edgecolor='#aaaaaa',hatch='////',label='Mixed sampled outcomes'))
        legends.extend([Line2D([],[],color='black',label='Post-fault projection'),Line2D([],[],color='#762a83',label='Fault-on projection'),Line2D([],[],marker='*',mfc='white',mec='black',ls='',ms=10,label='Stable root'),Line2D([],[],marker='D',mfc='none',mec='black',ls='',label='Saddle on section')])
        fig.legend(handles=legends,loc='lower center',ncol=4,fontsize=8,bbox_to_anchor=(.5,.025))
        fig.suptitle(rf'Voltage-adjusted case, k=0.80, $\tau_p={tau:g}$ s; both equilibria lie on the chord section')
        fig.text(.5,.005,'55 DAE-consistent initial states; 30/60/150-s classification. Actual trajectories are projections, not constrained to the section.',ha='center',fontsize=7)
        fig.tight_layout(rect=(0,.19,1,.93));save(fig,f'05_basin_k080_tau{suffix}')
        report['basins'][f'tau{tau:g}']=dict(points=55,counts=samples.category.value_counts().to_dict(),horizons=samples.horizon_s.value_counts().sort_index().to_dict())
        samples.to_csv(F/f'basin_k080_tau{suffix}_points.csv',index=False)
(F/'stronger_grid_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(events.to_string(index=False));print(json.dumps(report,indent=2));print('PASS: full roots, independent equations, unchanged capacitor/load/setpoints, event continuity, and chord sections')
