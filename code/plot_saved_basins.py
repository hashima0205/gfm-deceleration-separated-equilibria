"""Plot saved basins and optionally additional event trajectories; performs no integration."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[1]
BUNDLE=ROOT/'results/conditional_basins'
D=BUNDLE/'data'
COLORS={'R':'#9ecae1','N':'#fdae6b','P':'#bc80bd','U':'#d9d9d9'}
FAULT_COLOR='#762a83'
STYLES={.1:':',1.16:'--',1.18:'-.',1.2:'-'}
TAU_COLORS={.1:'#1f77b4',1.16:'#ff7f0e',1.18:'#2ca02c',1.2:'#d62728'}

def outcome(label):
    if label.startswith('recovery'):return 'R'
    return {'negative_running':'N','positive_running':'P'}.get(label,'U')

def save(fig,directory,name):
    for ext in ['pdf','png']:
        destination=directory/f'{name}.{ext}'
        temporary=directory/f'.{name}.{ext}.tmp'
        fig.savefig(temporary,format=ext,dpi=220,bbox_inches='tight')
        temporary.replace(destination)
    plt.close(fig)

def basin_figures(directory,meta,meshes,trajectories,faults=None):
    def panel(ax,tau,xlim,ylim,ylabel=False):
        samples,xy,labels,faces=meshes[tau]
        same=[];facecolors=[];mixed=[]
        for face in faces:
            values={labels[int(i)] for i in face}
            if len(values)==1:
                same.append(xy[face]);facecolors.append(COLORS[next(iter(values))])
            else:mixed.append(xy[face])
        if same:ax.add_collection(PolyCollection(same,facecolors=facecolors,edgecolors='none',alpha=.65,zorder=0))
        if mixed:ax.add_collection(PolyCollection(mixed,facecolors='#eeeeee',edgecolors='#bbbbbb',linewidths=.2,hatch='////',alpha=.55,zorder=1))
        binary=np.array([0. if c=='R' else 1. if c=='N' else np.nan for c in labels])
        mask=np.any(~np.isfinite(binary[faces]),axis=1)
        if np.any(binary==0) and np.any(binary==1) and not mask.all():
            triangulation=mtri.Triangulation(xy[:,0],xy[:,1],faces)
            triangulation.set_mask(mask)
            ax.tricontour(triangulation,np.nan_to_num(binary),levels=[.5],colors='#666666',linewidths=.7,linestyles='--',zorder=2)
        ax.scatter(xy[:,0],xy[:,1],s=4,c=[COLORS[c] for c in labels],edgecolors='#555555',linewidths=.2,zorder=3)
        if not len(samples):
            ax.text(.03,.97,'Basin not sampled',transform=ax.transAxes,va='top',fontsize=7,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
        elif 'wall_s' in samples.columns:
            ax.text(.03,.97,f'{len(samples)} points; 30/60/150 s',transform=ax.transAxes,va='top',fontsize=7,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
        d=trajectories[tau];x=d.theta_relative_rad.to_numpy();y=d.df_mean_Hz.to_numpy()
        ax.plot(x,y,color='#111111',lw=1,zorder=4)
        if faults is not None:
            fault=faults[tau]
            fx=fault.theta_relative_rad.to_numpy();fy=fault.df_mean_Hz.to_numpy()
            ax.plot(fx,fy,color=FAULT_COLOR,lw=1.5,zorder=6)
            ax.plot(fx[0],fy[0],'^',mfc='white',mec=FAULT_COLOR,ms=5,zorder=8)
            for fraction in [.35,.75]:
                a=int(fraction*(len(fx)-1));b=min(a+20,len(fx)-1)
                ax.annotate('',xy=(fx[b],fy[b]),xytext=(fx[a],fy[a]),
                            arrowprops=dict(arrowstyle='-|>',lw=.8,color=FAULT_COLOR,mutation_scale=8),zorder=7)
        for t in [1.,4.,10.,20.,45.,70.]:
            a=np.searchsorted(d.t_post_s,t);b=np.searchsorted(d.t_post_s,t+.6)
            if b>=len(x):continue
            if xlim[0]<x[a]<xlim[1] and ylim[0]<y[a]<ylim[1] and np.hypot((x[b]-x[a])/(xlim[1]-xlim[0]),(y[b]-y[a])/(ylim[1]-ylim[0]))>.008:
                ax.annotate('',xy=(x[b],y[b]),xytext=(x[a],y[a]),arrowprops=dict(arrowstyle='-|>',lw=.7,color='black',mutation_scale=8),zorder=5)
        ax.plot(x[0],y[0],'o',mfc='white',mec='black',ms=5,zorder=7)
        ax.plot(0,0,'*',mfc='white',mec='black',ms=9,zorder=8)
        ax.plot(meta['saddle_theta_relative_rad'],0,'D',mfc='none',mec='black',ms=5,zorder=8)
        ax.set_xlim(*xlim);ax.set_ylim(*ylim);ax.grid(alpha=.15)
        ax.set_title(r'$\tau_p$ = '+f'{tau:g} s')
        ax.set_xlabel(r'Remote mean internal angle $\theta$ (rad)')
        if ylabel:ax.set_ylabel(r'Remote mean frequency deviation $\Delta f$ (Hz)')
    legend=[Patch(facecolor=COLORS['R'],alpha=.65,label='Synchronization recovery'),
            Patch(facecolor=COLORS['N'],alpha=.65,label='Negative rotation'),
            Patch(facecolor='#eeeeee',hatch='////',edgecolor='#bbbbbb',label='Mixed sampled outcomes'),
            Line2D([],[],color='#666666',ls='--',label='Interpolated transition (guide)'),
            Line2D([],[],color='black',label='Post-fault trajectory'),
            Line2D([],[],marker='*',ls='',mfc='white',mec='black',ms=9,label='Stable equilibrium'),
            Line2D([],[],marker='D',ls='',mfc='none',mec='black',label='Saddle projection'),
            Line2D([],[],marker='o',ls='',mfc='white',mec='black',label='Clearing projection')]
    if faults is not None:
        legend += [Line2D([],[],color=FAULT_COLOR,lw=1.5,label='Fault-on trajectory'),
                   Line2D([],[],marker='^',ls='',mfc='white',mec=FAULT_COLOR,label='Fault application')]
    if any('U' in mesh[2] for mesh in meshes.values()):
        legend.append(Patch(facecolor=COLORS['U'],label='Unresolved / numerical failure'))
    n=len(meta['taus']);width=max(7.1,2.5*n) if n!=3 else 7.1
    limits=meta.get('plot_limits',[[-.4,1.6],[-.15,.05]])
    fig,axes=plt.subplots(1,n,figsize=(width,3.6),sharex=True,sharey=True,squeeze=False)
    for i,tau in enumerate(meta['taus']):panel(axes[0,i],tau,limits[0],limits[1],i==0)
    fig.legend(handles=legend,loc='lower center',ncol=3,fontsize=7,bbox_to_anchor=(.5,-.01))
    fig.suptitle('Conditional numerical resynchronization regions',fontsize=10)
    fig.tight_layout(rect=(0,.29 if faults is not None else .22,1,.96));save(fig,directory,'conditional_basins_comparison')
    for tau in meta['taus']:
        fig,ax=plt.subplots(figsize=(3.46,3.7));panel(ax,tau,limits[0],limits[1],True)
        fig.tight_layout();save(fig,directory,f'conditional_basin_tau{tau:.2f}')
    fig,axes=plt.subplots(2,n,figsize=(width,5.4),sharex='row',sharey='row',squeeze=False)
    clearing_x=(1.25,1.5);clearing_y=(-.025,-.01)
    if n!=3:
        starts=[d.iloc[0] for d in trajectories.values()]
        clearing_x=(min(1.25,min(d.theta_relative_rad for d in starts)-.04),
                    max(1.5,max(d.theta_relative_rad for d in starts)+.04))
        clearing_y=(min(-.025,min(d.df_mean_Hz for d in starts)-.005),
                    max(-.01,max(d.df_mean_Hz for d in starts)+.005))
    for i,tau in enumerate(meta['taus']):
        panel(axes[0,i],tau,(-.2,.15),(-.006,.006),i==0)
        panel(axes[1,i],tau,clearing_x,clearing_y,i==0)
    fig.suptitle('Local views: equilibrium neighborhood (top), clearing neighborhood (bottom)',fontsize=9)
    fig.tight_layout(rect=(0,0,1,.96));save(fig,directory,'conditional_basins_zoom')

def fault_figure(directory,faults):
    """Resolve the short fault interval in phase coordinates and time."""
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.3),constrained_layout=True)
    for tau in sorted(faults):
        style=STYLES.get(tau,':')
        d=faults[tau];label=rf'$\tau_p={tau:g}$ s'
        line=axes[0].plot(d.theta_relative_rad,d.df_mean_Hz,style,lw=1.3,label=label,
                          color=TAU_COLORS.get(tau))[0]
        color=line.get_color()
        axes[0].plot(d.theta_relative_rad.iloc[0],d.df_mean_Hz.iloc[0],'^',color=color,mfc='white',mec=color,ms=6)
        axes[0].plot(d.theta_relative_rad.iloc[-1],d.df_mean_Hz.iloc[-1],'o',color=color,mfc='white',mec=color,ms=5)
        axes[1].plot(d.t_post_s,d.df_mean_Hz,style,color=color,lw=1.3)
    axes[0].set_xlabel(r'Angle relative to post-fault equilibrium $\theta$ (rad)')
    axes[0].set_ylabel(r'Remote mean $\Delta f$ (Hz)')
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel('Time relative to clearing (s)')
    axes[1].set_ylabel(r'Remote mean $\Delta f$ (Hz)')
    axes[1].set_xticks([-.15,-.10,-.05,0.])
    axes[1].axvline(0,color='black',ls=':',lw=.7)
    for ax in axes:ax.grid(alpha=.2)
    axes[0].set_title('Fault-on phase trajectories')
    axes[1].set_title('Fault application at -0.15 s; clearing at 0 s')
    fig.suptitle('Fault interval: triangles = application, circles = clearing',fontsize=9)
    save(fig,directory,'fault_on_comparison')

def hidden_figure(directory):
    pairs=pd.read_csv(D/'hidden_state_comparison.csv')
    pairs=pairs[~pairs.classification_matches]
    fig,axes=plt.subplots(1,len(pairs),figsize=(7.1,2.8),squeeze=False)
    for ax,r in zip(axes[0],pairs.itertuples()):
        for case,label,style,color in [(r.equilibrium_case,'Equilibrium hidden states','-','#c55b16'),(r.clearing_case,'Clearing hidden states','--','#2166ac')]:
            d=pd.read_csv(D/f'hidden_trajectories/{case}.csv.gz')
            d=d[d.t_post_s<=100]
            ax.plot(d.t_post_s,d.theta_relative_rad/(2*np.pi),style,color=color,lw=1,label=label)
        ax.set_title(r'$\tau_p$ = '+f'{r.taup:g} s');ax.set_xlabel('Time from initialization (s)');ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Remote internal angle relative\nto equilibrium (turns)')
    axes[0,-1].legend(fontsize=6)
    fig.suptitle('Identical initial angle and frequency; different hidden states',fontsize=9)
    fig.tight_layout();save(fig,directory,'hidden_state_sensitivity')

def time_figures(directory,trajectories):
    fig,axes=plt.subplots(2,1,figsize=(5.4,4.4),sharex=True,constrained_layout=True)
    for tau in sorted(trajectories):
        style=STYLES.get(tau,':')
        d=trajectories[tau];t=d.t_post_s
        line=axes[0].plot(t,(d.theta_relative_rad-d.theta_relative_rad.iloc[0])/(2*np.pi),style,
                          color=TAU_COLORS.get(tau),label=rf'$\tau_p={tau:g}$ s')[0]
        axes[1].plot(t,d.df_mean_Hz,style,color=line.get_color())
        f,a=plt.subplots(2,2,figsize=(7.1,4.7),sharex=True,constrained_layout=True)
        for bus,ls in zip([1,2,3,7],['-','--','-.',':']):
            angle=d[f'b{bus}_delta_rad']
            line=a[0,0].plot(t,(angle-angle.iloc[0])/(2*np.pi),ls,label=f'Bus {bus}')[0]
            color=line.get_color()
            for ax,column in [(a[0,1],f'b{bus}_df_Hz'),(a[1,0],f'b{bus}_V_pu'),(a[1,1],f'b{bus}_Pcap_pu')]:ax.plot(t,d[column],ls,color=color)
            a[1,1].axhline(.5 if bus==1 else .4,color=color,ls=':',lw=.6)
        for ax,label in zip(a.flat,['Internal angle change [turns]','Internal frequency deviation [Hz]','Terminal voltage magnitude [pu]','Capacitor-side active power [pu]']):ax.set_ylabel(label);ax.grid(alpha=.25)
        for ax in a[1]:ax.set_xlabel('Time after clearing [s]')
        a[0,0].legend(ncol=2,fontsize=8);f.suptitle(rf'$\tau_p={tau:g}$ s, $T=0.15$ s')
        save(f,directory,f'each_gfm_tau{tau:.2f}')
    axes[0].set_ylabel('Remote mean internal\nangle change [turns]');axes[0].legend(fontsize=8)
    axes[1].set_ylabel('Remote mean internal\nfrequency deviation [Hz]');axes[1].set_xlabel('Time after clearing [s]')
    for ax in axes:ax.grid(alpha=.25)
    save(fig,directory,'saved_time_responses')

def add_case(tag,meta,meshes,trajectories,faults):
    """Add a freshly calculated trajectory without inventing a basin map."""
    import observations12 as observations
    from classify import measure
    row=pd.read_csv(ROOT/'data/cases.csv').drop_duplicates('case',keep='last').set_index('case').loc[tag]
    if row.retcode!='Success':raise ValueError(f'Incomplete case: {tag}')
    tau=float(row.taup)
    if tau in trajectories:raise ValueError(f'A trajectory for tau={tau} already exists')
    if row.tauq!=meta['tauq'] or row['T']!=meta['T'] or row.load_scale!=meta['load_scale']:
        raise ValueError('Additional case differs in tauq, fault duration, or load')
    if [row.P1,row.P2,row.P3,row.P7]!=meta['Pset']:raise ValueError('Different power references')
    audit=pd.read_csv(ROOT/'data/parameter_audit.csv')
    audit=audit[audit['case']==tag].drop_duplicates(['stage','bus'],keep='last')
    assert len(audit)==12 and set(audit.bus)=={1,2,3,7}
    assert (audit.taup==tau).all() and (audit.tauq==meta['tauq']).all()
    assert (audit.VC_F==1).all() and (audit.CC1_F==1).all()
    assert np.allclose(audit.Pset,audit.bus.map({1:row.P1,2:row.P2,3:row.P3,7:row.P7}),rtol=0,atol=0)
    roots=pd.read_csv(ROOT/'data/equilibria.csv')
    roots=roots[(roots['case']==tag)&(roots.stage=='post')].drop_duplicates('branch',keep='last')
    assert (roots.residual<1e-8).all() and (roots.max_real_eigen<0).any() and (roots.unstable_modes==1).any()
    d=observations.load(tag)
    d['theta_relative_rad']=d.theta_mean_rad-meta['origin_rad']
    post=d[d.stage.str.startswith('post')].drop_duplicates('t_post_s',keep='last').copy()
    fault=d[d.stage=='fault'].copy()
    assert np.isclose(post.t_post_s.iloc[-1],row.horizon)
    assert np.isclose(fault.t_post_s.iloc[0],-row['T']) and np.isclose(fault.t_post_s.iloc[-1],0.)
    assert np.allclose(fault[['theta_relative_rad','df_mean_Hz']].iloc[-1],
                       post[['theta_relative_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-8)
    stable=np.loadtxt(ROOT/'data/seeds'/f'{tag}_stable.csv',delimiter=',')
    saddle=np.loadtxt(ROOT/'data/seeds'/f'{tag}_saddle.csv',delimiter=',')
    assert np.allclose(stable,np.loadtxt(D/'references/section_stable.csv',delimiter=','),rtol=0,atol=1e-7)
    metrics=measure(tag,d,stable,saddle)
    print(f'Additional case: {tag}; {metrics["classification"]}; '
          f'{metrics["remote_internal_turns"]:.9g} post-fault turns; {row.horizon:g} s')
    trajectories[tau]=post
    faults[tau]=fault
    meshes[tau]=(pd.DataFrame(),np.empty((0,2)),[],np.empty((0,3),dtype=int))
    meta['taus']=sorted(trajectories)
    # Show the new physical orbit even if it leaves the existing map domain.
    coords=pd.concat([post,fault])
    meta['plot_limits']=[
        [min(-.4,float(coords.theta_relative_rad.min())-.05),max(1.6,float(coords.theta_relative_rad.max())+.05)],
        [min(-.15,float(coords.df_mean_Hz.min())-.01),max(.05,float(coords.df_mean_Hz.max())+.01)]]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'figures/saved_results/conditional_basins')
    parser.add_argument('--fault-trajectories',type=Path,default=ROOT/'data/fault_trajectories.csv.gz',
                        help='Optional verified fault-on export from export_fault_trajectories.py')
    parser.add_argument('--additional-case',help='Completed checkpoint case to compare; no basin map is inferred')
    parser.add_argument('--display-horizon',type=float,help='Common post-fault display horizon in seconds')
    parser.add_argument('--basin-samples',type=Path,help='Additional basin samples; paired plot_triangles_tauXX.csv required')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    meta=json.loads((D/'definition.json').read_text());meshes={};trajectories={}
    for tau in meta['taus']:
        samples=pd.read_csv(D/f'plot_samples_tau{tau:.2f}.csv')
        xy=samples[['theta_relative_rad','df_Hz']].to_numpy()
        faces=pd.read_csv(D/f'plot_triangles_tau{tau:.2f}.csv').to_numpy(dtype=int)
        meshes[tau]=(samples,xy,[outcome(c) for c in samples.classification],faces)
        trajectories[tau]=pd.read_csv(D/f'trajectories/fault_tau{tau:.2f}.csv.gz')
    faults=None
    if args.fault_trajectories.exists():
        fault_data=pd.read_csv(args.fault_trajectories)
        faults={}
        for tau in meta['taus']:
            d=fault_data[np.isclose(fault_data.taup,tau)].sort_values('t_post_s')
            if d.empty or not np.isclose(d.t_post_s.iloc[0],-meta['T']) or not np.isclose(d.t_post_s.iloc[-1],0):
                raise ValueError(f'Missing or incomplete fault interval for tau={tau}')
            if not np.isfinite(d[['theta_relative_rad','df_mean_Hz']].to_numpy()).all():
                raise ValueError(f'Nonfinite fault coordinates for tau={tau}')
            if not np.allclose(d[['theta_relative_rad','df_mean_Hz']].iloc[-1],
                               trajectories[tau][['theta_relative_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-7):
                raise ValueError(f'Fault and saved post-fault trajectories do not meet for tau={tau}')
            faults[tau]=d
    if args.additional_case:
        if faults is None:raise ValueError('Additional comparison requires the existing fault-on export')
        add_case(args.additional_case,meta,meshes,trajectories,faults)
    if args.basin_samples:
        samples=pd.read_csv(args.basin_samples)
        taus=samples.taup.unique()
        if len(taus)!=1 or float(taus[0]) not in trajectories:
            raise ValueError('Additional basin samples require one matching trajectory')
        tau=float(taus[0])
        xy=samples[['theta_relative_rad','df_Hz']].to_numpy()
        if not np.isfinite(xy).all():raise ValueError('Nonfinite basin coordinates')
        faces=pd.read_csv(args.basin_samples.with_name(f'plot_triangles_tau{tau:.2f}.csv')).to_numpy(dtype=int)
        if faces.min()<0 or faces.max()>=len(samples):raise ValueError('Invalid triangle indices')
        meshes[tau]=(samples,xy,[outcome(c) for c in samples.classification],faces)
    if args.display_horizon is not None:
        if args.display_horizon<=0:raise ValueError('Display horizon must be positive')
        trajectories={tau:d[d.t_post_s<=args.display_horizon].copy() for tau,d in trajectories.items()}
    if faults is not None: fault_figure(args.output,faults)
    basin_figures(args.output,meta,meshes,trajectories,faults)
    hidden_figure(args.output);time_figures(args.output,trajectories)
    print('Replotted saved data; no numerical integration:',args.output)
    print('Fault-on overlay:',args.fault_trajectories if faults is not None else 'not available (post-fault only)')

if __name__=='__main__':main()
