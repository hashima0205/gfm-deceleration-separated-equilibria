"""Replot immutable saved basin data and trajectories; performs no integration."""
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

def outcome(label):
    if label.startswith('recovery'):return 'R'
    return {'negative_running':'N','positive_running':'P'}.get(label,'U')

def save(fig,directory,name):
    for ext in ['pdf','png']:fig.savefig(directory/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

def basin_figures(directory,meta,meshes,trajectories):
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
        d=trajectories[tau];x=d.theta_relative_rad.to_numpy();y=d.df_mean_Hz.to_numpy()
        ax.plot(x,y,color='#111111',lw=1,zorder=4)
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
            Line2D([],[],color='black',label='Physical fault trajectory'),
            Line2D([],[],marker='*',ls='',mfc='white',mec='black',ms=9,label='Stable equilibrium'),
            Line2D([],[],marker='D',ls='',mfc='none',mec='black',label='Saddle projection'),
            Line2D([],[],marker='o',ls='',mfc='white',mec='black',label='Clearing projection')]
    fig,axes=plt.subplots(1,3,figsize=(7.1,3.6),sharex=True,sharey=True)
    for i,tau in enumerate(meta['taus']):panel(axes[i],tau,(-.4,1.6),(-.15,.05),i==0)
    fig.legend(handles=legend,loc='lower center',ncol=3,fontsize=7,bbox_to_anchor=(.5,-.01))
    fig.suptitle('Conditional numerical resynchronization regions',fontsize=10)
    fig.tight_layout(rect=(0,.22,1,.96));save(fig,directory,'conditional_basins_comparison')
    for tau in meta['taus']:
        fig,ax=plt.subplots(figsize=(3.46,3.7));panel(ax,tau,(-.4,1.6),(-.15,.05),True)
        fig.tight_layout();save(fig,directory,f'conditional_basin_tau{tau:.2f}')
    fig,axes=plt.subplots(2,3,figsize=(7.1,5.4),sharex='row',sharey='row')
    for i,tau in enumerate(meta['taus']):
        panel(axes[0,i],tau,(-.2,.15),(-.006,.006),i==0)
        panel(axes[1,i],tau,(1.25,1.5),(-.025,-.01),i==0)
    fig.suptitle('Local views: equilibrium neighborhood (top), clearing neighborhood (bottom)',fontsize=9)
    fig.tight_layout(rect=(0,0,1,.96));save(fig,directory,'conditional_basins_zoom')

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
    for tau,style in zip(sorted(trajectories),['--','-.','-']):
        d=trajectories[tau];t=d.t_post_s
        line=axes[0].plot(t,(d.theta_relative_rad-d.theta_relative_rad.iloc[0])/(2*np.pi),style,label=rf'$\tau_p={tau:g}$ s')[0]
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

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'figures/saved_results/conditional_basins')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    meta=json.loads((D/'definition.json').read_text());meshes={};trajectories={}
    for tau in meta['taus']:
        samples=pd.read_csv(D/f'plot_samples_tau{tau:.2f}.csv')
        xy=samples[['theta_relative_rad','df_Hz']].to_numpy()
        faces=pd.read_csv(D/f'plot_triangles_tau{tau:.2f}.csv').to_numpy(dtype=int)
        meshes[tau]=(samples,xy,[outcome(c) for c in samples.classification],faces)
        trajectories[tau]=pd.read_csv(D/f'trajectories/fault_tau{tau:.2f}.csv.gz')
    basin_figures(args.output,meta,meshes,trajectories)
    hidden_figure(args.output);time_figures(args.output,trajectories)
    print('Replotted saved data; no numerical integration:',args.output)

if __name__=='__main__':main()
