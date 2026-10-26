"""Plot actual saved samples; internal and terminal phases are kept distinct."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import observations12 as o

P=Path(__file__).resolve().parents[1]
F=P/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.linewidth':.7,'lines.linewidth':1.0,'pdf.fonttype':42})
cases=pd.read_csv(P/'data/cases.csv').drop_duplicates('case',keep='last')
fig,axes=plt.subplots(2,1,figsize=(5.4,4.4),sharex=True,constrained_layout=True)
count=0
for r in cases[cases.retcode=='Success'].itertuples(index=False):
    d=o.load(r.case)
    p=d[d.stage.str.startswith('post')].drop_duplicates('t_post_s',keep='last')
    t=p.t_post_s
    label=f'tau_p={r.taup:g} s ({r.horizon:g} s)'
    if not r.case.startswith('smoke_') or len(cases)==1:
        style={1.16:'--',1.18:'-.',1.2:'-'}.get(float(r.taup),':')
        line=axes[0].plot(t,(p.theta_mean_rad-p.theta_mean_rad.iloc[0])/(2*np.pi),style,label=label)[0]
        axes[1].plot(t,p.df_mean_Hz,style,color=line.get_color(),label=label)
    fig2,a=plt.subplots(2,2,figsize=(7.1,4.7),sharex=True,constrained_layout=True)
    for bus,style in zip([1,2,3,7],['-','--','-.',':']):
        angle=p[f'b{bus}_delta_rad']
        line=a[0,0].plot(t,(angle-angle.iloc[0])/(2*np.pi),style,label=f'Bus {bus}')[0]
        color=line.get_color()
        for ax,column in [(a[0,1],f'b{bus}_df_Hz'),(a[1,0],f'b{bus}_V_pu'),(a[1,1],f'b{bus}_Pcap_pu')]:
            ax.plot(t,p[column],style,color=color)
        a[1,1].axhline(getattr(r,f'P{bus}'),color=color,ls=':',lw=.6)
    for ax,label in zip(a.flat,['Internal angle change [turns]','Internal frequency deviation [Hz]','Terminal voltage magnitude [pu]','Capacitor-side active power [pu]']):
        ax.set_ylabel(label);ax.grid(alpha=.25)
    for ax in a[1]:ax.set_xlabel('Time after clearing [s]')
    a[0,0].legend(ncol=2,fontsize=8)
    fig2.suptitle(r.case)
    for suffix in ['pdf','png']:fig2.savefig(F/f'{r.case}.{suffix}',dpi=180)
    plt.close(fig2);count+=1
axes[0].set_ylabel('Remote mean internal\nangle change [turns]')
axes[1].set_ylabel('Remote mean internal\nfrequency deviation [Hz]')
axes[1].set_xlabel('Time after clearing [s]')
for ax in axes:ax.grid(alpha=.25)
if count:axes[0].legend(fontsize=8)
for suffix in ['pdf','png']:fig.savefig(F/f'tau_comparison.{suffix}',dpi=180)
plt.close(fig)
print(f'Created comparison and {count} individual-machine figure(s) in {F}')
