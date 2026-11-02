"""Draw pre/post equilibrium bus voltages and terminal phases; no simulation needed."""
from pathlib import Path
import ast
import re
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'docs/figures'
OUT=ROOT/'data/equilibrium_power_flow'
FIG.mkdir(parents=True,exist_ok=True)
OUT.mkdir(parents=True,exist_ok=True)
states=pd.read_csv(ROOT/'results/conditional_basins/data/states.csv')
meta=json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text())
source=(ROOT/'provenance/source_snapshot/studies/04_clearing/07_grid_strength/code/run.jl').read_text()
branches=ast.literal_eval(re.search(r'^const branches=(\[.*\])$',source,re.M).group(1))
assert len(branches)==11 and branches[6][:2]==(4,5)

def value(x,bus,symbol):
    return float(x[int(states.loc[(states.bus==bus)&(states.symbol==symbol),'index'].iloc[0])-1])

bus_rows=[];line_rows=[];balances=[]
for stage,file in [('pre','nf02_1_pre0.csv'),('post','nf02_1_post0.csv')]:
    x=np.loadtxt(ROOT/'data/seeds'/file,delimiter=',')
    U={b:1+0j if b==6 else complex(value(x,b,'busbar₊u_r'),value(x,b,'busbar₊u_i')) for b in range(1,9)}
    currents={b:0j for b in U}
    for j,(a,b,magnitude) in enumerate(branches,1):
        # All branches use 75 degrees; gridparams scales only the two grid ties.
        Z=(meta['k'] if j in [2,3] else 1)/magnitude*np.exp(1j*np.deg2rad(75))
        active=not(stage=='post' and (a,b)==(4,5))
        current=(U[a]-U[b])/Z if active else 0j
        currents[a]+=current;currents[b]-=current
        sa=U[a]*np.conj(current);sb=U[b]*np.conj(-current)
        assert abs((sa+sb).real-abs(current)**2*Z.real)<1e-10
        line_rows.append(dict(stage=stage,line=j,from_bus=a,to_bus=b,active=active,
                              R_pu=Z.real,X_pu=Z.imag,current_pu=abs(current),
                              P_from_pu=sa.real,P_to_pu=sb.real,loss_pu=(sa+sb).real))
    injections={b:U[b]*np.conj(currents[b]) for b in U}
    for j,b in enumerate([1,2,3,7]):
        ig=complex(value(x,b,f'inv{b}₊vsrc₊filter₊i_g_r'),value(x,b,f'inv{b}₊vsrc₊filter₊i_g_i'))
        vc=complex(value(x,b,f'inv{b}₊vsrc₊filter₊V_C_r'),value(x,b,f'inv{b}₊vsrc₊filter₊V_C_i'))
        assert abs(injections[b]-U[b]*np.conj(ig))<1e-8
        assert abs((vc*np.conj(ig)).real-meta['Pset'][j])<1e-8
    assert abs(injections[4])<1e-8 and abs(injections[8])<1e-8
    assert abs(injections[5].imag)<1e-8
    for b in U:
        bus_rows.append(dict(stage=stage,bus=b,V_pu=abs(U[b]),
                             terminal_angle_deg=float(np.angle(U[b],deg=True)),
                             P_injection_pu=injections[b].real,Q_injection_pu=injections[b].imag))
    load=-injections[5].real
    total=sum(injections[b].real for b in [1,2,3,7])
    loss=sum(row['loss_pu'] for row in line_rows if row['stage']==stage)
    assert abs(total+injections[6].real-load-loss)<1e-8
    balances.append(dict(stage=stage,GFM_terminal_P_pu=total,slack_net_P_pu=injections[6].real,
                         load_P_pu=load,line_loss_P_pu=loss,load_G_pu=load/abs(U[5])**2))
assert abs(balances[0]['load_G_pu']-balances[1]['load_G_pu'])<1e-8
buses=pd.DataFrame(bus_rows);lines=pd.DataFrame(line_rows);balance=pd.DataFrame(balances)
for name,table in [('buses',buses),('lines',lines),('balance',balance)]:table.to_csv(OUT/f'{name}.csv',index=False)

pos={1:(0,5),4:(3.5,5),6:(7,5),3:(0,1),5:(3.5,1),2:(0,-3),8:(7,-3),7:(3.5,-6)}
roles={1:'GFM',2:'GFM',3:'GFM',7:'GFM',4:'fault bus',5:'load',6:'infinite bus',8:'junction'}
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(12,9))
for ax,stage,title in zip(axes,['pre','post'],['Before fault: pre-fault equilibrium','After clearing: stable post-fault equilibrium']):
    table=buses[buses.stage==stage].set_index('bus')
    branch_table=lines[lines.stage==stage].set_index('line')
    for j,(a,b,magnitude) in enumerate(branches,1):
        pa=np.array(pos[a]);pb=np.array(pos[b]);opened=not branch_table.loc[j,'active']
        main=(a,b) in [(4,6),(5,6),(4,5)]
        ax.plot([pa[0],pb[0]],[pa[1],pb[1]],color='#b2182b' if opened else '#3a6c91' if main else '#a5adb5',
                lw=1.5 if main else 1.1,ls='--' if opened else '-',zorder=0)
        if opened:
            middle=(pa+pb)/2
            ax.plot(middle[0],middle[1],'x',color='#b2182b',ms=11,mew=2,zorder=2)
            ax.text(middle[0]-.3,middle[1],'OPEN',ha='right',va='center',color='#b2182b',fontsize=9,
                    bbox=dict(facecolor='white',edgecolor='none',pad=2),zorder=3)
        elif main:
            pfrom=float(branch_table.loc[j,'P_from_pu']);pto=float(branch_table.loc[j,'P_to_pu'])
            start,end=(pa,pb) if pfrom>0 else (pb,pa)
            sent,received=(pfrom,-pto) if pfrom>0 else (pto,-pfrom)
            arrow=FancyArrowPatch(start+.44*(end-start),start+.56*(end-start),
                                  arrowstyle='-|>',mutation_scale=12,lw=1.7,color='#24597f',zorder=2)
            ax.add_patch(arrow)
            text=f'{sent:.3f} → {received:.3f}'
            if (a,b)==(4,6):at=(5.25,6.12);ha='center'
            elif (a,b)==(4,5):at=(3.14,3.0);ha='right'
            else:at=(5.9,2.75);ha='left'
            ax.text(*at,text,ha=ha,va='center',fontsize=8,color='#174769',
                    bbox=dict(facecolor='white',edgecolor='none',pad=2),zorder=3)
    for b,(x,y) in pos.items():
        face='#eaf2fb' if b in [1,2,3,7] else '#fff0d9' if b==5 else '#e4f3e8' if b==6 else '#f5f5f5'
        ax.add_patch(FancyBboxPatch((x-1.45,y-.88),2.9,1.76,boxstyle='round,pad=0.035',
                                   facecolor=face,edgecolor='#546271',lw=.9,zorder=4))
        row=table.loc[b]
        ax.text(x,y+.55,f'Bus {b}: {roles[b]}',ha='center',va='center',weight='bold',fontsize=8.5,zorder=5)
        ax.text(x,y+.17,f'|U| = {row.V_pu:.3f} pu',ha='center',va='center',fontsize=8.5,zorder=5)
        ax.text(x,y-.18,f'∠U = {row.terminal_angle_deg:+.2f}°',ha='center',va='center',fontsize=8.5,zorder=5)
        ax.text(x,y-.55,f'P injection = {row.P_injection_pu:+.3f} pu',ha='center',va='center',fontsize=8,zorder=5)
    total=balance[balance.stage==stage].iloc[0]
    ax.text(3.5,-7.65,f'GFM terminals: {total.GFM_terminal_P_pu:.3f} pu    Grid net: +{total.slack_net_P_pu:.3f} pu\n'
            f'Load: {total.load_P_pu:.3f} pu    Line losses: {total.line_loss_P_pu:.3f} pu',
            ha='center',va='center',fontsize=9,linespacing=1.6,
            bbox=dict(facecolor='#f6f8fa',edgecolor='#d3d8df',boxstyle='round,pad=.6'))
    ax.set_title(title,fontsize=11,pad=20)
    ax.set(xlim=(-1.8,9.4),ylim=(-8.6,6.8),aspect='equal')
    ax.axis('off')
fig.suptitle('Bus-voltage magnitudes and terminal-voltage phases (Bus6 = 1 pu ∠ 0°)',fontsize=13,y=.985)
fig.text(.5,.043,'Blue arrows: selected active-power paths; labels are sending-end P → received P [pu].',ha='center',fontsize=9)
fig.text(.5,.022,'Bus P: net terminal injection (+ generation, − consumption). Post-fault panel is an equilibrium, not the clearing instant.',ha='center',fontsize=8.5)
fig.tight_layout(rect=(0,.06,1,.955),w_pad=2)
for ext in ['png','pdf','svg']:
    temporary=FIG/f'.equilibrium_network.{ext}.tmp'
    fig.savefig(temporary,format=ext,dpi=200,bbox_inches='tight')
    temporary.replace(FIG/f'equilibrium_network.{ext}')
fig.canvas.draw()
for ax,stage in zip(axes,['pre','post']):
    # Full-width single-panel images keep the bus labels legible in the README.
    bbox=ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.025,1.025)
    temporary=FIG/f'.equilibrium_network_{stage}.png.tmp'
    fig.savefig(temporary,format='png',dpi=200,bbox_inches=bbox)
    temporary.replace(FIG/f'equilibrium_network_{stage}.png')
plt.close(fig)
print(balance.to_string(index=False))
print('PASS: converter dispatch, nodal balances, constant load conductance, and total P balance')
print('Created:',FIG/'equilibrium_network.png')
