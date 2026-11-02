"""Plot only the taup=0.1 s map, separating terminal phase-branch changes."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from plot_saved_basins import save

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--section',choices=['fixed','chord'],default='fixed')
parser.add_argument('--event-overlay',action='store_true',help='Overlay the taup=0.1 physical event in full and local views')
parser.add_argument('--angle-reference',choices=['infinite-bus','equilibrium'],default='infinite-bus',
                    help='Display phase relative to Bus6 (default), or reproduce the old equilibrium-relative axis')
args=parser.parse_args()
chord=args.section=='chord'
if args.event_overlay and not chord:parser.error('--event-overlay currently requires --section chord')
D=ROOT/'basin_recalculation'/('tau010_chord' if chord else 'tau010_adaptive')
F=ROOT/'figures'/('tau010_chord_basin' if chord else 'tau010_adaptive_basin')
F.mkdir(parents=True,exist_ok=True)
meta=json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text())
states=pd.read_csv(ROOT/'results/conditional_basins/data/states.csv')
stable=np.loadtxt(ROOT/'data/seeds/tau010_T015_h150_stable.csv',delimiter=',')
saddle=np.loadtxt(ROOT/'data/seeds/tau010_T015_h150_saddle.csv',delimiter=',')
indices=[int(states.loc[(states.bus==b)&(states.symbol==f'inv{b}₊droop₊δ'),'index'].iloc[0])-1 for b in [1,2,3,7]]
stable_phase=float(np.mean(stable[np.array(indices[1:])]))
assert abs(stable_phase-meta['origin_rad'])<1e-10
# Bus6 is fixed at 1+0j. Keep unwrapped angles and translate every plotted object equally.
offset=stable_phase if args.angle_reference=='infinite-bus' else 0.
def display_angle(relative):return relative+offset
xlabel=(r'Remote mean internal angle relative to Bus6 $\bar{\delta}-\delta_6$ (rad)'
        if args.angle_reference=='infinite-bus' else r'Remote mean internal angle relative to equilibrium $\theta$ (rad)')
samples=pd.read_csv(D/'plot_samples_tau0.10.csv')
assert len(samples)==35 and np.allclose(samples.taup,.1)
rows=[]
initial=pd.read_csv(D/'initial_states.csv',float_precision='round_trip').set_index('case')
for r in samples.itertuples(index=False):
    end=np.loadtxt(D/f'{r.case}_end.csv',delimiter=',')
    branches=np.rint((end[indices]-stable[indices])/(2*np.pi)).astype(int)
    if chord and np.max(abs(initial.loc[r.case].to_numpy()-saddle))<1e-9:
        label='saddle_initial'
    elif r.classification.startswith('recovery'):
        assert np.max(abs(np.angle(np.exp(1j*(end[indices]-stable[indices])))))<.002
        label='original_branch' if (branches==0).all() else 'shifted_branch'
    else:
        label='other'
    rows.append(dict(case=r.case,recapture_type=label,
                     **{f'bus{b}_terminal_branch':int(k) for b,k in zip([1,2,3,7],branches)}))
branches=pd.DataFrame(rows)
branches.to_csv(D/'phase_branches.csv',index=False)
samples=samples.merge(branches,on='case',validate='one_to_one',sort=False)
assert list(samples['case'])==list(branches['case'])
xy=samples[['theta_relative_rad','df_Hz']].to_numpy(copy=True)
xy[:,0]=display_angle(xy[:,0])
samples['theta_infinite_bus_rad']=samples.theta_relative_rad+stable_phase
samples.to_csv(D/'plot_samples_with_angle_references.csv',index=False)
faces=pd.read_csv(D/'plot_triangles_tau0.10.csv').to_numpy(dtype=int)
colors={'original_branch':'#9ecae1','shifted_branch':'#fdbb84','other':'#d9d9d9','saddle_initial':'#bdbdbd'}
event=None
if args.event_overlay:
    import observations12 as observations
    case='tau010_T015_h150'
    checkpoint=ROOT/'checkpoints'/case
    assert np.array_equal(np.loadtxt(checkpoint/'post150_parameters.csv',delimiter=','),
                          np.loadtxt(D/'post_parameters.csv',delimiter=','))
    event=observations.load(case)
    event['theta_relative_rad']=event.theta_mean_rad-meta['origin_rad']
    event['theta_infinite_bus_rad']=event.theta_mean_rad
    saddle_theta=float(np.mean(saddle[np.array(indices[1:])]-stable[np.array(indices[1:])]))
    event['chord_alpha']=event.theta_relative_rad/saddle_theta
    fault=event[event.stage=='fault'].copy()
    post=event[event.stage.str.startswith('post')].drop_duplicates('t_post_s',keep='last').copy()
    assert np.isclose(fault.t_post_s.iloc[0],-.15) and np.isclose(fault.t_post_s.iloc[-1],0.)
    assert np.allclose(fault[['theta_relative_rad','df_mean_Hz']].iloc[-1],
                       post[['theta_relative_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-8)
    event[event.stage!='pre'][['t_absolute_s','t_post_s','stage','theta_relative_rad','theta_infinite_bus_rad','df_mean_Hz','chord_alpha']].to_csv(D/'event_projection.csv',index=False)
    assert np.allclose(event.theta_relative_rad+stable_phase,event.theta_infinite_bus_rad,rtol=0,atol=1e-12)
    dynamic=states.loc[states.kind=='differential','index'].to_numpy(dtype=int)-1
    pf=[int(states.loc[(states.bus==b)&(states.symbol==f'inv{b}₊droop₊Pfilt'),'index'].iloc[0])-1 for b in [2,3,7]]
    params=pd.read_csv(D/'parameters.csv').set_index('bus')
    checks=[]
    for name,file in [('fault_application','fault_initial.csv'),('clearing_before_switch','fault_end.csv'),
                      ('clearing_after_switch','post030_initial.csv'),('post150','post150_end.csv')]:
        x=np.loadtxt(checkpoint/file,delimiter=',')
        theta=float(np.mean(x[np.array(indices[1:])]-stable[np.array(indices[1:])]))
        frequency=float(np.mean(params.loc[[2,3,7],'Kp'].to_numpy()*(params.loc[[2,3,7],'Pset'].to_numpy()-x[pf])/(2*np.pi)))
        alpha=theta/saddle_theta
        section=stable+alpha*(saddle-stable)
        section[pf]-=frequency*2*np.pi/params.loc[[2,3,7],'Kp'].to_numpy()
        difference=x[dynamic]-section[dynamic]
        checks.append(dict(point=name,theta_relative_rad=theta,theta_infinite_bus_rad=theta+stable_phase,df_mean_Hz=frequency,chord_alpha=alpha,
                           off_section_differential_components=int(np.sum(abs(difference)>1e-8)),
                           on_chord_differential_plane=bool(np.allclose(x[dynamic],section[dynamic],rtol=0,atol=1e-8))))
    before=np.loadtxt(checkpoint/'fault_end.csv',delimiter=',')
    after=np.loadtxt(checkpoint/'post030_initial.csv',delimiter=',')
    assert np.array_equal(before[dynamic],after[dynamic])
    pd.DataFrame(checks).to_csv(D/'event_section_checks.csv',index=False)
    print(pd.DataFrame(checks).to_string(index=False))
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def draw_basin(ax):
    for face in faces:
        labels=set(samples.recapture_type.iloc[face])
        if len(labels)==1:
            ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(labels))],edgecolors='none',alpha=.7))
        else:
            ax.add_collection(PolyCollection([xy[face]],facecolors='#eeeeee',edgecolors='#bbbbbb',linewidths=.3,hatch='////',alpha=.8))
    ax.scatter(xy[:,0],xy[:,1],s=17,c=[colors[c] for c in samples.recapture_type],edgecolors='#444444',linewidths=.5,zorder=3)
    ax.plot(display_angle(0),0,'*',mfc='white',mec='black',ms=12,zorder=7)
    ax.plot(display_angle(meta['saddle_theta_relative_rad']),0,'D',mfc='white',mec='black',ms=6,zorder=7)
    ax.grid(alpha=.15)

def draw_event(ax):
    ax.plot(display_angle(fault.theta_relative_rad),fault.df_mean_Hz,color='#762a83',lw=1.5,zorder=5)
    ax.plot(display_angle(post.theta_relative_rad),post.df_mean_Hz,color='black',lw=1.,zorder=5)
    ax.plot(display_angle(fault.theta_relative_rad.iloc[0]),fault.df_mean_Hz.iloc[0],'^',color='#762a83',mfc='white',ms=6,zorder=8)
    ax.plot(display_angle(post.theta_relative_rad.iloc[0]),post.df_mean_Hz.iloc[0],'o',color='black',mfc='white',ms=6,zorder=8)
    for frame,times,dt,color in [(fault,[-.11,-.05],.015,'#762a83'),(post,[.3,1.,3.,8.],.15,'black')]:
        for time in times:
            a=np.searchsorted(frame.t_post_s,time);b=np.searchsorted(frame.t_post_s,time+dt)
            if b>=len(frame):continue
            ax.annotate('',xy=(display_angle(frame.theta_relative_rad.iloc[b]),frame.df_mean_Hz.iloc[b]),
                        xytext=(display_angle(frame.theta_relative_rad.iloc[a]),frame.df_mean_Hz.iloc[a]),
                        arrowprops=dict(arrowstyle='-|>',lw=.7,color=color,mutation_scale=8),zorder=6)

if args.event_overlay:
    fig,(wide,ax)=plt.subplots(1,2,figsize=(10.5,4.8))
    for axis in [wide,ax]:draw_basin(axis);draw_event(axis)
    wide.set(xlim=(xy[:,0].min()-.02,max(display_angle(1.45),.08) if args.angle_reference=='infinite-bus' else 1.45),ylim=(-.16,.06),
             xlabel=xlabel,ylabel=r'Remote mean $\Delta f$ (Hz)',
             title='Full event projection')
    if args.angle_reference=='infinite-bus':
        wide.axvline(0,color='#777777',ls=':',lw=.8)
        wide.text(0,.057,'Bus6: 0 rad',ha='right',va='top',fontsize=7,color='#555555')
    wide.text(.47,.96,'No basin samples outside\nthe colored rectangle',transform=wide.transAxes,va='top',fontsize=8)
else:
    fig,ax=plt.subplots(figsize=(6.2,4.7));draw_basin(ax)
ax.annotate('Saddle equilibrium\n(on this section)' if chord else 'Saddle projection\n(outside this initial-state slice)',
            xy=(display_angle(meta['saddle_theta_relative_rad']),0),xytext=(display_angle(-.10),.035) if chord else (display_angle(.30),.022),
            fontsize=8,arrowprops=dict(arrowstyle='->',lw=.7))
ax.set(xlim=(xy[:,0].min()-.005,xy[:,0].max()+.005) if chord else (display_angle(-.4),display_angle(1.6)),
       ylim=(-.053,.053) if chord else (-.15,.05),xlabel=xlabel,
       ylabel=r'Remote mean frequency deviation $\Delta f$ (Hz)')
ax.set_title('Local section' if args.event_overlay else ('Equilibrium-chord section' if chord else 'Fixed-hidden-state section')+r', $\tau_p=0.1$ s')
counts=samples.recapture_type.value_counts()
legend=[Patch(facecolor=colors['original_branch'],label=f'Original phase branch ({counts.get("original_branch",0)} points)'),
        Patch(facecolor=colors['shifted_branch'],label=f'Changed phase branch ({counts.get("shifted_branch",0)} points)'),
        Patch(facecolor='#eeeeee',hatch='////',edgecolor='#bbbbbb',label='Mixed sampled branch outcomes'),
        Line2D([],[],marker='*',ls='',mfc='white',mec='black',ms=10,label='Stable equilibrium'),
        Line2D([],[],marker='D',ls='',mfc='white',mec='black',label='Saddle equilibrium' if chord else 'Saddle projection')]
if counts.get('other',0):legend.append(Patch(facecolor=colors['other'],label=f'Unresolved / other ({counts["other"]} points)'))
if args.event_overlay:
    legend += [Line2D([],[],color='#762a83',label='Fault-on projection'),Line2D([],[],color='black',label='Post-fault projection'),
               Line2D([],[],marker='^',ls='',color='#762a83',mfc='white',label='Fault application'),
               Line2D([],[],marker='o',ls='',color='black',mfc='white',label='Clearing')]
    fig.suptitle(r'$\tau_p=0.1$ s: physical event and equilibrium-chord section',fontsize=11)
fig.legend(handles=legend,loc='lower center',ncol=3 if args.event_overlay else 2,fontsize=8,bbox_to_anchor=(.5,.035))
fig.text(.5,.008,'Event curves are coordinate projections; their full states can lie off the section.' if args.event_overlay else
         '35 sampled initial states; stopping at 30/60/150 s. Colors do not certify triangle interiors.',ha='center',fontsize=7)
fig.tight_layout(rect=(0,.25 if args.event_overlay else .22,1,.94 if args.event_overlay else 1))
save(fig,F,'conditional_basin_tau0.10')
held=np.array(meta['held_differential_indices'])-1
count=int(np.sum(abs(saddle[held]-stable[held])>1e-8))
report=dict(points=len(samples),recapture_counts=counts.to_dict(),
            section=args.section,
            angle_reference=args.angle_reference,angle_translation_rad=offset,
            stable_phase_infinite_bus_rad=stable_phase,
            saddle_phase_infinite_bus_rad=stable_phase+meta['saddle_theta_relative_rad'],
            saddle_held_states_differing=count,
            note='Branch labels use unwrapped terminal angles relative to the reference equilibrium. '+
                 ('The chord section contains both full equilibria.' if chord else 'The saddle is outside the fixed-hidden-state section.'))
(D/'phase_branch_summary.json').write_text(json.dumps(report,indent=2)+'\n')
(F/'README.md').write_text(f'''# τp=0.1秒：35点の有限時間回復領域

このフォルダは当該条件の領域図（PDF・PNG）のみを保持します。

- 青：参照平衡点と同じ内部位相の枝へ回復、{counts.get('original_branch',0)}点。
- 橙：異なる内部位相の枝へ回復、{counts.get('shifted_branch',0)}点。
- 灰：未確定・その他、{counts.get('other',0)}点。鞍点そのものはひし形で表示。
- 斜線：三角形の頂点で終端位相の枝が異なる。内部の分類は未計算。
- 星：事故後安定平衡点。ひし形：{'断面上にある全状態の鞍点' if chord else '全状態の鞍点を2軸へ投影した位置'}。

回復は2π周期の基準です。青と橙を合わせて無スリップの安全領域とは呼びません。判定時間は30・60・150秒。三角形内部全点や、その後の長時間挙動は未確認です。

横軸の表示基準は{'無限大母線Bus6（位相0 rad）' if args.angle_reference=='infinite-bus' else '事故後安定平衡点'}です。内部位相の遠方3台平均を表示します。Bus6基準の安定点は{stable_phase:.9f} rad、鞍点は{stable_phase+meta['saddle_theta_relative_rad']:.9f} radです。保存済みのtheta_relative_radは安定点基準のまま保ち、表示時に基準角を加えています。位相を折り返す処理は行いません。

{'52微分状態を安定平衡点から鞍点へ結ぶ方向と、遠方3台の共通周波数方向で断面を定義します。残り46状態も位置に応じて変わります。14代数変数は整合化します。α=0・Δf=0で安定点、α=1・Δf=0で鞍点そのものです。' if chord else '固定した46微分状態のうち27状態が鞍点では異なるため、ひし形の位置にこの断面の初期値を作っても全状態の鞍点にはなりません。'}

再作図：リポジトリ直下で `python code/plot_tau010_basin.py --section {args.section} --angle-reference {args.angle_reference}`。
''')
if args.event_overlay:
    with (F/'README.md').open('a') as f:
        f.write('\n事故中・事故後の軌道を追加した表示です。左は全体、右は計算した断面の局所拡大。紫線は事故中、黒線は事故後、三角が投入、白丸が除去です。軌道は同じθ・Δfでの座標投影であり、全状態が断面上にあることを意味しません。着色していない場所の回復領域は未計算です。\n\n再現時は上のコマンドに `--event-overlay` を追加してください。\n')
print(json.dumps(report,indent=2))
print('Created only taup=0.1 s basin figure:',F)
