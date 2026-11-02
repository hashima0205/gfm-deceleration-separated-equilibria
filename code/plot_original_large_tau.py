"""Restore the original large-tau basin/fault figures with Bus6 angle reference."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

R=Path(__file__).resolve().parents[1]

def restore_large_tau():
    D=R/'results/conditional_basins/data';F=R/'figures/01_original';F.mkdir(exist_ok=True,parents=True)
    meta=json.loads((D/'definition.json').read_text());origin=meta['origin_rad']
    fault=pd.read_csv(R/'data/fault_trajectories.csv.gz')
    colors={'R':'#9ecae1','N':'#fdae6b','P':'#bc80bd','U':'#d9d9d9'}
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    data={};counts=[]
    for tau in meta['taus']:
        s=pd.read_csv(D/f'plot_samples_tau{tau:.2f}.csv')
        faces=pd.read_csv(D/f'plot_triangles_tau{tau:.2f}.csv').to_numpy(int)
        xy=s[['theta_relative_rad','df_Hz']].to_numpy();xy[:,0]+=origin
        labels=np.array(['R' if c.startswith('recovery') else {'negative_running':'N','positive_running':'P'}.get(c,'U') for c in s.classification])
        p=pd.read_csv(D/f'trajectories/fault_tau{tau:.2f}.csv.gz')
        f=fault[np.isclose(fault.taup,tau)].sort_values('t_post_s')
        assert len(s)==88 and len(f)>1
        assert np.allclose(f[['theta_relative_rad','df_mean_Hz']].iloc[-1],p[['theta_relative_rad','df_mean_Hz']].iloc[0],rtol=0,atol=1e-7)
        assert np.isclose(f.t_post_s.iloc[0],-.15) and np.isclose(f.t_post_s.iloc[-1],0)
        data[tau]=(s,xy,labels,faces,p,f)
        counts.append({'taup':tau,'basin_points':len(s),'fault_on_samples':len(f),'post_horizon_s':p.t_post_s.iloc[-1]})
    def panel(ax,tau,xlim,ylim):
        s,xy,labels,faces,p,f=data[tau]
        for face in faces:
            cc=set(labels[face]);same=len(cc)==1
            ax.add_collection(PolyCollection([xy[face]],facecolors=colors[next(iter(cc))] if same else '#eeeeee',edgecolors='none' if same else '#aaaaaa',lw=.25,alpha=.7))
        ax.scatter(xy[:,0],xy[:,1],c=[colors[c] for c in labels],s=8,edgecolors='#555555',lw=.2)
        ax.plot(p.theta_relative_rad+origin,p.df_mean_Hz,color='#222222',lw=1)
        ax.plot(f.theta_relative_rad+origin,f.df_mean_Hz,color='#762a83',lw=1.7)
        ax.plot(f.theta_relative_rad.iloc[0]+origin,f.df_mean_Hz.iloc[0],'^',mfc='white',mec='#762a83',ms=6)
        ax.plot(p.theta_relative_rad.iloc[0]+origin,p.df_mean_Hz.iloc[0],'o',mfc='white',mec='black',ms=5)
        ax.plot(origin,0,'*',mfc='white',mec='black',ms=10)
        ax.plot(origin+meta['saddle_theta_relative_rad'],0,'D',mfc='none',mec='black',ms=5)
        ax.set(xlim=xlim,ylim=ylim,xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title=rf'$\tau_p={tau:g}$ s; 88 basin points')
        ax.grid(alpha=.2)
    handles=[Patch(facecolor=colors['R'],label='Recovery'),Patch(facecolor=colors['N'],label='Negative rotation'),Patch(facecolor='#eeeeee',edgecolor='#aaaaaa',label='Mixed sampled outcomes'),
             Line2D([],[],color='#222222',label='Post-fault projection'),Line2D([],[],color='#762a83',label='Fault-on projection'),
             Line2D([],[],marker='*',mfc='white',mec='black',ls='',ms=10,label='Stable equilibrium'),Line2D([],[],marker='D',mfc='none',mec='black',ls='',label='Saddle projection'),
             Line2D([],[],marker='o',mfc='white',mec='black',ls='',label='Clearing projection')]
    def save(fig,name):
        for ext in ['png','pdf']:fig.savefig(F/f'{name}.{ext}',dpi=180,bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.5))
    for ax,tau in zip(axes,meta['taus']):panel(ax,tau,(origin-.4,origin+1.6),(-.15,.05))
    fig.suptitle('Original setting: conditional recovery regions and fault-on / post-fault trajectories')
    fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=8,bbox_to_anchor=(.5,.025))
    fig.text(.5,.005,'46 other initial differential states fixed at the stable equilibrium; trajectories may lie outside the section.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.18,1,.94));save(fig,'04_basins_large_tau')
    fig,axes=plt.subplots(2,3,figsize=(13,7))
    for col,tau in enumerate(meta['taus']):
        panel(axes[0,col],tau,(origin-.2,origin+.15),(-.006,.006))
        panel(axes[1,col],tau,(origin+1.25,origin+1.5),(-.025,-.01))
    fig.suptitle('Original setting: equilibrium neighborhood (top), clearing neighborhood (bottom)')
    fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=8,bbox_to_anchor=(.5,.005))
    fig.tight_layout(rect=(0,.10,1,.95));save(fig,'05_basins_large_tau_zoom')
    fig,axes=plt.subplots(1,2,figsize=(9,3.7),constrained_layout=True)
    for tau in meta['taus']:
        f=data[tau][-1];label=rf'$\tau_p={tau:g}$ s'
        axes[0].plot(f.theta_relative_rad+origin,f.df_mean_Hz,label=label)
        axes[1].plot(f.t_post_s,f.df_mean_Hz,label=label)
    axes[0].set(xlabel='Mean internal phase relative to Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title='Fault-on trajectory')
    axes[1].set(xlabel='Time relative to clearing [s]',ylabel='Mean frequency deviation [Hz]',title='Fault-on frequency')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    save(fig,'06_fault_on_large_tau')
    pd.DataFrame(counts).to_csv(F/'large_tau_calculation_status.csv',index=False)
    update_readmes()
    print('PASS: restored 264 basin points and three verified fault-on trajectories; Bus6 angle reference')

def gain_current_section():
    path=R/'figures/02_voltage_adjusted/B100_Kp01_tau_comparison_validation.json'
    if not path.exists():return ''
    q=json.loads(path.read_text())
    text="""## ドループ1%と2%の比較（同じB=1.0、k）

1%・2%は共通電力ベースに対する値で、インバータ定格ベースとの一致は仮定していない。全4台のKpだけを7.53982237から3.76991118 rad/s/puへ半減した。事故前・事故後の平衡点は共通。事故操作も同一で、300秒の応答を比較した。

![1%と2%の事故応答](12_trajectories_B100_Kp_comparison.png)

τp=0.5秒では両ゲインで回復。τp=1.0秒では2%が減速方向に持続回転する一方、1%は回復した。

![1%での引き込み領域と事故軌道](09_basins_B100_Kp01_tau050_tau100.png)

![フィルタ電力状態を共通座標としたゲイン比較](11_basins_B100_Kp_state_comparison.png)

位相・周波数を同じ格子にすると、Kpが異なるため初期Pfは一致しない。図11では縦軸を共通のPset−Pfに変更した。格子の範囲・密度はこの座標では異なるので、同一の全66状態を持つ点の組も別途照合した。事故軌道は断面への投影である。

| τp [s] | 1%：参照枝回復 | 1%：別枝回復 | 1%：負方向回転 | 1%：未確定 | 全状態が一致する比較点数 |
|---|---:|---:|---:|---:|---:|
"""
    for tau in [.5,1.]:
        c=q['gain_comparison'][str(tau)];n=c['mp01_counts']
        text+=f"| {tau:g} | {n.get('same',0)} | {n.get('shift',0)} | {n.get('negative',0)} | {n.get('unknown',0)} | {c['matched_full_state_points']} |\n"
    text+="""
各τ83点（鞍点そのもの1点を含む）。未確定は300秒以内に判定基準を満たさなかった点で、脱調とは区別する。点数は領域の面積・体積を表さない。

同じ全66状態の48点ずつを比較したところ、τp=1.0秒で2%が負方向回転した25点は、1%で7点回復、16点負方向回転、2点未確定になった。τp=0.5秒で2%が負方向回転した17点は、1%で3点が別位相枝へ回復、14点が未確定。したがって、改善例は確認できるが、未確定点を含めた領域全体の包含関係までは確定していない。

電圧一定・非飽和の単純化では、τ δ̈+δ̇=Kp(Pset−P(δ))。小振幅の固有角周波数は√(Kp P′/τ)、減衰比は1/(2√(τ Kp P′))。Kpを下げる場合とτを上げる場合はともに応答を遅くするが、減衰比への作用は逆向き。この近似だけから多次元の引き込み領域の単調性までは決められない。

[監査結果](B100_Kp01_tau_comparison_validation.json) ／ [全状態が一致する点の分類比較](B100_Kp01_matched_full_state_comparison.csv) ／ [事故応答PDF](12_trajectories_B100_Kp_comparison.pdf)

再計算：`julia --project=provenance code/run_B100_Kp01.jl`（全条件を1プロセスで実行）。監査・再作図：`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/check_plot_B100_Kp01.py`。生データ：`basin_recalculation/B100_Kp01_tau_comparison/`。

"""
    return text

def b100_current_section():
    path=R/'figures/02_voltage_adjusted/B100_tau_comparison_validation.json'
    if not path.exists():return ''
    q=json.loads(path.read_text())
    events=pd.read_csv(R/'figures/02_voltage_adjusted/B100_tau_comparison_events.csv')
    names={'recovery':'回復','recovery_after_slip':'スリップ後回復','negative_running':'減速方向の持続回転','unresolved':'未確定'}
    text=f'''<!-- B100-current:start -->
## B=1.0、ドループ2%でのτp=0.5秒と1.0秒

k={q['k']:.10f}、補償B=1.0 pu。事故前の負荷端電圧{q['pre_V5_pu']:.5f} pu、事故後安定点{q['post_eq_V5_pu']:.5f} pu。対応する定常枝の折返しk≈0.68134116の0.01%手前を使用。事故前電圧は1 puに近いが、系統条件はなお平衡点消失の限界に近い。

![B=1.0の引き込み領域と事故軌道：τ=0.5秒／1.0秒](07_basins_B100_tau050_tau100.png)

[領域図PDF](07_basins_B100_tau050_tau100.pdf) ／ [事故応答の時系列](08_trajectories_B100_tau050_tau100.png) ／ [系統図](equilibrium_network_B100.png) ／ [分類・監査要約](B100_tau_comparison_validation.json)

| τp [s] | 回復：参照枝 | 回復：異なる位相枝 | 負方向回転 | 未確定 | 今回の実事故 |
|---|---:|---:|---:|---:|---|
'''
    for tau in [.5,1.]:
        c=q['basins'][str(tau)]['counts'];event=events[np.isclose(events.taup,tau)].iloc[0]
        text+=f"| {tau:g} | {c.get('same',0)} | {c.get('shift',0)} | {c.get('negative',0)} | {c.get('unknown',0)} | {names.get(event.classification,event.classification)} |\n"
    text+='''
各τに鞍点そのものの初期値を別に1点含む。上段は広域格子、下段は平衡点・鞍点の近傍を細分した表示。両τで全66状態が同一の83初期値を使用した。52微分状態は安定点と鞍点を通る弦で変化させ、遠方3台の周波数を共通にずらす。広域には弦の外挿を含む。代数整合化後は全状態を積分し、30→60→150→300秒の段階分類を行った。

斜線は頂点の分類が混在する三角形で、境界の厳密な位置を表さない。分類点は非一様に配置しており、点数比を引き込み領域の面積比とは解釈しない。実事故軌道は断面外の投影。事故軌道が回復色の上に見えることだけでは、回復の可否は決まらない。

事故操作はBus4の20 puシャント事故0.15秒、除去時に線路4–5開放。P指令[0.5,0.4,0.4,0.4] pu、負荷G=2.2682301152 pu、τq=0.02秒を共通にした。事故前の負荷Pは約2.25281 pu、補償Qは約0.99320 pu。Bの妥当性を電圧だけで判断した試験ではなく、今回は電圧を1 pu近傍にした設定でτを比較した。

再計算：`JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=provenance code/run_B100_tau_comparison.jl`。監査・再作図：`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/check_plot_B100_tau_comparison.py`。生データは`basin_recalculation/B100_tau_comparison/`。各条件でJuliaを再起動せず、1プロセスで処理する。

## 以前の比較条件（B≈0.55618）
<!-- B100-current:end -->
'''
    text=text.replace('## 以前の比較条件（B≈0.55618）',gain_current_section()+'## 以前の比較条件（B≈0.55618）')
    root=R/'figures/README.md'
    if root.exists():
        s=root.read_text();old='| [02_voltage_adjusted](02_voltage_adjusted/README.md) | 常時接続コンデンサで事故前V5=0.95000 pu | 補償条件のk比較、k=0.80の引き込み領域 |'
        new='| [02_voltage_adjusted](02_voltage_adjusted/README.md) | 今回：B=1.0、k=0.681273、事故前V5≈0.997 pu（以前の条件も併記） | τp=0.5／1.0秒の引き込み領域と事故軌道 |'
        s=s.replace(old,new)
        note='最新の比較：[B=1.0・τp=0.5／1.0秒の領域図](02_voltage_adjusted/07_basins_B100_tau050_tau100.png)。\n\n'
        if '最新の比較：' not in s:s=s.replace('図の位相・周波数は遠方3台',note+'図の位相・周波数は遠方3台',1)
        if gain_current_section():
            s=s.replace(note,'最新の比較：[ドループ1%と2%の事故応答](02_voltage_adjusted/12_trajectories_B100_Kp_comparison.png) ／ [共通状態座標の引き込み領域](02_voltage_adjusted/11_basins_B100_Kp_state_comparison.png)。\n\n')
        root.write_text(s)
    return text

def update_readmes():
    original="""

## 大きいτの検討と計算状況

- [04_basins_large_tau](04_basins_large_tau.png)（[PDF](04_basins_large_tau.pdf)）：τp=1.16・1.18・1.20秒の引き込み領域と事故中・事故後軌道。各88点、計264点。
- [05_basins_large_tau_zoom](05_basins_large_tau_zoom.png)（[PDF](05_basins_large_tau_zoom.pdf)）：平衡点近傍と事故除去点近傍の拡大。
- [06_fault_on_large_tau](06_fault_on_large_tau.png)（[PDF](06_fault_on_large_tau.pdf)）：事故中0.15秒の軌道と周波数を単独表示。

| τp [s] | 平衡点基準の条件付き引き込み領域 | 事故中・事故後軌道 |
|---|---|---|
| 0.1 | 計算済み：両平衡点を通る35点 | 計算済み |
| 1.16・1.18・1.20 | 計算済み：各88点 | 計算済み |

大きいτの断面は残り46状態を安定点に固定した断面。鞍点は投影であり、この断面上の点ではない。両平衡点を通るτp=0.1秒の断面とは定義が異なる。
"""
    adjusted="""

## 引き込み領域の計算状況（元のk=1.0299）

| τp [s] | 事故後安定平衡点 | 平衡点基準の引き込み領域 | 事故中軌道 | 事故後軌道 |
|---|---|---|---|---|
| 0.1 | 再探索でも未発見 | 基準安定点がなく未計算 | 計算済み：0.15秒 | 計算済み：150秒 |
| 1.2 | 再探索でも未発見 | 基準安定点がなく未計算 | 計算済み：0.15秒 | 計算済み：150秒 |

対象は代表例（事故前V5=0.95 pu、B=0.55618062 pu）。平衡点の位置はτpに依存しないため、定常方程式の探索は両τで共通。安定性と引き込み領域はτに依存する。
"""
    status=R/'basin_recalculation/capacitor_equilibrium_search/summary.json'
    if status.exists():
        q=json.loads(status.read_text())
        adjusted+=f"""
- [03_postfault_equilibrium_search](03_postfault_equilibrium_search.png)（[PDF](03_postfault_equilibrium_search.pdf)）：元の2平衡点の枝の追跡と、現在のBにおける多初期値探索。

既知の全状態モデルの6平衡点と照合した定常方程式（最大残差{q['maximum_validation_residual']:.2e}）を使用。{q['seed_count']}初期値から探索し、現在のBでは平衡点候補{q['distinct_reduced_roots']}個。元の安定点と鞍点につながる枝の局所的な折返しはB={q.get('local_fold_B_pu',float('nan')):.9f} pu。現在のBはこれより大幅に大きい。

全平衡点の不存在を証明したものではない。回復先が確認できないため、領域を回復色で塗る計算は実施していない。探索データは`basin_recalculation/capacitor_equilibrium_search/`、要約はこのフォルダの`equilibrium_search_status.json`。現条件で引き込み領域の計算が完了したという扱いにはしない。
"""
    strength=R/'figures/02_voltage_adjusted/stronger_grid_validation.json'
    if strength.exists():
        q=json.loads(strength.read_text())
        events=pd.read_csv(R/'basin_recalculation/capacitor_k_strength/events.csv')
        adjusted+="""

## kを下げて系統を強めた比較

コンデンサB=0.55618062 pu、負荷G=2.2682301152 pu、GFMのP指令と事故操作を維持し、線路4–6・5–6のR/Xだけを変更。事故前電圧を再び0.95 puへ合わせ直す操作はしていない。

| k | 事故前V5 [pu] | 事故後安定点V5 [pu] | τp=0.1秒の実事故 | τp=1.2秒の実事故 |
|---:|---:|---:|---|---|
"""
        labels={'recovery':'回復','recovery_after_slip':'回復（スリップあり）','negative_running':'減速方向の持続回転','unresolved':'未確定'}
        for k in sorted(events.k.unique(),reverse=True):
            d=events[events.k==k];lo=d[np.isclose(d.taup,.1)].iloc[0];hi=d[np.isclose(d.taup,1.2)].iloc[0]
            adjusted+=f"| {k:g} | {lo.pre_V5_pu:.5f} | {lo.post_eq_V5_pu:.5f} | {labels.get(lo.classification,lo.classification)} | {labels.get(hi.classification,hi.classification)} |\n"
        adjusted+="""

- [04_stronger_grid_trajectories](04_stronger_grid_trajectories.png)（[PDF](04_stronger_grid_trajectories.pdf)）：k=0.80・0.75・0.70の実事故応答。事故中から事故後150秒まで計算。

全状態モデルで各k・両τの事故前安定点・事故後安定点・鞍点を確認。独立の定常方程式とも照合した。元のk=1.0299で平衡点未発見という上の結果とは条件を区別する。

定常枝の折返しはk≈0.807313。k=0.80はこの境界より約0.9%小さく、なお限界に近い設定である。事故前電圧0.95 pu以上でも安定平衡点と減速脱調が共存する例にはなったが、広い運転範囲で起こることを示したものではない。
"""
        for tau,suffix in [(.1,'010'),(1.2,'120')]:
            key=f'tau{tau:g}'
            if key in q.get('basins',{}):
                name=f'05_basin_k080_tau{suffix}'
                adjusted+=f"\n- [{name}]({name}.png)（[PDF]({name}.pdf)）：k=0.80、τp={tau:g}秒の引き込み領域と事故中・事故後軌道。55点。\n"
            else:
                adjusted+=f"\nτp={tau:g}秒の引き込み領域は計算中。\n"
        if q.get('basins'):
            adjusted+="""

両τで全66状態が同一の初期値55点を使用。52微分状態を安定点と鞍点の弦で変化させ、遠方3台の周波数を共通にずらした後、14代数変数を整合化する。両平衡点そのものを含む。事故軌道はこの断面外への投影であり、投影が回復色の上にあってもその色から実事故の回復は判定しない。分類は30秒、未確定のみ60・150秒へ延長する有限時間の格子計算。

再計算：`julia --project=provenance code/run_capacitor_k_strength.jl`（1プロセス）。監査・再作図：`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/check_plot_capacitor_k_strength.py`。生データは`basin_recalculation/capacitor_k_strength/`、数値表と監査要約はこのフォルダに保存。
"""
    fastpath=R/'figures/02_voltage_adjusted/fast_tau_search_validation.json'
    if fastpath.exists():
        fast=json.loads(fastpath.read_text())
        adjusted+="""

## τp=0.1秒の追加探索

- [06_fast_tau_search](06_fast_tau_search.png)（[PDF](06_fast_tau_search.pdf)）：kを折返しへ近づける探索、補償シャントBの変更、事故シャント強度の変更を集約。
- [探索した25条件の数値](fast_tau_search_results.csv) ／ [検証要約](fast_tau_search_validation.json)

固定G・P指令のまま、B=0.55618でk=0.804～0.80731の6条件、B=0.8・1.0・1.2で各折返しの1%・0.1%・0.01%手前の9条件、事故シャント強度2・5・10・40・80 puの追加10条件を確認した。事故時間はすべて0.15秒。事故前負荷端電圧は0.9548～1.0159 pu、各条件の事故前・事故後平衡点は全状態モデルで局所安定だった。

25条件すべて回復。必要に応じて30→60→150→300秒へ延長したが、この範囲ではτp=0.1秒の持続する減速脱調は未発見。一般的な不存在を示した結果ではない。元のk=0.80・τp=1.2秒の減速脱調例は引き続き成立する。

再現：`code/prepare_fast_tau_search.py`で条件を準備し、`basin_recalculation/fast_tau_search/finish.flag`を置いてから`julia --project=provenance code/run_fast_tau_search.jl`を実行する。モデル構築は1回。finish.flagがない場合は条件CSVへの追加を待つため、同じプロセスに追加条件を投入できる。監査・再作図は`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/check_plot_fast_tau_search.py`。
"""
        if fast.get('negative_example') is not None:
            adjusted+="\n[τp=0.1秒の減速脱調候補](07_fast_tau_example.png)を参照。\n"
    for folder,heading,addition in [('01_original','\n\n## 大きいτの検討と計算状況',original),('02_voltage_adjusted','\n\n## 引き込み領域の計算状況（元のk=1.0299）',adjusted)]:
        path=R/'figures'/folder/'README.md'
        if path.exists():
            raw=path.read_text()
            start='<!-- B100-current:start -->';end='<!-- B100-current:end -->'
            if start in raw:
                a=raw.index(start);b=raw.index(end,a)+len(end);raw=raw[:a]+raw[b:]
            text=raw.split('\n\n## 引き込み領域の計算状況')[0].split(heading)[0].split('\n\n## 計算済み／未計算')[0]
            text=text+addition
            current=b100_current_section() if folder=='02_voltage_adjusted' else ''
            if current:
                title,rest=text.split('\n',1);text=title+'\n\n'+current+'\n'+rest.lstrip('\n')
            path.write_text(text)

if __name__=='__main__':restore_large_tau()
