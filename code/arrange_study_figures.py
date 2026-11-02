"""Curated original/voltage-adjusted figures; archive old generated layouts."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

R=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--original-only',action='store_true')
a=parser.parse_args()
F=R/'figures';O=F/'01_original';V=F/'02_voltage_adjusted'
ARCHIVE=R/'basin_recalculation/figure_archive_before_case_layout'
for folder in [O,V]:folder.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def save(fig,folder,name):
    for ext in ['png','pdf']:fig.savefig(folder/f'{name}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def trajectories(folder,adjusted):
    D=R/'basin_recalculation/capacitor_probe_predictor'
    fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
    for tau,suffix,color in [(.1,'010','#0072b2'),(1.2,'120','#d55e00')]:
        case=('v095' if adjusted else 'v090')+'_tau'+suffix
        d=pd.read_csv(D/f'{case}.csv');f=pd.read_csv(D/f'{case}_fault.csv')
        label=rf'$\tau_p={tau:g}$ s'
        axes[0].plot(d.theta_infinite_bus_rad,d.df_mean_Hz,color=color,lw=1,label=label)
        axes[0].plot(f.theta_infinite_bus_rad,f.df_mean_Hz,color=color,lw=1.5,ls='--')
        axes[0].plot(d.theta_infinite_bus_rad.iloc[0],d.df_mean_Hz.iloc[0],'o',color=color,ms=4)
        time=np.r_[f.t_absolute_s.to_numpy()-1.15,d.t_post_s.to_numpy()]
        phase=np.r_[f.theta_infinite_bus_rad.to_numpy(),d.theta_infinite_bus_rad.to_numpy()]
        freq=np.r_[f.df_mean_Hz.to_numpy(),d.df_mean_Hz.to_numpy()]
        axes[1].plot(time,phase,color=color,lw=1,label=label)
        axes[2].plot(time,freq,color=color,lw=1,label=label)
    axes[0].set(xlabel='Mean internal phase / Bus6 [rad]',ylabel='Mean frequency deviation [Hz]',title='Phase-plane projection')
    axes[1].set(xlabel='Time after clearing [s]',ylabel='Mean internal phase / Bus6 [rad]',title='Internal phase; 150-s observation',xlim=(-.15,150))
    axes[2].set(xlabel='Time after clearing [s]',ylabel='Mean frequency deviation [Hz]',title='Frequency; 150-s observation',xlim=(-.15,150))
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8)
    for ax in axes[1:]:ax.axvline(0,color='black',ls=':',lw=.6)
    fig.suptitle('Voltage adjusted: pre V5=0.950 pu, B=0.55618 pu' if adjusted else 'Original: pre V5=0.90364 pu, B=0')
    save(fig,folder,'02_trajectories_tau010_tau120' if adjusted else '03_trajectories_tau010_tau120')

from plot_original_large_tau import restore_large_tau, update_readmes

trajectories(O,False)
restore_large_tau()
for ext in ['png','pdf']:
    src=F/'tau010_chord_basin'/f'conditional_basin_tau0.10.{ext}'
    if not src.exists():src=ARCHIVE/'figures/tau010_chord_basin'/src.name
    assert src.exists()
    shutil.copyfile(src,O/f'01_equilibria_section_tau010.{ext}')
if a.original_only:raise SystemExit(0)
trajectories(V,True)
latest=R/'basin_recalculation/B100_tau_comparison'
if all((latest/f'basin_tau{suffix}.csv').exists() and len(pd.read_csv(latest/f'basin_tau{suffix}.csv'))==83 for suffix in ['050','100']):
    import subprocess,sys
    subprocess.run([sys.executable,str(R/'code/check_plot_B100_tau_comparison.py')],check=True)
gain_latest=R/'basin_recalculation/B100_Kp01_tau_comparison'
if all((gain_latest/f'basin_tau{suffix}.csv').exists() and len(pd.read_csv(gain_latest/f'basin_tau{suffix}.csv'))==83 for suffix in ['050','100']):
    import subprocess,sys
    subprocess.run([sys.executable,str(R/'code/check_plot_B100_Kp01.py')],check=True)
if (R/"basin_recalculation/capacitor_equilibrium_search/summary.json").exists():
    from plot_capacitor_equilibrium_search import main as plot_equilibrium_search
    plot_equilibrium_search()
# Reversible archive: keep only two case folders and their common index active.
for item in list(F.iterdir()):
    if item.name in ['01_original','02_voltage_adjusted','README.md']:continue
    dest=ARCHIVE/'figures'/item.name;dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        serial=2
        while (ARCHIVE/'regenerated'/str(serial)/item.name).exists():serial+=1
        dest=ARCHIVE/'regenerated'/str(serial)/item.name
        dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.move(str(item),str(dest))
for name in ['TRAJECTORY_VOLTAGE_PROBES.md','CAPACITOR_PROBE.md']:
    src=R/'docs'/name;dest=ARCHIVE/'notes'/name
    if src.exists():dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(src),str(dest))
for name in ['tau010_trajectory_neighborhood.png','voltage_load_probe.png','capacitor_responses.png']:
    src=R/'docs/figures'/name;dest=ARCHIVE/'notes/figures'/name
    if src.exists():dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(src),str(dest))
(O/'README.md').write_text('''# 元の設定

B=0、事故前Bus5電圧0.90364 pu、負荷G=2.2682301152 pu。GFMのP指令[0.5,0.4,0.4,0.4] pu、τq=0.02秒。Bus4の20 puシャント事故0.15秒後、除去と同時に線路4–5開放。

- [01_equilibria_section_tau010](01_equilibria_section_tau010.png)（[PDF](01_equilibria_section_tau010.pdf)）：安定点と鞍点の両方を通る断面＋事故中・事故後軌道の投影。τp=0.1秒。断面外の軌道を投影している。
- [03_trajectories_tau010_tau120](03_trajectories_tau010_tau120.png)（[PDF](03_trajectories_tau010_tau120.pdf)）：τp=0.1秒は回復、1.2秒は減速方向の持続回転。破線は事故中、丸は除去時。

各図はPNG/PDF。位相は無限大母線Bus6基準、折り返しなし。領域は有限時間の格子分類で、塗りの全点を保証しない。
''')
(V/'README.md').write_text('''# コンデンサによる電圧調整後

元のk=1.0299の代表例：事故前Bus5電圧0.95000 pu、並列B=0.55618062 pu。コンデンサは事故前・中・後すべて接続。負荷G・GFMのP指令・線路・事故操作は元の設定と共通。事故前の負荷Pは1.85217→2.04707 pu、コンデンサ供給Qは0.50195 pu。

- [02_trajectories_tau010_tau120](02_trajectories_tau010_tau120.png)（[PDF](02_trajectories_tau010_tau120.pdf)）：代表例だけの位相平面・位相時系列・周波数時系列。τp=0.1秒・1.2秒を比較。破線は事故中、丸は除去時。

元のk=1.0299では事故後安定平衡点が未確認で、平衡点基準の引き込み領域は未計算。kを下げた比較は下の別節に記す。電圧だけでなく負荷Pも増えており、電圧だけの効果は切り分けていない。

各図はPNG/PDF。位相はBus6基準。

''')
(F/'README.md').write_text('''# 検討別の図

| フォルダ | 条件 | 主に見る図 |
|---|---|---|
| [01_original](01_original/README.md) | 元の設定、コンデンサなし、事故前V5=0.90364 pu | 引き込み領域＋事故軌道、τの異なる実事故応答 |
| [02_voltage_adjusted](02_voltage_adjusted/README.md) | 常時接続コンデンサで事故前V5=0.95000 pu | 補償条件のk比較、k=0.80の引き込み領域 |

図の位相・周波数は遠方3台（Bus2・3・7）の平均内部量。位相はBus6基準。引き込み領域と事故中・事故後の実軌道を中心に表示する。

再作図：`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/arrange_study_figures.py`。

trajectory_neighborhoodは全点が同じ分類で境界情報が得られなかったため一覧から外した。図と数値は`basin_recalculation/figure_archive_trajectory_neighborhood/`、生計算は元ディレクトリに保持している。

従来の重複図・途中の比較図と旧レポートは`basin_recalculation/figure_archive_before_case_layout/`へ移動した。生の計算データは保持している。今後の図はこの2分類で整理する。
''')
update_readmes()
from plot_case_networks import plot_case_networks
plot_case_networks()
print('PASS: two curated figure folders; legacy figures archived without deleting numerical data')
