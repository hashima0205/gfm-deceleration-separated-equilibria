"""Pre/post voltage and P/Q flow diagrams for the two curated study settings."""
from pathlib import Path
import re
R=Path(__file__).resolve().parents[1]

def update_network_readme():
    p=R/'figures/README.md';text=p.read_text();start='<!-- network-panels:start -->';end='<!-- network-panels:end -->'
    block='''<!-- network-panels:start -->

## 事故前後の電圧・潮流

各バスに端子電圧の大きさ、Bus6基準の端子位相、系統へのP/Q注入を表示。P/Qとも正は供給、負は消費。矢印は主要線路のPの向きで、P/Qの数値はその向きに沿う送電端→受電端の値（符号付き）。Qが負なら矢印とは逆向きの供給を表す。

事故後パネルは事故後系統の安定平衡点。事故除去直後ではなく、τp=1.2秒で減速脱調する応答の到達点でもない。平衡点の位置は両τで共通。

### 01：元の設定

k=1.0299、補償シャントなし、負荷G=2.2682301152 pu。

![01：事故前と事故後安定点の電圧・P/Q潮流](01_original/equilibrium_network.png)

[PDF](01_original/equilibrium_network.pdf) ／ [バス数値](01_original/network_data/buses.csv) ／ [全線路潮流](01_original/network_data/lines.csv)

### 02：電圧補償あり・k=0.80

B=0.55618062 pu、k=0.80、GとP指令は01と共通。事故前の負荷端電圧は約0.955 pu、事故後安定点では約0.848 pu。元のk=1.0299・B=0.55618062の条件は安定な事故後平衡点が未確認のため、この事故後パネルには使わない。

![02：事故前と事故後安定点の電圧・P/Q潮流](02_voltage_adjusted/equilibrium_network.png)

[PDF](02_voltage_adjusted/equilibrium_network.pdf) ／ [バス数値](02_voltage_adjusted/network_data/buses.csv) ／ [全線路潮流](02_voltage_adjusted/network_data/lines.csv)

再作図：`MPLCONFIGDIR=/tmp/gfm-matplotlib .venv/bin/python code/plot_case_networks.py`。
<!-- network-panels:end -->'''
    if (R/'figures/02_voltage_adjusted/equilibrium_network_B100.png').exists():
        block=block.replace('<!-- network-panels:end -->',"""
### 02：今回のB=1.0・τ比較

B=1.0 pu、k=0.68127303（対応する折返しk≈0.68134116）。事故前負荷端電圧0.99659 pu、事故後安定点0.85708 pu。τp=0.5秒と1.0秒の平衡点は共通。

![02 B=1.0：事故前と事故後安定点の電圧・P/Q潮流](02_voltage_adjusted/equilibrium_network_B100.png)

[PDF](02_voltage_adjusted/equilibrium_network_B100.pdf) ／ [バス数値](02_voltage_adjusted/network_data_B100/buses.csv) ／ [全線路潮流](02_voltage_adjusted/network_data_B100/lines.csv)
<!-- network-panels:end -->""")
    if start in text:
        a=text.index(start);b=text.index(end,a)+len(end);text=text[:a]+block+text[b:]
    else:text=text.rstrip()+'\n\n'+block+'\n'
    p.write_text(text)

def plot_case_networks():
    source=(R/'code/plot_equilibrium_network.py').read_text()
    cases=[
        ('01_original',1.0299,0.,R/'data/seeds/nf02_1_pre0.csv',R/'data/seeds/nf02_1_post0.csv','01 Original: k=1.0299, B=0'),
        ('02_voltage_adjusted',.8,.5561806200995012,R/'basin_recalculation/capacitor_k_strength/k080_pre.csv',R/'basin_recalculation/capacitor_k_strength/k080_stable.csv','02 Voltage adjusted: k=0.80, B=0.55618')]
    latest=R/'basin_recalculation/B100_tau_comparison'
    if (latest/'stable.csv').exists():
        cases.append(('02_voltage_adjusted',.6812730259441179,1.,latest/'pre.csv',latest/'stable.csv','02 B=1.0 comparison: k=0.68127303'))
    for folder,k,B,pre,post,label in cases:
        code=source.replace("FIG=ROOT/'docs/figures'",f"FIG=ROOT/'figures/{folder}'")
        code=code.replace("OUT=ROOT/'data/equilibrium_power_flow'","OUT=FIG/'network_data'")
        code=code.replace("meta=json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text())",f"meta=json.loads((ROOT/'results/conditional_basins/data/definition.json').read_text());meta['k']={k!r};CAP_B={B!r}")
        code=code.replace("for stage,file in [('pre','nf02_1_pre0.csv'),('post','nf02_1_post0.csv')]:",f"for stage,file in [('pre',Path({str(pre)!r})),('post',Path({str(post)!r}))]:")
        code=code.replace("x=np.loadtxt(ROOT/'data/seeds'/file,delimiter=',')","x=np.loadtxt(file,delimiter=',')")
        code=code.replace("P_from_pu=sa.real,P_to_pu=sb.real,loss_pu=(sa+sb).real)","P_from_pu=sa.real,P_to_pu=sb.real,Q_from_pu=sa.imag,Q_to_pu=sb.imag,loss_pu=(sa+sb).real,Q_loss_pu=(sa+sb).imag)")
        code=code.replace("assert abs(injections[5].imag)<1e-8","assert abs(injections[5].imag-CAP_B*abs(U[5])**2)<1e-8")
        code=code.replace("load_G_pu=load/abs(U[5])**2)","load_G_pu=load/abs(U[5])**2,cap_Q_pu=injections[5].imag,cap_B_pu=injections[5].imag/abs(U[5])**2)")
        code=code.replace("roles={1:'GFM'","roles={1:'GFM'")
        code=code.replace("text=f'{sent:.3f} → {received:.3f}'", "qfrom=float(branch_table.loc[j,'Q_from_pu']);qto=float(branch_table.loc[j,'Q_to_pu'])\n            qsent,qreceived=(qfrom,-qto) if pfrom>0 else (qto,-qfrom)\n            text=f'P {sent:+.3f} → {received:+.3f}\\nQ {qsent:+.3f} → {qreceived:+.3f}'")
        code=code.replace("f'P injection = {row.P_injection_pu:+.3f} pu'","f'P {row.P_injection_pu:+.3f}  Q {row.Q_injection_pu:+.3f}'")
        code=code.replace("f'Load: {total.load_P_pu:.3f} pu    Line losses: {total.line_loss_P_pu:.3f} pu'","f'Load: {total.load_P_pu:.3f} pu    Line losses: {total.line_loss_P_pu:.3f} pu\\n'\n            f'Capacitor Q: {total.cap_Q_pu:.3f} pu'")
        code=code.replace("fig.suptitle('Bus-voltage magnitudes and terminal-voltage phases (Bus6 = 1 pu ∠ 0°)'",f"fig.suptitle({(label+'; Bus6 = 1 pu ∠ 0°')!r}")
        code=code.replace("Blue arrows: selected active-power paths; labels are sending-end P → received P [pu].","Arrows follow P; labels: signed sending-end → received P/Q [pu] along that direction.")
        code=code.replace("Bus P: net terminal injection (+ generation, − consumption). Post-fault panel is an equilibrium, not the clearing instant.","Bus P/Q: net injections (+ supply, − consumption). Post-fault panel: stable equilibrium, not clearing instant.")
        shunt="""    if CAP_B>0:
        ax.plot([4.95,5.45,5.45],[1.,1.,.53],color='#7b4c92',lw=1.4,zorder=2)
        ax.plot([5.17,5.73],[.53,.53],color='#7b4c92',lw=2,zorder=2)
        ax.plot([5.17,5.73],[.29,.29],color='#7b4c92',lw=2,zorder=2)
        ax.plot([5.45,5.45],[.29,-.20],color='#7b4c92',lw=1.4,zorder=2)
        for yy,ww in [(-.20,.24),(-.34,.16),(-.46,.08)]:
            ax.plot([5.45-ww,5.45+ww],[yy,yy],color='#7b4c92',lw=1.1,zorder=2)
        ax.text(5.95,.34,f'Shunt B={CAP_B:.3f} pu',color='#7b4c92',fontsize=8,va='center')
"""
        code=code.replace('    for b,(x,y) in pos.items():',shunt+'    for b,(x,y) in pos.items():')
        if B==1.:
            code=code.replace("OUT=FIG/'network_data'","OUT=FIG/'network_data_B100'")
            code=code.replace('equilibrium_network','equilibrium_network_B100')
        exec(compile(code,str(R/'code/plot_equilibrium_network.py'),'exec'),{'__file__':str(R/'code/plot_equilibrium_network.py'),'__name__':'case_network_renderer'})
    update_network_readme()
    print('PASS: both cases written to figures, network panels embedded in figures/README.md')
if __name__=='__main__':plot_case_networks()
