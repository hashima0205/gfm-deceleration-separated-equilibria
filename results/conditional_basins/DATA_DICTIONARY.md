# 列名・単位・平均方法

遠方群はBus2・3・7、全台はBus1・2・3・7です。平均は等重みです。内部位相δは固定60Hz同期座標系における[rad]で、Bus1基準ではありません。端子電圧位相とは区別します。

## 座標と時刻

| 列・値 | 定義・単位 |
|---|---|
| theta_relative_rad | 遠方平均内部角−事故後安定平衡点の遠方平均内部角 [rad] |
| definition.json の origin_rad | 上式の基準角、−1.5548103219232832 rad |
| df_Hz | 初期状態として指定した遠方平均内部周波数偏差 [Hz] |
| df_mean_Hz | 各時刻の遠方平均内部周波数偏差 [Hz]。内部角速度ωから求める |
| t_absolute_s | 事故投入前からの絶対時刻 [s]。投入1.0、除去1.15 |
| t_post_s（実事故時系列） | 絶対時刻−1.15 [s] |
| t_post_s（隠れ状態比較） | 人為的な初期状態からの経過時間 [s] |
| remote_internal_turns | 終端と開始時の遠方平均内部角差/(2π) [回転] |

領域図はtheta_relative_rad、通常の時系列位相変化はその先頭値を引いて2πで割ったものです。隠れ状態図はtheta_relative_rad/(2π)。この3つの表示を混同しないでください。

## jobs.csv／classification.csv／plot_samples_tau*.csv

| 列 | 定義 |
|---|---|
| case | 元研究での計算ID |
| taup | 全台の有効電力測定フィルタ時定数 [s]。τqは全台.02秒 |
| section | equilibrium：46状態の初期値を平衡点に置く。clearing：対応する実事故除去状態に置く |
| kind | coarse、refine等の点を追加した目的。strictは積分精度比較 |
| horizon_s | 指定観測時間 [s] |
| completed_horizon_s | 保存結果が実際に到達した時間 [s] |
| level | A：reltol=1e−9,abstol=1e−11,dtmax=.005秒。B：1e−10,1e−12,.0025秒 |
| step_s | 通常の保存間隔 [s]。最初20msには.0001秒の追加保存 |
| worker | 元研究の実行分担番号。物理パラメータではない |
| classification | recovery、negative_running等。有限観測の分類 |
| max_algebraic_residual | 保存点での代数方程式の最大絶対残差。モデルのスケーリングのまま |
| max_remote_angle_spread_rad | 遠方各機の内部角の最大機間差 [rad] |
| max_remote_df_spread_Hz | 遠方各機の内部周波数偏差の最大機間差 [Hz] |
| *_state_error_inf_mixed_units | 安定根との最大成分誤差。成分は混合単位であり単一の物理距離ではない |
| cc_integral_output_error_max_pu | 電流制御PI積分状態の誤差に積分ゲインを掛けた出力差 [pu] |
| vc_integral_output_error_max_pu | 電圧制御PI積分状態の誤差に積分ゲインを掛けた出力差 [pu] |
| tail_voltage_min_pu,max_pu | 末尾評価区間の端子電圧範囲 [pu] |
| persistent_low_voltage | 元の判定では末尾区間中ずっと最大電圧が.2 pu未満の機がある場合True |

背景図はequilibriumかつlevel=Aの88点/τを使います。plot_triangles_tau*.csvのvertex0/1/2は、対応するplot_samples_tau*.csvの**0始まりの行番号**です。読み込み後の行順を変えないでください。

## initial_states.csv／references／states.csv

initial_states.csvはcaseとx1～x66の列を持ち、各計算の保存初期状態を示します。52微分状態と14代数変数を含みます。添字・接続バス・変数名・微分/代数の区別はstates.csvに対応します。すべてが同じ単位ではなく、各変数のモデル定義に従います。

代数変数は積分前に整合化し、52微分状態は保持します。definition.jsonのheld_differential_indicesは初期値をそろえる46微分状態の1始まりの添字です。遠方内部角3変数とPフィルタ3変数を一様にずらして指定する2座標を作ります。46状態は初期化後には動きます。

section_stable/saddle.csvは事故後の全状態、clearing_*.csvは事故除去時の全状態、parameters_*.csvは当該τでの事故後パラメータベクトルです。これらのベクトル添字は固定したJuliaモデルの順序に依存します。

## trajectories／hidden_trajectories

実事故時系列は元の保存点のうち次の観測列を選んだものです。追加の時間間引きはありません。制御器の全観測量や52状態の全時間履歴を同梱したものではありません。

| 各機の列 | 単位・測定点 |
|---|---|
| bN_delta_rad | 内部位相 [rad] |
| bN_omega_rad_s | 内部角速度 [rad/s]。60 Hzの基準を含む |
| bN_df_Hz | 内部角速度/(2π)−60 [Hz]。端子位相の微分ではない |
| bN_U_r_pu,U_i_pu | バス端子電圧の実部・虚部 [pu] |
| bN_V_pu | 上記複素端子電圧の絶対値 [pu] |
| bN_Pcap_pu | Re(Vc conj(Ig))、コンデンサ側の有効電力 [pu]。指令に対応 |
| bN_Pterminal_pu | Re(U conj(Ig))、grid側フィルタ後の端子有効電力 [pu] |

hidden_trajectoriesは同じ初期座標で分類が異なった3組6本の平均内部角・周波数だけを保存します。全て1000秒、図の表示は最初100秒です。対応はhidden_state_comparison.csvにあります。

## 診断表

machine_metrics.csvは各機の分類と内部/端子位相変化、末尾周波数・電圧等を保存します。末尾周波数平均は台形積分による時間平均。terminal_phase_reliableは元の保存点で電圧が.1 pu以上かによる信頼性判定で、保存点間の位相誤差まで除外するものではありません。

initialization_audit.csvは元研究での代数整合化・座標誤差・微分状態保持の監査です。numerical_comparison.csvは同じ保存刻みでの積分精度A/B比較。anchor_reproduction.csvは実事故軌道と領域計算の照合です。共通保存時刻での差と、時刻補間を含む差を区別してください。
