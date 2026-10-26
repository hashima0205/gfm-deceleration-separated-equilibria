# 保存済みの引き込み領域・時系列・各機応答

研究32の保存結果を整理したフォルダです。**全台τp=1.16／1.18／1.20秒、全台τq=.02秒、事故時間T=.15秒、k=1.02990**で、リポジトリの主ケースと同じ条件です。図を見る・再作図するために、新しいシミュレーションは必要ありません。

## 図を選ぶ

| 確認したいこと | PDF | PNG |
|---|---|---|
| 3つの時定数で背景と実事故軌道を比較 | [全体比較](figures/conditional_basins_comparison.pdf) | [PNG](figures/conditional_basins_comparison.png) |
| 平衡点・鞍点と事故除去点の近傍 | [拡大比較](figures/conditional_basins_zoom.pdf) | [PNG](figures/conditional_basins_zoom.png) |
| τp=1.16秒の領域 | [単独図](figures/conditional_basin_tau1.16.pdf) | [PNG](figures/conditional_basin_tau1.16.png) |
| τp=1.18秒の領域 | [単独図](figures/conditional_basin_tau1.18.pdf) | [PNG](figures/conditional_basin_tau1.18.png) |
| τp=1.20秒の領域 | [単独図](figures/conditional_basin_tau1.20.pdf) | [PNG](figures/conditional_basin_tau1.20.png) |
| 実事故後の平均位相・周波数 | [時系列](figures/saved_time_responses.pdf) | [PNG](figures/saved_time_responses.png) |
| 各GFMの位相・周波数・電圧・P（1.16秒） | [各機応答](figures/each_gfm_tau1.16.pdf) | [PNG](figures/each_gfm_tau1.16.png) |
| 各GFMの位相・周波数・電圧・P（1.18秒） | [各機応答](figures/each_gfm_tau1.18.pdf) | [PNG](figures/each_gfm_tau1.18.png) |
| 各GFMの位相・周波数・電圧・P（1.20秒） | [各機応答](figures/each_gfm_tau1.20.pdf) | [PNG](figures/each_gfm_tau1.20.png) |
| 同じ角度・周波数でも残りの状態で結果が変わるか | [隠れ状態比較](figures/hidden_state_sensitivity.pdf) | [PNG](figures/hidden_state_sensitivity.png) |

領域図の横軸は遠方Bus2・3・7の平均内部角θ[rad]（事故後安定平衡点を0）、縦軸は平均内部周波数偏差Δf[Hz]（60Hzを0）。星は安定点、ひし形は鞍点の投影、白丸は実事故除去点の投影です。青が有限時間回復、橙が負方向回転、斜線が分類混在、破線が補間の目安です。

時系列図の横軸は事故除去後時間[s]、位相の縦軸は除去時からの変化[回転]です。隠れ状態比較では横軸は人為的な初期化からの時間、位相は安定平衡点基準の[回転]です。位相の原点と単位の違いに注意してください。

## 何点を計算した結果か

保存結果は281条件です。背景として描く初期状態は各τ88点、計264点。残りは実事故除去状態3点、隠れ状態比較12点、厳しい積分精度の比較2点です。

| τp [s] | 背景の計算点 | 回復 | 負方向回転 |
|---:|---:|---:|---:|
| 1.16 | 88 | 43 | 45 |
| 1.18 | 88 | 38 | 50 |
| 1.20 | 88 | 32 | 56 |

背景は通常300秒、追加確認点は1000秒で判定されています。各点の評価時間はdata/jobs.csvとclassification.csvを参照してください。この割合を全状態空間での安定確率と解釈することはできません。

初期値では、遠方3台の内部角・Pフィルタ状態を一様にずらし、残り46微分状態を安定平衡点の値にそろえます。14代数変数を整合化した後、全52微分状態を積分します。46状態を積分中も固定しているわけではありません。実事故の状態はこの断面上にあるとは限らず、同じ2座標で残りの状態を変更した比較も含めています。

## データを探す

- [data/classification.csv](data/classification.csv)：全281条件の分類・回転数・観測時間。
- [data/jobs.csv](data/jobs.csv)：各計算点のτ、座標、初期値断面、精度、保存刻み。
- [data/initial_states.csv](data/initial_states.csv)：各条件の初期66変数。添字は[data/states.csv](data/states.csv)に対応。
- [data/machine_metrics.csv](data/machine_metrics.csv)：平均だけでなく各GFM別の分類・周波数・位相・電圧・電流など。
- [data/outcomes.jsonl](data/outcomes.jsonl)：各条件の保存済み結果JSONを1行1条件にまとめたもの。
- data/plot_samples_tau*.csv と plot_triangles_tau*.csv：背景の計算点と、描画に使う三角形。
- data/trajectories/：実事故後の選択した観測量。保存時刻は間引いていません。
- data/hidden_trajectories/：異なる結果となった3組・計6本の1000秒平均位相／周波数時系列。
- data/references/：事故後安定点、鞍点、各τのパラメータ、事故除去時の全状態。
- [data/numerical_comparison.csv](data/numerical_comparison.csv)：元研究での積分精度A/B比較。
- [provenance/recalculation_check.json](provenance/recalculation_check.json)：今回の共有用入口で再計算した代表2点の比較。

詳細な列名・単位は[DATA_DICTIONARY.md](DATA_DICTIONARY.md)、初期状態の作り方は[計算方法](../../docs/SAVED_BASINS.md)を参照してください。全281点の生の高次元時系列は同梱していません。全初期状態とパラメータから再計算できます。

## 再作図

リポジトリ直下で次を実行します。

```bash
python code/check_saved_basins.py
python code/plot_saved_basins.py
```

出力先はfigures/saved_results/conditional_basins/。保存済みの三角形を使用し、分類点を増やしたり境界を再探索したりしません。初期値と分類表を照合し、時系列は既存保存点を使います。

## 出典と確認範囲

元データは研究32で計算済みです。繰り返し記録された時間延長ジョブは最後の設定にまとめ、元の284行のマニフェストもprovenance/original_jobs.csvとして保持しています。出典ファイルのSHA-256はprovenance/sources.json、公開する保存データのSHA-256はSHA256SUMSにあります。

公開準備時は全データの整合性と再作図を確認し、代表の回復点・負回転点のみ追加で再計算しました。全281点を再積分したという意味ではありません。provenance/original_*.mdは元研究の記録で、その中の実行手順は元の完全版フォルダ用です。現在の実行手順はこのリポジトリのREADME・docs/SAVED_BASINS.mdを使ってください。
