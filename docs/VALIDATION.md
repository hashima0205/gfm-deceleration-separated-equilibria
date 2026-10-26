# 共有用コードの再計算確認

2026-10-27、共有用フォルダの新規チェックポイントから、事故前平衡点・事故中・事故除去後まで計算した。元の事故除去状態・時系列を新規計算の初期値として流用していない。使用したモデルの共通snapshotは元研究31と全ファイルのSHA-256が同一である。

独立した新規リポジトリ用に構成を変更した後も、再び新規チェックポイントからvalidate_delivery.csvの4条件（smoke＋demoの3条件）を実行し、下表と全状態照合が一致した。以前のフォルダ構成での生成結果はローカルのvalidation_previous_layout/へ保全し、この確認には使用していない。

環境：Julia 1.12.5、リポジトリ共通のprovenance/Project.toml・Manifest.toml、Python環境はcode/requirements.txtに対応。OPENBLAS_NUM_THREADS=1、JULIA_NUM_THREADS=1、GKSwstype=100。Rodas5P、reltol=1e−9、abstol=1e−11、dtmax=.005秒。保存刻みはREADME記載のとおり。

実行コマンド（リポジトリ直下）：

```bash
julia --history-file=no --project=provenance code/run.jl smoke.csv
python code/check_reference.py --smoke
julia --history-file=no --project=provenance code/run.jl demo.csv
python code/check_reference.py --all
python code/summarize.py
python code/plot_results.py
```

| 全台τp [s] | 事故除去後観測 [s] | 新規内部角変化 [回転] | 分類 | 元値との差 [回転] |
|---:|---:|---:|---|---:|
| 1.20（smoke） | 30 | −0.258413397670 | 未判定 | 終端全状態で照合 |
| 1.16 | 300 | −0.216354594663 | 回復 | 0 |
| 1.18 | 300 | −41.666656611154 | 負方向回転継続 | 0 |
| 1.20 | 1000 | −164.399278475694 | 負方向回転継続 | 0 |

各ケースの事故前1.0秒、事故除去1.15秒、評価終端について、元の全66状態との最大絶対差は今回の環境では0だった。異なる環境でのビット単位一致を要求するものではなく、check_reference.pyは明示した絶対・相対許容差で比較する。状態成分の単位は混在するため、最大差を物理的な単一距離とは解釈しない。

全台・事故前/中/後のτp/τq、正のP指令、VC_F=CC1_F=1を確認。各条件で安定な事故後根と、不安定モード1個の鞍点の残差・固有値を再検証した。平均角間隔0.13944473586354644 radも保存全状態から再計算して一致した。集計・作図まで正常終了し、3条件の比較図と各機の角度・周波数・電圧・P応答図を確認した。

実行表と終端比較の小さな確認データをtests/reference/sharing_reproduction_summary.csv、sharing_endpoint_validation.csvに保存した。計算ログ・全時系列・生成図はローカルのlogs/、checkpoints/、figures/に保持し、Gitの対象からは除外している。

この共有確認では精度Bを再実行していない。元研究でτp=1.20秒・T=.15秒を精度Bで1000秒再計算した結果はoriginal_numerical_comparison.csvにある。共有確認は別ソルバ・EMT/PWM・実機の追加検証ではない。
