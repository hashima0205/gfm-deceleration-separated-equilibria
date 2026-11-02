# 元の設定

B=0、事故前Bus5電圧0.90364 pu、負荷G=2.2682301152 pu。GFMのP指令[0.5,0.4,0.4,0.4] pu、τq=0.02秒。Bus4の20 puシャント事故0.15秒後、除去と同時に線路4–5開放。

- [01_equilibria_section_tau010](01_equilibria_section_tau010.png)（[PDF](01_equilibria_section_tau010.pdf)）：安定点と鞍点の両方を通る断面＋事故中・事故後軌道の投影。τp=0.1秒。断面外の軌道を投影している。
- [03_trajectories_tau010_tau120](03_trajectories_tau010_tau120.png)（[PDF](03_trajectories_tau010_tau120.pdf)）：τp=0.1秒は回復、1.2秒は減速方向の持続回転。破線は事故中、丸は除去時。

各図はPNG/PDF。位相は無限大母線Bus6基準、折り返しなし。領域は有限時間の格子分類で、塗りの全点を保証しない。


## 大きいτの検討と計算状況

- [04_basins_large_tau](04_basins_large_tau.png)（[PDF](04_basins_large_tau.pdf)）：τp=1.16・1.18・1.20秒の引き込み領域と事故中・事故後軌道。各88点、計264点。
- [05_basins_large_tau_zoom](05_basins_large_tau_zoom.png)（[PDF](05_basins_large_tau_zoom.pdf)）：平衡点近傍と事故除去点近傍の拡大。
- [06_fault_on_large_tau](06_fault_on_large_tau.png)（[PDF](06_fault_on_large_tau.pdf)）：事故中0.15秒の軌道と周波数を単独表示。

| τp [s] | 平衡点基準の条件付き引き込み領域 | 事故中・事故後軌道 |
|---|---|---|
| 0.1 | 計算済み：両平衡点を通る35点 | 計算済み |
| 1.16・1.18・1.20 | 計算済み：各88点 | 計算済み |

大きいτの断面は残り46状態を安定点に固定した断面。鞍点は投影であり、この断面上の点ではない。両平衡点を通るτp=0.1秒の断面とは定義が異なる。
