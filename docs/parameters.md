<!-- このファイルは scripts/export_docs.py が生成します。直接編集しないでください。 -->

# 使用した文献値

モデルが使う物理定数の一覧。すべて `adhesion/literature.py` で出典と一体で
定義されており、コード中に数値を直書きしている箇所はない。

内訳はCostan et al. (2024) から 18 件、Klamecka et al. (2015) から 2 件、Rassam et al. (2015) から 2 件。

| 記号 | 量 | 値 | 出典 | 備考 |
|---|---|---:|---|---|
| `l` | 細胞の全長 | 1.80 ± 0.25 µm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | スフェロシリンダー（円柱＋両端半球）として扱ったときの端から端まで。 |
| `r` | 細胞の半径 | 0.40 ± 0.10 µm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `2M` | アドヘシン対が架橋できる全長 | 8.00 nm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 片方の細胞表面から、ナノボディと抗原の両方を経て、相手の細胞表面までの全長。2細胞をまたぐ総額なので、両細胞ぶんを足してはいけない。 |
| `F_b` | 単一アドヘシン結合の破断力 | 16.1 ± 0.4 pN | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 光ピンセットによる測定（n = 14 細胞）。本文の値は 16.1 ± 0.4 pN、要約表では 16.2 と表記され系統誤差 1.8 pN が併記されている。本モデルは本文の値を採る。この値は anti-EPEA / EPEA ペアの測定である（EPEA は4アミノ酸のペプチドタグ）。原著の表が期待値として挙げる 28–45 pN は anti-GFP / GFP ペアの値なので、両者は別の分子ペアであり食い違いではない。原著自身も『パラメータはアドヘシンのペアごとに、とくに F_b は変わる』と註記している。すなわち破断力はペアの選択で変えられる設計変数である。 |
| `D` | 外膜内でのアドヘシンの拡散係数 | 0.36 ± 0.06 µm²/s | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | FRAP による測定。外膜タンパクの拡散係数はタンパクによって2桁にわたり広がる。Spector et al. (2010) によれば OmpF は 0.006 µm²/s で直径 100 nm の領域に閉じ込められ、BtuB は 0.05–0.10、TonB による繋ぎ止めを外した変異体では0.27 µm²/s に達する。Rassam et al. (2015) の BtuB は 0.013、Cir は 0.019。本値 0.36 はこの範囲の上端にあたるが外れ値ではない。ただし可動割合が0.37 にとどまるため、全分子が自由に動くわけではない。 |
| `M_f` | 拡散可能なアドヘシンの割合 | 0.37 ± 0.04 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 残り63%は不動。外膜タンパクの島構造と整合的。 |
| `f_p` | プラスミドコピー数の非対称に対する補正係数 | 0.50 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | ナノボディ側は中コピー数プラスミド（20–30 コピー）、抗原側は低コピー数プラスミド（10–12 コピー）に載っている。抗原の数が約半分で結合数の律速になるため、著者らは総分子数を 0.5 倍して結合ペア数を見積もっている。「幾何的に届く分子のうち結合する割合」ではないことに注意。 |
| `t_0.5` | 結合が形成される半減時間 | 0.10 s | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 測定は装置の時間分解能（1 s）に阻まれ、上限 t_0.5 < 1 s しか決まっていない。原著は「典型的な on-rate は 100 ms のオーダーと予想される」と述べており、ここではその予想値を採る。桁の見積もりとして扱うこと。 |
| `k_off` | 荷重なしでの自発的な解離速度 | 0.00 1/s | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 荷重速度に対する破断力の傾きが有意に 0 と異ならず、Costan et al. は「k_off はゼロとみなすのが最良の推定であり、アドヘシンは自発的には外れない」と結論している。Glass & Riedel-Kruse (2018) も自分たちが作った凝集体を「essentially static（本質的に静的）」と呼び、動的で可逆な接着を得るには競合阻害剤の分泌や intimin 特異的プロテアーゼを別途組み込む必要がある、と述べている。したがって接着は自発的には不可逆である。ただし力学的には可逆で、光ピンセットでは結合と解離を繰り返し起こせる。 |
| `N_100` | アドヘシン数（100 ng/mL aTc 誘導時） | 7,300 ± 1,800 分子/細胞 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | anti-GFP / GFP ペアでの測定。 |
| `N_300` | アドヘシン数（300 ng/mL aTc、飽和誘導時） | 15,300 ± 4,100 分子/細胞 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | anti-GFP / GFP ペアでの測定（n = 20 細胞）。 |
| `σ_300` | 表面密度（飽和誘導時） | 3,400 ± 1,200 分子/µm² | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 参考までに、外膜タンパク全種の総密度は 10⁵ /µm² のオーダーとされるが、この比較値は一次文献を確認できていない（UNVERIFIED_ESTIMATES 参照）。 |
| `a_max` | 最大産生速度 | 9,900 ± 2,100 分子/細胞/h | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `b_eff` | アドヘシンの実効減衰速度 | 0.78 ± 0.12 1/h | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 大半は細胞分裂による希釈（0.82 ± 0.30 /h）。分解自体は 0.05 ± 0.01 /h と遅い。 |
| `n_∥` | 平行に並んだ2細胞間の結合数 | 130.0 ± 30.0 対 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 100 ng/mL 誘導時。モデルの検証対象。 |
| `n_⊥` | 直交した2細胞間の結合数 | 7.00 ± 2.00 対 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 100 ng/mL 誘導時。接触面積が最小になる配置。モデルの検証対象。 |
| `A_con,min` | 接触面積が細胞表面積に占める割合（下限） | 0.00 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `A_con,max` | 接触面積が細胞表面積に占める割合（上限） | 0.04 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 配置に依存して 0.002–0.036 の範囲で変わると報告されている。 |
| `d_island` | 外膜タンパクが作る島の直径 | 0.50 µm | [Rassam et al. (2015)](https://doi.org/10.1038/nature14461) | 外膜タンパクは数百分子からなる島に凝集し、島の内部に拡散が閉じ込められる。極半球の測地距離は 0.628 µm なので、島ひとつが極をほぼ覆う大きさにあたる。これより細かく局在を絞れるかは自明でない（仮定 A11）。 |
| `D_native` | 内在性外膜タンパク（BtuB）の拡散係数 | 0.01 µm²/s | [Rassam et al. (2015)](https://doi.org/10.1038/nature14461) | 単一分子追跡による。Cir では 0.019 µm²/s。コリシン標識した BtuB・Cir は3 分たっても FRAP による蛍光回復を示さず、事実上不動であった。Costan et al. が合成アドヘシンで測った 0.36 µm²/s とは約 20 倍の開きがある。 |
| `Δx` | 結合解離のエネルギー障壁までの距離 | 1.59 nm | [Klamecka et al. (2015)](https://doi.org/10.1088/1478-3975/12/5/056009) | 単一末端でアンカーした野生型 GFP の値。eGFP では 1.36 nm、両末端アンカーでは 2.02 nm となり、引く向きによって約 50 % 変わる。細胞表面での引かれ方は AFM の配置と異なるため、この幅は無視できない。 |
| `F_AFM` | anti-GFP ナノボディと野生型 GFP の破断力（AFM、範囲の上端） | 45.0 pN | [Klamecka et al. (2015)](https://doi.org/10.1088/1478-3975/12/5/056009) | AFM による直接測定。野生型 GFP で 28–45 pN、eGFP で 41–56 pN。引く速度 300 nm/s–10 µm/s の範囲での最頻破断力。Costan et al. が期待値として挙げる 28–45 pN はこの野生型の値にあたる。本モデルは接着に anti-EPEA ペアの 16.1 pN を使っており、この値は比較の参照としてのみ保持する。 |

## 値の性格について

上の表は2種類の値を含む。

- **モデルの入力**: 細胞寸法、架橋長、破断力、拡散係数、発現量など。
  これらはモデルに与える。
- **検証の対象**: 報告された結合数（`n_∥`, `n_⊥`）と接触面積割合（`A_con`）。
  これらはモデルに与えず、モデルが独立に再現できるかを確かめるために使う。
  検証結果は `tests/test_validation.py` に固定してある。

## 引用文献

- Costan, Stefana A. and Ryan, Paul M. and Kim, Honesty and Wolgemuth, Charles W. and Riedel-Kruse, Ingmar H. (2024). Biophysical characterization of synthetic adhesins for predicting and tuning engineered living material properties. *Matter **7**, 2125-2143*. <https://doi.org/10.1016/j.matt.2024.03.019>
- Glass, David S. and Riedel-Kruse, Ingmar H. (2018). A Synthetic Bacterial Cell-Cell Adhesion Toolbox for Programming Multicellular Morphologies and Patterns. *Cell **174**, 649-658*. <https://doi.org/10.1016/j.cell.2018.06.041>
- Lim, Hoong Chuin and Bernhardt, Thomas G. (2019). A PopZ-linked apical recruitment assay for studying protein-protein interactions in the bacterial cell envelope. *Molecular Microbiology **112**, 1757-1768*. <https://doi.org/10.1111/mmi.14391>
- Lin, Da-Wei and Liu, Yang and Lee, Yue-Qi and Yang, Po-Jiun and Ho, Chia-Tse and Hong, Jui-Chung and Hsiao, Jye-Chian and Others (2021). Construction of intracellular asymmetry and asymmetric division in Escherichia coli. *Nature Communications **12**, 888*. <https://doi.org/10.1038/s41467-021-21135-1>
- Chou, Tsung-Wei and Huang, Ing-Jen and Hsu, Ya-Chu and Wu, Szu-Hua and Hong, Jui-Chung and Ko, Yu-Hao and Su, Ting-Ran and Jiang, Jie-Hong Roland and Huang, Hsiao-Chun (2025). Construction and phenotypic classification of synthetic dual-pole Escherichia coli cells. *Communications Biology **8***. <https://doi.org/10.1038/s42003-025-08495-w>
- Klamecka, Kamila and Severin, Philip M. and Milles, Lukas F. and Gaub, Hermann E. and Leonhardt, Heinrich (2015). Energy profile of nanobody-GFP complex under force. *Physical Biology **12**, 056009*. <https://doi.org/10.1088/1478-3975/12/5/056009>
- Rassam, Patrice and Copeland, Nicholas A. and Others (2015). Supramolecular assemblies underpin turnover of outer membrane proteins in bacteria. *Nature **523**, 333-336*. <https://doi.org/10.1038/nature14461>

---
生成日: 2026-08-31
