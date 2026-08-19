<!-- このファイルは scripts/export_docs.py が生成します。直接編集しないでください。 -->

# 使用した文献値

モデルが使う物理定数の一覧。すべて `adhesion/literature.py` で出典と一体で
定義されており、コード中に数値を直書きしている箇所はない。

内訳はCostan et al. (2024) から 15 件、Tirado et al. (1979) から 1 件。

| 記号 | 量 | 値 | 出典 | 備考 |
|---|---|---:|---|---|
| `l` | 細胞の全長 | 1.80 ± 0.25 µm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | スフェロシリンダー（円柱＋両端半球）として扱ったときの端から端まで。 |
| `r` | 細胞の半径 | 0.40 ± 0.10 µm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `2M` | アドヘシン対が架橋できる全長 | 8.00 nm | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 片方の細胞表面から、ナノボディと抗原の両方を経て、相手の細胞表面までの全長。2細胞をまたぐ総額なので、両細胞ぶんを足してはいけない。 |
| `F_b` | 単一アドヘシン結合の破断力 | 16.1 ± 0.4 pN | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 光ピンセットによる測定。荷重速度 5–100 pN/s の範囲で値は変わらなかった。 |
| `D` | 外膜内でのアドヘシンの拡散係数 | 0.36 ± 0.06 µm²/s | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `M_f` | 拡散可能なアドヘシンの割合 | 0.37 ± 0.04 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 残り63%は不動。外膜タンパクの島構造と整合的。 |
| `N_100` | アドヘシン数（100 ng/mL aTc 誘導時） | 7,300 ± 1,800 分子/細胞 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `N_300` | アドヘシン数（300 ng/mL aTc、飽和誘導時） | 15,300 ± 4,100 分子/細胞 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `σ_300` | 表面密度（飽和誘導時） | 3,400 ± 1,200 分子/µm² | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 外膜タンパク全種の総密度（10⁵ /µm² のオーダー）の数%にすぎない。 |
| `a_max` | 最大産生速度 | 9,900 ± 2,100 分子/細胞/h | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `b_eff` | アドヘシンの実効減衰速度 | 0.78 ± 0.12 1/h | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 大半は細胞分裂による希釈（0.82 ± 0.30 /h）。分解自体は 0.05 ± 0.01 /h と遅い。 |
| `n_∥` | 平行に並んだ2細胞間の結合数 | 130.0 ± 30.0 対 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 100 ng/mL 誘導時。モデルの検証対象。 |
| `n_⊥` | 直交した2細胞間の結合数 | 7.00 ± 2.00 対 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 100 ng/mL 誘導時。接触面積が最小になる配置。モデルの検証対象。 |
| `A_con,min` | 接触面積が細胞表面積に占める割合（下限） | 0.00 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) |  |
| `A_con,max` | 接触面積が細胞表面積に占める割合（上限） | 0.04 | [Costan et al. (2024)](https://doi.org/10.1016/j.matt.2024.03.019) | 配置に依存して 0.002–0.036 の範囲で変わると報告されている。 |
| `η` | 培地の動粘性（水で近似） | 1.00 mPa·s | [Tirado et al. (1979)](https://doi.org/10.1063/1.438613) | Tirado らの抵抗係数式と組で使うため、便宜上ここに置く。 |

## 値の性格について

上の表は2種類の値を含む。

- **モデルの入力**: 細胞寸法、架橋長、破断力、拡散係数、発現量など。
  これらはモデルに与える。
- **検証の対象**: 報告された結合数（`n_∥`, `n_⊥`）と接触面積割合（`A_con`）。
  これらはモデルに与えず、モデルが独立に再現できるかを確かめるために使う。
  検証結果は `tests/test_validation.py` に固定してある。

## 引用文献

- Costan, Sara A. and Kim, Honesty and Others (2024). Biophysical characterization of synthetic adhesins for predicting and tuning engineered living material properties. *Matter **7**, 2836-2854*. <https://doi.org/10.1016/j.matt.2024.03.019>
- Glass, David S. and Riedel-Kruse, Ingmar H. (2018). A Synthetic Bacterial Cell-Cell Adhesion Toolbox for Programming Multicellular Morphologies and Patterns. *Cell **174**, 649-658*. <https://doi.org/10.1016/j.cell.2018.06.041>
- Rassam, Patrice and Copeland, Nicholas A. and Others (2015). Supramolecular assemblies underpin turnover of outer membrane proteins in bacteria. *Nature **523**, 333-336*. <https://doi.org/10.1038/nature14461>
- Tirado, Maria M. and Garcia de la Torre, Jose (1979). Translational friction coefficients of rigid, symmetric top macromolecules. Application to circular cylinders. *The Journal of Chemical Physics **71**, 2581-2587*. <https://doi.org/10.1063/1.438613>

---
生成日: 2026-08-20
