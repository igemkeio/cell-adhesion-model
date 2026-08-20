# 細胞接着モデル

極局在させたアドヘシンで大腸菌どうしを接着させたとき、できあがる多細胞構造が
**一列に伸びる鎖になるのか、横に固まる束になるのか**を予測する粗視化モデル。

物理定数はすべて Costan et al., *Matter* (2024) の実測値に固定してあり、
調整パラメータを持たない。自由に動かせるのは**発現量**と**局在**の2つだけである。

## 主要な結果

**実測されている局在（一様分布）は、鎖ではなく束を予測する。**

| 局在（N = 7,300） | 横並び | 濃い極 | 薄い極 | 接合強度 | 構造 |
|---|---:|---:|---:|---:|---|
| 一様（実測値） | 100 | 8 | 8 | 1,607 pN | **束** |
| 極 80 % 対称二極 | 48 | 15 | 15 | 765 pN | 束 |
| 極 95 % 対称二極 | 25 | 18 | 18 | 409 pN | 束 |
| 極 95 % 二極 + 先端 0.20 µm | 9 | 75 | 75 | 1,209 pN | **鎖** |
| 極 95 % 二極 + 先端 0.10 µm | 8 | 265 | 265 | 4,262 pN | **鎖** |
| 極 95 % **単極** + 先端 0.10 µm | 8 | 529 | 0 | 8,525 pN | **二量体どまり** |

要点は3つある。

**1. 極に集めるだけでは足りず、極の先端に絞る必要がある。**
極への局在比を 0.95 まで上げても、半球全体に広がったままでは束のままである。
接触パッチの半径は 56.6 nm しかなく、極半球（測地距離で 0.63 µm）の大部分は
接触に寄与しない。効くのは先端集中幅 *w* で、パッチ内の分子数は *w*⁻² で増える。

| 先端集中幅 *w* | 鎖になるのに必要な極局在比 |
|---|---:|
| 集中なし（極半球に一様） | 事実上達成不能 |
| 0.30 µm | > 0.819 |
| 0.25 µm | > 0.762 |
| 0.20 µm | > 0.679 |
| 0.15 µm | > 0.554 |
| 0.10 µm | > 0.444（わずかな偏りで足りる） |

**2. 発現量では形は変わらない。**
極接触も側面接触も分子数に比例するため、比を取ると分子数が消える。
すなわち**発現量は接合の強さを決め、局在は形を決める。両者は干渉しない。**
プロモーターを強くしても束が鎖になることはない。

**3. 両極に集める必要がある。**
鎖が伸びるには1細胞が2つの接続を持たねばならず、その律速は分子が少ない側の極である。
単極局在では接合1つあたりは2倍強くなる代わりに、相手を1つしか持てず二量体で止まる。
外膜タンパクに自然に生じる偏りは古い極への片側性の蓄積（Rassam et al. 2015）であり、
まさにこの単極型にあたる。**放置して出てくる局在は鎖には使えない。**

![主要結果](results/figures/main_results.png)

## 用語の図解

本文と図に出てくる用語のうち、初見で分かりにくく、かつ結論を左右するものを図にまとめた。
用語の一覧は [docs/glossary.md](docs/glossary.md) にある。

![用語の図解](docs/figures/glossary_diagram.png)

とくに **A**（測地距離）、**B**（接触パッチの半径がどこから来るか）、
**D**（*f*pole = 0.444 が一様分布であること）は誤解しやすい。

## 検証

モデルには調整パラメータがない。したがって以下の一致は較正の結果ではなく予測である。

| 量 | このモデル | Costan et al. (2024) | 一致 |
|---|---:|---:|---|
| 横並び接触の結合数 | 100 対 | 130 ± 30 対 | 1 SD 以内 |
| 極接触の結合数 | 8 対 | 7 ± 2 対 | 1 SD 以内 |
| 接触面積割合 *A*con | 0.0023 – 0.0274 | 0.002 – 0.036 | レンジ一致 |

とくに *A*con が独立に一致した点が重要である。原著はこれを配置ごとの経験値として
与えているが、本モデルは細胞形状から幾何的に導出している。**同じ量が二つの独立な
経路から得られた**ことになる。

検証は `tests/test_validation.py` にテストとして固定してあり、将来モデルを変更して
文献との整合が崩れた場合はテストが失敗する。

## 使った文献値と置いた仮定

科学的な再現性のため、この2つは明示的に分離して管理している。

- **[docs/glossary.md](docs/glossary.md)** — 用語集。先端集中幅・接触パッチ・極局在比など、
  本文と図に出てくる用語の定義、記号、コード上の名前、代表値の対応表。**最初にここを読むとよい。**
- **[docs/parameters.md](docs/parameters.md)** — 使用した文献値の一覧。
  値・単位・不確かさ・出典が対応表になっている。
- **[docs/assumptions.md](docs/assumptions.md)** — モデルが置いている 10 個の仮定。
  各仮定について、そう置いた根拠と、**もし破れていたら結果がどう変わるか**を併記。

どちらも `adhesion/literature.py` と `adhesion/assumptions.py` から自動生成される
（`python3 scripts/export_docs.py`）。**コード中に出典のない数値を直書きしている
箇所はなく**、ドキュメントと実装が食い違うことはない。

```python
# adhesion/literature.py — 数値は必ず出典と一体で定義する
RUPTURE_FORCE = LiteratureValue(
    symbol="F_b",
    value=16.1,
    uncertainty=0.4,
    unit="pN",
    quantity="単一アドヘシン結合の破断力",
    source="costan2024",
    note="光ピンセットによる測定。荷重速度 5–100 pN/s の範囲で値は変わらなかった。",
)
```

仮定も同様に第一級のオブジェクトとして扱い、各モジュールが自分の依存する仮定を
宣言する。存在しない ID を参照すると読み込み時に失敗するため、仮定を変更すれば
参照側も必ず追随する。

```python
# adhesion/model.py
MODEL_ASSUMPTIONS = uses("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A10")
```

## パラメータ

| 分類 | パラメータ | 備考 |
|---|---|---|
| **自由 ①** | 発現量 *N*（0 – 15,300 分子/細胞） | aTc 誘導で実験的に制御。原著が誘導濃度との対応を実測 |
| **自由 ②** | 局在（`f_pole`, `tip_width_um`, `polar_split`） | 設計で変えたい量。実測の出発点は一様分布 |
| 固定 | 細胞寸法、架橋長、破断力、拡散係数ほか | すべて実測値（[docs/parameters.md](docs/parameters.md)） |
| 従属 | 結合数、接触面積、表面密度、接合強度 | 上記から計算される |

**結合の数は自由パラメータではない。** 幾何と発現量から決まる従属量であり、
1結合あたりの破断力は実測値に固定されるため、接着の強さの軸には自由パラメータが
1つも残らない。接合強度は単純な積になる。

```
接合強度 = 結合数 × 16.1 pN
```

## 震盪培養での挙動

散逸率 ε = 0.01 – 1 W/kg のとき、せん断速度は γ̇ = √(ε/ν) = 100 – 1000 s⁻¹。
このときコルモゴロフ長は 32 – 100 µm で細胞（1.8 µm）よりはるかに大きく、
細胞は乱流ではなく局所的な層流せん断の中にいる。

震盪の主な役割は**接合を壊すことではなく、細胞を出会わせること**である。
せん断による衝突頻度はブラウン運動の 27 – 267 倍に達する一方、
γ̇ = 300 s⁻¹ で 10 細胞の鎖にかかる張力は 12.3 pN にすぎず、
これを結合数で割ると単結合の破断力 16.1 pN には遠く届かない。

したがって律速は力学ではなく、**衝突したときに極どうしが向き合う確率**である。

## 構成

```
adhesion/
  literature.py    文献値。値・単位・不確かさ・出典を一体で保持
  assumptions.py   モデルが置く仮定。根拠と、破れた場合の影響つき
  geometry.py      接触幾何。表面の求積と、架橋しうる分子の数え上げ
  model.py         発現量と局在から接着と形態を予測
  shear.py         震盪培養での衝突頻度とせん断荷重
scripts/
  validate.py      文献値との照合表を出力
  run_scan.py      局在の掃引 → CSV
  make_figures.py  主要結果の図
  make_glossary_figure.py  用語の図解
  export_docs.py   docs/ をコードから生成
tests/             60 件。うち検証テストが文献との整合を固定する
docs/              自動生成されるパラメータ表と仮定一覧
references/        BibTeX
results/           図とデータ
```

## 実行方法

```bash
pip install -r requirements.txt

python3 -m unittest discover -s tests -t .   # 60 件
python3 scripts/validate.py                   # 文献値との照合
python3 scripts/run_scan.py                   # 掃引 → results/data/
python3 scripts/make_figures.py               # 図 → results/figures/
python3 scripts/make_glossary_figure.py       # 用語の図解 → docs/figures/
python3 scripts/export_docs.py                # docs/ を再生成
```

単発の予測は次のように書ける。

```python
from adhesion.model import predict, UNIFORM_F_POLE

predict(7300, UNIFORM_F_POLE).morphology                  # 束（実測の局在）
predict(7300, 0.95, tip_width_um=0.10).morphology         # 鎖
predict(7300, 0.95, 1.0, tip_width_um=0.10).morphology    # 二量体どまり
```

## モデルの適用範囲

このモデルは**接触した2細胞がどの配置を取るかという相対的な有利不利**を予測する。
以下は扱っていない（仮定 A7）。

- 衝突するまでの待ち時間、結合の寿命、解離速度定数
- 凝集にかかる時間、鎖長分布の時間発展
- 接触後の細胞の回転による再配置

速度論を入れるには解離速度定数と2次元親和性が必要だが、いずれも本系では
測定されていない。

## 残された問い

モデル側で決着していない自由度はほぼない。残るのは生物学の問いが1つである。

> **外膜提示されたアドヘシンを、対称二極に、極先端 0.2 µm 以内へ濃縮できるか。**

できれば鎖ができ、できなければ束になる。現時点で判明しているのは次の4点である。

- この構成での実測分布は**一様**である（Costan et al. 2024）
- 外膜に自然に生じる偏りは**単極**（古い極への蓄積）で、鎖には使えない（Rassam et al. 2015）
- PopZ・DivIVA といった既知の極局在因子は細胞質・内膜に作用するもので、外膜タンパクには届かない
- ただし拡散は十分速く（*D* = 0.36 µm²/s なら細胞長を約 4.5 秒で横断）、**トラップさえあれば局在は数秒で成立する**

最後の点が唯一の希望であり、問題は輸送ではなく外膜側に留め具が存在しないことにある。

## 引用文献

- Costan, S. A., Kim, H. et al. (2024). Biophysical characterization of synthetic
  adhesins for predicting and tuning engineered living material properties.
  *Matter* **7**, 2836–2854. <https://doi.org/10.1016/j.matt.2024.03.019>
- Glass, D. S. & Riedel-Kruse, I. H. (2018). A Synthetic Bacterial Cell-Cell Adhesion
  Toolbox for Programming Multicellular Morphologies and Patterns.
  *Cell* **174**, 649–658. <https://doi.org/10.1016/j.cell.2018.06.041>
- Rassam, P., Copeland, N. A. et al. (2015). Supramolecular assemblies underpin
  turnover of outer membrane proteins in bacteria.
  *Nature* **523**, 333–336. <https://doi.org/10.1038/nature14461>
- Tirado, M. M. & García de la Torre, J. (1979). Translational friction coefficients
  of rigid, symmetric top macromolecules. Application to circular cylinders.
  *J. Chem. Phys.* **71**, 2581–2587. <https://doi.org/10.1063/1.438613>

BibTeX は [references/references.bib](references/references.bib) にある。

## ライセンス

MIT License（[LICENSE](LICENSE)）
