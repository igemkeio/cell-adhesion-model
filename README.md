Claude codeを使ってモデリングして検証しています．全文目を通していますが，おかしなところがあれば大森まで，
termの使い方など所々不自然なところありますが，用語集なども作りましたのでぜひそちらを参照ください


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

さらに原著の Supplemental Note には *A*con の導出が式の形で載っているため、
出力だけでなく**式そのものを突き合わせられる**。端 - 端接触について、

| | 式 | 値 |
|---|---|---:|
| 原著 (S32) | *A*ee = 2π*rM* | 0.01005 µm² |
| 本モデル | patch = π(2*M*)*r* | 0.01005 µm² |

この2つは**恒等的に同じ式**である。本モデルは原著の式を参照せず細胞形状から
独立に導出しているので、**同じ量が二つの独立な経路から得られた**ことになる。

原著が 2D 投影の近似と自ら断っている配置（直交・斜め）では、本モデルは曲面上で
直接数え上げるぶん小さい値を返す。詳細は
**[docs/source_comparison.md](docs/source_comparison.md)** にまとめた。
この差は結論には影響しない（極接触と側面接触の比は本モデルで 11.3、原著の値では 16 で、
どちらでも「局在させなければ束になる」は変わらない）。

検証は `tests/test_validation.py` にテストとして固定してあり、将来モデルを変更して
文献との整合が崩れた場合はテストが失敗する。

## 使った文献値と置いた仮定

科学的な再現性のため、この2つは明示的に分離して管理している。

- **[docs/project_design.md](docs/project_design.md)** — 本モデルが対象にしている設計（「One Direction」）と、
  接着機構（ペリプラズムでの CoilR–CoilE 係留）の説明。**最初にここを読むとよい。**
- **[docs/model_scope_audit.md](docs/model_scope_audit.md)** — モデルが設計資料に沿っているかの監査。
  一致している点と、ずれている点を列挙している。
- **[docs/localization_feasibility.md](docs/localization_feasibility.md)** — 極局在は実現できるかの文献調査。
- **[docs/glossary.md](docs/glossary.md)** — 用語集。先端集中幅・接触パッチ・極局在比など、
  本文と図に出てくる用語の定義、記号、コード上の名前、代表値の対応表。**最初にここを読むとよい。**
- **[docs/source_comparison.md](docs/source_comparison.md)** — 原著の導出式との照合。
  どこが一致し、どこがなぜ違うかを式のレベルで記録している。
- **[docs/parameters.md](docs/parameters.md)** — 使用した文献値の一覧。
  値・単位・不確かさ・出典が対応表になっている。
- **[docs/assumptions.md](docs/assumptions.md)** — モデルが置いている 13 個の仮定。
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
    note="光ピンセットによる測定（n = 14 細胞）。要約表では 16.2、系統誤差 1.8 pN。",
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

なお破断力 16.1 pN は **anti-EPEA / EPEA** ペアの測定値である。原著は3つの
Nb-Ag ペアを使い分けており、分子数・拡散・分布は anti-GFP / GFP、破断力は
anti-EPEA / EPEA で測っている。原著自身が「パラメータはペアごとに、とくに
*F*b は変わる」と註記しているとおり、**破断力はペアの選択で変えられる設計変数**
である。参考までに anti-GFP / GFP では AFM で 28–45 pN が報告されている
（Klamecka et al. 2015）。

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
  kinetics.py      荷重下での結合寿命（別ペアの値による桁の見積もり）
scripts/
  validate.py      文献値との照合表を出力
  run_scan.py      局在の掃引 → CSV
  make_figures.py  主要結果の図
  make_glossary_figure.py  用語の図解
  export_docs.py   docs/ をコードから生成
tests/             80 件。うち検証テストが文献との整合を式レベルで固定する
docs/              自動生成されるパラメータ表と仮定一覧
references/        BibTeX
results/           図とデータ
```

## 実行方法

```bash
pip install -r requirements.txt

python3 -m unittest discover -s tests -t .   # 80 件
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

- 衝突するまでの待ち時間、凝集にかかる時間、鎖長分布の時間発展
- 接触後の細胞の回転による再配置
- 接触中の拡散によるアドヘシンの集積（仮定 A12）

本系のアドヘシンペアについて解離速度定数も2次元親和性も測られていないため、
速度論を入れても値を埋められない。

ただし別のペア（anti-GFP / GFP）については解離速度定数
*k*off = 1.45×10⁻⁴ s⁻¹ と障壁位置 Δ*x* = 1.59 nm が測られており、
`adhesion.kinetics` で桁の目安を出せる。それによれば

- 無荷重での結合寿命は **約 1.9 時間**
- 荷重が効き始めるスケールは *kT*/Δ*x* ≈ **2.6 pN**
- γ̇ = 300 s⁻¹ で 10 細胞の鎖にかかる 1 結合あたりの荷重は 0.046 pN であり、
  寿命はほとんど変わらない（無荷重の 98 %）

**時間を入れても力学が律速にならないという結論は変わらない見込みである。**

## 設計との対応と、残された問い

本モデルが対象にしているのは「One Direction」（iGEM Keio-Japan 2026）の設計である。
詳細は [docs/project_design.md](docs/project_design.md)。要点だけ記すと、

- 外膜の Intimin は、**ペリプラズムでの CoilR–CoilE 結合**によって内膜側の
  極マーカー（PopZ）へ係留される。PopZ も DivIVA も細胞質・内膜側の因子で
  そのままでは外膜タンパクに届かないが、この係留がその隔たりを橋渡しする
- 2つの極に **Ag と Nb を分けて**提示するので、Ag 極は Nb 極としか結合できない。
  したがって鎖は**頭-尾に揃い、向きを持つ**
- 本モデルはこれを `distinct_adhesins_per_pole=True`（既定）として扱い、
  形態を `DIRECTED_CHAIN`（鎖・向きあり）として出力する

なお **Ag/Nb を極ごとに分けても結合数は変わらない**。等量発現なら接合あたりの
結合数は `min(n_Ag, n_Nb) = N·f_pole/2` で、同種を両極に置いた場合と一致する。
変わるのはどの極とどの極が結合できるかであり、そこから向きが生まれる。

### 本モデルが答えたこと

- 一様分布では束になる（極-極 8 対に対し側面 100 対）
- 向きのある鎖には**両極**が必要（片極だけでは二量体で止まる）
- 極に集めるだけでは足りず、**先端に絞る**必要がある
- 発現量では形は変わらない。発現量は強さを、局在は形を決める
- 震盪培養のせん断は接合を壊さない

### 残された問い

**モデルの外にあるもの**

- 係留機構が実際にどれだけ先端に絞れるか（= `tip_width_um` の実測値）
- PopZ 極と DivIVA 極に、それぞれ何分子の Intimin を係留できるか
- 外膜タンパクの島構造（直径 0.5 µm）が係留にどう影響するか（仮定 A11）。
  ただし係留は拡散捕捉ではなく直接の繋ぎ止めなので、島の制約を回避できる可能性がある

**モデル自体に足りていないもの**（[docs/model_scope_audit.md](docs/model_scope_audit.md)）

- **ブラウン運動**。資料は「多価結合とブラウン運動を考慮したシミュレーション」を
  指定しているが、本モデルは仮定 A7 で時間を持たない。方向性接着の様子を
  可視化するにはこれが要る
- **架橋長**。2*M* = 8 nm は Costan et al. の構築の値であり、この設計の Intimin は
  LysM を CoilE に置換し Spacer を持つ。到達距離は変わるはずである
- **誘導系**。本モデルは aTc 濃度に紐づいているが、資料が使うのは IPTG と
  プロモーター強度である

## 文献の扱い方

**値は全文を確認した文献からのみ取る。**検索結果の要約や、他の論文が引用している
孫引きの数値をモデルに入れてはならない。この規則は
`tests/test_literature.py::test_values_only_come_from_fully_verified_sources`
で強制してあり、未確認の文献から値を取るとテストが失敗する。

| 文献 | 全文確認 | 供給している値 |
|---|---|---:|
| Costan et al. (2024) ＋ Supplement | ✅ | 18 |
| Rassam et al. (2015) | ✅ | 2 |
| Klamecka et al. (2015) | ✅ | 2 |
| Glass & Riedel-Kruse (2018) | ⚠️ 未入手 | 0（背景の言及のみ） |

Klamecka については照合の結果、**当初 3 値を取っていたうち 2 値が誤りだった**。
解離速度定数 *k*off と解離定数 *K*D は Klamecka 自身の測定ではなく、
Kirchhofer et al. (2010) および Kubala et al. (2010) からの引用値を
Bell–Evans 回帰に固定して使っているだけだった。両論文とも未入手のため、
この 2 値はモデルから削除した（`adhesion/kinetics.py` では呼び出し側が
明示的に渡す形にしてある）。

### 出典を確認できていない見積もり

- **外膜タンパクの総面密度（10⁵ /µm² のオーダー）** — 極先端への濃縮が物理的に
  可能かを論じる比較対象に使っているが、一次文献を確認していない。
  仮定 A11 の議論の強さはこの値に依存する。`literature.UNVERIFIED_ESTIMATES` に
  明示してある

### 全文未確認のまま書いた文書

- **[docs/localization_feasibility.md](docs/localization_feasibility.md)** —
  検索結果の要約のみに基づく。**モデルには取り込んでいない。**
  IcsA の極局在、アミダーゼ欠損による鎖形成、intimin の一様分布などは
  すべてこの状態にある

## 読むべき文献

現時点で入手・全文確認が必要なもの。優先度順。

| 優先 | 文献 | 何のために |
|---|---|---|
| 高 | **Lim & Bernhardt (2019) Mol. Microbiol.**（POLAR 系） | 本設計の中核であるペリプラズム CoilR–CoilE 係留の原典。係留がどれだけ極に絞れるか（`tip_width_um`）の唯一の手掛かり |
| 高 | **Lin et al. (2021) Nat. Commun.**（PopZ） | 発現量で単極/両極が切り替わる条件。実験①の設計根拠であり、`polar_split` の実測値に相当する |
| 高 | **Glass & Riedel-Kruse (2018) Cell** | アドヘシン・ツールボックスの原典。鎖ができた前例があるか、intimin 構築の詳細（架橋長 8 nm の内訳） |
| 中 | **Rassam の島構造を裏づける外膜密度の一次文献** | 上記 UNVERIFIED_ESTIMATES を解消するため |
| 中 | **Kirchhofer et al. (2010) Nat. Struct. Mol. Biol.** | 削除した *k*off と *K*D の原典。結合寿命の議論を復活させたい場合に必要 |
| 低 | **Leo et al. (2016) JBC** / **Oberhettinger et al. (2012) PLoS ONE** | LysM 除去と CoilE 付与の妥当性（実験⑤の根拠） |
| 低 | Salema et al. (2013)、Spector et al. (2010) | Costan が期待値の比較に使っている文献。モデルは変わらない |

## 引用文献

- Costan, S. A., Ryan, P. M., Kim, H., Wolgemuth, C. W. & Riedel-Kruse, I. H. (2024). Biophysical characterization of synthetic
  adhesins for predicting and tuning engineered living material properties.
  *Matter* **7**, 2125–2143. <https://doi.org/10.1016/j.matt.2024.03.019>
- Glass, D. S. & Riedel-Kruse, I. H. (2018). A Synthetic Bacterial Cell-Cell Adhesion
  Toolbox for Programming Multicellular Morphologies and Patterns.
  *Cell* **174**, 649–658. <https://doi.org/10.1016/j.cell.2018.06.041>
- Klamecka, K., Severin, P. M., Milles, L. F., Gaub, H. E. & Leonhardt, H. (2015).
  Energy profile of nanobody-GFP complex under force.
  *Phys. Biol.* **12**, 056009. <https://doi.org/10.1088/1478-3975/12/5/056009>
- Rassam, P., Copeland, N. A. et al. (2015). Supramolecular assemblies underpin
  turnover of outer membrane proteins in bacteria.
  *Nature* **523**, 333–336. <https://doi.org/10.1038/nature14461>
- Tirado, M. M. & García de la Torre, J. (1979). Translational friction coefficients
  of rigid, symmetric top macromolecules. Application to circular cylinders.
  *J. Chem. Phys.* **71**, 2581–2587. <https://doi.org/10.1063/1.438613>

BibTeX は [references/references.bib](references/references.bib) にある。

## ライセンス

MIT License（[LICENSE](LICENSE)）
