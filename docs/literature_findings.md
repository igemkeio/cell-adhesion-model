# 全文照合で分かったこと

一次文献8本の全文を読んで照合した結果。設計に直接効くものを先に置く。

照合済み: Costan 2024（本文＋Supplement）／ Rassam 2015 ／ Klamecka 2015 ／
Lim &amp; Bernhardt 2019 ／ Lin 2021 ／ Glass &amp; Riedel-Kruse 2018 ／
Chou 2025 ／ Kirchhofer 2010 ／ Spector 2010 ／ Leo 2015 ／ Salema 2013

---

## 1. 実験③（第2極づくり）は、すでに解かれている

**Chou et al. (2025) *Commun. Biol.*** — [doi:10.1038/s42003-025-08495-w](https://doi.org/10.1038/s42003-025-08495-w)

> By engineering a **PodJ fragment as the second scaffold** and integrating
> **SpmX, a negative regulator of PodJ condensation, as a direct chimera with PopZ**,
> we effectively confined **PopZ and PodJ to opposite poles** within *E. coli*.

**第2極を作るという課題は、この論文が解決している。**しかも Lin et al. (2021) と
著者が重なる（Hong, J.-C.）ので、資料が引いている Lin 2021 の続報にあたる。

資料の設計との違い：

| | 資料の計画（実験③） | Chou 2025 |
|---|---|---|
| 第2の極マーカー | DivIVA | **PodJ 断片** |
| 排除機構 | PopZ に分解機能を付ける等、**これから開発** | **SpmX を PopZ に直接キメラ化**（PodJ 凝集の負の制御因子） |
| 結果 | 未達 | **PopZ と PodJ が反対の極に確定** |

資料 p10 では **SpmXΔC** を「PopZ 結合アダプター」として使う計画になっている。
つまりチームは SpmX を既に知っているが、Chou 2025 は同じ分子を**別の役割**
（PodJ の排除因子）で使って第2極を作っている。**実験③の設計を見直す価値がある。**

## 2. DivIVA を第2極に使う計画には、文献上の障害がある

**Lin et al. (2021)** — [doi:10.1038/s41467-021-21135-1](https://doi.org/10.1038/s41467-021-21135-1)

> Both mRFP-PopZ and DivIVA-sfGFP maintained **unipolar and bipolar localization,
> respectively**, when the two were co-expressed ... **PopZ does not have a
> significant effect on prohibiting DivIVA for membrane targeting**

- **PopZ は単極、DivIVA はもともと両極**である
- **PopZ は DivIVA を排除しない**

資料 p14 の註記（「Lin 2021 では PopZ は DivIVA を排除しないため、排除機構の付与が
要点」）は正しい。ただし **DivIVA が本来両極である**以上、排除機構を付けても
DivIVA は PopZ 極にも残りうる。Chou 2025 が DivIVA ではなく PodJ を選んだのは
この問題を回避するためと思われる。

あわせて確認できた実験①の根拠：

- 中程度の Anderson プロモーター **J23116** が極性の指標（極の鋭さ）で最高スコア
- **強いプロモーターにすると単極から両極へ切り替わる**

## 3. CoilE/CoilR は成立する。ただし試されたのはβバレルではない

**Lim &amp; Bernhardt (2019)** — [doi:10.1111/mmi.14391](https://doi.org/10.1111/mmi.14391)

POLAR は **PopZ-Linked Apical Recruitment** の略。本設計の中核をなす論文である。

**成立が確認できたこと**

> H3H4-GFP–TM-CoilE was found to **specifically recruit periplasmic CoilR-mScar
> to the cell poles**

内膜に固定した融合体がペリプラズムの相手を極へ引き寄せる、という機構そのものは
実証済みである。CoilE/CoilR は Moll et al. (2001) と Doh et al. (2018) の設計ペプチド。

**ただし3つの但し書きがある**

**(a) 向きが資料と逆である。**
Lim &amp; Bernhardt では **CoilE が内膜側**、**CoilR がペリプラズムの相手**。
資料 p9 は「内膜に固定した **CoilR** と、外膜 Intimin から伸びる **CoilE**」で、
**逆になっている**。ヘテロ二量体なのでどちらでも働く可能性はあるが、
原典と向きが違うことは意識しておくべきである。

**(b) 外膜側で試されたのはリポタンパクであって、βバレル外膜タンパクではない。**
外膜との相互作用を試した実験（Fig. 5）で使われたのは **LpoP**、すなわち外膜
**リポタンパク**である。本文に **"β-barrel" という語は一度も現れない**。

Intimin はβバレル外膜タンパクで、外膜に埋め込まれ、島構造に閉じ込められている
（Rassam 2015）。**リポタンパクを極へ引き寄せられることは、βバレルを引き寄せ
られることを意味しない。**これは仮定 A11 とは別の、新しい未検証点である。

**(c) 極リクルートが害になる場合がある。**
活性のあるペプチドグリカン合成酵素を極へ引き寄せると細胞に有害で、
外膜のブレブ形成が観察されている。Intimin では同じ問題は起きないと思われるが、
「極へ何かを引き寄せると細胞が困ることがある」という前例として記憶しておく価値がある。

**そして `tip_width_um` は、この論文からは決まらない。**
局在の評価はデモグラフ（細胞長に沿った蛍光の集団プロット）で行われており、
**極の広がりを nm で測った値はない**。モデルの最重要パラメータは依然として未定である。

## 4. 本プロジェクトは、原典が提案した将来方向そのものである

**Glass &amp; Riedel-Kruse (2018)** — [doi:10.1016/j.cell.2018.06.041](https://doi.org/10.1016/j.cell.2018.06.041)

> adhesins distributed **isotropically** over the outer membrane, but sub-cellular
> localization ... would allow spatial symmetry breaking and the production of
> **linear chains or sheets**

**アドヘシン・ツールボックスの原典が、「局在させれば線状の鎖が作れる」と
将来の方向として明示している。**本プロジェクトはそれを実行するものにあたる。

これで**アドヘシンが一様分布であることは3つの独立な経路から確認された**
（Costan の実測 p = 0.56、EPEC の免疫標識、この記述）。仮定 A10 は堅い。

構築の詳細も確認できた。

- 由来は **EHEC O157:H7 の intimin N末端**（reverse autotransporter）
- 構成は シグナルペプチド ＋ **LysM**（ペプチドグリカン結合）＋ **βバレル**
- その C 末端に接着ドメイン（ナノボディまたは抗原）を融合
- 誘導は **ATc（pTet/TetR）** または **アラビノース（pBAD/AraC）**

資料の実験⑤（LysM を除いて CoilE を付ける）は、この LysM をペリプラズム側の
CoilE に置き換えるということであり、構築上は理にかなっている。

### せん断についての、本モデルと食い違う記述

> Anecdotally, **mechanical agitation reversibly disrupts** the essentially static
> aggregates engineered in this work

**機械的な撹拌が凝集を可逆的に壊す**、と報告されている。本モデルは
「震盪のせん断は接合を壊さない（1結合あたり荷重が破断力の 1/56）」と結論しており、
一見矛盾する。

矛盾しない読み方はある。凝集塊が崩れることと、個々の結合が切れることは別である。
接触あたりの結合数が少なければ、結合が切れなくても細胞は引き離されうる。
また Costan は k_off をゼロとみなしており、**結合そのものは自発的には外れない**。

とはいえ **anecdotally（逸話的に）と断られた観察であり、定量データではない。**
本モデルの結論を覆すものではないが、**震盪条件を振る実験で確かめるべき点**として
記録しておく。

## 5. 削除した k_off と K_D は、まだ復活できない

**Kirchhofer et al. (2010)** — [doi:10.1038/nsmb.1727](https://doi.org/10.1038/nsmb.1727)

全文7ページを確認したが、**k_off = 1.45×10⁻⁴ /s も K_D = 1.4 nM も本文にない**。
本文は「binding constants (Kd) ... (**Supplementary Fig. 2b and Supplementary
Table 1**)」と述べており、値は補足資料にある。

**補足資料を入手するまで、この2値は使えない。**

## モデルへの反映

| 発見 | モデルへの影響 |
|---|---|
| アドヘシンの一様分布が3経路で確認 | 仮定 A10 が強化された。変更なし |
| CoilE/CoilR の極リクルートは実証済み | 設計の前提が支持された。変更なし |
| **βバレル外膜タンパクでは未実証** | **新しい未検証点。A11 とは別** |
| `tip_width_um` は文献から決まらない | 最重要パラメータは依然として未定 |
| DivIVA は本来両極 | 実験③の難度が上がる |
| PodJ + SpmX で第2極は達成済み | 実験③の設計を見直す材料 |
| 撹拌が凝集を壊すという逸話的観察 | せん断の結論に対する反証候補。要実験確認 |

---
照合日: 2026-08-21。上記はすべて一次文献の全文を読んで確認した内容である。
