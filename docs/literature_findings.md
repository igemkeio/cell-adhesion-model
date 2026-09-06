# 全文照合で分かったこと

一次文献 **11 本**の全文を読んで照合した結果。設計に直接効くものを先に置く。

| # | 文献 | 全文照合 |
|---|---|---|
| 1 | Costan et al. (2024) *Matter*（本文＋Supplemental Note 1–6） | ✅ |
| 2 | Rassam et al. (2015) *Nature*（本文＋Extended Data） | ✅ |
| 3 | Klamecka et al. (2015) *Phys. Biol.* | ✅ |
| 4 | Lim &amp; Bernhardt (2019) *Mol. Microbiol.* | ✅ |
| 5 | Lin et al. (2021) *Nat. Commun.* | ✅ |
| 6 | Glass &amp; Riedel-Kruse (2018) *Cell* | ✅ |
| 7 | Chou et al. (2025) *Commun. Biol.* | ✅ |
| 8 | Kirchhofer et al. (2010) *Nat. Struct. Mol. Biol.* | ✅ 本文のみ（補足は未入手） |
| 9 | Spector et al. (2010) *Biophys. J.* | ✅ |
| 10 | Leo et al. (2015) *Mol. Microbiol.* | ✅ |
| 11 | Salema et al. (2013) *PLoS One* | ✅ |

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

## 5. 実験⑤（LysM 除去）には、見落としやすい副作用がある

**Leo et al. (2015) *Mol. Microbiol.*** — [doi:10.1111/mmi.12840](https://doi.org/10.1111/mmi.12840)

Intimin のペリプラズムドメイン（LysM を含む）が何をしているかを調べた論文。
実験⑤はこの部分を CoilE に置き換える計画なので、直接関係する。

**除去してよさそうな理由**

> the periplasmic domain of Intimin does bind to peptidoglycan both *in vitro* and
> *in vivo*, but **only under acidic conditions**. ... dissociation constant of 0.8 µM

ペプチドグリカン結合は**酸性条件でしか起きない**。通常の培養（中性付近）では
LysM は結合していないと考えられるので、除去しても失うものは少ない。
著者らも「消化管を通過する際の機械的・化学的ストレスへの耐性に役立つのだろう」
という感染過程の文脈で解釈している。

**しかし見落としやすい副作用がある**

> In contrast to previous reports, we demonstrate that the periplasmic domain of
> Intimin **mediates dimerisation**.

**ペリプラズムドメインは二量体化を担っている。**しかもこれは LysM を持つ
inverse autotransporter に共通の性質だとしている。

したがって LysM を CoilE に置き換えると、

- Intimin どうしの**自己二量体化がなくなる**
- 代わりに CoilE–CoilR という**ヘテロ二量体**で内膜アンカーに繋がる

設計の意図から見ればこれは望ましい置き換えである（自己会合を、極への係留に
差し替えることになる）。ただし**野生型と会合状態が変わる**ことは意識しておくべきで、
表面提示の効率や接着挙動が変わる可能性がある。資料の実験⑤が「膜に出ることの確認」を
含んでいるのは妥当である。

## 6. 拡散係数の「食い違い」の正体が分かった

**Spector et al. (2010) *Biophys. J.*** — [doi:10.1016/j.bpj.2010.10.029](https://doi.org/10.1016/j.bpj.2010.10.029)

以前「Costan の 0.36 µm²/s と Rassam の 0.013 µm²/s が約20倍食い違う」と書いたが、
**外膜タンパクの拡散係数はもともと2桁にわたって広がる**ことが分かった。

| タンパク | *D* (µm²/s) | 備考 |
|---|---:|---|
| OmpF | **0.006 ± 0.002** | **直径 100 nm の領域に閉じ込められている** |
| BtuB（抗体標識） | 0.05 ± 0.01 | |
| BtuB（コリシン標識） | 0.10 ± 0.02 | |
| BtuB（TonB box 変異体） | **0.27 ± 0.06** | 5倍速い。TonB が繋ぎ止めていた |
| Rassam の BtuB / Cir | 0.013 / 0.019 | |
| Costan の合成アドヘシン | 0.36 | 可動割合 0.37 |

**繋ぎ止めが外れた BtuB は 0.27 µm²/s に達する。**したがって Costan の 0.36 は
範囲の上端ではあるが、外れ値ではない。「食い違い」ではなく**タンパク依存の広がり**
として理解すべきである。Costan が期待値として挙げる 0.05–0.6 µm²/s も
この論文の値に由来する。

### A11（島構造による下限）の前提が変わる

より重要なのは、**OmpF が直径 100 nm の領域に閉じ込められている**という点である。

Rassam が報告した島は直径 0.5 µm で、これが「先端集中度 *w* ≤ 0.20 µm は
島より細かいので到達できないかもしれない」という仮定 A11 の根拠だった。
しかし **100 nm スケールの閉じ込めも実在する**のなら、閉じ込めの単位は
タンパクによって変わり、0.5 µm が下限とは限らない。

**A11 の懸念は弱まる。**ただし「100 nm に閉じ込められる」ことと「100 nm に
集められる」ことは別なので、否定されたわけではない。

## 7. 表示系は K-12 で確立している

**Salema et al. (2013) *PLoS One*** — [doi:10.1371/journal.pone.0075126](https://doi.org/10.1371/journal.pone.0075126)

> the native β-domains of **EhaA autotransporter and intimin** ... are effective
> systems for the display of immune libraries of single domain Abs ... on the
> surface of ***E. coli* K-12 cells**

**EHEC O157:H7 由来の intimin βドメインを、K-12 で表示している。**
K-12 は O抗原を持たないため、以前述べた「立体障害の問題は K-12 なら生じない」
という見立てと整合する。表示系そのものは確立していると考えてよい。

なお Costan が期待値 6,000–8,000 分子/細胞の出典としてこの論文を挙げているが、
本文の検索では該当する数値を特定できなかった（定量は Western blot と既知量の
標準品との比較で行われており、値は図中にあると思われる）。**この期待値は
未確認のままである。**

## 8. 接着は自発的には不可逆である — 束から鎖への組み替えはない

複数の文献が一致してこれを支持する。

**Costan et al. (2024) Supplemental Note 4**

> the slope is not significantly different than 0, and due to instrument resolution,
> ultimately the best estimate is that **k_off should be assumed zero, i.e.,
> adhesins do not spontaneously disengage**

**Glass &amp; Riedel-Kruse (2018)**

> mechanical agitation reversibly disrupts the **essentially static** aggregates
> engineered in this work, but genetically encoding an excreted competitive
> inhibitor ... or intimin-specific protease **could control dynamic, reversible
> adhesion**

同論文はさらに「本研究が扱ったのは主に定常状態の凝集体であって、その背後にある
**動力学ではない**」と述べている。

したがって：

| | 可逆か |
|---|---|
| **自発的な解離** | **しない**（k_off ≈ 0） |
| **力による解離** | **する**。光ピンセットで結合・解離を繰り返せる |
| **凝集体の撹拌による崩壊** | する（ただし逸話的な観察） |
| **動的・可逆な接着** | **別途作り込みが必要**（競合阻害剤やプロテアーゼ） |

なお Costan は繰り返し結合・解離させたときにアドヘシンが傷むかも調べており、
接触時間 1 s で 26 回引いた実験では最初の10回で 1.1 ± 0.2 本、最後の10回で
0.5 ± 0.2 本と、**繰り返しにつれて形成される結合が減る**傾向を報告している
（原著は「統計的有意には達しない、t検定 p = 0.045」と評価している）。

### モデルの結論の読み方が変わる

本モデルは「一様分布では側面接触が 100 対、極接触が 8 対だから束になる」と
予測してきた。不可逆であることを踏まえると、これは**単なる平衡の傾きではない**。

**最初に成立した接触がそのまま固定される。**細胞が束から鎖へ組み替わる経路はない。
したがって形態を決めるのは「どちらが熱力学的に有利か」ではなく、
**「どちらの接触が先に成立するか」**である。

これは以前の解析と噛み合う。せん断下の接触時間（γ̇ = 300 で 3.3 ms）は結合形成の
時間（t₀.₅ ≈ 100 ms）より短く、その間に熱運動で向きを直せる角度は 4.4° で
接触パッチの角度許容（約 8°）に届かない。**向きを直す余地がないまま、
先に触れた配置で固定される。**

### 設計上の含意

方向性のある組み立てを狙うなら、**可逆性はむしろ欲しい性質**かもしれない。
不可逆だと誤った配置が固定されるが、可逆なら正しい配置へ焼きなませる。

Glass はその手段として**分泌型の競合阻害剤**と **intimin 特異的プロテアーゼ**を
挙げている。本モデルは可逆性を扱っていないので、この方向を評価するには
仮定 A7 を外して速度論を入れる必要がある。

## 9. 削除した k_off と K_D は、まだ復活できない

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
| **OmpF は直径 100 nm の領域に閉じ込められる** | **A11 の懸念が弱まった。0.5 µm が下限とは限らない** |
| 外膜タンパクの拡散は 0.006–0.27 µm²/s に広がる | 「20倍の食い違い」という以前の記述を撤回。タンパク依存の広がりである |
| Intimin ペリプラズムドメインは二量体化を担う | 実験⑤で LysM を除くと会合状態が変わる。要注意 |
| LysM の PG 結合は酸性条件のみ | 通常培養では除去の影響が小さい見込み |
| 表示系は K-12 で確立 | O抗原がないため立体障害の懸念は小さい |

---
照合日: 2026-08-21。上記はすべて一次文献の全文を読んで確認した内容である。
