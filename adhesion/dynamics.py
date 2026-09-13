"""細胞集団が時間とともに衝突し、鎖や束を作っていく過程をシミュレートする。

`adhesion.model` は「すでに接触した2細胞がどちらの配置を取るか」という相対的な
有利不利しか答えない（仮定 A7）。ここではその制約を外し、多数の細胞が z=0 平面上を
ブラウン運動で動き回り、実際に出会い、不可逆に結合し、鎖・束へと育っていく過程を
陽に時間発展させる。

設計の要点
----------
1. **平面拘束・厳密な幾何**：細胞は z=0 平面内で並進・回転する（2次元）。ただし
   接触判定はスフェロシリンダー（線分＋半径）の厳密な最近接距離計算であり、
   `adhesion.geometry` の3次元の求積を単純化したものではない。
2. **結合数は集団共通の定数**：同じ集団内の細胞はすべて同じ発現量・局在
   （`n_total`, `f_pole`, `tip_width_um`, `polar_split`）を持つので、「極どうしが
   接触したときに何本結合するか」「側面どうしが接触したときに何本結合するか」は
   `adhesion.model.predict` が1度返す2つの定数（`polar_bond_count`,
   `side_bond_count`）で決まる。動的シミュレーションが新たに扱うのは「いつ・どの
   組み合わせで接触するか」という幾何と時間の部分だけであり、結合数そのものの
   計算式は静的モデルと共有する。
3. **結合の解離は既定で0、ただしパラメータ化されている（`k_off_per_s`）**：
   既定値は Costan et al. (2024) の結論（k_off ≈ 0、`SPONTANEOUS_OFF_RATE`）を
   そのまま延長し、一度成立した結合は解離しない。したがって形態は熱力学的な
   有利不利ではなく、「どの接触が先に成立するか」という経路依存の過程になる
   （`docs/literature_findings.md` の該当箇所を参照）。ただし原著のこの結論は
   荷重速度 5–100 pN/s という狭い範囲で傾きが検出できなかったことに基づく
   検出限界の主張であり（Supplemental Note 4 原文で確認済み）、この結合
   （anti-EPEA/EPEA）は他の測定に使われている例外的に安定な anti-GFP/GFP
   ペアとは異なる、より低親和性のペアである。短時間のシミュレーション（鎖形成・
   震盪デモ）ではこの不確かさは結論を左右しないが、24時間の沈降モードでは
   無視できない可能性がある（A24）。`k_off_per_s` を正の値にすると、結合は
   `_apply_spontaneous_dissociation` により毎ステップ確率的に解離しうる。
4. **結合形成は速度論的**：`literature.BOND_FORMATION_HALF_TIME`（≈0.1秒）を
   ベルヌーイ過程の特徴時間として使う。接触が長く続く静置条件ではほぼ確実に
   結合が成立し、震盪条件で接触時間が短くなる場合には結合が成立しないまま
   離れることもありうる（`shear.contact_duration_s` と比較するための土台になる）。
5. **結合は有限剛性の柔軟なジョイント**：細胞は剛体クラスタにまとめず、各細胞が
   独立にブラウン運動する（自由排水近似、A16）。結合は2種の位置ベース拘束を持つ
   ジョイントとして扱う：(a) 距離拘束（硬め）で付着点を架橋長に保ち、(b) 曲げ拘束
   （柔らかめ、A19）で相対角を成立時の値へゆるく戻す。曲げ剛性が有限なので接合は
   熱運動でしなり、鎖は半柔軟、束は乱れた配置を取る。結合は成立した瞬間の接触
   角度をそのまま自然長とするので、斜めの接触やT字接合も生じ、平行や一直線を
   強制しない。剛性ゼロだと結合が伸びきり細胞が貫入するため（初版の破綻）、
   有限の曲げ剛性が要る。この剛性は実測値ではなく、鎖・束の柔らかさを決める
   実効パラメータである。排除体積は硬い拘束として貫入を防ぐ（結合が架橋長を
   保つことと貫入しないことは `tests/test_dynamics.py` で固定している）。
6. **震盪培養（`shear_rate_per_s > 0`）**：既定（0）は静置・拡散のみ。正の値を
   渡すと、`shear.py` と同じ物理量（せん断速度）で平面クエット流を課す。
   流れ方向（x）は周期境界、勾配方向（y）は壁とし（A21）、細胞は決定論的な
   移流とブラウン運動を両方受ける。回転には Jeffery (1922) の軌道式を等価
   アスペクト比の回転楕円体近似で適用し（A20）、流れに対して向きを変えながら
   タンブリングする。結合形成の判定式自体は静置条件と同じ（4）だが、震盪では
   接触が長続きしないため、同じ式が自然に「出会っても結合しないことがある」
   という帰結を生む——`contact_duration_s` を別途参照する必要はなく、実際の
   軌道から接触の持続時間が emergent に決まる。せん断荷重が破断力を超えた
   極結合（鎖の背骨）は切れる（A22）。文献の結論（震盪培養の典型的な速度では
   接合を壊さない）どおり、通常の shear_rate_per_s では発火しないはずである。
7. **沈降（`well_depth_um` を指定）**：Glass & Riedel-Kruse (2018) の実際の
   凝集アッセイは、震盪ではなく「短時間のボルテックス混合 → 深底96穴プレートで
   静置 → 上層の光学密度で判定」という手順である（要旨レベルで確認）。震盪
   モードの「もし振盪し続けたら」という仮想条件とは対照的に、こちらが実際の
   アッセイに対応する。大きい凝集塊ほど Stokes 則で速く沈み（クラスタの
   等価球半径は A16 と同じ n^(1/3) スケーリング）、鉛直方向に離れすぎた細胞
   どうしは結合できない。`supernatant_fraction` が、実際に測定される「上層に
   残る割合」に対応する。密度差の値は要旨レベルの参照であり全文照合していない
   ため、A18・A19 と同じ性格の実効パラメータとして扱う（A23）。ボルテックスに
   よる初期混合そのものは、既存のランダムな初期配置がすでに表しているとみなし、
   別途シミュレートしない。
8. **能動運動（`swim_speed_um_s` を指定、既定は無効）**：既定では水平面内の
   運動は `_brownian_step` の受動拡散だけである。`swim_speed_um_s > 0` に
   すると、結合していない自由な細胞は run（直進遊泳）と tumble（向きを
   変える）をポアソン過程的に切り替える run-and-tumble 運動をする（A25）。
   結合した細胞は能動推進を止める（クラスタとの能動-弾性的な相互作用は範囲外
   とする簡略化）。速度・run/tumble時間・旋回角度の既定値は Berg & Brown
   (1972) の野生株 AW405 の値（全文照合済み）。

   **重要——この機能は既定で無効であり、そうしているのには理由がある。**
   当初は「凝集は速く進むので受動拡散だけでは説明しにくい」という動機で
   追加したが、検証したところこの動機は成り立たなかった。スモルコフスキーの
   会合速度で見積もると、実際のアッセイの典型的な菌体濃度（OD600 0.1–1、
   すなわち ~10⁸–10⁹ 細胞/mL）では受動拡散だけで凝集の半減時間は ~4–40 分
   であり、1時間以内の凝集は拡散だけで十分に説明がつく（`_A_EFF_UM` と
   `D_TRANS0_UM2_S` から計算できる）。加えて、(a) 実際のアッセイ条件——定常期
   まで培養し、深底プレートで静置（酸素律速になりうる）、合成アドヘシンを
   大量発現（外膜への負荷）——では、指数増殖期の活発な AW405 ほど運動性が
   保たれるかは疑わしい。(b) アッセイの読み出しは沈降（上清の光学密度）で
   あり、遊泳はむしろ沈降に抗する向きに働く——つまり能動運動は測定される
   凝集を速めるとは限らず、遅らせる可能性すらある。したがって運動性モードは
   「もし細胞が活発に泳いだら会合はどれだけ速まりうるか」を探る任意の
   探索ツールであって、実際のアッセイの既定条件ではない（A25）。

このモジュールは matplotlib に依存しない（可視化は `scripts/simulate_chain_formation.py`
が担う）。

依存する仮定: A7（の時間発展への拡張）、A15、A16、A17、A18、A19、A20、A21、A22、
A23、A24、A25。詳細は `adhesion.assumptions` を参照。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from . import literature as lit
from . import model
from . import shear
from .assumptions import uses

DYNAMICS_ASSUMPTIONS = uses(
    "A5", "A7", "A8", "A9", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23",
    "A24", "A25",
)

# --------------------------------------------------------------------------
# 幾何定数（`adhesion.model` と同じ細胞寸法を使う）
# --------------------------------------------------------------------------
CELL_LENGTH_UM = model.CELL_LENGTH_UM
CELL_RADIUS_UM = model.CELL_RADIUS_UM
HALF_CYLINDER_UM = 0.5 * (CELL_LENGTH_UM - 2.0 * CELL_RADIUS_UM)
BRIDGE_UM = model.BRIDGE_UM
MIN_BONDS_FOR_CONTACT = model.MIN_BONDS_FOR_CONTACT

BondKind = Literal["polar", "lateral"]

# --------------------------------------------------------------------------
# 単体の拡散係数。細胞と同体積の球に対する Stokes–Einstein で見積もる。
# `shear.py` の `reorientation_during_contact_deg` が既定値として使っている
# 回転拡散係数 0.889 rad^2/s と同じ近似（同一の等価球体積）であり、独立に
# 再計算しても一致することを `tests/test_dynamics.py` で確認している。
# --------------------------------------------------------------------------


def _equivalent_sphere_radius_um() -> float:
    """細胞（スフェロシリンダー）と同じ体積を持つ球の半径（µm）。"""

    r = CELL_RADIUS_UM
    cylinder_length = CELL_LENGTH_UM - 2.0 * r
    volume_um3 = (4.0 / 3.0) * math.pi * r**3 + math.pi * r**2 * cylinder_length
    return (3.0 * volume_um3 / (4.0 * math.pi)) ** (1.0 / 3.0)


_A_EFF_UM = _equivalent_sphere_radius_um()
_A_EFF_M = _A_EFF_UM * 1e-6

#: 単体の並進拡散係数（µm^2/s）。kT / (6 π η a)。
D_TRANS0_UM2_S = (
    shear.THERMAL_ENERGY_J / (6.0 * math.pi * shear.DYNAMIC_VISCOSITY_PA_S * _A_EFF_M)
) * 1e12

#: 単体の回転拡散係数（rad^2/s）。kT / (8 π η a^3)。
D_ROT0_RAD2_S = shear.THERMAL_ENERGY_J / (
    8.0 * math.pi * shear.DYNAMIC_VISCOSITY_PA_S * _A_EFF_M**3
)

# --------------------------------------------------------------------------
# 沈降（A23）。Glass & Riedel-Kruse (2018) の実際の凝集アッセイは、震盪ではなく
# 「短時間のボルテックス混合 → 深底96穴プレートで24時間静置 → 上層のOD600測定」
# という手順である（要旨レベルで確認、全文未照合）。大きい凝集塊ほど重力で速く
# 沈み、測定領域（上層）から抜けることで初めて「凝集した」と判定される。
# したがって、この実際のアッセイに即した時間発展では、水平面内の運動に加えて
# 鉛直方向の沈降を持つ必要がある。
#
# 密度差・E. coli の浮遊密度は、複数の報告（培地の浸透圧に依存し
# 1.08–1.16 g/cm^3 の範囲）を要旨レベルで参照したものであり、本リポジトリでは
# 全文照合していない。A18（アビディティ基準点）・A19（ジョイント剛性）と同じ
# 性格の、桁の見積もりとしての実効パラメータとして扱う。
# --------------------------------------------------------------------------
#: E. coli の浮遊密度（g/cm^3）。培地の浸透圧に依存する（要旨レベル参照、A23）。
CELL_DENSITY_G_CM3 = 1.10
#: 培地（LB・PBS等）の密度（g/cm^3）。既存の `shear.DYNAMIC_VISCOSITY_PA_S` と
#: 同じく水で近似する。
MEDIUM_DENSITY_G_CM3 = 1.00
_DENSITY_DIFFERENCE_KG_M3 = (CELL_DENSITY_G_CM3 - MEDIUM_DENSITY_G_CM3) * 1000.0
_GRAVITY_M_S2 = 9.81


def _settling_velocity_um_s(cluster_size: int) -> float:
    """クラスタの沈降速度（µm/s、下向きを正）。Stokes 則。

    等価球半径は A16 の並進・回転拡散係数と同じ n^(1/3) スケーリングを使う。
    重力による沈降は熱ゆらぎではなく力の釣り合いで決まる決定論的な量なので、
    A16（水平面内の運動は非剛体・自由排水）とは独立に、クラスタ全体としての
    大きさで決まるとする。1細胞では 24 時間で高々数 µm しか沈まないが、
    十分大きな凝集塊（およそ20細胞以上）は深底プレートの液柱を24時間で
    沈みきる——これが実際のアッセイで「大きい塊だけが上清から抜ける」という
    選別が起こる理由である。
    """

    radius_um = _A_EFF_UM * cluster_size ** (1.0 / 3.0)
    radius_m = radius_um * 1e-6
    velocity_m_s = (
        (2.0 / 9.0) * _DENSITY_DIFFERENCE_KG_M3 * _GRAVITY_M_S2 * radius_m**2
        / shear.DYNAMIC_VISCOSITY_PA_S
    )
    return velocity_m_s * 1e6

# --------------------------------------------------------------------------
# 数値法のパラメータ（物理定数ではなく、拘束緩和の技術的な選択）
# --------------------------------------------------------------------------
#: 各ステップで拘束（結合・排除体積・壁）を緩和する反復回数。
N_RELAXATION_ITERS = 10
#: 結合の距離拘束の強さ（1反復あたりの補正率, 0-1）。硬め＝結合はほぼ架橋長を保つ。
BOND_DISTANCE_STIFFNESS = 0.7
#: 結合の曲げ拘束の強さ（1反復あたりの補正率, 0-1）。柔らかめ＝接合は自然長の
#: 相対角へゆるく戻るだけで、熱運動で角度が揺らぐ。これが0だと初版のように
#: 結合が伸びきって細胞が貫入する（＝剛性ゼロの破綻）。1に近いほど剛体に近づく。
#: 実測値ではなく、鎖・束の柔らかさを決める実効パラメータ（A19）。
BOND_BENDING_STIFFNESS = 0.3
#: 排除体積（貫入）を押し返す強さ（1反復あたりの補正率, 0-1）。硬め。
EXCLUSION_STIFFNESS = 1.0
#: 結合候補として調べる表面間ギャップの許容幅（µm）。架橋長 8 nm 自体は細胞
#: スケールでは無視できるほど小さいので、離散時間の残差を吸収するための
#: 数値的な捕捉半径として少し広めに取る。
BOND_CAPTURE_GAP_UM = 0.05
#: 極とみなす、線分パラメータ s（0=尾側端, 1=頭側端）からの許容幅。
POLE_PARAMETER_TOLERANCE = 0.03
#: 半球を「先端寄り半分」と「赤道寄り半分」に等面積で分ける極角のコサイン
#: （cos 60° = 0.5）。線分＋半径というスフェロシリンダーの表現では、半球上の
#: どこに接触しても最近接軸上点は同じ端点に丸め込まれてしまう（s=0 または 1）。
#: s だけで極接触と判定すると、実際には赤道寄り（＝実質的に側面と地続きの領域）
#: での接触まで極接触として過大評価してしまう。そこで、接触方向と軸方向との
#: なす角がこの閾値より小さい場合だけを「先端寄り」＝極接触とみなし、それ以外は
#: 側面接触として扱う。
POLE_APEX_COS_THRESHOLD = 0.5
#: 結合数から実効オンレートへの換算の基準点。`BOND_FORMATION_HALF_TIME` は
#: Costan et al. (2024) の光ピンセット実験（30 ng/mL 誘導、接触あたり0-2結合と
#: いう少数の多重度）での見かけの値であり、そのままでは「弱い」接触
#: （結合数が一桁程度）でもオンレートが速すぎ、幾何的にどちらの接触型が先に
#: 起こりやすいか（棒状細胞の衝突頻度）だけで局在の効果が埋もれてしまう。
#: そこで、この基準結合数を跨いで実効オンレートが大きく変わるように、
#: `BOND_FORMATION_HALF_TIME` をこの結合数に対する値とみなして規格化する。
#: 一次文献に基準点そのものの実測はないため、桁の見積もりとして置いた自由
#: パラメータである（母集団密度・視野サイズなどと同じ性格）。
REFERENCE_BOND_COUNT = 40.0


def _unit_vector(theta: float) -> np.ndarray:
    return np.array([math.cos(theta), math.sin(theta)])


@dataclass(frozen=True)
class DynamicsConfig:
    """時間発展シミュレーションの設定。

    局在パラメータ（`f_pole`, `tip_width_um`, `polar_split`, `n_total`,
    `distinct_adhesins_per_pole`）は `adhesion.model.predict` とまったく同じ意味で、
    既定値も同じにしてある。新規に導入するのは、母集団の規模・視野・時間刻みなど、
    「時間と多数の細胞」を扱うために避けられないパラメータだけである。

    母集団密度・視野サイズ・総シミュレーション時間は、Costan et al. (2024) にも
    他の一次文献にも較正対象となる経時データが存在しないため、「動画として
    見やすい」という基準で選んだ自由パラメータである（本文・Supplemental
    Information を確認済み）。実測の凝集速度データが得られ次第、差し替えるべき値。
    """

    n_cells: int = 24
    box_size_um: float = 20.0

    n_total: float = float(lit.ADHESINS_AT_MID_INDUCTION)
    f_pole: float = 0.95
    polar_split: float = 0.5
    tip_width_um: float | None = 0.10
    distinct_adhesins_per_pole: bool = True

    dt_s: float = 0.02
    total_time_s: float = 120.0
    frame_stride: int = 25

    #: 震盪培養のせん断速度（1/s）。0（既定）なら静置・拡散のみ。>0 にすると
    #: x 方向を周期境界とした平面クエット流（A21）を課し、細胞は shear.py と
    #: 同じ物理量（せん断速度）に応じて流され・回転し・出会う。ただし
    #: Glass & Riedel-Kruse (2018) の実際の凝集アッセイは震盪中には行われない
    #: （育成培養が 300 rpm、凝集アッセイ自体は短いボルテックス混合のあと
    #: 静置）。したがってこのモードは、実際に報告された凝集アッセイの条件では
    #: なく、「もし凝集の間も振盪し続けたら」という仮想的な条件である。
    shear_rate_per_s: float = 0.0

    #: 深底プレートでの沈降を有効にする鉛直方向の液柱の深さ（µm）。None（既定）
    #: なら沈降を計算しない（水平面内の運動のみ）。正の値を渡すと、Glass &
    #: Riedel-Kruse (2018) の実際の凝集アッセイ（震盪ではなく、短時間の
    #: ボルテックス混合 → 深底96穴プレートで静置 → 上層の光学密度で判定）に
    #: 対応する時間発展になる（A23）。深底96穴プレートの物理的な深さは製品
    #: 仕様でおよそ 30–44 mm だが、実際の充填量（液柱の高さ）は原著に明記
    #: されていないため、この値自体も見やすさ・計算量との兼ね合いで選ぶ
    #: 自由パラメータである。ボルテックスによる初期混合そのものは、既存の
    #: 一様ランダムな初期配置（`_initial_state`）がすでに表しているとみなす
    #: （別途シミュレートしない）。
    well_depth_um: float | None = None
    #: 上清として測定される領域の割合。原著は上部 ~25% を測定に使う。
    #: `well_depth_um` が None のときは使われない。
    supernatant_fraction: float = 0.25

    #: 結合の自発的な解離速度定数（1/s）。既定 0 は Costan et al. (2024) の
    #: 結論（k_off はゼロとみなすのが最良の推定）に合わせた不可逆凝集（A15）。
    #: ただしこの結論は、荷重速度 5–100 pN/s という狭い範囲で傾きが検出でき
    #: なかったことに基づく「検出限界」の主張であり、k_off が数学的に厳密な
    #: ゼロだと証明したものではない（Supplemental Note 4 原文で確認済み）。
    #: この結合（anti-EPEA / EPEA）は、Costan et al. が他の測定に使っている
    #: 例外的に安定な anti-GFP / GFP ペア（市販品「GFP-Trap」、Kd ~ 1 pM）とは
    #: 別の、より低親和性（ナノモル）のペアである。他のナノボディ-抗原ペアの
    #: SPR 測定では k_off ~ 10⁻³–10⁻² /s（半減期にして数十秒〜十数分）が
    #: 珍しくない。短時間のシミュレーション（鎖形成・震盪デモ）ではこの不確かさ
    #: は結論を左右しないが、24時間の沈降モードでは無視できない可能性がある
    #: （A24）。非零の値を渡すと、結合は毎ステップ確率的に解離しうる。
    k_off_per_s: float = 0.0

    #: 遊泳速度（µm/s）。0（既定）なら能動運動なし（水平面内は受動拡散のみ）。
    #: 正の値を渡すと、結合していない自由な細胞は run-and-tumble 運動をする
    #: （A25）。既定を 0（無効）にしているのは意図的である——実際のアッセイの
    #: 典型濃度では受動拡散だけで凝集は数分〜数十分で進み、能動運動を持ち込む
    #: 必要がないうえ、アッセイ条件（定常期・静置・アドヘシン負荷）で運動性が
    #: 保たれるかも疑わしく、遊泳は沈降の読み出しに抗する向きにも働くためである
    #: （詳細は A25 とモジュール docstring の 8）。有効にする場合の目安は
    #: `lit.SWIM_SPEED`（Berg & Brown 1972、野生株 AW405、全文照合済み）。
    swim_speed_um_s: float = 0.0
    #: run 相の平均継続時間（秒）。指数分布を仮定する（ポアソン過程——原著
    #: 自身が自己相関関数からこれを確認している）。既定値は Berg & Brown
    #: (1972) の野生株 AW405 の値（全文照合済み、`lit.MEAN_RUN_TIME`）。
    mean_run_time_s: float = float(lit.MEAN_RUN_TIME)
    #: tumble 相の平均継続時間（秒）。既定値は同上（`lit.MEAN_TWIDDLE_TIME`。
    #: 原著の用語は twiddle だが、本モデルは後年の標準的な呼び方に合わせる）。
    mean_tumble_time_s: float = float(lit.MEAN_TWIDDLE_TIME)
    #: tumble で向きを変える角度の平均（度）。既定値は同上（`lit.TUMBLE_TURN_ANGLE`）
    #: で、3次元追跡での測定値を本モデルの平面内の運動にそのまま転用している。
    #: 実際の分布は0°に偏らない非一様な形（原著 Fig. 3）だが、ここでは平均・
    #: 標準偏差が一致するガンマ分布で近似する（A25）。
    tumble_turn_angle_deg: float = float(lit.TUMBLE_TURN_ANGLE)
    #: tumble で向きを変える角度の標準偏差（度）。既定値は Berg & Brown (1972)
    #: の野生株 AW405 の値。ガンマ分布の形状パラメータを決めるのに使う。
    tumble_turn_angle_std_deg: float = float(lit.TUMBLE_TURN_ANGLE.uncertainty or 36.0)

    seed: int = 0

    def __post_init__(self) -> None:
        if self.polar_split != 0.5:
            raise ValueError(
                "動的シミュレーションは今のところ polar_split = 0.5（対称二極）"
                "のみに対応する。非対称分配は極ごとの結合数定数が異なり、"
                "頭・尾のラベルと結合数の対応づけを別途行う必要がある。"
            )
        if self.n_cells < 2:
            raise ValueError("n_cells は2以上である必要があります")
        if self.box_size_um <= CELL_LENGTH_UM:
            raise ValueError("box_size_um は細胞長より十分大きい必要があります")
        if self.shear_rate_per_s < 0.0:
            raise ValueError("shear_rate_per_s は非負である必要があります")
        if self.well_depth_um is not None and self.well_depth_um <= 0.0:
            raise ValueError("well_depth_um は正である必要があります")
        if not 0.0 < self.supernatant_fraction <= 1.0:
            raise ValueError("supernatant_fraction は 0 より大きく 1 以下である必要があります")
        if self.k_off_per_s < 0.0:
            raise ValueError("k_off_per_s は非負である必要があります")
        if self.swim_speed_um_s < 0.0:
            raise ValueError("swim_speed_um_s は非負である必要があります")
        if self.mean_run_time_s <= 0.0:
            raise ValueError("mean_run_time_s は正である必要があります")
        if self.mean_tumble_time_s <= 0.0:
            raise ValueError("mean_tumble_time_s は正である必要があります")
        if not 0.0 <= self.tumble_turn_angle_deg <= 180.0:
            raise ValueError("tumble_turn_angle_deg は0以上180以下である必要があります")

    @property
    def motile(self) -> bool:
        """能動的な run-and-tumble 運動を計算するか（A25）。"""

        return self.swim_speed_um_s > 0.0

    @property
    def periodic_x(self) -> bool:
        """せん断流の流れ方向（x）を周期境界にするか（A21）。"""

        return self.shear_rate_per_s > 0.0

    @property
    def settling_enabled(self) -> bool:
        """沈降（A23）を計算するか。"""

        return self.well_depth_um is not None


@dataclass(frozen=True)
class Bond:
    """成立した結合1本。

    axial_i / axial_j  結合成立時の、細胞の中心軸に沿った付着点の位置（µm、
                        中心からの符号付き距離）。頭側端は +HALF_CYLINDER_UM、
                        尾側端は -HALF_CYLINDER_UM に固定する。
    delta_i / delta_j   付着点の半径方向が細胞自身の向きから見てどの角度に
                        あるか（rad）。頭側極は 0、尾側極は π、側面はそのときの
                        実際の接触方向を保持する。細胞が回転すればこの角度も
                        一緒に回転する。
    rest_angle          結合成立時の相対角 (angle_j - angle_i)。曲げ拘束が
                        ゆるくこの角度へ戻そうとする自然長。斜めに触れて結合
                        すれば斜めが自然長になり、平行を強制しない。
    """

    i: int
    j: int
    kind: BondKind
    axial_i: float
    delta_i: float
    axial_j: float
    delta_j: float
    rest_angle: float = 0.0


@dataclass
class SimulationState:
    positions: np.ndarray  # (N, 2) µm
    angles: np.ndarray  # (N,) rad
    head_free: np.ndarray  # (N,) bool。頭側極（+側）がまだ空いているか
    tail_free: np.ndarray  # (N,) bool。尾側極（-側）がまだ空いているか
    bonds: list[Bond] = field(default_factory=list)
    #: 沈降の深さ（µm、下向きが正）。A23。沈降を計算しない設定では常に0のまま。
    depth_um: np.ndarray | None = None
    #: run-and-tumble 運動（A25）の相。True なら run（直進遊泳）、False なら
    #: tumble（向きを変える）。運動性を計算しない設定では使われない。
    running: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.depth_um is None:
            self.depth_um = np.zeros(self.positions.shape[0])
        if self.running is None:
            self.running = np.ones(self.positions.shape[0], dtype=bool)


@dataclass(frozen=True)
class Frame:
    time_s: float
    positions: np.ndarray
    angles: np.ndarray
    bonds: list[Bond]
    #: 沈降の深さ（µm）。沈降を計算しない設定では常に0。A23。
    depth_um: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass(frozen=True)
class ClusterInfo:
    """接続成分1つの分類。"""

    cells: list[int]
    kind: Literal["monomer", "chain", "bundle"]
    directed: bool


@dataclass(frozen=True)
class SimulationResult:
    config: DynamicsConfig
    prediction: model.AdhesionPrediction
    polar_bond_count: float
    side_bond_count: float
    frames: list[Frame]


def attachment_point(
    position: np.ndarray, theta: float, axial: float, delta: float
) -> np.ndarray:
    """付着点の実座標（µm）。中心軸上の位置 + 半径方向へのオフセット。"""

    axis = _unit_vector(theta)
    radial = _unit_vector(theta + delta)
    return position + axial * axis + CELL_RADIUS_UM * radial


def _closest_points_on_segments(
    a0: np.ndarray, a1: np.ndarray, b0: np.ndarray, b1: np.ndarray
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """2線分の最近接点と、それぞれの線分パラメータ（0-1）を返す。

    スフェロシリンダーは線分と半径のミンコフスキー和なので、この最近接距離から
    半径2つ分を引けば、球冠部分も含めて表面間ギャップが厳密に求まる。
    """

    d1 = a1 - a0
    d2 = b1 - b0
    r = a0 - b0
    a = float(d1 @ d1)
    e = float(d2 @ d2)
    f = float(d2 @ r)
    eps = 1e-12

    if a <= eps and e <= eps:
        s = t = 0.0
    elif a <= eps:
        s = 0.0
        t = float(np.clip(f / e, 0.0, 1.0))
    else:
        c = float(d1 @ r)
        if e <= eps:
            t = 0.0
            s = float(np.clip(-c / a, 0.0, 1.0))
        else:
            b = float(d1 @ d2)
            denom = a * e - b * b
            s = float(np.clip((b * f - c * e) / denom, 0.0, 1.0)) if denom > eps else 0.0
            t = (b * s + f) / e
            if t < 0.0:
                t = 0.0
                s = float(np.clip(-c / a, 0.0, 1.0))
            elif t > 1.0:
                t = 1.0
                s = float(np.clip((b - c) / a, 0.0, 1.0))

    p1 = a0 + s * d1
    p2 = b0 + t * d2
    return p1, p2, s, t


def _segment_endpoints(position: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    axis = _unit_vector(theta)
    return position - HALF_CYLINDER_UM * axis, position + HALF_CYLINDER_UM * axis


def _classify_contact_region(
    s: float, theta: float, own_axis_point: np.ndarray, partner_point: np.ndarray
) -> Literal["tail", "head", "side"]:
    """接触点が細胞のどの領域にあるかを判定する。

    最近接軸上点が線分の端に丸め込まれている（s ≈ 0 または 1）だけでは、まだ
    極接触とは言えない。丸め込みは半球上のどこに接触しても起こるからである。
    接触方向（軸上端点から相手側への向き）と、その端の外向き軸方向との
    なす角が小さい（半球の先端寄り半分に入っている）場合だけを極接触とする。
    """

    if POLE_PARAMETER_TOLERANCE < s < 1.0 - POLE_PARAMETER_TOLERANCE:
        return "side"

    is_head_end = s >= 1.0 - POLE_PARAMETER_TOLERANCE
    outward = partner_point - own_axis_point
    norm = float(np.linalg.norm(outward))
    if norm < 1e-9:
        return "head" if is_head_end else "tail"
    outward = outward / norm
    axis_direction = _unit_vector(theta) if is_head_end else -_unit_vector(theta)
    cos_angle = float(outward @ axis_direction)
    if cos_angle < POLE_APEX_COS_THRESHOLD:
        return "side"
    return "head" if is_head_end else "tail"


def _initial_state(config: DynamicsConfig, rng: np.random.Generator) -> SimulationState:
    n = config.n_cells
    positions = np.zeros((n, 2))
    angles = rng.uniform(0.0, 2.0 * math.pi, size=n)
    margin = HALF_CYLINDER_UM + CELL_RADIUS_UM
    low, high = margin, config.box_size_um - margin

    placed: list[np.ndarray] = []
    placed_angles: list[float] = []
    for k in range(n):
        for _attempt in range(500):
            candidate = rng.uniform(low, high, size=2)
            ok = True
            for other_pos, other_theta in zip(placed, placed_angles):
                a0, a1 = _segment_endpoints(candidate, angles[k])
                b0, b1 = _segment_endpoints(other_pos, other_theta)
                p1, p2, _, _ = _closest_points_on_segments(a0, a1, b0, b1)
                if float(np.linalg.norm(p1 - p2)) < 2.2 * CELL_RADIUS_UM:
                    ok = False
                    break
            if ok:
                positions[k] = candidate
                placed.append(candidate)
                placed_angles.append(angles[k])
                break
        else:
            # 500回試して置き場所が見つからない場合は、密度が高すぎる。
            # そのまま置いて、初回の拘束緩和で押し出させる。
            positions[k] = rng.uniform(low, high, size=2)
            placed.append(positions[k])
            placed_angles.append(angles[k])

    head_free = np.ones(n, dtype=bool)
    tail_free = np.ones(n, dtype=bool)
    return SimulationState(positions, angles, head_free, tail_free, [])


def _surface_gap(state: SimulationState, i: int, j: int) -> float:
    """細胞 i, j の表面間ギャップ（µm）。負なら貫入。"""

    a0, a1 = _segment_endpoints(state.positions[i], state.angles[i])
    b0, b1 = _segment_endpoints(state.positions[j], state.angles[j])
    p1, p2, _, _ = _closest_points_on_segments(a0, a1, b0, b1)
    return float(np.linalg.norm(p2 - p1)) - 2.0 * CELL_RADIUS_UM


def _wrap_to_pi(angle: float) -> float:
    return math.atan2(math.sin(angle), math.cos(angle))


# --------------------------------------------------------------------------
# 震盪培養（せん断流）。既定（shear_rate_per_s = 0）では一切呼ばれない。
# --------------------------------------------------------------------------
#: スフェロシリンダーを近似する回転楕円体のアスペクト比（全長 / 直径）。
#: Jeffery 軌道の式（A20）にはこの1つの形状パラメータしか現れない。
_ASPECT_RATIO = CELL_LENGTH_UM / (2.0 * CELL_RADIUS_UM)


def _minimum_image_dx(dx: float, box_size_um: float) -> float:
    """周期境界での最小像変位（µm）。"""

    return dx - box_size_um * round(dx / box_size_um)


def _image_of_j(
    pos_i: np.ndarray, pos_j: np.ndarray, box_size_um: float, periodic_x: bool
) -> np.ndarray:
    """周期境界が有効なとき、j の位置を i に最も近い像に置き換える。

    x 方向だけを周期境界にする（y は壁で仕切られた平面クエット流、A21）。
    `periodic_x=False` では pos_j をそのまま返すので、静置条件の経路は
    このモジュールが増える前とまったく同じ計算になる。
    """

    if not periodic_x:
        return pos_j
    dx = _minimum_image_dx(float(pos_j[0] - pos_i[0]), box_size_um)
    return np.array([pos_i[0] + dx, pos_j[1]])


def _wrap_periodic_x(state: SimulationState, box_size_um: float) -> None:
    state.positions[:, 0] = np.mod(state.positions[:, 0], box_size_um)


def _shear_flow_velocity_x(
    y_um: np.ndarray, box_size_um: float, shear_rate_per_s: float
) -> np.ndarray:
    """平面クエット流の x 方向流速（µm/s）。箱の中央で 0、上下の壁でそれぞれ逆符号。

    `shear.py` の一様せん断（A9）をそのまま速度場として実現したもの。
    """

    return shear_rate_per_s * (y_um - 0.5 * box_size_um)


def _jeffery_rotation_rate(theta: np.ndarray, shear_rate_per_s: float) -> np.ndarray:
    """単純せん断中の細長い粒子の回転角速度（rad/s）。Jeffery (1922) の軌道式。

    厳密には回転楕円体に対する式であり、スフェロシリンダーを等価アスペクト比
    （全長/直径）の回転楕円体で近似する（A20）。流れ方向（θ=0）付近で回転が
    遅く、垂直付近で速く回転する「タンブリング」を再現する。
    """

    r = _ASPECT_RATIO
    return (shear_rate_per_s / (r * r + 1.0)) * (
        np.cos(theta) ** 2 + r * r * np.sin(theta) ** 2
    )


def _apply_shear_drift(state: SimulationState, config: DynamicsConfig, dt_s: float) -> None:
    """せん断流による決定論的な移流・回転を加える（ブラウン項に上乗せ）。"""

    if config.shear_rate_per_s <= 0.0:
        return
    state.positions[:, 0] += (
        _shear_flow_velocity_x(state.positions[:, 1], config.box_size_um, config.shear_rate_per_s)
        * dt_s
    )
    state.angles += _jeffery_rotation_rate(state.angles, config.shear_rate_per_s) * dt_s


#: せん断中の接触持続時間 `shear.contact_duration_s` (≈ 1/γ̇) に対して、
#: dt をこの倍率以下に抑える（震盪モードのみ）。大きすぎる dt は「出会っても
#: 短時間しか触れ合わない」という物理を潰し、結合形成確率を過大評価する。
_SHEAR_DT_SAFETY_FACTOR = 0.1


def _safe_dt_for_shear(config: DynamicsConfig) -> float:
    """接触の持続時間を取りこぼさない上限まで、必要なら dt を自動的に縮める。

    震盪では接触時間が `shear.contact_duration_s(shear_rate) = 1/γ̇` 程度しかない
    （`docs/literature_findings.md` 参照）。dt がこれに対して大きすぎると、
    実際には一瞬しか触れ合わない接触を1ステップ丸ごと「接触中」として
    ベルヌーイ試行してしまい、結合確率を過大評価する。震盪が出会いを増やす
    だけでなく短命にもする、という定性的な効果（A18 の docstring 参照）を
    保つために、dt をこの持続時間に対して十分小さく取る。
    """

    if config.shear_rate_per_s <= 0.0:
        return config.dt_s
    safe_dt = _SHEAR_DT_SAFETY_FACTOR * shear.contact_duration_s(config.shear_rate_per_s)
    return min(config.dt_s, safe_dt)


#: 遊泳中の1ステップの移動量が細胞半径のこの割合を超えないよう dt を抑える
#: （運動性モードのみ）。大きすぎる dt は、接触判定（`BOND_CAPTURE_GAP_UM`
#: スケール）を素通りしたり、細胞どうしがすり抜けたりする原因になる。
_MOTILITY_DT_SAFETY_FACTOR = 0.2


def _safe_dt_for_motility(config: DynamicsConfig) -> float:
    """遊泳による移動量を取りこぼさない上限まで、必要なら dt を自動的に縮める。"""

    if not config.motile:
        return config.dt_s
    safe_dt = _MOTILITY_DT_SAFETY_FACTOR * CELL_RADIUS_UM / config.swim_speed_um_s
    return min(config.dt_s, safe_dt)


#: 回転に対する実効慣性（µm²）。細長い棒の中心まわりの慣性モーメント L²/12。
#: 拘束を付着点・接触点で作用させたとき、並進と回転へどう配分するかを決める。
_ROT_INERTIA_UM2 = CELL_LENGTH_UM**2 / 12.0
_INV_ROT = 1.0 / _ROT_INERTIA_UM2


def _cross2(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def _solve_distance(
    state: SimulationState,
    i: int,
    j: int,
    point_i: np.ndarray,
    point_j: np.ndarray,
    ref_i: np.ndarray,
    ref_j: np.ndarray,
    target_um: float,
    stiffness: float,
    repulsive_only: bool,
) -> None:
    """2つの物体固定点の距離を target に近づける位置ベース拘束。

    補正を点の位置で作用させ、剛体の並進と回転（慣性 `_ROT_INERTIA_UM2`）に
    配分する。これにより、単なる平行移動では解けない交差した棒どうしの貫入も、
    回転して離れることで解消できる。`repulsive_only=True` なら縮める向き
    （貫入の解消）だけに働く。

    `ref_i` / `ref_j` はてこの腕を測る基準点（通常は各細胞の中心）。周期境界
    （A21）で j の最近接像を使って point_j を計算した場合、`ref_j` にも同じ
    像の中心を渡す。実際の書き込みは常に `state.positions[j]` へ行うので、
    像を使ったことによる不整合は生じない。
    """

    d = point_j - point_i
    distance = float(np.linalg.norm(d))
    if distance < 1e-9:
        return
    error = distance - target_um
    if repulsive_only and error >= 0.0:
        return
    normal = d / distance
    r_i = point_i - ref_i
    r_j = point_j - ref_j
    cross_i = _cross2(r_i, normal)
    cross_j = _cross2(r_j, normal)
    weight = 2.0 + _INV_ROT * (cross_i * cross_i + cross_j * cross_j)
    lam = stiffness * error / weight
    state.positions[i] += lam * normal
    state.angles[i] += _INV_ROT * lam * cross_i
    state.positions[j] -= lam * normal
    state.angles[j] -= _INV_ROT * lam * cross_j


def _brownian_step(state: SimulationState, dt_s: float, rng: np.random.Generator) -> None:
    """各細胞を独立にブラウン運動させる（自由排水近似・非剛体、A16）。

    結合した細胞も剛体としてまとめて動かすのではなく、それぞれが単体の拡散係数で
    動く。結合は下の `_resolve_constraints` が扱う柔軟なジョイントであり、
    クラスタは内部自由度（曲げ・伸び）を持つ半柔軟な集合体としてふるまう。
    """

    n = state.positions.shape[0]
    state.positions += math.sqrt(2.0 * D_TRANS0_UM2_S * dt_s) * rng.normal(size=(n, 2))
    state.angles += math.sqrt(2.0 * D_ROT0_RAD2_S * dt_s) * rng.normal(size=n)


def _resolve_constraints(
    state: SimulationState,
    candidate_pairs: list[tuple[int, int]],
    box_size_um: float,
    iterations: int,
    periodic_x: bool = False,
) -> None:
    """結合・排除体積・壁の拘束を、位置ベースの反復緩和で満たす。

    結合は2種の拘束を持つ有限剛性のジョイントである。
    (1) 距離拘束（硬め）：付着点どうしを架橋長に保つ。
    (2) 曲げ拘束（柔らかめ）：相対角を成立時の値へゆるく戻す。剛性が有限なので
        接合は熱運動でしなり、鎖は半柔軟、束は乱れた配置を取る。剛性ゼロだと
        初版のように結合が伸びきり細胞が貫入するので、有限の曲げ剛性が要る。
    排除体積は硬い拘束として貫入を押し返す。

    `periodic_x=True`（震盪モード、A21）では、j 側の位置を i に最も近い周期像へ
    置き換えてから幾何を計算する（`_image_of_j`）。書き込みは常に実座標
    `state.positions[j]` へ行うので、像の選び方によらず結果は一貫する。
    `periodic_x=False` では常に j の実座標がそのまま使われ、静置条件の経路は
    震盪モードを導入する前とまったく同じ計算になる。
    """

    margin = HALF_CYLINDER_UM + CELL_RADIUS_UM
    for _ in range(iterations):
        # (1) 結合の距離拘束（付着点で作用させ、並進＋回転で満たす）。
        for bond in state.bonds:
            pos_j = _image_of_j(
                state.positions[bond.i], state.positions[bond.j], box_size_um, periodic_x
            )
            p_i = attachment_point(
                state.positions[bond.i], state.angles[bond.i], bond.axial_i, bond.delta_i
            )
            p_j = attachment_point(pos_j, state.angles[bond.j], bond.axial_j, bond.delta_j)
            _solve_distance(
                state, bond.i, bond.j, p_i, p_j,
                state.positions[bond.i], pos_j,
                BRIDGE_UM, BOND_DISTANCE_STIFFNESS, False,
            )

            # (2) 曲げ拘束（柔らかめ）：相対角を成立時の値へゆるく戻す。
            err = _wrap_to_pi(
                float(state.angles[bond.j] - state.angles[bond.i]) - bond.rest_angle
            )
            twist = 0.5 * BOND_BENDING_STIFFNESS * err
            state.angles[bond.i] += twist
            state.angles[bond.j] -= twist

        # 排除体積（接触点で作用させ、交差した棒どうしが回転して離れる）。
        for i, j in candidate_pairs:
            pos_j = _image_of_j(state.positions[i], state.positions[j], box_size_um, periodic_x)
            a0, a1 = _segment_endpoints(state.positions[i], state.angles[i])
            b0, b1 = _segment_endpoints(pos_j, state.angles[j])
            p1, p2, _, _ = _closest_points_on_segments(a0, a1, b0, b1)
            if float(np.linalg.norm(p2 - p1)) - 2.0 * CELL_RADIUS_UM >= 0.0:
                continue
            _solve_distance(
                state, i, j, p1, p2,
                state.positions[i], pos_j,
                2.0 * CELL_RADIUS_UM, EXCLUSION_STIFFNESS, True,
            )

        # 壁。周期境界にした軸（x）は壁で挟まない。
        if not periodic_x:
            np.clip(state.positions[:, 0], margin, box_size_um - margin, out=state.positions[:, 0])
        np.clip(state.positions[:, 1], margin, box_size_um - margin, out=state.positions[:, 1])


def _broad_phase_pairs(
    state: SimulationState, cutoff_um: float, box_size_um: float = 0.0, periodic_x: bool = False
) -> list[tuple[int, int]]:
    n = state.positions.shape[0]
    diff = state.positions[:, None, :] - state.positions[None, :, :]
    if periodic_x:
        diff[..., 0] -= box_size_um * np.round(diff[..., 0] / box_size_um)
    dist = np.linalg.norm(diff, axis=-1)
    rows, cols = np.triu_indices(n, k=1)
    close = dist[rows, cols] <= cutoff_um
    return list(zip(rows[close].tolist(), cols[close].tolist()))


def _already_bonded(state: SimulationState, i: int, j: int) -> bool:
    return any((bond.i == i and bond.j == j) or (bond.i == j and bond.j == i) for bond in state.bonds)


def _binding_probability(bond_count: float, dt_s: float) -> float:
    """時間刻み dt の間にその接触で結合が成立する確率。

    `BOND_FORMATION_HALF_TIME` は Costan et al. (2024) の光ピンセット実験（30 ng/mL
    誘導、接触あたり0-2結合という少数の分子多重度）での見かけの半減時間である。
    ここでは、独立な結合形成の試行が結合しうる分子の数だけ並列に起こるという
    多価結合のアビディティの標準的な描像に従い、実効オンレートを bond_count に
    比例させる：多数の分子が接触面にある配置ほど、速く結合が「ロック」される。

    比例係数は `REFERENCE_BOND_COUNT` で規格化してある。これがないと、幾何的な
    接触の残留時間（2次元の線分＋半径で衝突を近似したことに伴う副作用）だけで
    どちらの接触型が先に結合するかが決まってしまい、局在パラメータ（f_pole,
    tip_width_um）の効果が結合数の大小に正しく反映されなくなる。
    """

    rate_per_s = (
        (bond_count / REFERENCE_BOND_COUNT)
        * math.log(2.0)
        / float(lit.BOND_FORMATION_HALF_TIME)
    )
    return 1.0 - math.exp(-rate_per_s * dt_s)


def _bond_decision(
    state: SimulationState,
    i: int,
    j: int,
    config: DynamicsConfig,
    polar_bond_count: float,
    side_bond_count: float,
    dt_s: float,
    rng: np.random.Generator,
) -> tuple[BondKind, str, str, np.ndarray, np.ndarray, float, float] | None:
    """細胞 i, j が今この瞬間に結合するかを判定する（状態は変えない）。

    結合する場合は (kind, region_i, region_j, p1, p2, s, t) を返す。しないなら
    None。幾何（どちらの領域が触れているか）と、結合数に応じた確率（A18）で
    決める。p1, p2, s, t（最近接点と線分パラメータ）は `_record_bond` がそのまま
    使うので、周期境界（A21）で j の像をどちらに取ったかも含めて再計算せず
    引き継がれ、常に一貫する。沈降（A23）が有効なとき、鉛直方向に離れすぎた
    細胞どうしは物理的に接触できないため、その時点で判定を打ち切る。
    """

    if config.settling_enabled:
        if abs(float(state.depth_um[i] - state.depth_um[j])) > VERTICAL_CONTACT_GAP_UM:
            return None

    pos_j = _image_of_j(state.positions[i], state.positions[j], config.box_size_um, config.periodic_x)
    a0, a1 = _segment_endpoints(state.positions[i], state.angles[i])
    b0, b1 = _segment_endpoints(pos_j, state.angles[j])
    p1, p2, s, t = _closest_points_on_segments(a0, a1, b0, b1)
    gap = float(np.linalg.norm(p1 - p2)) - 2.0 * CELL_RADIUS_UM
    if gap > BOND_CAPTURE_GAP_UM:
        return None

    region_i = _classify_contact_region(s, state.angles[i], p1, p2)
    region_j = _classify_contact_region(t, state.angles[j], p2, p1)

    if region_i == "side" or region_j == "side":
        # 少なくとも一方が側面である接触は、すべて側面接触として扱う（A17）。
        if side_bond_count < MIN_BONDS_FOR_CONTACT:
            return None
        if rng.random() >= _binding_probability(side_bond_count, dt_s):
            return None
        return ("lateral", "side", "side", p1, p2, s, t)

    # ここに来るのは region_i, region_j がともに極（head/tail）の場合。
    if config.distinct_adhesins_per_pole and region_i == region_j:
        # Ag 極どうし・Nb 極どうしは結合しない（設計上の方向性）。
        return None

    slot_i_free = state.head_free[i] if region_i == "head" else state.tail_free[i]
    slot_j_free = state.head_free[j] if region_j == "head" else state.tail_free[j]
    if not (slot_i_free and slot_j_free):
        return None

    if polar_bond_count < MIN_BONDS_FOR_CONTACT:
        return None
    if rng.random() >= _binding_probability(polar_bond_count, dt_s):
        return None
    return ("polar", region_i, region_j, p1, p2, s, t)


def _record_bond(
    state: SimulationState,
    i: int,
    j: int,
    kind: BondKind,
    region_i: str,
    region_j: str,
    p1: np.ndarray,
    p2: np.ndarray,
    s: float,
    t: float,
) -> bool:
    """いまの接触配置のまま結合を記録する（スナップしない）。

    付着点は現在の接触点、自然長の相対角は現在の相対角にする。したがって斜めに
    触れて結合すれば斜めが自然長になり、平行や一直線を強制しない。以後は
    `_resolve_constraints` の柔軟なジョイントがこの配置をゆるく保つ。
    p1, p2, s, t は `_bond_decision` が計算したものをそのまま受け取る
    （周期境界での j の像を再計算せずに済む）。
    True を返せば結合が成立した。
    """

    # 付着点を現在の接触からローカル座標で求める。
    if region_i == "head":
        axial_i, delta_i = HALF_CYLINDER_UM, 0.0
    elif region_i == "tail":
        axial_i, delta_i = -HALF_CYLINDER_UM, math.pi
    else:
        normal_i = p2 - p1
        axial_i = -HALF_CYLINDER_UM + s * (2.0 * HALF_CYLINDER_UM)
        delta_i = math.atan2(normal_i[1], normal_i[0]) - float(state.angles[i])

    if region_j == "head":
        axial_j, delta_j = HALF_CYLINDER_UM, 0.0
    elif region_j == "tail":
        axial_j, delta_j = -HALF_CYLINDER_UM, math.pi
    else:
        normal_j = p1 - p2
        axial_j = -HALF_CYLINDER_UM + t * (2.0 * HALF_CYLINDER_UM)
        delta_j = math.atan2(normal_j[1], normal_j[0]) - float(state.angles[j])

    rest_angle = _wrap_to_pi(float(state.angles[j] - state.angles[i]))
    state.bonds.append(Bond(i, j, kind, axial_i, delta_i, axial_j, delta_j, rest_angle))

    if kind == "polar":
        if region_i == "head":
            state.head_free[i] = False
        else:
            state.tail_free[i] = False
        if region_j == "head":
            state.head_free[j] = False
        else:
            state.tail_free[j] = False
    return True


def _try_form_bond(
    state: SimulationState,
    i: int,
    j: int,
    config: DynamicsConfig,
    polar_bond_count: float,
    side_bond_count: float,
    dt_s: float,
    rng: np.random.Generator,
) -> bool:
    """i, j の結合を試みる。成立すれば現在の配置のまま記録して True を返す。"""

    if _already_bonded(state, i, j):
        return False
    decision = _bond_decision(
        state, i, j, config, polar_bond_count, side_bond_count, dt_s, rng
    )
    if decision is None:
        return False
    kind, region_i, region_j, p1, p2, s, t = decision
    return _record_bond(state, i, j, kind, region_i, region_j, p1, p2, s, t)


#: 鉛直方向にこの距離（µm）以上離れたクラスタどうしは、物理的にもう接触
#: できないとみなし、結合を試みない（A23）。細胞直径のオーダーに取ってある。
VERTICAL_CONTACT_GAP_UM = 2.0 * CELL_RADIUS_UM


def _apply_settling(state: SimulationState, config: DynamicsConfig, dt_s: float) -> None:
    """クラスタごとの沈降（A23）を1ステップぶん進める。

    `well_depth_um` が None（既定）なら何もしない。液柱の底に達したクラスタは
    それ以上沈まない（`well_depth_um` でクリップする）。
    """

    if not config.settling_enabled:
        return
    for members in connected_components(config.n_cells, state.bonds):
        velocity = _settling_velocity_um_s(len(members))
        idx = np.array(members)
        state.depth_um[idx] += velocity * dt_s
    np.clip(state.depth_um, 0.0, config.well_depth_um, out=state.depth_um)


def _bonded_cell_ids(bonds: list[Bond]) -> set[int]:
    ids: set[int] = set()
    for bond in bonds:
        ids.add(bond.i)
        ids.add(bond.j)
    return ids


def _apply_active_motility(
    state: SimulationState, config: DynamicsConfig, dt_s: float, rng: np.random.Generator
) -> None:
    """結合していない自由な細胞の run-and-tumble 運動を1ステップぶん進める（A25）。

    `swim_speed_um_s = 0`（既定）では何もせず、水平面内の運動は
    `_brownian_step` の受動拡散だけになる（この関数を追加する前と同じ）。

    正の値では、各自由な細胞は run（直進遊泳）と tumble（向きを変える）の
    2相を、それぞれ指数分布する持続時間でポアソン過程的に切り替える
    （結合形成・解離と同じ様式のベルヌーイ試行——Berg & Brown 1972 が自己
    相関関数からこの統計量そのものをポアソン過程だと確認している）。run 中は
    現在の向きに沿って `swim_speed_um_s * dt` だけ直進し、tumble が終わる
    瞬間に向きをランダムに変える。方向転換角度の大きさは、平均
    `tumble_turn_angle_deg`・標準偏差 `tumble_turn_angle_std_deg` に一致する
    ガンマ分布からサンプルする——原著 Fig. 3 が示す「0°に偏らない非一様な
    分布」を、一様分布より現実に近い形で近似するため（A25）。符号（左右
    どちらに曲がるか）は独立に一様ランダムとする。受動的な並進・回転拡散
    （`_brownian_step`）はこれに関係なく常時かかる——遊泳中も熱ゆらぎ自体は
    消えないが、遊泳速度がずっと大きいため無視できるほど小さい寄与になる。

    結合した細胞（クラスタの一部になった細胞）は能動推進を止める——鞭毛の
    推力とクラスタの粘性抵抗・結合拘束との相互作用（能動-弾性的な鎖の座屈や
    回転など）は本モデルの範囲外とする単純化である。
    """

    if not config.motile:
        return
    bonded = _bonded_cell_ids(state.bonds)
    run_switch_prob = 1.0 - math.exp(-dt_s / config.mean_run_time_s)
    tumble_switch_prob = 1.0 - math.exp(-dt_s / config.mean_tumble_time_s)
    mean_turn_deg = config.tumble_turn_angle_deg
    std_turn_deg = config.tumble_turn_angle_std_deg
    gamma_shape = (mean_turn_deg / std_turn_deg) ** 2
    gamma_scale_deg = std_turn_deg**2 / mean_turn_deg

    for k in range(state.positions.shape[0]):
        if k in bonded:
            continue
        if state.running[k]:
            axis = _unit_vector(state.angles[k])
            state.positions[k] += config.swim_speed_um_s * dt_s * axis
            if rng.random() < run_switch_prob:
                state.running[k] = False
        else:
            if rng.random() < tumble_switch_prob:
                turn = math.radians(rng.gamma(gamma_shape, gamma_scale_deg))
                if rng.random() < 0.5:
                    turn = -turn
                state.angles[k] += turn
                state.running[k] = True


def _free_pole_slot_for_axial(state: SimulationState, cell: int, axial: float) -> None:
    """極結合が切れたとき、その極スロットを再び空ける。"""

    if axial > 0.0:
        state.head_free[cell] = True
    else:
        state.tail_free[cell] = True


def _apply_shear_rupture(
    state: SimulationState, n_cells: int, shear_rate_per_s: float
) -> None:
    """震盪下でのせん断荷重が破断力を超える極結合（鎖の背骨）を切る。

    既存の静的な解析（`shear.chain_tension_pn`、仮定 A8 の等荷重分配）を、
    いまできているクラスタにそのまま適用する（A22）。判定は極結合（鎖の背骨）
    だけを対象にする。`chain_tension_pn` は細長い剛体棒（鎖）を仮定した
    slender-body 近似の式であり、束の側面結合が実際にどう荷重を分担するかは
    この式の対象外だからである。文献の結論（震盪培養の典型的な速度では
    接合を壊さない）どおり、通常はここが発火しないはずである
    （`tests/test_dynamics.py` で確認している）。
    """

    if shear_rate_per_s <= 0.0 or not state.bonds:
        return
    components = connected_components(n_cells, state.bonds)
    to_remove: set[int] = set()
    for cells in components:
        if len(cells) < 2:
            continue
        cell_set = set(cells)
        polar_indices = [
            idx
            for idx, bond in enumerate(state.bonds)
            if bond.kind == "polar" and bond.i in cell_set and bond.j in cell_set
        ]
        if not polar_indices:
            continue
        tension_pn = shear.chain_tension_pn(len(cells), shear_rate_per_s)
        load_per_bond_pn = tension_pn / len(polar_indices)
        if load_per_bond_pn > model.RUPTURE_FORCE_PN:
            to_remove.update(polar_indices)

    for idx in sorted(to_remove, reverse=True):
        removed = state.bonds.pop(idx)
        _free_pole_slot_for_axial(state, removed.i, removed.axial_i)
        _free_pole_slot_for_axial(state, removed.j, removed.axial_j)


def _apply_spontaneous_dissociation(
    state: SimulationState, k_off_per_s: float, dt_s: float, rng: np.random.Generator
) -> None:
    """結合の自発的な解離（A15・A24）を1ステップぶん進める。

    `k_off_per_s = 0`（既定）ではまったく呼ばれず、結合は不可逆（Costan et al.
    の best estimate どおり）。0 より大きい値を渡すと、各結合が独立に確率
    `1 - exp(-k_off_per_s * dt)` で解離する（結合形成と対称なベルヌーイ過程）。
    この既定値の生物学的な妥当性については `DynamicsConfig.k_off_per_s` の
    docstring と仮定 A24 を参照——狭い荷重速度範囲での検出限界に基づく主張で
    あり、この結合（anti-EPEA/EPEA）は Costan et al. が他の測定に使っている
    例外的に安定な anti-GFP/GFP ペアとは異なる、より低親和性のペアである。
    """

    if k_off_per_s <= 0.0 or not state.bonds:
        return
    p_break = 1.0 - math.exp(-k_off_per_s * dt_s)
    survivors: list[Bond] = []
    for bond in state.bonds:
        if rng.random() < p_break:
            if bond.kind == "polar":
                _free_pole_slot_for_axial(state, bond.i, bond.axial_i)
                _free_pole_slot_for_axial(state, bond.j, bond.axial_j)
            continue
        survivors.append(bond)
    state.bonds = survivors


def run(config: DynamicsConfig) -> SimulationResult:
    """細胞集団の時間発展シミュレーションを実行する。

    各ステップは次の段からなる。
    1. 各細胞を独立にブラウン運動させる（`_brownian_step`、非剛体・自由排水）。
    2. 運動性モード（`swim_speed_um_s > 0`）では、結合していない自由な細胞に
       run-and-tumble 運動を上乗せする（`_apply_active_motility`、A25）。
    3. 震盪モード（`shear_rate_per_s > 0`）では、平面クエット流による決定論的な
       移流・Jeffery 回転を上乗せする（`_apply_shear_drift`、A20・A21）。
    4. 沈降モード（`well_depth_um` を指定）では、クラスタごとの重力沈降を
       進める（`_apply_settling`、A23）。鉛直方向に離れすぎた細胞どうしは
       結合判定を打ち切る（`VERTICAL_CONTACT_GAP_UM`）。
    5. 結合（柔軟なジョイント）・排除体積・壁の拘束を反復緩和で満たす
       （`_resolve_constraints`）。震盪モードでは流れ方向（x）を周期境界にする。
    6. 接触した細胞間で結合判定を行い、成立すれば現在の接触配置のまま記録する。
       既定（`k_off_per_s = 0`）では一度できた結合は解離しない（A15、
       k_off ≈ 0）——ただし震盪モードでは、せん断荷重が破断力を超えた極結合
       だけは切れうる（A22）。`k_off_per_s` を正にすると、すべての結合が
       自発的にも確率的に解離しうる（A24、`_apply_spontaneous_dissociation`）。

    震盪モードと沈降モードは独立に有効化できるが、Glass & Riedel-Kruse (2018)
    の実際の凝集アッセイに最も近いのは**沈降モードのみ**（震盪なし）である。
    実際のアッセイは、育成培養こそ震盪するが、凝集そのものは短時間の
    ボルテックス混合のあと静置で進む。長時間（24時間）の沈降モードでは、
    k_off ≈ 0 の生物学的な妥当性に限界があるため（A24）、`k_off_per_s` に
    非零の値を指定することを検討するとよい。
    """

    rng = np.random.default_rng(config.seed)

    prediction = model.predict(
        config.n_total,
        config.f_pole,
        config.polar_split,
        config.tip_width_um,
        distinct_adhesins_per_pole=config.distinct_adhesins_per_pole,
    )
    polar_bond_count = prediction.bonds_heavy_pole
    side_bond_count = prediction.bonds_side_by_side

    state = _initial_state(config, rng)

    dt_s = min(_safe_dt_for_shear(config), _safe_dt_for_motility(config))
    n_steps = int(round(config.total_time_s / dt_s))
    broad_cutoff = CELL_LENGTH_UM + 0.5
    periodic_x = config.periodic_x
    relaxation_iters = N_RELAXATION_ITERS

    frames = [
        Frame(
            0.0, state.positions.copy(), state.angles.copy(), list(state.bonds),
            state.depth_um.copy(),
        )
    ]

    for step in range(1, n_steps + 1):
        _brownian_step(state, dt_s, rng)
        _apply_active_motility(state, config, dt_s, rng)
        _apply_shear_drift(state, config, dt_s)
        if periodic_x:
            _wrap_periodic_x(state, config.box_size_um)
        _apply_settling(state, config, dt_s)

        candidate_pairs = _broad_phase_pairs(
            state, broad_cutoff, config.box_size_um, periodic_x
        )
        _resolve_constraints(
            state, candidate_pairs, config.box_size_um, relaxation_iters, periodic_x
        )

        for i, j in candidate_pairs:
            if _already_bonded(state, i, j):
                continue
            decision = _bond_decision(
                state, i, j, config, polar_bond_count, side_bond_count, dt_s, rng
            )
            if decision is None:
                continue
            kind, region_i, region_j, p1, p2, s, t = decision
            _record_bond(state, i, j, kind, region_i, region_j, p1, p2, s, t)

        _apply_shear_rupture(state, config.n_cells, config.shear_rate_per_s)
        _apply_spontaneous_dissociation(state, config.k_off_per_s, dt_s, rng)

        if step % config.frame_stride == 0 or step == n_steps:
            frames.append(
                Frame(
                    step * dt_s,
                    state.positions.copy(),
                    state.angles.copy(),
                    list(state.bonds),
                    state.depth_um.copy(),
                )
            )

    return SimulationResult(config, prediction, polar_bond_count, side_bond_count, frames)


def connected_components(n_cells: int, bonds: list[Bond]) -> list[list[int]]:
    """結合グラフの連結成分を返す（Union-Find）。"""

    parent = list(range(n_cells))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for bond in bonds:
        union(bond.i, bond.j)

    groups: dict[int, list[int]] = {}
    for k in range(n_cells):
        groups.setdefault(find(k), []).append(k)
    return list(groups.values())


def classify_clusters(
    n_cells: int, bonds: list[Bond], distinct_adhesins_per_pole: bool
) -> list[ClusterInfo]:
    """各連結成分を単量体・鎖・束に分類する。"""

    components = connected_components(n_cells, bonds)
    infos = []
    for cells in components:
        if len(cells) == 1:
            infos.append(ClusterInfo(cells, "monomer", directed=False))
            continue
        cell_set = set(cells)
        in_cluster = [b for b in bonds if b.i in cell_set and b.j in cell_set]
        has_lateral = any(b.kind == "lateral" for b in in_cluster)
        kind = "bundle" if has_lateral else "chain"
        infos.append(ClusterInfo(cells, kind, directed=distinct_adhesins_per_pole))
    return infos


def supernatant_fraction(frame: Frame, well_depth_um: float, top_fraction: float) -> float:
    """上層（`top_fraction` 割合）に残っている細胞の割合。

    Glass & Riedel-Kruse (2018) の実際の凝集アッセイが測っている量
    （上層の光学密度）に対応する。大きい凝集塊ほど速く沈むので、この値が
    下がっていくことが「凝集が進んだ」ことを意味する。
    """

    if frame.positions.shape[0] == 0:
        return 1.0
    threshold = well_depth_um * top_fraction
    return float(np.mean(frame.depth_um <= threshold))


__all__ = [
    "D_TRANS0_UM2_S",
    "D_ROT0_RAD2_S",
    "Bond",
    "ClusterInfo",
    "DynamicsConfig",
    "Frame",
    "SimulationResult",
    "SimulationState",
    "attachment_point",
    "classify_clusters",
    "connected_components",
    "run",
    "supernatant_fraction",
]
