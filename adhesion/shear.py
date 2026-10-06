"""震盪培養で細胞にかかる流体力を見積もる。

震盪の役割は2つある。細胞どうしを出会わせること（衝突頻度）と、
できた接合を引き剥がそうとすること（せん断荷重）である。本モデルの範囲では、
前者が支配的で後者は効かない。

依存する仮定: A9（層流せん断）、A8（荷重の均等分配）。
"""

from __future__ import annotations

import math

from . import literature as lit
from .assumptions import uses

SHEAR_ASSUMPTIONS = uses("A8", "A9")

#: 水の動粘度（m²/s）。20-37 度で 1e-6 のオーダー。
KINEMATIC_VISCOSITY_M2_S = 1e-6

#: 培地の粘度（Pa·s）。水の標準物性値で近似する（特定の文献に由来する値ではない）。
DYNAMIC_VISCOSITY_PA_S = 1e-3

#: 室温での熱エネルギー（J）。
THERMAL_ENERGY_J = 4.11e-21

_CELL_LENGTH_M = float(lit.CELL_LENGTH) * 1e-6
_CELL_DIAMETER_M = 2 * float(lit.CELL_RADIUS) * 1e-6


def shear_rate(dissipation_w_per_kg: float) -> float:
    """乱流の散逸率からコルモゴロフ・スケールのせん断速度（1/s）を求める。

    γ̇ = √(ε/ν)。フラスコ震盪の典型値 ε = 0.01–1 W/kg で γ̇ = 100–1000 /s。
    """

    if dissipation_w_per_kg < 0:
        raise ValueError("散逸率は非負である必要があります")
    return math.sqrt(dissipation_w_per_kg / KINEMATIC_VISCOSITY_M2_S)


def kolmogorov_length_um(dissipation_w_per_kg: float) -> float:
    """最小渦のスケール（µm）。η = (ν³/ε)^(1/4)。

    これが細胞より十分大きければ、細胞は渦の内側にいて一様なせん断を受ける（A9）。
    ε = 0.01–1 W/kg で 32–100 µm となり、細胞（1.8 µm）よりはるかに大きい。
    """

    if dissipation_w_per_kg <= 0:
        raise ValueError("散逸率は正である必要があります")
    length_m = (KINEMATIC_VISCOSITY_M2_S**3 / dissipation_w_per_kg) ** 0.25
    return length_m * 1e6


def chain_tension_pn(
    n_cells: int,
    shear_rate_per_s: float,
    cell_length_um: float | None = None,
    cell_radius_um: float | None = None,
) -> float:
    """せん断中の細胞鎖の中央にかかる張力（pN）。

    cell_length_um / cell_radius_um を省略すると文献値の細胞寸法を使う。

    伸長方位に置かれた細長い剛体棒に対する slender-body 近似:

        T = π η γ̇ L² / (8 ln(L/d))

    鎖長の2乗で効くため、長い鎖ほど不利になる。
    """

    if n_cells < 2:
        raise ValueError("鎖は2細胞以上である必要があります")
    if shear_rate_per_s < 0:
        raise ValueError("せん断速度は非負である必要があります")

    cell_length_m = _CELL_LENGTH_M if cell_length_um is None else cell_length_um * 1e-6
    cell_diameter_m = (
        _CELL_DIAMETER_M if cell_radius_um is None else 2.0 * cell_radius_um * 1e-6
    )
    chain_length_m = n_cells * cell_length_m
    slenderness = chain_length_m / cell_diameter_m
    tension_n = (
        math.pi
        * DYNAMIC_VISCOSITY_PA_S
        * shear_rate_per_s
        * chain_length_m**2
        / (8 * math.log(slenderness))
    )
    return tension_n * 1e12


def load_per_bond_pn(n_cells: int, shear_rate_per_s: float, n_bonds: float) -> float:
    """1結合あたりの荷重（pN）。荷重は均等に分配されるとする（A8）。"""

    if n_bonds <= 0:
        return math.inf
    return chain_tension_pn(n_cells, shear_rate_per_s) / n_bonds


def contact_duration_s(shear_rate_per_s: float) -> float:
    """せん断中の2細胞が擦れ違うまでの接触持続時間（s）。

    せん断場では相対速度が勾配で決まるため、接触は概ね 1/γ̇ の時間で終わる。
    これを結合形成にかかる時間と比べると、その衝突で結合ができるかが分かる。
    """

    if shear_rate_per_s <= 0:
        return math.inf
    return 1.0 / shear_rate_per_s


def reorientation_during_contact_deg(
    shear_rate_per_s: float,
    rotational_diffusion_rad2_s: float = 0.889,
) -> float:
    """接触している間に熱運動で向きを変えられる角度（度）。

    √(2 D_rot t) を接触時間 t = 1/γ̇ について評価する。回転拡散係数の既定値は
    細胞と同体積の球に対する kT/(8πηa³) で、桁の見積もりである。
    """

    return math.degrees(
        math.sqrt(2 * rotational_diffusion_rad2_s * contact_duration_s(shear_rate_per_s))
    )


def thermal_force_scale_pn(barrier_position_nm: float = 1.59) -> float:
    """熱ゆらぎが結合に与える力のスケール kT/Δx（pN）。

    結合1本にかかる力がこれより十分小さければ、その力は熱ゆらぎに埋もれる。
    既定の障壁位置は anti-GFP / GFP ペアの測定値（Klamecka et al. 2015）。
    """

    return THERMAL_ENERGY_J / (barrier_position_nm * 1e-9) * 1e12


def collision_rate_per_s(
    shear_rate_per_s: float,
    cell_concentration_per_ml: float,
) -> float:
    """1細胞が単位時間に他細胞と衝突する回数（1/s）。

    Smoluchowski のせん断（orthokinetic）カーネル k = (4/3) γ̇ (a₁+a₂)³ を使う。

    このカーネルは本来球どうしの衝突に対するもので、a₁+a₂ は2球の半径の和である。
    桿菌に厳密に対応する量はないため、ここでは長軸半分（0.9 µm）と短軸半径
    （0.4 µm）の和を衝突断面の代表長さとして採る。棒状粒子の向きによって実際の
    衝突断面は変わるので、この値は桁の見積もりであって精密値ではない。
    短軸側だけを採れば約 0.5 倍、長軸側だけなら約 1.7 倍に振れる。
    """

    collision_radius_m = (0.5 * _CELL_LENGTH_M) + (0.5 * _CELL_DIAMETER_M)
    kernel_m3_per_s = (4.0 / 3.0) * shear_rate_per_s * collision_radius_m**3
    concentration_per_m3 = cell_concentration_per_ml * 1e6
    return kernel_m3_per_s * concentration_per_m3


def brownian_collision_rate_per_s(cell_concentration_per_ml: float) -> float:
    """ブラウン運動（perikinetic）による衝突頻度（1/s）。比較の基準に使う。

    k = 8 k_B T / (3 η)。せん断項との比を取ると、震盪が衝突をどれだけ加速するかが
    分かる。
    """

    kernel_m3_per_s = 8 * THERMAL_ENERGY_J / (3 * DYNAMIC_VISCOSITY_PA_S)
    concentration_per_m3 = cell_concentration_per_ml * 1e6
    return kernel_m3_per_s * concentration_per_m3


def shear_to_brownian_ratio(shear_rate_per_s: float) -> float:
    """せん断による衝突がブラウン運動の何倍かを返す（濃度によらない）。"""

    if shear_rate_per_s <= 0:
        return 0.0
    return collision_rate_per_s(shear_rate_per_s, 1.0) / brownian_collision_rate_per_s(
        1.0
    )
