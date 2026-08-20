"""荷重下での結合寿命を Bell モデルで見積もる。

本体のモデル（`adhesion.model`）は時間を持たず、接触した2細胞がどの配置を取るか
という相対的な有利不利だけを扱う（仮定 A7）。このモジュールはそこに含めない
補助的な見積もりで、「できた接合がどれくらい持つか」の桁を確かめるために使う。

重要な注意
----------
解離速度定数と障壁位置は anti-GFP / GFP ペアの測定値であり、破断力 16.1 pN を
測ったペア（anti-EPEA / EPEA）とは異なる。したがってここで得られる寿命は
**ナノボディ–抗原結合としての桁の目安**であって、本系の実測値ではない。
異なるペアの値を混ぜていることを承知のうえで使うこと。

Bell モデル
-----------
一定荷重 F のもとでの結合寿命は

    tau(F) = (1 / k_off) * exp(-F * dx / kT)

障壁位置 dx だけ荷重方向に伸びると解離障壁が F*dx だけ下がる、という描像による。
荷重が kT/dx（本系では約 2.6 pN）を大きく下回るあいだ、寿命はほとんど変わらない。
"""

from __future__ import annotations

import math

from . import literature as lit
from .assumptions import uses

KINETICS_ASSUMPTIONS = uses("A5", "A8")

#: 室温での熱エネルギー（pN·nm）。
THERMAL_ENERGY_PN_NM = 4.11

OFF_RATE_PER_S = float(lit.BOND_OFF_RATE)
BARRIER_POSITION_NM = float(lit.BOND_BARRIER_POSITION)

#: 荷重の効きが顕著になる力のスケール（pN）。kT / dx。
FORCE_SCALE_PN = THERMAL_ENERGY_PN_NM / BARRIER_POSITION_NM


def bond_lifetime_s(force_pn: float = 0.0) -> float:
    """一定荷重 F のもとでの結合1本の平均寿命（秒）。"""

    if force_pn < 0:
        raise ValueError("荷重は非負である必要があります")
    return math.exp(-force_pn * BARRIER_POSITION_NM / THERMAL_ENERGY_PN_NM) / OFF_RATE_PER_S


def lifetime_reduction(force_pn: float) -> float:
    """無荷重に対する寿命の比 tau(F)/tau(0)。1 に近ければ荷重は効いていない。"""

    return bond_lifetime_s(force_pn) / bond_lifetime_s(0.0)


def rupture_force_pn(loading_rate_pn_per_s: float) -> float:
    """動的力分光で観測される最頻破断力（pN）。Bell–Evans。

        F* = (kT / dx) * ln(r * dx / (kT * k_off))

    破断力が荷重速度に依存することを示すための関数であり、モデル本体では使わない。
    荷重速度が低いほど破断力は小さくなる。
    """

    if loading_rate_pn_per_s <= 0:
        raise ValueError("荷重速度は正である必要があります")
    argument = (
        loading_rate_pn_per_s * BARRIER_POSITION_NM
        / (THERMAL_ENERGY_PN_NM * OFF_RATE_PER_S)
    )
    if argument <= 1.0:
        return 0.0
    return FORCE_SCALE_PN * math.log(argument)


__all__ = [
    "BARRIER_POSITION_NM",
    "FORCE_SCALE_PN",
    "OFF_RATE_PER_S",
    "THERMAL_ENERGY_PN_NM",
    "bond_lifetime_s",
    "lifetime_reduction",
    "rupture_force_pn",
]
