"""荷重下での結合寿命を Bell モデルで見積もる。

本体のモデル（`adhesion.model`）は時間を持たず、接触した2細胞がどの配置を取るか
という相対的な有利不利だけを扱う（仮定 A7）。このモジュールはそこに含めない
補助的な見積もりで、「できた接合がどれくらい持つか」の桁を確かめるために使う。

重要な注意
----------
障壁位置 Δx は Klamecka et al. (2015) が AFM で自ら測った値であり、全文照合済み。
ただし **解離速度定数 k_off はこのモジュールに固定値として持たない**。

理由は出典にある。Klamecka らは Bell–Evans 回帰の際に k_off を自分で測らず、
Kirchhofer et al. (2010) の文献値 1.45e-4 /s を固定して使っている。したがって
k_off を Klamecka の測定値として引くのは誤りであり、Kirchhofer の全文を確認する
までは使えない。呼び出し側が明示的に渡す形にしてある。

さらに Costan et al. (2024) は本系のアドヘシンについて「k_off はゼロとみなすのが
最良の推定であり、アドヘシンは自発的には外れない」と結論している。したがって
以下の寿命の見積もりは、別のペア（anti-GFP / GFP）を借りた場合に何が言えるかを
示すものであって、本系の予測ではない。

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

BARRIER_POSITION_NM = float(lit.BOND_BARRIER_POSITION)

#: 参考値。Klamecka et al. が回帰に用いた Kirchhofer et al. (2010) の文献値。
#: **全文未確認のため既定値にはしない。** 使う場合は呼び出し側が明示的に渡し、
#: 出典を自分で確認すること。
UNVERIFIED_OFF_RATE_PER_S = 1.45e-4

#: 荷重の効きが顕著になる力のスケール（pN）。kT / dx。
FORCE_SCALE_PN = THERMAL_ENERGY_PN_NM / BARRIER_POSITION_NM


def bond_lifetime_s(force_pn: float, off_rate_per_s: float) -> float:
    """一定荷重 F のもとでの結合1本の平均寿命（秒）。

    off_rate_per_s は呼び出し側が明示的に渡す。本系については Costan et al. が
    ゼロとみなすのが最良の推定としているため、既定値は置かない。
    別のペアの値を借りる場合は `UNVERIFIED_OFF_RATE_PER_S` を出発点にできるが、
    その出典（Kirchhofer et al. 2010）を自分で確認すること。
    """

    if force_pn < 0:
        raise ValueError("荷重は非負である必要があります")
    if off_rate_per_s <= 0:
        raise ValueError("解離速度定数は正である必要があります")
    return math.exp(-force_pn * BARRIER_POSITION_NM / THERMAL_ENERGY_PN_NM) / off_rate_per_s


def lifetime_reduction(force_pn: float) -> float:
    """無荷重に対する寿命の比 tau(F)/tau(0)。1 に近ければ荷重は効いていない。

    k_off は比を取ると打ち消えるので、この量だけは未確認の値に依存しない。
    """

    if force_pn < 0:
        raise ValueError("荷重は非負である必要があります")
    return math.exp(-force_pn * BARRIER_POSITION_NM / THERMAL_ENERGY_PN_NM)


def rupture_force_pn(
    loading_rate_pn_per_s: float, off_rate_per_s: float
) -> float:
    """動的力分光で観測される最頻破断力（pN）。Bell–Evans。

        F* = (kT / dx) * ln(r * dx / (kT * k_off))

    破断力が荷重速度に依存することを示すための関数であり、モデル本体では使わない。
    荷重速度が低いほど破断力は小さくなる。
    """

    if loading_rate_pn_per_s <= 0:
        raise ValueError("荷重速度は正である必要があります")
    argument = (
        loading_rate_pn_per_s * BARRIER_POSITION_NM
        / (THERMAL_ENERGY_PN_NM * off_rate_per_s)
    )
    if argument <= 1.0:
        return 0.0
    return FORCE_SCALE_PN * math.log(argument)


__all__ = [
    "BARRIER_POSITION_NM",
    "FORCE_SCALE_PN",
    "off_rate_per_s",
    "THERMAL_ENERGY_PN_NM",
    "bond_lifetime_s",
    "lifetime_reduction",
    "rupture_force_pn",
]
