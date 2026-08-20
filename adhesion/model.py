"""発現量と局在から、2細胞が作る接着とその形態を予測する。

自由パラメータは2つだけである。

    n_total   細胞あたりのアドヘシン数。aTc 誘導で実験的に制御できる
    局在      f_pole / tip_width_um / polar_split の3つ組

結合の数は自由パラメータではない。幾何と n_total から決まる従属量であり、
1結合あたりの破断力は実測値に固定されるため、接着の強さの軸には自由パラメータが
残らない。

依存する仮定: A1–A8, A10。詳細は `adhesion.assumptions` を参照。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from . import literature as lit
from .assumptions import uses
from .geometry import (
    POSE_END_TO_END,
    POSE_SIDE_BY_SIDE,
    SurfaceQuadrature,
    adhesin_counts,
    contactable_adhesins,
    polar_area_fraction,
    spherocylinder_surface,
)

CELL_LENGTH_UM = float(lit.CELL_LENGTH)
CELL_RADIUS_UM = float(lit.CELL_RADIUS)
BRIDGE_UM = float(lit.BRIDGE_LENGTH) / 1000.0
RUPTURE_FORCE_PN = float(lit.RUPTURE_FORCE)

#: プラスミドコピー数の非対称に対する補正係数（A4）。抗原側が律速になる。
FUNCTIONAL_FRACTION = float(lit.COPY_NUMBER_CORRECTION)

#: 局在を操作しない場合の分布（A10）。Costan et al. (2024) の実測は一様。
UNIFORM_F_POLE = polar_area_fraction(CELL_LENGTH_UM, CELL_RADIUS_UM)

#: 結合が1本も成立しないとみなす下限。
MIN_BONDS_FOR_CONTACT = 1.0

MODEL_ASSUMPTIONS = uses("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A10")

_SURFACE_CACHE: dict[int, SurfaceQuadrature] = {}


def cell_surface(n_points: int = 20000) -> SurfaceQuadrature:
    """細胞表面の求積を返す（同じ解像度なら使い回す）。"""
    if n_points not in _SURFACE_CACHE:
        _SURFACE_CACHE[n_points] = spherocylinder_surface(
            CELL_LENGTH_UM, CELL_RADIUS_UM, n_points
        )
    return _SURFACE_CACHE[n_points]


class Morphology(Enum):
    """接着した細胞が作る多細胞構造。"""

    NO_ADHESION = "接着しない"
    BUNDLE = "束"
    DIMER = "二量体どまり"
    CHAIN = "鎖"

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True)
class AdhesionPrediction:
    """1つの発現量・局在条件についての予測。"""

    n_total: float
    f_pole: float
    polar_split: float
    tip_width_um: float | None

    bonds_side_by_side: float
    bonds_heavy_pole: float
    bonds_light_pole: float
    peak_density_per_um2: float

    @property
    def morphology(self) -> Morphology:
        """できあがる多細胞構造。

        鎖が伸びるには1細胞が2つの接続を持つ必要があり、その律速は分子が少ない
        側の極である。片極が空なら二量体で組み立てが止まる。
        """

        if max(self.bonds_side_by_side, self.bonds_heavy_pole) < MIN_BONDS_FOR_CONTACT:
            return Morphology.NO_ADHESION
        if self.bonds_side_by_side >= self.bonds_heavy_pole:
            return Morphology.BUNDLE
        if self.bonds_light_pole < MIN_BONDS_FOR_CONTACT:
            return Morphology.DIMER
        return Morphology.CHAIN

    @property
    def junction_force_pn(self) -> float:
        """優勢な接触様式での接合強度。結合数 × 単結合の破断力（A8）。"""
        bonds = max(self.bonds_side_by_side, self.bonds_heavy_pole)
        return bonds * RUPTURE_FORCE_PN

    @property
    def bonds_per_junction(self) -> float:
        return max(self.bonds_side_by_side, self.bonds_heavy_pole)


def _pole_contact_bonds(
    surface: SurfaceQuadrature,
    molecules_on_pole: float,
    tip_width_um: float | None,
) -> float:
    """指定した数の分子を載せた極どうしが接触したときの結合数。

    2細胞は自由に向きを変えられるので、互いに同じ側の極（濃い極どうし、または
    薄い極どうし）を向かい合わせられる。そこで両極に同数を置いた分布を作り、
    先端-先端の接触を評価する。
    """

    if molecules_on_pole <= 0:
        return 0.0
    counts = adhesin_counts(
        surface,
        n_total=2.0 * molecules_on_pole,
        f_pole=1.0,
        polar_split=0.5,
        tip_width_um=tip_width_um,
    )
    accessible = contactable_adhesins(
        surface, counts, POSE_END_TO_END, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
    )
    return FUNCTIONAL_FRACTION * accessible


def _peak_density(
    n_total: float,
    f_pole: float,
    polar_split: float,
    tip_width_um: float | None,
) -> float:
    """極先端でのアドヘシン面密度（分子/µm²）。

    外膜のタンパク収容量（全種合計で 10⁵ /µm² のオーダー）と比較して、
    局在が物理的に成立しうるかを見るために使う。
    """

    on_heavy_pole = n_total * f_pole * max(polar_split, 1.0 - polar_split)
    if tip_width_um is None:
        return on_heavy_pole / (2 * math.pi * CELL_RADIUS_UM**2)
    return on_heavy_pole / (2 * math.pi * tip_width_um**2)


def predict(
    n_total: float = float(lit.ADHESINS_AT_MID_INDUCTION),
    f_pole: float = UNIFORM_F_POLE,
    polar_split: float = 0.5,
    tip_width_um: float | None = None,
    n_points: int = 20000,
) -> AdhesionPrediction:
    """接着の予測を返す。

    n_total       細胞あたりのアドヘシン数
    f_pole        両極に置く割合。既定は一様分布（A10）
    polar_split   極分子の左右分配。0.5 で対称、1.0 で単極
    tip_width_um  極先端へのガウス集中幅（A6）。None なら半球に一様
    """

    if not 0.5 <= polar_split <= 1.0:
        raise ValueError("polar_split は 0.5（対称）から 1.0（単極）で指定します")

    surface = cell_surface(n_points)

    counts = adhesin_counts(surface, n_total, f_pole, polar_split, tip_width_um)
    side = FUNCTIONAL_FRACTION * contactable_adhesins(
        surface, counts, POSE_SIDE_BY_SIDE, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
    )

    on_heavy = n_total * f_pole * polar_split
    on_light = n_total * f_pole * (1.0 - polar_split)

    return AdhesionPrediction(
        n_total=n_total,
        f_pole=f_pole,
        polar_split=polar_split,
        tip_width_um=tip_width_um,
        bonds_side_by_side=side,
        bonds_heavy_pole=_pole_contact_bonds(surface, on_heavy, tip_width_um),
        bonds_light_pole=_pole_contact_bonds(surface, on_light, tip_width_um),
        peak_density_per_um2=_peak_density(n_total, f_pole, polar_split, tip_width_um),
    )


def morphology_threshold(
    tip_width_um: float | None,
    polar_split: float = 0.5,
    n_points: int = 20000,
    tolerance: float = 1e-4,
) -> float | None:
    """鎖が束に打ち勝つのに必要な f_pole を二分探索で求める。

    この閾値は発現量に依存しない。極接触も側面接触も n_total に比例するため、
    比を取ると n_total が消えるからである。すなわち発現量は接合の強さを決め、
    局在は形を決める。両者は干渉しない。

    到達可能な f_pole では鎖にならない場合は None を返す。
    """

    probe_n = float(lit.ADHESINS_AT_MID_INDUCTION)

    def pole_wins(f_pole: float) -> bool:
        result = predict(probe_n, f_pole, polar_split, tip_width_um, n_points)
        return result.bonds_heavy_pole > result.bonds_side_by_side

    low, high = UNIFORM_F_POLE, 1.0 - 1e-9
    if not pole_wins(high):
        return None
    if pole_wins(low):
        return UNIFORM_F_POLE
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if pole_wins(middle):
            high = middle
        else:
            low = middle
    return high
