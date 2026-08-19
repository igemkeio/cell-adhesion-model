"""2細胞が接触したときの幾何を扱う。

細胞はスフェロシリンダー（円柱の両端に半球）として表し、表面を等面積の求積点で
覆う。ある配置において「相手の表面まで架橋長以内、かつ法線が向き合っている」
求積点を数えることで、結合しうるアドヘシンの数を得る。

依存する仮定: A1（剛体）、A2（架橋長は2細胞間の総額）、A3（法線の向き）、
A5（1分子1結合）。詳細は `adhesion.assumptions` を参照。

解析式との対応
--------------
2つの球面が接するとき、接点から横方向に距離 r 離れた位置での隙間は h ≈ r²/R。
これが架橋長 ℓ に等しくなるのが r_c = √(ℓR) で、ここが接触パッチの縁になる。
したがって

    極-極の接触パッチ面積 = π ℓ R           （ℓ に比例）
    側面-側面の接触帯面積 = 2 √(ℓR) (l - 2R)  （√ℓ に比例）

側面接触のほうが1桁大きい。局在させない限り細胞が束になるのはこのためである。
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

Region = Literal[-1, 0, 1]  # -1: 左の半球, 0: 円柱側面, +1: 右の半球

# A3。法線がこの角度以内で向き合っていれば結合しうるとする。
FACING_ANGLE_DEG = 60.0


@dataclass(frozen=True)
class SurfaceQuadrature:
    """スフェロシリンダー表面の等面積求積。

    points        表面上の点（細胞中心を原点とするローカル座標、µm）
    normals       各点での外向き単位法線
    areas         各点が代表する面積（µm²）。合計は全表面積になる
    regions       -1: 左半球, 0: 円柱側面, +1: 右半球
    tip_distances 属する極の先端からの測地距離（µm）。側面では 0
    """

    points: np.ndarray
    normals: np.ndarray
    areas: np.ndarray
    regions: np.ndarray
    tip_distances: np.ndarray

    @property
    def total_area(self) -> float:
        return float(self.areas.sum())


def _validate(length_um: float, radius_um: float) -> float:
    """半円柱部の長さを返す。形状として成立しない寸法は弾く。"""
    if radius_um <= 0:
        raise ValueError("radius_um は正である必要があります")
    if length_um < 2 * radius_um:
        raise ValueError("length_um は 2 * radius_um 以上である必要があります")
    return 0.5 * (length_um - 2 * radius_um)


def total_surface_area(length_um: float, radius_um: float) -> float:
    """スフェロシリンダーの全表面積（µm²）。"""
    half_cylinder = _validate(length_um, radius_um)
    return 4 * math.pi * radius_um**2 + 4 * math.pi * radius_um * half_cylinder


def polar_area_fraction(length_um: float, radius_um: float) -> float:
    """両端の半球が全表面積に占める割合。

    局在させない一様分布は f_pole がこの値に等しい状態にあたる。
    l = 1.80 µm, r = 0.40 µm では約 0.444。
    """
    half_cylinder = _validate(length_um, radius_um)
    cap = 4 * math.pi * radius_um**2
    side = 4 * math.pi * radius_um * half_cylinder
    return cap / (cap + side)


def contact_patch_area(radius_um: float, bridge_um: float) -> float:
    """極-極接触で架橋しうる領域の面積（µm²）。π ℓ R。"""
    return math.pi * bridge_um * radius_um


def contact_strip_area(length_um: float, radius_um: float, bridge_um: float) -> float:
    """側面-側面接触で架橋しうる領域の面積（µm²）。2 √(ℓR) (l - 2R)。"""
    half_width = math.sqrt(bridge_um * radius_um)
    return 2 * half_width * (length_um - 2 * radius_um)


def _golden_sequence(count: int) -> np.ndarray:
    """黄金比による低食い違い列。方位角を偏りなく散らすために使う。"""
    return np.mod(np.arange(count, dtype=float) * ((math.sqrt(5.0) - 1.0) / 2.0), 1.0)


def spherocylinder_surface(
    length_um: float,
    radius_um: float,
    n_points: int = 20000,
) -> SurfaceQuadrature:
    """スフェロシリンダー表面を等面積の求積点で覆う。

    円柱側面は軸方向に等間隔、方位角は黄金比列で散らす。半球部は cos(極角) を
    等間隔に取ることで等面積になる（アルキメデスの定理）。
    """

    half_cylinder = _validate(length_um, radius_um)
    cap_area = 4 * math.pi * radius_um**2
    side_area = 4 * math.pi * radius_um * half_cylinder
    total = cap_area + side_area

    n_side = max(200, int(round(n_points * side_area / total)))
    n_cap_each = max(200, int(round(n_points * cap_area / total)) // 2)

    # --- 円柱側面 ---
    index = np.arange(n_side, dtype=float)
    axial = -half_cylinder + 2 * half_cylinder * (index + 0.5) / n_side
    azimuth = 2 * math.pi * _golden_sequence(n_side)
    side_normals = np.column_stack(
        (np.zeros(n_side), np.cos(azimuth), np.sin(azimuth))
    )
    side_points = np.column_stack(
        (axial, radius_um * np.cos(azimuth), radius_um * np.sin(azimuth))
    )

    # --- 両端の半球。cos(theta) を等間隔にすると等面積になる ---
    cap_index = np.arange(n_cap_each, dtype=float)
    cos_theta = (cap_index + 0.5) / n_cap_each
    cap_azimuth = 2 * math.pi * _golden_sequence(n_cap_each)
    sin_theta = np.sqrt(1.0 - cos_theta**2)

    right_normals = np.column_stack(
        (cos_theta, sin_theta * np.cos(cap_azimuth), sin_theta * np.sin(cap_azimuth))
    )
    left_normals = right_normals * np.array([-1.0, 1.0, 1.0])
    right_points = np.array([half_cylinder, 0.0, 0.0]) + radius_um * right_normals
    left_points = np.array([-half_cylinder, 0.0, 0.0]) + radius_um * left_normals

    points = np.vstack((side_points, left_points, right_points))
    normals = np.vstack((side_normals, left_normals, right_normals))
    areas = np.concatenate(
        (
            np.full(n_side, side_area / n_side),
            np.full(n_cap_each, 0.5 * cap_area / n_cap_each),
            np.full(n_cap_each, 0.5 * cap_area / n_cap_each),
        )
    )
    regions = np.concatenate(
        (
            np.zeros(n_side, dtype=np.int8),
            np.full(n_cap_each, -1, dtype=np.int8),
            np.full(n_cap_each, 1, dtype=np.int8),
        )
    )
    # 極先端からの測地距離。半球上では R * theta。
    tip = radius_um * np.arccos(np.clip(cos_theta, -1.0, 1.0))
    tip_distances = np.concatenate((np.zeros(n_side), tip, tip))

    return SurfaceQuadrature(points, normals, areas, regions, tip_distances)


def adhesin_counts(
    surface: SurfaceQuadrature,
    n_total: float,
    f_pole: float,
    polar_split: float = 0.5,
    tip_width_um: float | None = None,
) -> np.ndarray:
    """各求積点が代表するアドヘシン分子数を返す。

    f_pole        両極（半球部）に置く分子の割合。一様分布は polar_area_fraction に等しい
    polar_split   極分子のうち左の極に置く割合。0.5 で対称、1.0 で単極
    tip_width_um  極先端からの測地距離に対するガウス幅（A6）。None なら半球に一様

    合計は必ず n_total に一致する。
    """

    if not 0.0 <= f_pole <= 1.0:
        raise ValueError("f_pole は 0 以上 1 以下である必要があります")
    if not 0.0 <= polar_split <= 1.0:
        raise ValueError("polar_split は 0 以上 1 以下である必要があります")
    if n_total < 0:
        raise ValueError("n_total は非負である必要があります")
    if tip_width_um is not None and tip_width_um <= 0:
        raise ValueError("tip_width_um は正である必要があります")

    counts = np.zeros(len(surface.points), dtype=float)
    targets = {
        -1: n_total * f_pole * polar_split,
        0: n_total * (1.0 - f_pole),
        1: n_total * f_pole * (1.0 - polar_split),
    }
    for region, target in targets.items():
        mask = surface.regions == region
        if target <= 0 or not np.any(mask):
            continue
        weight = surface.areas[mask].astype(float)
        if region != 0 and tip_width_um is not None:
            weight = weight * np.exp(
                -0.5 * (surface.tip_distances[mask] / tip_width_um) ** 2
            )
        counts[mask] = target * weight / weight.sum()
    return counts


def _gap_and_normal_to_partner(
    points: np.ndarray,
    partner_centre: np.ndarray,
    partner_axis: np.ndarray,
    half_cylinder: float,
    radius_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """各点から相手細胞の表面までの最短距離と、その位置での相手の外向き法線。"""

    relative = points - partner_centre
    axial = np.clip(relative @ partner_axis, -half_cylinder, half_cylinder)
    nearest_axis_point = partner_centre + axial[:, None] * partner_axis
    radial = points - nearest_axis_point
    distance = np.linalg.norm(radial, axis=1)
    outward = radial / np.maximum(distance, 1e-12)[:, None]
    return distance - radius_um, outward


@dataclass(frozen=True)
class Pose:
    """2細胞の相対配置。

    relative_angle_deg  2細胞の長軸のなす角。0 なら平行
    contact_normal_deg  細胞1から見た接触方向。0 なら長軸方向（先端-先端）、
                        90 なら長軸に垂直（横並び）
    """

    relative_angle_deg: float
    contact_normal_deg: float

    @property
    def description(self) -> str:
        """配置の呼び名。長軸のなす角と接触方向の両方で決まる。

        接触方向だけでは平行な横並びと直交を区別できないので、両方を見る。
        """

        axes_angle = self.relative_angle_deg % 180.0
        normal_angle = self.contact_normal_deg % 180.0
        axes_parallel = axes_angle <= 30.0 or axes_angle >= 150.0
        contact_is_axial = normal_angle <= 30.0 or normal_angle >= 150.0

        if axes_parallel:
            return "先端-先端" if contact_is_axial else "平行・横並び"
        if 60.0 <= axes_angle <= 120.0:
            return "直交"
        return "斜め"


POSE_END_TO_END = Pose(relative_angle_deg=0.0, contact_normal_deg=0.0)
POSE_SIDE_BY_SIDE = Pose(relative_angle_deg=0.0, contact_normal_deg=90.0)
POSE_PERPENDICULAR = Pose(relative_angle_deg=90.0, contact_normal_deg=90.0)


def contactable_adhesins(
    surface: SurfaceQuadrature,
    counts: np.ndarray,
    pose: Pose,
    length_um: float,
    radius_um: float,
    bridge_um: float,
    facing_angle_deg: float = FACING_ANGLE_DEG,
) -> float:
    """指定した配置で、架橋しうるアドヘシンの数を返す。

    両細胞は同じ分布を持つとし、それぞれの側で接触可能な分子を数えて少ない方を
    採用する。これにより1分子1結合の制約（A5）を満たす。
    """

    half_cylinder = _validate(length_um, radius_um)
    theta = math.radians(pose.relative_angle_deg)
    alpha = math.radians(pose.contact_normal_deg)

    axis_1 = np.array([1.0, 0.0, 0.0])
    axis_2 = np.array([math.cos(theta), math.sin(theta), 0.0])
    normal = np.array([math.cos(alpha), math.sin(alpha), 0.0])

    # 2細胞がちょうど触れる位置に配置する
    extent_1 = half_cylinder * abs(float(axis_1 @ normal)) + radius_um
    extent_2 = half_cylinder * abs(float(axis_2 @ normal)) + radius_um
    centre_1 = np.zeros(3)
    centre_2 = normal * (extent_1 + extent_2)

    rotation = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0],
            [math.sin(theta), math.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    facing_cosine = math.cos(math.radians(facing_angle_deg))

    def count_on_side(
        points: np.ndarray,
        normals: np.ndarray,
        partner_centre: np.ndarray,
        partner_axis: np.ndarray,
    ) -> float:
        gap, partner_normal = _gap_and_normal_to_partner(
            points, partner_centre, partner_axis, half_cylinder, radius_um
        )
        facing = np.einsum("ij,ij->i", normals, partner_normal) <= -facing_cosine
        reachable = (gap >= -1e-7) & (gap <= bridge_um) & facing
        return float(counts[reachable].sum())

    from_cell_1 = count_on_side(
        surface.points + centre_1, surface.normals, centre_2, axis_2
    )
    from_cell_2 = count_on_side(
        surface.points @ rotation.T + centre_2,
        surface.normals @ rotation.T,
        centre_1,
        axis_1,
    )
    return min(from_cell_1, from_cell_2)


def contact_area_fraction(
    surface: SurfaceQuadrature,
    pose: Pose,
    length_um: float,
    radius_um: float,
    bridge_um: float,
    facing_angle_deg: float = FACING_ANGLE_DEG,
) -> float:
    """接触領域が全表面積に占める割合。

    Costan et al. (2024) が A_con として報告している量に対応する。
    分子の分布によらない純粋な幾何量なので、一様分布で計算する。
    """

    uniform = adhesin_counts(
        surface,
        n_total=1.0,
        f_pole=polar_area_fraction(length_um, radius_um),
    )
    return contactable_adhesins(
        surface, uniform, pose, length_um, radius_um, bridge_um, facing_angle_deg
    )
