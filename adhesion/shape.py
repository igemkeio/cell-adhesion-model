"""細胞の形（スフェロシリンダー）を1つの値にまとめる。

細胞は円柱の両端に半球を付けた形で、独立な寸法は2つしかない。

    radius_um   半径 r。円柱の断面であり、両端の半球の曲率半径でもある
    length_um   全長 l（端から端まで）。円柱部分の長さは l - 2r

半球の半径は円柱の半径と同じ（接合部で表面が滑らかにつながる）ので、
**半球の曲率 κ = 1/r は円柱の太さと切り離せない**。「円柱部分の長さ」と
「半球の曲率」を独立に振りたいときは、全長 l ではなく円柱長 l_c = l - 2r を
固定して r を動かす（`CellShape.from_cylinder`）。全長 l を固定したまま r を
動かすと、円柱長も同時に変わってしまう。

細胞の形に依存する量（表面積・体積・等価球半径・アスペクト比など）はここに
集め、`adhesion.model` と `adhesion.dynamics` は同じ `CellShape` を共有する。
これにより、形を変えたときに片方だけが古い寸法のままになる事故を防ぐ。
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class CellShape:
    """スフェロシリンダー形の細胞の寸法。

    length_um  全長 l（µm）。円柱部分の長さ l_c = l - 2r に両端の半球の直径を足したもの
    radius_um  半径 r（µm）。半球の曲率半径でもある（曲率は 1/r）

    l = 2r のとき円柱部分が 0 になり、細胞は球になる。
    """

    length_um: float
    radius_um: float

    def __post_init__(self) -> None:
        if not self.radius_um > 0:
            raise ValueError("radius_um は正である必要があります")
        if self.length_um < 2 * self.radius_um:
            raise ValueError("length_um は 2 * radius_um 以上である必要があります")

    @classmethod
    def from_cylinder(cls, cylinder_length_um: float, radius_um: float) -> "CellShape":
        """円柱部分の長さと半径から作る。両者を独立に振りたいときはこちら。"""

        if cylinder_length_um < 0:
            raise ValueError("cylinder_length_um は非負である必要があります")
        return cls(length_um=cylinder_length_um + 2.0 * radius_um, radius_um=radius_um)

    # --- 寸法 -------------------------------------------------------------
    @property
    def cylinder_length_um(self) -> float:
        """円柱部分の長さ l_c = l - 2r（µm）。"""
        return self.length_um - 2.0 * self.radius_um

    @property
    def half_cylinder_um(self) -> float:
        """円柱部分の半分の長さ（µm）。細胞中心から半球の付け根までの距離。"""
        return 0.5 * (self.length_um - 2.0 * self.radius_um)

    @property
    def curvature_per_um(self) -> float:
        """半球の曲率 κ = 1/r（1/µm）。"""
        return 1.0 / self.radius_um

    @property
    def aspect_ratio(self) -> float:
        """全長 / 直径。球で 1。"""
        return self.length_um / (2.0 * self.radius_um)

    # --- 面積・体積 ---------------------------------------------------------
    @property
    def surface_area_um2(self) -> float:
        """全表面積（µm²）。半球2つ（= 球1つ）と円柱側面の和。"""
        r = self.radius_um
        return 4.0 * math.pi * r**2 + 2.0 * math.pi * r * self.cylinder_length_um

    @property
    def cap_area_um2(self) -> float:
        """両端の半球の面積の合計（µm²）。"""
        return 4.0 * math.pi * self.radius_um**2

    @property
    def side_area_um2(self) -> float:
        """円柱側面の面積（µm²）。"""
        return 2.0 * math.pi * self.radius_um * self.cylinder_length_um

    @property
    def polar_area_fraction(self) -> float:
        """両端の半球が全表面積に占める割合。局在させない一様分布の f_pole にあたる。"""
        return self.cap_area_um2 / self.surface_area_um2

    @property
    def volume_um3(self) -> float:
        """体積（µm³）。"""
        r = self.radius_um
        return (4.0 / 3.0) * math.pi * r**3 + math.pi * r**2 * self.cylinder_length_um

    @property
    def equivalent_sphere_radius_um(self) -> float:
        """同じ体積を持つ球の半径（µm）。拡散係数と沈降速度の見積もりに使う。"""

        r = self.radius_um
        cylinder_length = self.length_um - 2.0 * r
        volume_um3 = (4.0 / 3.0) * math.pi * r**3 + math.pi * r**2 * cylinder_length
        return (3.0 * volume_um3 / (4.0 * math.pi)) ** (1.0 / 3.0)

    @property
    def planar_footprint_um2(self) -> float:
        """2次元に射影した面積（µm²）。長方形（l_c × 2r）と円（πr²）の和。

        平面拘束のシミュレーションで、箱に対する細胞の占有率を見積もるのに使う。
        """
        r = self.radius_um
        return self.cylinder_length_um * 2.0 * r + math.pi * r**2

    @property
    def rotational_inertia_um2(self) -> float:
        """拘束の配分に使う回転の実効慣性（µm²）。細い棒の L²/12。"""
        return self.length_um**2 / 12.0

    def describe(self) -> str:
        return (
            f"l={self.length_um:.3g} µm（円柱 {self.cylinder_length_um:.3g} µm）, "
            f"r={self.radius_um:.3g} µm（曲率 {self.curvature_per_um:.3g} /µm）"
        )
