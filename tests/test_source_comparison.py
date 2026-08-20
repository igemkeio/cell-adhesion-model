"""原著の導出式との照合を固定する。

Costan et al. (2024) の Supplemental Note 5 は接触面積を式の形で与えている。
本モデルは同じ量を独立に導出しているので、式レベルで突き合わせられる。
ここが崩れたら、幾何の実装が原著から乖離したことを意味する。

詳細な比較は docs/source_comparison.md を参照。
"""

import math
import unittest

from adhesion.geometry import (
    POSE_END_TO_END,
    POSE_SIDE_BY_SIDE,
    contact_area_fraction,
    contact_patch_area,
    contact_strip_area,
    total_surface_area,
)
from adhesion.model import BRIDGE_UM, CELL_LENGTH_UM, CELL_RADIUS_UM, cell_surface

#: 原著の記法。M は結合長の半分で、2M が2細胞をまたぐ全長にあたる。
HALF_BRIDGE_UM = BRIDGE_UM / 2
CYLINDER_LENGTH_UM = CELL_LENGTH_UM - 2 * CELL_RADIUS_UM


class EndToEndFormulaTests(unittest.TestCase):
    """端 - 端接触。原著が厳密に計算しており、式が一致するはずの配置。"""

    def test_patch_area_equals_the_published_formula(self):
        """本モデルの π(2M)r が、原著 (S32) の 2πrM と一致すること。"""
        published = 2 * math.pi * CELL_RADIUS_UM * HALF_BRIDGE_UM
        ours = contact_patch_area(CELL_RADIUS_UM, BRIDGE_UM)
        self.assertAlmostEqual(ours, published, places=12)

    def test_patch_area_matches_the_published_value(self):
        """原著が (S32) で示した数値 0.010 µm² と一致すること。"""
        self.assertAlmostEqual(
            contact_patch_area(CELL_RADIUS_UM, BRIDGE_UM), 0.010, delta=0.0005
        )

    def test_fractional_area_matches_the_published_value(self):
        """原著 (S33) の Â_ee = M/l と同程度になること。"""
        published = HALF_BRIDGE_UM / CELL_LENGTH_UM
        ours = contact_area_fraction(
            cell_surface(60000), POSE_END_TO_END, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        )
        self.assertAlmostEqual(ours, published, delta=0.0005)


class ParallelContactTests(unittest.TestCase):
    """平行接触。原著の式は隙間の扱いが異なるので、厳密解と照合する。"""

    def test_strip_area_matches_exact_cylinder_geometry(self):
        """放物線近似の 2√(2Mr)·L が、厳密な円筒幾何と一致すること。

        平行に接する2円柱では隙間の合計が 2r(1 − cos θ) なので、
        これが 2M に等しくなる θ から半幅 rθ を求める。
        """
        theta = math.acos(1 - HALF_BRIDGE_UM / CELL_RADIUS_UM)
        exact = 2 * CELL_RADIUS_UM * theta * CYLINDER_LENGTH_UM
        ours = contact_strip_area(CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM)
        self.assertAlmostEqual(ours / exact, 1.0, delta=0.002)

    def test_our_value_is_smaller_than_the_published_one(self):
        """原著の記載値 0.165 µm² より小さいこと。隙間の基準が異なるため。"""
        ours = contact_strip_area(CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM)
        self.assertLess(ours, 0.165)
        self.assertGreater(ours, 0.5 * 0.165)


class ConsistencyTests(unittest.TestCase):
    """全配置で同一の基準を使っていること。"""

    def test_all_contact_areas_use_the_same_gap_criterion(self):
        """端 - 端と平行の面積比が、同一基準から導かれる値になること。

        どちらも隙間 2M を基準にすると、比は
            2√(2Mr)(l−2r) / (π·2M·r) = 2(l−2r) / (π√(2Mr))
        になる。
        """
        expected = 2 * CYLINDER_LENGTH_UM / (
            math.pi * math.sqrt(BRIDGE_UM * CELL_RADIUS_UM)
        )
        ratio = contact_strip_area(
            CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        ) / contact_patch_area(CELL_RADIUS_UM, BRIDGE_UM)
        self.assertAlmostEqual(ratio, expected, places=9)
        self.assertAlmostEqual(ratio, 11.3, delta=0.2)

    def test_conclusion_is_robust_to_the_published_areas(self):
        """原著の値を使っても「局在なしでは束」という結論が変わらないこと。"""
        published_ratio = 0.036 / 0.002
        ours = contact_area_fraction(
            cell_surface(60000), POSE_SIDE_BY_SIDE, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        ) / contact_area_fraction(
            cell_surface(60000), POSE_END_TO_END, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        )
        self.assertGreater(published_ratio, 5.0)
        self.assertGreater(ours, 5.0)

    def test_total_area_matches_the_published_expression(self):
        """全表面積が原著の 2πrl と一致すること。"""
        published = 2 * math.pi * CELL_RADIUS_UM * CELL_LENGTH_UM
        self.assertAlmostEqual(
            total_surface_area(CELL_LENGTH_UM, CELL_RADIUS_UM), published, places=12
        )


if __name__ == "__main__":
    unittest.main()
