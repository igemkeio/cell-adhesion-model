"""Costan et al. (2024) の報告値との照合。

このモデルの信頼性はここに集約される。調整パラメータを一切使わずに、
独立に測定された量を再現できるかどうかを固定する。将来モデルを変更したときに
文献との整合が崩れたら、ここで検出される。
"""

import unittest

from adhesion import literature as lit
from adhesion.geometry import (
    POSE_END_TO_END,
    POSE_PERPENDICULAR,
    POSE_SIDE_BY_SIDE,
    contact_area_fraction,
    total_surface_area,
)
from adhesion.model import (
    BRIDGE_UM,
    CELL_LENGTH_UM,
    CELL_RADIUS_UM,
    UNIFORM_F_POLE,
    Morphology,
    cell_surface,
    predict,
)

N_MID = float(lit.ADHESINS_AT_MID_INDUCTION)


class ReportedBondCountTests(unittest.TestCase):
    """報告された2細胞間の結合数を再現するか。"""

    def setUp(self) -> None:
        self.prediction = predict(N_MID, UNIFORM_F_POLE)

    def test_parallel_contact_matches_within_two_sigma(self):
        reported = lit.REPORTED_BONDS_PARALLEL
        deviation = abs(self.prediction.bonds_side_by_side - reported.value)
        self.assertLess(deviation, 2 * reported.uncertainty)

    def test_pole_contact_matches_within_three_sigma(self):
        reported = lit.REPORTED_BONDS_PERPENDICULAR
        deviation = abs(self.prediction.bonds_heavy_pole - reported.value)
        self.assertLess(deviation, 3 * reported.uncertainty)

    def test_side_contact_dominates_pole_contact_for_uniform_display(self):
        """一様分布では側面接触が1桁優勢になる。"""
        ratio = self.prediction.bonds_side_by_side / self.prediction.bonds_heavy_pole
        self.assertGreater(ratio, 5.0)


class ContactAreaTests(unittest.TestCase):
    """報告された接触面積割合 A_con のレンジに収まるか。"""

    def setUp(self) -> None:
        self.surface = cell_surface()

    def _fraction(self, pose) -> float:
        return contact_area_fraction(
            self.surface, pose, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        )

    def test_all_poses_fall_inside_reported_range(self):
        low = lit.REPORTED_CONTACT_FRACTION_MIN.value
        high = lit.REPORTED_CONTACT_FRACTION_MAX.value
        for pose in (POSE_END_TO_END, POSE_SIDE_BY_SIDE, POSE_PERPENDICULAR):
            with self.subTest(pose=pose.description):
                fraction = self._fraction(pose)
                self.assertGreaterEqual(fraction, low * 0.9)
                self.assertLessEqual(fraction, high * 1.1)

    def test_smallest_contact_is_near_the_reported_lower_bound(self):
        self.assertAlmostEqual(
            self._fraction(POSE_END_TO_END),
            lit.REPORTED_CONTACT_FRACTION_MIN.value,
            delta=0.002,
        )

    def test_side_by_side_is_the_largest_contact(self):
        side = self._fraction(POSE_SIDE_BY_SIDE)
        for pose in (POSE_END_TO_END, POSE_PERPENDICULAR):
            self.assertGreater(side, self._fraction(pose))


class MeasuredLocalizationTests(unittest.TestCase):
    """実測された分布が何を予測するか。"""

    def test_uniform_display_is_the_measured_baseline(self):
        """A10。極と中央で有意差がないという報告に対応する。"""
        from adhesion.geometry import polar_area_fraction

        self.assertAlmostEqual(
            UNIFORM_F_POLE,
            polar_area_fraction(CELL_LENGTH_UM, CELL_RADIUS_UM),
            places=12,
        )

    def test_measured_localization_predicts_bundles_not_chains(self):
        self.assertEqual(predict(N_MID, UNIFORM_F_POLE).morphology, Morphology.BUNDLE)

    def test_measured_surface_density_is_reproduced(self):
        """一様分布での面密度が報告値と一致するか。"""
        n_max = float(lit.ADHESINS_AT_MAX_INDUCTION)
        area = total_surface_area(CELL_LENGTH_UM, CELL_RADIUS_UM)
        reported = lit.SURFACE_DENSITY_AT_MAX
        self.assertAlmostEqual(n_max / area, reported.value, delta=reported.uncertainty)


if __name__ == "__main__":
    unittest.main()
