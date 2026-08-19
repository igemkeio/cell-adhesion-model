"""幾何計算のテスト。

求積が正しく面積・分子数を保存しているか、解析式と数値解が一致するかを確認する。
"""

import math
import unittest

import numpy as np

from adhesion.geometry import (
    POSE_END_TO_END,
    POSE_PERPENDICULAR,
    POSE_SIDE_BY_SIDE,
    adhesin_counts,
    contact_patch_area,
    contact_strip_area,
    contactable_adhesins,
    polar_area_fraction,
    spherocylinder_surface,
    total_surface_area,
)

LENGTH, RADIUS, BRIDGE = 1.80, 0.40, 0.008


class QuadratureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = spherocylinder_surface(LENGTH, RADIUS, 20000)

    def test_quadrature_reproduces_analytic_surface_area(self):
        self.assertAlmostEqual(
            self.surface.total_area, total_surface_area(LENGTH, RADIUS), places=9
        )

    def test_normals_are_unit_vectors(self):
        norms = np.linalg.norm(self.surface.normals, axis=1)
        self.assertTrue(np.allclose(norms, 1.0, atol=1e-12))

    def test_points_lie_on_the_surface(self):
        """各点が実際にスフェロシリンダー表面上にあること。"""
        half_cylinder = 0.5 * (LENGTH - 2 * RADIUS)
        axial = np.clip(self.surface.points[:, 0], -half_cylinder, half_cylinder)
        nearest = np.column_stack((axial, np.zeros(len(axial)), np.zeros(len(axial))))
        distance = np.linalg.norm(self.surface.points - nearest, axis=1)
        self.assertTrue(np.allclose(distance, RADIUS, atol=1e-12))

    def test_polar_region_area_matches_analytic_fraction(self):
        polar = self.surface.areas[self.surface.regions != 0].sum()
        self.assertAlmostEqual(
            polar / self.surface.total_area,
            polar_area_fraction(LENGTH, RADIUS),
            places=9,
        )

    def test_tip_distance_is_zero_at_the_apex_and_quarter_circumference_at_equator(self):
        caps = self.surface.regions != 0
        self.assertAlmostEqual(float(self.surface.tip_distances[caps].min()), 0.0, delta=0.01)
        self.assertAlmostEqual(
            float(self.surface.tip_distances[caps].max()),
            math.pi * RADIUS / 2,
            delta=0.01,
        )

    def test_invalid_geometry_is_rejected(self):
        with self.assertRaises(ValueError):
            spherocylinder_surface(0.5, 0.4, 1000)
        with self.assertRaises(ValueError):
            spherocylinder_surface(1.8, -0.1, 1000)


class PoseTests(unittest.TestCase):
    def test_named_poses_have_distinct_descriptions(self):
        """呼び名が衝突すると、結果を辞書に集めたとき静かに上書きされる。"""
        poses = (POSE_END_TO_END, POSE_SIDE_BY_SIDE, POSE_PERPENDICULAR)
        names = [pose.description for pose in poses]
        self.assertEqual(len(set(names)), len(names))

    def test_description_distinguishes_parallel_from_perpendicular_axes(self):
        self.assertEqual(POSE_SIDE_BY_SIDE.description, "平行・横並び")
        self.assertEqual(POSE_PERPENDICULAR.description, "直交")
        self.assertEqual(POSE_END_TO_END.description, "先端-先端")


class DistributionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.surface = spherocylinder_surface(LENGTH, RADIUS, 20000)

    def test_molecule_count_is_conserved(self):
        for f_pole in (0.0, 0.3, polar_area_fraction(LENGTH, RADIUS), 0.95, 1.0):
            for width in (None, 0.10):
                with self.subTest(f_pole=f_pole, width=width):
                    counts = adhesin_counts(self.surface, 7300, f_pole, 0.5, width)
                    self.assertAlmostEqual(float(counts.sum()), 7300.0, places=6)

    def test_polar_split_directs_molecules_to_the_correct_pole(self):
        counts = adhesin_counts(self.surface, 1000, 1.0, polar_split=1.0)
        left = counts[self.surface.regions == -1].sum()
        right = counts[self.surface.regions == 1].sum()
        self.assertAlmostEqual(float(left), 1000.0, places=6)
        self.assertAlmostEqual(float(right), 0.0, places=6)

    def test_uniform_fraction_produces_uniform_surface_density(self):
        f_uniform = polar_area_fraction(LENGTH, RADIUS)
        counts = adhesin_counts(self.surface, 7300, f_uniform)
        density = counts / self.surface.areas
        self.assertLess(float(np.std(density) / np.mean(density)), 1e-9)

    def test_tip_focusing_concentrates_molecules_near_the_apex(self):
        wide = adhesin_counts(self.surface, 7300, 0.95, 0.5, None)
        narrow = adhesin_counts(self.surface, 7300, 0.95, 0.5, 0.10)
        near_tip = self.surface.tip_distances < 0.1
        self.assertGreater(float(narrow[near_tip].sum()), float(wide[near_tip].sum()))

    def test_invalid_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            adhesin_counts(self.surface, 100, f_pole=1.5)
        with self.assertRaises(ValueError):
            adhesin_counts(self.surface, -1, f_pole=0.5)
        with self.assertRaises(ValueError):
            adhesin_counts(self.surface, 100, 0.5, tip_width_um=0.0)


class AnalyticAgreementTests(unittest.TestCase):
    """解析式と数値求積が一致すること。"""

    def setUp(self) -> None:
        self.surface = spherocylinder_surface(LENGTH, RADIUS, 60000)

    def test_pole_contact_matches_analytic_patch_area(self):
        """極接触の分子数 ≈ 面密度 × π ℓ R。"""
        f_uniform = polar_area_fraction(LENGTH, RADIUS)
        counts = adhesin_counts(self.surface, 7300, f_uniform)
        numeric = contactable_adhesins(
            self.surface, counts, POSE_END_TO_END, LENGTH, RADIUS, BRIDGE
        )
        density = 7300 / total_surface_area(LENGTH, RADIUS)
        analytic = density * contact_patch_area(RADIUS, BRIDGE)
        self.assertAlmostEqual(numeric / analytic, 1.0, delta=0.15)

    def test_side_contact_is_close_to_analytic_strip_area(self):
        """側面接触は解析式より少し大きい。半球の一部も接触に寄与するため。"""
        f_uniform = polar_area_fraction(LENGTH, RADIUS)
        counts = adhesin_counts(self.surface, 7300, f_uniform)
        numeric = contactable_adhesins(
            self.surface, counts, POSE_SIDE_BY_SIDE, LENGTH, RADIUS, BRIDGE
        )
        density = 7300 / total_surface_area(LENGTH, RADIUS)
        analytic = density * contact_strip_area(LENGTH, RADIUS, BRIDGE)
        self.assertGreater(numeric, analytic)
        self.assertLess(numeric / analytic, 1.3)

    def test_side_contact_area_is_about_eleven_times_the_pole_patch(self):
        ratio = contact_strip_area(LENGTH, RADIUS, BRIDGE) / contact_patch_area(
            RADIUS, BRIDGE
        )
        self.assertAlmostEqual(ratio, 11.3, delta=0.2)

    def test_patch_area_scales_linearly_with_bridge_length(self):
        self.assertAlmostEqual(
            contact_patch_area(RADIUS, 2 * BRIDGE) / contact_patch_area(RADIUS, BRIDGE),
            2.0,
            places=9,
        )

    def test_strip_area_scales_with_square_root_of_bridge_length(self):
        ratio = contact_strip_area(LENGTH, RADIUS, 4 * BRIDGE) / contact_strip_area(
            LENGTH, RADIUS, BRIDGE
        )
        self.assertAlmostEqual(ratio, 2.0, places=9)

    def test_quadrature_is_converged(self):
        """求積の解像度を変えても結果が変わらないこと。"""
        f_uniform = polar_area_fraction(LENGTH, RADIUS)
        values = []
        for n_points in (20000, 60000):
            surface = spherocylinder_surface(LENGTH, RADIUS, n_points)
            counts = adhesin_counts(surface, 7300, f_uniform)
            values.append(
                contactable_adhesins(
                    surface, counts, POSE_SIDE_BY_SIDE, LENGTH, RADIUS, BRIDGE
                )
            )
        self.assertAlmostEqual(values[0] / values[1], 1.0, delta=0.05)


if __name__ == "__main__":
    unittest.main()
