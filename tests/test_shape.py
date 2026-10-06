"""細胞の形（円柱の長さ・半球の曲率）をパラメータにしたときのテスト。

既定の形状（文献値）では従来とまったく同じ計算になること、形を変えたときに
静的モデルと時間発展モデルの両方が同じ寸法を使い、解析的に分かる量が正しく動く
ことを確かめる。
"""

from __future__ import annotations

import math
import unittest
import warnings

import numpy as np

from adhesion import dynamics as dyn
from adhesion import literature as lit
from adhesion import model, shear
from adhesion.geometry import polar_area_fraction, total_surface_area
from adhesion.shape import CellShape

N_MID = float(lit.ADHESINS_AT_MID_INDUCTION)


class CellShapeTests(unittest.TestCase):
    def test_default_shape_is_the_literature_value(self):
        self.assertEqual(model.DEFAULT_SHAPE.length_um, 1.80)
        self.assertEqual(model.DEFAULT_SHAPE.radius_um, 0.40)
        self.assertAlmostEqual(model.DEFAULT_SHAPE.cylinder_length_um, 1.0)

    def test_from_cylinder_sets_cylinder_and_radius_independently(self):
        """円柱長を固定して半径だけ動かすと、全長は 2Δr だけ変わる。"""

        thin = CellShape.from_cylinder(1.0, 0.25)
        fat = CellShape.from_cylinder(1.0, 0.60)
        self.assertAlmostEqual(thin.cylinder_length_um, fat.cylinder_length_um)
        self.assertAlmostEqual(fat.length_um - thin.length_um, 2 * (0.60 - 0.25))
        self.assertAlmostEqual(fat.curvature_per_um, 1 / 0.60)

    def test_areas_agree_with_the_geometry_module(self):
        for cylinder, radius in [(1.0, 0.4), (0.3, 0.25), (3.0, 0.6)]:
            shape = CellShape.from_cylinder(cylinder, radius)
            self.assertAlmostEqual(
                shape.surface_area_um2,
                total_surface_area(shape.length_um, shape.radius_um),
            )
            self.assertAlmostEqual(
                shape.polar_area_fraction,
                polar_area_fraction(shape.length_um, shape.radius_um),
            )

    def test_zero_cylinder_is_a_sphere(self):
        sphere = CellShape.from_cylinder(0.0, 0.5)
        self.assertAlmostEqual(sphere.volume_um3, 4.0 / 3.0 * math.pi * 0.5**3)
        self.assertAlmostEqual(sphere.surface_area_um2, 4.0 * math.pi * 0.5**2)
        self.assertAlmostEqual(sphere.polar_area_fraction, 1.0)
        self.assertAlmostEqual(sphere.equivalent_sphere_radius_um, 0.5)

    def test_rejects_impossible_dimensions(self):
        with self.assertRaises(ValueError):
            CellShape(length_um=0.5, radius_um=0.4)  # 全長が直径より短い
        with self.assertRaises(ValueError):
            CellShape(length_um=1.8, radius_um=0.0)
        with self.assertRaises(ValueError):
            CellShape.from_cylinder(-0.1, 0.4)


class StaticModelShapeTests(unittest.TestCase):
    def test_omitting_shape_means_the_default_shape(self):
        a = model.predict(N_MID, 0.95, 0.5, 0.10)
        b = model.predict(N_MID, 0.95, 0.5, 0.10, shape=model.DEFAULT_SHAPE)
        self.assertEqual(a.bonds_heavy_pole, b.bonds_heavy_pole)
        self.assertEqual(a.bonds_side_by_side, b.bonds_side_by_side)

    def test_uniform_f_pole_depends_on_the_shape(self):
        self.assertEqual(model.uniform_f_pole(), model.UNIFORM_F_POLE)
        longer = CellShape.from_cylinder(3.0, 0.4)
        self.assertLess(model.uniform_f_pole(longer), model.UNIFORM_F_POLE)

    def test_polar_bonds_grow_with_cap_radius_and_ignore_cylinder_length(self):
        """先端に集中させた分子の極接触は、接触パッチ π ℓ R に比例し、円柱長に依らない。"""

        def polar(cylinder, radius):
            shape = CellShape.from_cylinder(cylinder, radius)
            return model.predict(N_MID, 0.95, 0.5, 0.10, shape=shape).bonds_heavy_pole

        self.assertGreater(polar(1.0, 0.6), polar(1.0, 0.4))
        self.assertGreater(polar(1.0, 0.4), polar(1.0, 0.25))
        self.assertAlmostEqual(polar(0.3, 0.4) / polar(3.0, 0.4), 1.0, delta=0.03)

    def test_side_bonds_follow_the_strip_formula(self):
        """f_pole を固定した側面結合は √(ℓR) の接触帯で決まり、半径とともに減る。

        側面の分子密度は 1/(R l_c)、接触帯の面積は √(ℓR) l_c なので、
        円柱長は消えて S ∝ R^(-1/2) になる。
        """

        def side(cylinder, radius):
            shape = CellShape.from_cylinder(cylinder, radius)
            return model.predict(N_MID, 0.95, 0.5, 0.10, shape=shape).bonds_side_by_side

        self.assertAlmostEqual(side(0.5, 0.4) / side(3.0, 0.4), 1.0, delta=0.03)
        expected = math.sqrt(0.25 / 0.6)
        self.assertAlmostEqual(side(1.0, 0.6) / side(1.0, 0.25), expected, delta=0.06)

    def test_uniform_side_bonds_grow_with_cylinder_length(self):
        def side(cylinder):
            shape = CellShape.from_cylinder(cylinder, 0.4)
            result = model.predict(N_MID, model.uniform_f_pole(shape), 0.5, None, shape=shape)
            return result.bonds_side_by_side

        self.assertLess(side(0.3), side(1.0))
        self.assertLess(side(1.0), side(3.0))

    def test_morphology_threshold_uses_the_shape_specific_uniform_baseline(self):
        long_cell = CellShape.from_cylinder(3.0, 0.4)
        threshold = model.morphology_threshold(0.20, shape=long_cell)
        self.assertIsNotNone(threshold)
        self.assertGreater(threshold, model.uniform_f_pole(long_cell))


class DynamicsShapeTests(unittest.TestCase):
    def test_diffusion_depends_only_on_volume(self):
        """等価球近似（A16, A26）では、同体積なら形が違っても拡散係数は同じ。"""

        base = CellShape.from_cylinder(1.0, 0.4)
        # 同体積で半径を 0.3 にしたときの円柱長を求める
        r = 0.3
        cylinder = (base.volume_um3 - 4.0 / 3.0 * math.pi * r**3) / (math.pi * r**2)
        other = CellShape.from_cylinder(cylinder, r)
        self.assertAlmostEqual(other.volume_um3, base.volume_um3)
        self.assertAlmostEqual(
            dyn.diffusion_coefficients(base)[0], dyn.diffusion_coefficients(other)[0]
        )
        # 大きい細胞ほど遅い
        big = CellShape.from_cylinder(3.0, 0.4)
        self.assertLess(
            dyn.diffusion_coefficients(big)[0], dyn.diffusion_coefficients(base)[0]
        )

    def test_default_diffusion_matches_module_constants(self):
        d_trans, d_rot = dyn.diffusion_coefficients(dyn.DEFAULT_SHAPE)
        self.assertEqual(d_trans, dyn.D_TRANS0_UM2_S)
        self.assertEqual(d_rot, dyn.D_ROT0_RAD2_S)

    def test_config_rejects_a_cylinder_that_is_too_short(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(shape=CellShape.from_cylinder(0.0, 0.4))
        dyn.DynamicsConfig(shape=CellShape.from_cylinder(dyn.MIN_CYLINDER_LENGTH_UM, 0.4))

    def test_bond_counts_come_from_the_configured_shape(self):
        shape = CellShape.from_cylinder(2.0, 0.55)
        config = dyn.DynamicsConfig(
            shape=shape, n_cells=4, box_size_um=12.0, total_time_s=0.2,
            f_pole=0.95, tip_width_um=0.10,
        )
        result = dyn.run(config)
        expected = model.predict(
            config.n_total, 0.95, 0.5, 0.10, shape=shape,
            distinct_adhesins_per_pole=config.distinct_adhesins_per_pole,
        )
        self.assertEqual(result.polar_bond_count, expected.bonds_heavy_pole)
        self.assertEqual(result.side_bond_count, expected.bonds_side_by_side)

    def test_bond_counts_override_replaces_the_computed_counts(self):
        config = dyn.DynamicsConfig(
            shape=CellShape.from_cylinder(2.0, 0.55), n_cells=4, box_size_um=12.0,
            total_time_s=0.2, bond_counts_override=(123.0, 4.5),
        )
        result = dyn.run(config)
        self.assertEqual((result.polar_bond_count, result.side_bond_count), (123.0, 4.5))
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(bond_counts_override=(-1.0, 2.0))

    def test_diffusion_scale_scales_the_rotational_random_walk(self):
        """回転拡散の倍率を 0.25 にすると、向きの揺らぎの分散がおよそ 1/4 になる。"""

        warnings.simplefilter("ignore", dyn.CoarseTimeStepWarning)  # dt に依らない性質なので粗い dt で速く

        def angle_variance(scale: float) -> float:
            changes = []
            for seed in range(60):
                config = dyn.DynamicsConfig(
                    n_cells=2, box_size_um=60.0, total_time_s=2.0, seed=seed,
                    frame_stride=10**9, diffusion_scale=(1.0, scale), dt_s=0.02,
                )
                frames = dyn.run(config).frames
                changes.extend(frames[-1].angles - frames[0].angles)
            return float(np.var(changes))

        ratio = angle_variance(0.25) / angle_variance(1.0)
        self.assertAlmostEqual(ratio, 0.25, delta=0.1)
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(diffusion_scale=(1.0, 0.0))

    def test_invariants_hold_for_non_default_shapes(self):
        """形を変えても、貫入せず、結合は架橋長付近に留まり、箱から出ない。"""

        for cylinder, radius in [(3.0, 0.3), (0.2, 0.7)]:
            with self.subTest(cylinder=cylinder, radius=radius):
                shape = CellShape.from_cylinder(cylinder, radius)
                config = dyn.DynamicsConfig(
                    shape=shape, n_cells=12, box_size_um=16.0, total_time_s=12.0,
                    seed=1, f_pole=0.95, tip_width_um=0.10, frame_stride=10**9,
                )
                last = dyn.run(config).frames[-1]
                n = config.n_cells
                worst = min(
                    dyn._surface_gap(last, i, j, shape)
                    for i in range(n) for j in range(i + 1, n)
                )
                self.assertGreater(worst, -0.08)
                self.assertTrue(np.all(last.positions >= 0.0))
                self.assertTrue(np.all(last.positions <= config.box_size_um))
                self.assertTrue(last.bonds, "この設定では結合が少なくとも1本できるはず")
                lengths = [
                    float(
                        np.linalg.norm(
                            dyn.attachment_point(
                                last.positions[b.j], last.angles[b.j], b.axial_j, b.delta_j,
                                shape.radius_um,
                            )
                            - dyn.attachment_point(
                                last.positions[b.i], last.angles[b.i], b.axial_i, b.delta_i,
                                shape.radius_um,
                            )
                        )
                    )
                    for b in last.bonds
                ]
                # 既存テスト（r = 0.4 µm で 0.06 µm）と同じ比率で、半径に比例させる。
                # 緩和の残差は細胞の大きさ（てこの腕の長さ）に比例するため。バグ時は
                # 細胞長スケール（1 µm 以上）まで伸びるので、これでも十分に検出できる。
                self.assertLess(float(np.median(lengths)), 0.2 * shape.radius_um)

    def test_attachment_points_sit_on_the_cell_surface(self):
        """付着点は中心軸から半径ぶんだけ離れた位置にある（形に依らず）。"""

        for radius in (0.25, 0.4, 0.7):
            point = dyn.attachment_point(
                np.zeros(2), 0.0, 0.3, math.pi / 2, radius
            )
            self.assertAlmostEqual(float(point[0]), 0.3)
            self.assertAlmostEqual(float(point[1]), radius)


class TimeStepTests(unittest.TestCase):
    """既定の時間刻みが、掃引する全ての形で収束した値であることの目安を固定する。"""

    #: docs/shape_scan.md の掃引範囲。円柱長だけを変える（r = 0.4）掃引と、半径だけを変える
    #: （l_c = 1.0）掃引。両方を同時に変えた形（格子）は掃引していない。
    SWEEP_SHAPES = [CellShape.from_cylinder(c, 0.4) for c in (0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0)] + [
        CellShape.from_cylinder(1.0, r) for r in (0.25, 0.3, 0.4, 0.5, 0.6, 0.8)
    ]

    def test_default_time_step_is_the_converged_value(self):
        config = dyn.DynamicsConfig()
        self.assertEqual(config.dt_s, 0.0025)
        # フレームは 0.5 s ごと
        self.assertAlmostEqual(config.frame_stride * config.dt_s, 0.5)

    def test_default_time_step_is_fine_enough_for_every_swept_shape(self):
        for shape in self.SWEEP_SHAPES:
            with self.subTest(shape=shape.describe()):
                ratio = dyn.step_displacement_ratio(shape, dyn.DynamicsConfig().dt_s)
                self.assertLessEqual(ratio, dyn.STEP_DISPLACEMENT_LIMIT)

    def test_thin_and_long_cells_are_too_coarse_even_at_the_default(self):
        """細くて長い細胞（r = 0.25 µm かつ l_c ≥ 2 µm）は、既定の dt でも粗い。

        掃引していない領域。`DynamicsConfig` が警告するので、dt を小さくして使うこと。
        """

        shape = CellShape.from_cylinder(2.0, 0.25)
        self.assertGreater(dyn.step_displacement_ratio(shape, dyn.DynamicsConfig().dt_s),
                           dyn.STEP_DISPLACEMENT_LIMIT)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dyn.DynamicsConfig(shape=shape)
        self.assertTrue(any(issubclass(w.category, dyn.CoarseTimeStepWarning) for w in caught))

    def test_the_old_default_was_too_coarse_for_the_literature_shape(self):
        """旧既定の 0.02 s は、側面結合を過大に数えた領域（δ ≳ 0.7）に入っていた。"""

        ratio = dyn.step_displacement_ratio(dyn.DEFAULT_SHAPE, 0.02)
        self.assertGreater(ratio, dyn.STEP_DISPLACEMENT_LIMIT)

    def test_displacement_ratio_scales_with_the_square_root_of_dt(self):
        a = dyn.step_displacement_ratio(dyn.DEFAULT_SHAPE, 0.01)
        b = dyn.step_displacement_ratio(dyn.DEFAULT_SHAPE, 0.04)
        self.assertAlmostEqual(b / a, 2.0)

    def test_thinner_cells_need_a_finer_time_step(self):
        thin = dyn.step_displacement_ratio(CellShape.from_cylinder(1.0, 0.25), 0.0025)
        fat = dyn.step_displacement_ratio(CellShape.from_cylinder(1.0, 0.6), 0.0025)
        self.assertGreater(thin, fat)

    def test_a_coarse_time_step_warns_and_the_default_does_not(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dyn.DynamicsConfig()
        self.assertEqual([w for w in caught if issubclass(w.category, dyn.CoarseTimeStepWarning)], [])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dyn.DynamicsConfig(dt_s=0.02)
        self.assertTrue(any(issubclass(w.category, dyn.CoarseTimeStepWarning) for w in caught))
        # 細すぎて既定でも粗い形は警告される
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dyn.DynamicsConfig(shape=CellShape.from_cylinder(1.0, 0.2))
        self.assertTrue(any(issubclass(w.category, dyn.CoarseTimeStepWarning) for w in caught))

    def test_rejects_a_non_positive_time_step(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(dt_s=0.0)


class ShearShapeTests(unittest.TestCase):
    def test_chain_tension_defaults_to_the_literature_cell(self):
        explicit = shear.chain_tension_pn(10, 300.0, 1.80, 0.40)
        self.assertAlmostEqual(shear.chain_tension_pn(10, 300.0), explicit)

    def test_longer_cells_carry_more_tension(self):
        self.assertGreater(
            shear.chain_tension_pn(10, 300.0, 3.6, 0.40),
            shear.chain_tension_pn(10, 300.0, 1.8, 0.40),
        )


if __name__ == "__main__":
    unittest.main()
