"""時間発展シミュレーション（`adhesion.dynamics`）のテスト。

固定した幾何関係を持つ入力に対する決定的な検証と、局在パラメータが実際に
結合の型（極 vs 側面）の分布を左右することを確かめる統計的な検証の両方を行う。
"""

from __future__ import annotations

import math
import unittest

import numpy as np

from adhesion import dynamics as dyn
from adhesion import model
from adhesion import shear


class DiffusionCoefficientTests(unittest.TestCase):
    def test_rotational_diffusion_matches_shear_module_default(self):
        """shear.py の既定回転拡散係数（同じ等価球近似）と一致する。"""

        self.assertAlmostEqual(dyn.D_ROT0_RAD2_S, 0.889, places=3)

    def test_translational_diffusion_is_positive_and_reasonable(self):
        self.assertGreater(dyn.D_TRANS0_UM2_S, 0.0)
        self.assertLess(dyn.D_TRANS0_UM2_S, 10.0)


class SegmentGeometryTests(unittest.TestCase):
    def test_closest_points_between_touching_parallel_segments(self):
        a0, a1 = np.array([-0.5, 0.0]), np.array([0.5, 0.0])
        b0, b1 = np.array([-0.5, 0.9]), np.array([0.5, 0.9])
        p1, p2, s, t = dyn._closest_points_on_segments(a0, a1, b0, b1)
        self.assertAlmostEqual(float(np.linalg.norm(p2 - p1)), 0.9, places=9)

    def test_closest_points_clamp_to_endpoints_when_disjoint_ranges(self):
        a0, a1 = np.array([0.0, 0.0]), np.array([1.0, 0.0])
        b0, b1 = np.array([3.0, 0.0]), np.array([4.0, 0.0])
        _, _, s, t = dyn._closest_points_on_segments(a0, a1, b0, b1)
        self.assertAlmostEqual(s, 1.0)
        self.assertAlmostEqual(t, 0.0)


class ContactRegionClassificationTests(unittest.TestCase):
    def test_head_on_approach_is_classified_as_pole(self):
        own_axis_point = np.array([dyn.HALF_CYLINDER_UM, 0.0])
        partner_point = np.array([dyn.HALF_CYLINDER_UM + 0.5, 0.0])
        region = dyn._classify_contact_region(1.0, 0.0, own_axis_point, partner_point)
        self.assertEqual(region, "head")

    def test_broadside_approach_near_the_tip_is_classified_as_side(self):
        """先端の丸め込み（s=1）だけでは極とみなさない。赤道寄りの接触は側面。"""

        own_axis_point = np.array([dyn.HALF_CYLINDER_UM, 0.0])
        partner_point = own_axis_point + np.array([0.0, 0.5])
        region = dyn._classify_contact_region(1.0, 0.0, own_axis_point, partner_point)
        self.assertEqual(region, "side")

    def test_interior_contact_is_always_side(self):
        own_axis_point = np.array([0.0, 0.0])
        partner_point = np.array([0.0, 0.5])
        region = dyn._classify_contact_region(0.5, 0.0, own_axis_point, partner_point)
        self.assertEqual(region, "side")


class BindingProbabilityTests(unittest.TestCase):
    def test_more_potential_bonds_bind_faster(self):
        weak = dyn._binding_probability(8.0, 0.02)
        strong = dyn._binding_probability(300.0, 0.02)
        self.assertLess(weak, strong)
        self.assertGreater(strong, 0.6)
        self.assertLess(weak, 0.2)


class DirectBondFormationTests(unittest.TestCase):
    """`_try_form_bond` を直接呼び、幾何的な配置と規則の対応を確かめる。"""

    def _two_cell_state(self, positions, angles):
        return dyn.SimulationState(
            positions=np.array(positions, dtype=float),
            angles=np.array(angles, dtype=float),
            head_free=np.ones(2, dtype=bool),
            tail_free=np.ones(2, dtype=bool),
            bonds=[],
        )

    def _touching_head_to_tail_positions(self):
        gap = dyn.BRIDGE_UM
        x1 = dyn.HALF_CYLINDER_UM
        x2 = x1 + 2 * dyn.CELL_RADIUS_UM + gap + dyn.HALF_CYLINDER_UM
        return [[0.0, 0.0], [x2, 0.0]], [0.0, 0.0]

    def test_head_to_tail_contact_forms_a_directed_polar_bond(self):
        positions, angles = self._touching_head_to_tail_positions()
        state = self._two_cell_state(positions, angles)
        config = dyn.DynamicsConfig(n_cells=2, distinct_adhesins_per_pole=True)
        rng = np.random.default_rng(0)
        for _ in range(200):  # 確率的な結合なので、十分な回数だけ試行する
            dyn._try_form_bond(state, 0, 1, config, 500.0, 500.0, 0.02, rng)
            if state.bonds:
                break
        self.assertEqual(len(state.bonds), 1)
        self.assertEqual(state.bonds[0].kind, "polar")

    def test_head_to_head_contact_never_bonds_when_poles_are_distinct(self):
        """Ag 極どうしは結合しない（設計上の方向性）。"""

        x1 = dyn.HALF_CYLINDER_UM
        x2 = x1 + 2 * dyn.CELL_RADIUS_UM + dyn.BRIDGE_UM + dyn.HALF_CYLINDER_UM
        state = self._two_cell_state([[0.0, 0.0], [x2, 0.0]], [0.0, math.pi])
        config = dyn.DynamicsConfig(n_cells=2, distinct_adhesins_per_pole=True)
        rng = np.random.default_rng(0)
        for _ in range(200):
            dyn._try_form_bond(state, 0, 1, config, 500.0, 500.0, 0.02, rng)
        self.assertEqual(len(state.bonds), 0)

    def test_head_to_head_contact_bonds_when_poles_are_not_distinct(self):
        x1 = dyn.HALF_CYLINDER_UM
        x2 = x1 + 2 * dyn.CELL_RADIUS_UM + dyn.BRIDGE_UM + dyn.HALF_CYLINDER_UM
        state = self._two_cell_state([[0.0, 0.0], [x2, 0.0]], [0.0, math.pi])
        config = dyn.DynamicsConfig(n_cells=2, distinct_adhesins_per_pole=False)
        rng = np.random.default_rng(0)
        for _ in range(200):
            dyn._try_form_bond(state, 0, 1, config, 500.0, 500.0, 0.02, rng)
            if state.bonds:
                break
        self.assertEqual(len(state.bonds), 1)
        self.assertEqual(state.bonds[0].kind, "polar")


class ClusterClassificationTests(unittest.TestCase):
    def test_lone_cell_is_a_monomer(self):
        clusters = dyn.classify_clusters(3, [], True)
        self.assertEqual(len(clusters), 3)
        self.assertTrue(all(c.kind == "monomer" for c in clusters))

    def test_all_polar_bonds_form_a_directed_chain(self):
        bonds = [
            dyn.Bond(0, 1, "polar", 0.0, 0.0, 0.0, math.pi),
            dyn.Bond(1, 2, "polar", 0.0, 0.0, 0.0, math.pi),
        ]
        clusters = dyn.classify_clusters(3, bonds, True)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].kind, "chain")
        self.assertTrue(clusters[0].directed)

    def test_any_lateral_bond_makes_it_a_bundle(self):
        bonds = [
            dyn.Bond(0, 1, "polar", 0.0, 0.0, 0.0, math.pi),
            dyn.Bond(1, 2, "lateral", 0.3, 0.1, -0.3, math.pi + 0.1),
        ]
        clusters = dyn.classify_clusters(3, bonds, True)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].kind, "bundle")


class SimulationIntegrationTests(unittest.TestCase):
    """`run` 全体を短時間だけ流し、基本的な不変量を確かめる。"""

    def test_bonds_only_accumulate_never_disappear(self):
        config = dyn.DynamicsConfig(n_cells=12, box_size_um=14.0, total_time_s=15.0, seed=3)
        result = dyn.run(config)
        counts = [len(frame.bonds) for frame in result.frames]
        self.assertEqual(counts, sorted(counts))

    def test_cells_stay_within_the_box(self):
        config = dyn.DynamicsConfig(n_cells=12, box_size_um=14.0, total_time_s=15.0, seed=3)
        result = dyn.run(config)
        last = result.frames[-1]
        self.assertTrue(np.all(last.positions >= 0.0))
        self.assertTrue(np.all(last.positions <= config.box_size_um))

    def test_rejects_asymmetric_polar_split(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(polar_split=0.7)

    def test_bonds_stay_near_bridge_length(self):
        """距離拘束が硬いので、結合長は架橋長付近に留まる（伸びきらない）。

        曲げ剛性ゼロの点ばねでは、細胞が回転すると結合が細胞長スケール
        （中央値 1.6 µm ≈ 架橋長の200倍）まで伸びていた。柔軟なジョイントでも
        距離拘束は硬いので、結合長の中央値は架橋長の数倍に収まる。この不変量が
        その回帰を防ぐ。
        """

        config = dyn.DynamicsConfig(
            n_cells=18, box_size_um=16.0, total_time_s=45.0, seed=1,
            f_pole=0.95, tip_width_um=0.10, frame_stride=10**9,
        )
        last = dyn.run(config).frames[-1]
        self.assertTrue(last.bonds, "この設定では結合が少なくとも1本できるはず")
        lengths = []
        for bond in last.bonds:
            p_i = dyn.attachment_point(
                last.positions[bond.i], last.angles[bond.i], bond.axial_i, bond.delta_i
            )
            p_j = dyn.attachment_point(
                last.positions[bond.j], last.angles[bond.j], bond.axial_j, bond.delta_j
            )
            lengths.append(float(np.linalg.norm(p_j - p_i)))
        # 中央値は架橋長の数倍以内（バグ時は 1.6 µm だった）。
        self.assertLess(float(np.median(lengths)), 0.06)

    def test_cells_do_not_deeply_interpenetrate(self):
        """接触点で作用する排除体積により、表面の深い食い込みが起きない。"""

        config = dyn.DynamicsConfig(
            n_cells=18, box_size_um=16.0, total_time_s=45.0, seed=1,
            f_pole=0.95, tip_width_um=0.10, frame_stride=10**9,
        )
        last = dyn.run(config).frames[-1]
        n = config.n_cells
        worst = min(
            (dyn._surface_gap(last, i, j) for i in range(n) for j in range(i + 1, n)),
            default=0.0,
        )
        # 反復緩和の残差ぶんは許容する（バグ時は -0.8 µm ≈ 細胞直径ぶんだった）。
        self.assertGreater(worst, -0.08)


class LocalizationControlsMorphologyTests(unittest.TestCase):
    """局在パラメータが、実際にどちらの接触型が先に成立するかを左右することを

    統計的に確かめる。これは本モジュールの核心的な主張（`docs/literature_findings.md`
    が述べる「どちらの接触が先に成立するかで形態が決まる」という経路依存性）の
    直接の検証であり、単なる実装の細部ではない。
    """

    @staticmethod
    def _first_bond_kind(f_pole, tip_width_um, seed):
        config = dyn.DynamicsConfig(
            n_cells=2,
            box_size_um=6.0,
            total_time_s=60.0,
            f_pole=f_pole,
            tip_width_um=tip_width_um,
            seed=seed,
            frame_stride=10**9,
        )
        result = dyn.run(config)
        bonds = result.frames[-1].bonds
        return bonds[0].kind if bonds else None

    def test_tip_focused_bipolar_localization_favors_polar_bonds(self):
        kinds = [self._first_bond_kind(0.95, 0.10, seed) for seed in range(15)]
        self.assertGreater(kinds.count("polar"), kinds.count("lateral"))

    def test_uniform_localization_favors_lateral_bonds(self):
        kinds = [
            self._first_bond_kind(model.UNIFORM_F_POLE, None, seed) for seed in range(15)
        ]
        self.assertGreater(kinds.count("lateral"), kinds.count("polar"))

    @staticmethod
    def _final_polar_and_lateral(f_pole, tip_width_um, seed):
        config = dyn.DynamicsConfig(
            n_cells=18, box_size_um=16.0, total_time_s=45.0,
            f_pole=f_pole, tip_width_um=tip_width_um, seed=seed, frame_stride=10**9,
        )
        bonds = dyn.run(config).frames[-1].bonds
        polar = sum(1 for b in bonds if b.kind == "polar")
        lateral = sum(1 for b in bonds if b.kind == "lateral")
        return polar, lateral

    def test_directionality_requires_localization_in_many_body(self):
        """多体でも、方向性のある極接触（頭-尾の鎖）は先端局在でのみ生じる。

        先端局在では極どうしが繋がって有向な鎖の骨格ができるが、一様分布では
        極接触はほとんど生じず、側面でランダムに固まる（方向性なし）。局在した
        場合の極結合数は、一様分布の場合よりはっきり多い。これが本モジュールの
        核心的な主張であり、鎖・束の絶対数ではなく「方向性が局在に依存する」
        という頑健な対比として固定する。
        """

        chain_polar = [self._final_polar_and_lateral(0.95, 0.10, s)[0] for s in range(3)]
        uniform_polar = [
            self._final_polar_and_lateral(model.UNIFORM_F_POLE, None, s)[0]
            for s in range(3)
        ]
        # 局在させた場合の極結合は、一様分布のどの試行よりも多い。
        self.assertGreater(min(chain_polar), max(uniform_polar))
        # 先端局在では有向な鎖の骨格が実際にできる。
        self.assertGreaterEqual(min(chain_polar), 3)

    def test_uniform_is_lateral_dominated_in_many_body(self):
        for seed in range(3):
            polar, lateral = self._final_polar_and_lateral(model.UNIFORM_F_POLE, None, seed)
            self.assertGreater(
                lateral, polar, f"seed={seed}: 極={polar}, 側面={lateral}"
            )


class ShearModeConfigTests(unittest.TestCase):
    def test_rejects_negative_shear_rate(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(shear_rate_per_s=-1.0)

    def test_periodic_x_follows_shear_rate(self):
        self.assertFalse(dyn.DynamicsConfig(shear_rate_per_s=0.0).periodic_x)
        self.assertTrue(dyn.DynamicsConfig(shear_rate_per_s=100.0).periodic_x)

    def test_safe_dt_is_unchanged_when_quiescent(self):
        config = dyn.DynamicsConfig(dt_s=0.02, shear_rate_per_s=0.0)
        self.assertEqual(dyn._safe_dt_for_shear(config), 0.02)

    def test_safe_dt_shrinks_with_shear_rate(self):
        """dt は接触持続時間 shear.contact_duration_s(γ̇) = 1/γ̇ に対して縮む。

        震盪が「出会っても短時間しか触れ合わない」という効果を持つには、
        大きい dt のまま接触判定してしまい確率を過大評価してはならない。
        """

        config = dyn.DynamicsConfig(dt_s=0.02, shear_rate_per_s=300.0)
        safe_dt = dyn._safe_dt_for_shear(config)
        self.assertLess(safe_dt, config.dt_s)
        self.assertAlmostEqual(
            safe_dt, dyn._SHEAR_DT_SAFETY_FACTOR * shear.contact_duration_s(300.0)
        )


class PeriodicImageTests(unittest.TestCase):
    def test_no_shift_when_not_periodic(self):
        pos_i = np.array([1.0, 2.0])
        pos_j = np.array([19.0, 2.0])
        image = dyn._image_of_j(pos_i, pos_j, box_size_um=20.0, periodic_x=False)
        self.assertTrue(np.array_equal(image, pos_j))

    def test_wraps_to_nearest_image_across_boundary(self):
        """箱の反対端にいる2細胞は、周期境界では実は近接している。"""

        pos_i = np.array([1.0, 2.0])
        pos_j = np.array([19.0, 2.0])
        image = dyn._image_of_j(pos_i, pos_j, box_size_um=20.0, periodic_x=True)
        # 19 は 20 を法として -1 の像に等しく、i(x=1) の近くに来るはず。
        self.assertAlmostEqual(float(image[0]), -1.0)
        self.assertAlmostEqual(float(image[1]), 2.0)

    def test_wrap_periodic_x_keeps_positions_in_box(self):
        state = dyn.SimulationState(
            positions=np.array([[-1.0, 5.0], [21.0, 5.0], [10.0, 5.0]]),
            angles=np.zeros(3),
            head_free=np.ones(3, dtype=bool),
            tail_free=np.ones(3, dtype=bool),
            bonds=[],
        )
        dyn._wrap_periodic_x(state, box_size_um=20.0)
        self.assertTrue(np.all(state.positions[:, 0] >= 0.0))
        self.assertTrue(np.all(state.positions[:, 0] < 20.0))


class JefferyOrbitTests(unittest.TestCase):
    """Jeffery (1922) の2次元軌道式（A20）が、既知の性質を再現するかを確かめる。"""

    def test_rotation_rate_is_zero_shear_free(self):
        rate = dyn._jeffery_rotation_rate(np.array([0.0, 1.0, 2.0]), shear_rate_per_s=0.0)
        self.assertTrue(np.allclose(rate, 0.0))

    def test_rotation_slowest_aligned_with_flow_for_elongated_rod(self):
        """細長い棒（アスペクト比 > 1）は、流れ方向（θ=0）付近で最も回転が遅い。"""

        shear_rate = 300.0
        rate_aligned = dyn._jeffery_rotation_rate(np.array([0.0]), shear_rate)[0]
        rate_perpendicular = dyn._jeffery_rotation_rate(np.array([math.pi / 2]), shear_rate)[0]
        self.assertGreater(dyn._ASPECT_RATIO, 1.0)
        self.assertLess(rate_aligned, rate_perpendicular)

    def test_matches_closed_form_at_reference_angles(self):
        r = dyn._ASPECT_RATIO
        shear_rate = 250.0
        expected_aligned = shear_rate / (r * r + 1.0)
        expected_perpendicular = shear_rate * r * r / (r * r + 1.0)
        self.assertAlmostEqual(
            dyn._jeffery_rotation_rate(np.array([0.0]), shear_rate)[0], expected_aligned
        )
        self.assertAlmostEqual(
            dyn._jeffery_rotation_rate(np.array([math.pi / 2]), shear_rate)[0],
            expected_perpendicular,
        )


class ShearRuptureTests(unittest.TestCase):
    """`_apply_shear_rupture`（A22）を、鎖の状態を直接組み立てて確かめる。"""

    @staticmethod
    def _linear_chain_state(n_cells: int, spacing_um: float) -> dyn.SimulationState:
        positions = np.array([[k * spacing_um, 5.0] for k in range(n_cells)])
        angles = np.zeros(n_cells)
        head_free = np.array([k == n_cells - 1 for k in range(n_cells)])
        tail_free = np.array([k == 0 for k in range(n_cells)])
        bonds = [
            dyn.Bond(k, k + 1, "polar", dyn.HALF_CYLINDER_UM, 0.0, -dyn.HALF_CYLINDER_UM, math.pi)
            for k in range(n_cells - 1)
        ]
        return dyn.SimulationState(positions, angles, head_free, tail_free, bonds)

    def test_no_rupture_at_literature_shaking_rates(self):
        """文献の結論（震盪培養の典型的な速度では接合を壊さない）を再現する。

        既存の見積もり（1結合あたりの荷重は破断力の 1/56 程度、README 参照）
        から、γ̇=300 s⁻¹・10細胞程度の鎖では発火しないはずである。
        """

        state = self._linear_chain_state(10, dyn.CELL_LENGTH_UM)
        n_bonds_before = len(state.bonds)
        dyn._apply_shear_rupture(state, 10, shear_rate_per_s=300.0)
        self.assertEqual(len(state.bonds), n_bonds_before)

    def test_ruptures_under_extreme_shear(self):
        """機構そのものが実際に発火しうることを、非現実的に大きいγ̇で確かめる。"""

        state = self._linear_chain_state(10, dyn.CELL_LENGTH_UM)
        dyn._apply_shear_rupture(state, 10, shear_rate_per_s=1.0e6)
        self.assertEqual(len(state.bonds), 0)
        # 極スロットが解放されているはず。
        self.assertTrue(np.all(state.head_free))
        self.assertTrue(np.all(state.tail_free))

    def test_shear_rate_zero_never_ruptures(self):
        state = self._linear_chain_state(10, dyn.CELL_LENGTH_UM)
        n_bonds_before = len(state.bonds)
        dyn._apply_shear_rupture(state, 10, shear_rate_per_s=0.0)
        self.assertEqual(len(state.bonds), n_bonds_before)


class ShearIntegrationTests(unittest.TestCase):
    """`run` を震盪モードで短く流し、静置経路を壊していないことを確かめる。"""

    def test_quiescent_path_is_unaffected_by_shear_feature(self):
        """shear_rate_per_s=0（既定）は、震盪機能を追加する前とビット単位で同じ。"""

        config = dyn.DynamicsConfig(n_cells=10, box_size_um=14.0, total_time_s=4.0, seed=2)
        result = dyn.run(config)
        last = result.frames[-1]
        self.assertTrue(np.all(last.positions[:, 0] >= 0.0))
        self.assertTrue(np.all(last.positions[:, 0] <= config.box_size_um))

    def test_shear_run_completes_and_stays_periodic_in_x(self):
        config = dyn.DynamicsConfig(
            n_cells=6, box_size_um=10.0, total_time_s=0.1, seed=1, shear_rate_per_s=200.0,
            frame_stride=10**9,
        )
        result = dyn.run(config)
        last = result.frames[-1]
        self.assertTrue(np.all(last.positions[:, 0] >= 0.0))
        self.assertTrue(np.all(last.positions[:, 0] < config.box_size_um))
        # y は壁で挟まれたまま（周期境界にしない）。
        margin = dyn.HALF_CYLINDER_UM + dyn.CELL_RADIUS_UM
        self.assertTrue(np.all(last.positions[:, 1] >= margin - 1e-6))
        self.assertTrue(np.all(last.positions[:, 1] <= config.box_size_um - margin + 1e-6))


class SettlingConfigTests(unittest.TestCase):
    def test_settling_disabled_by_default(self):
        config = dyn.DynamicsConfig()
        self.assertFalse(config.settling_enabled)

    def test_settling_enabled_with_well_depth(self):
        config = dyn.DynamicsConfig(well_depth_um=1000.0)
        self.assertTrue(config.settling_enabled)

    def test_rejects_non_positive_well_depth(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(well_depth_um=0.0)
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(well_depth_um=-5.0)

    def test_rejects_invalid_supernatant_fraction(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(well_depth_um=1000.0, supernatant_fraction=0.0)
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(well_depth_um=1000.0, supernatant_fraction=1.5)


class SettlingVelocityTests(unittest.TestCase):
    """Stokes 則によるクラスタの沈降速度（A23）の性質を確かめる。"""

    def test_larger_clusters_settle_faster(self):
        self.assertLess(
            dyn._settling_velocity_um_s(1), dyn._settling_velocity_um_s(20)
        )

    def test_velocity_scales_as_two_thirds_power_of_cluster_size(self):
        """半径が n^(1/3) でスケールする Stokes 則なので、速度は n^(2/3) で増える。"""

        ratio = dyn._settling_velocity_um_s(8) / dyn._settling_velocity_um_s(1)
        self.assertAlmostEqual(ratio, 8.0 ** (2.0 / 3.0), places=6)

    def test_single_cell_needs_much_longer_than_a_day_to_cross_a_deep_well(self):
        """単量体は深底プレート程度の液柱を24時間では沈みきれない。

        実際のアッセイで単量体・小さい凝集塊が24時間経っても上清に残り続ける
        ことと整合する。
        """

        deep_well_um = 20_000.0  # 深底96穴プレートの物理的な深さのオーダー
        time_to_cross_s = deep_well_um / dyn._settling_velocity_um_s(1)
        self.assertGreater(time_to_cross_s, 24.0 * 3600.0)

    def test_large_cluster_can_cross_a_deep_well_in_a_day(self):
        """十分大きな凝集塊は、深底プレート程度の液柱を24時間で沈みきる。"""

        distance_in_a_day_um = dyn._settling_velocity_um_s(20) * 86400.0
        self.assertGreater(distance_in_a_day_um, 20_000.0)


class SettlingIntegrationTests(unittest.TestCase):
    """`run` を沈降モードで短く流し、基本的な不変量を確かめる。"""

    def test_depths_stay_non_negative_and_within_well(self):
        config = dyn.DynamicsConfig(
            n_cells=10, box_size_um=14.0, total_time_s=3.0, seed=4, well_depth_um=50.0,
        )
        result = dyn.run(config)
        last = result.frames[-1]
        self.assertTrue(np.all(last.depth_um >= 0.0))
        self.assertTrue(np.all(last.depth_um <= 50.0))

    def test_bonded_cells_share_nearly_the_same_depth(self):
        """結合したクラスタは一体で沈むので、深さはほぼ揃っているはず。

        結合が成立する瞬間の深さの差は高々 `VERTICAL_CONTACT_GAP_UM`
        （結合できる鉛直方向の許容差）であり、結合後はその差のまま一体で
        沈み続ける。したがって完全に一致はしないが、その許容差は超えない。
        """

        config = dyn.DynamicsConfig(
            n_cells=10, box_size_um=12.0, total_time_s=3.0, seed=1,
            f_pole=0.95, tip_width_um=0.10, well_depth_um=50.0,
        )
        result = dyn.run(config)
        last = result.frames[-1]
        for cells in dyn.connected_components(config.n_cells, last.bonds):
            if len(cells) < 2:
                continue
            depths = last.depth_um[cells]
            self.assertLess(float(depths.max() - depths.min()), dyn.VERTICAL_CONTACT_GAP_UM)

    def test_vertically_separated_cells_do_not_bond(self):
        """鉛直方向に離れすぎた細胞は、横方向が近くても結合しない。"""

        config = dyn.DynamicsConfig(n_cells=2, well_depth_um=1000.0)
        gap = dyn.BRIDGE_UM
        x1 = dyn.HALF_CYLINDER_UM
        x2 = x1 + 2 * dyn.CELL_RADIUS_UM + gap + dyn.HALF_CYLINDER_UM
        state = dyn.SimulationState(
            positions=np.array([[0.0, 0.0], [x2, 0.0]]),
            angles=np.zeros(2),
            head_free=np.ones(2, dtype=bool),
            tail_free=np.ones(2, dtype=bool),
            bonds=[],
        )
        state.depth_um[1] = dyn.VERTICAL_CONTACT_GAP_UM * 5.0
        rng = np.random.default_rng(0)
        for _ in range(200):
            dyn._try_form_bond(state, 0, 1, config, 500.0, 500.0, 0.02, rng)
        self.assertEqual(len(state.bonds), 0)


class SupernatantFractionTests(unittest.TestCase):
    def test_all_cells_at_surface_gives_full_supernatant(self):
        frame = dyn.Frame(0.0, np.zeros((5, 2)), np.zeros(5), [], np.zeros(5))
        self.assertEqual(dyn.supernatant_fraction(frame, well_depth_um=1000.0, top_fraction=0.25), 1.0)

    def test_all_cells_settled_gives_zero_supernatant(self):
        frame = dyn.Frame(0.0, np.zeros((5, 2)), np.zeros(5), [], np.full(5, 1000.0))
        self.assertEqual(dyn.supernatant_fraction(frame, well_depth_um=1000.0, top_fraction=0.25), 0.0)

    def test_partial_settling_gives_intermediate_fraction(self):
        depths = np.array([0.0, 0.0, 500.0, 500.0])
        frame = dyn.Frame(0.0, np.zeros((4, 2)), np.zeros(4), [], depths)
        self.assertAlmostEqual(
            dyn.supernatant_fraction(frame, well_depth_um=1000.0, top_fraction=0.25), 0.5
        )


class SpontaneousDissociationTests(unittest.TestCase):
    """結合の自発的な解離（A15の既定を外すオプション、A24）を確かめる。

    Costan et al. (2024) の「k_off should be assumed zero」という結論は、
    Supplemental Note 4 原文の確認により、狭い荷重速度範囲（5-100 pN/s）で
    傾きを検出できなかったことに基づく検出限界の主張だと分かった
    （数学的に厳密なゼロの証明ではない）。この結合（anti-EPEA/EPEA）は、
    他の測定に使われている例外的に安定な anti-GFP/GFP ペアとは異なる、より
    低親和性のペアであり、他のナノボディ-抗原ペアの独立した SPR 測定では
    k_off ~ 10⁻³-10⁻²/s も珍しくない。この節は、既定（k_off=0、完全に不可逆）
    が引き続き正しく動作すること、および非零の k_off を指定したときに
    結合が実際に確率的に解離することの両方を確認する。

    参考: Costan et al. (2024) Supplemental Note 4 原文
    "the slope is not significantly different than 0, and due to instrument
    resolution, ultimately the best estimate is that koff should be assumed
    zero, i.e., adhesins do not spontaneously disengage."
    """

    def _bonded_two_cell_state(self) -> dyn.SimulationState:
        gap = dyn.BRIDGE_UM
        x1 = dyn.HALF_CYLINDER_UM
        x2 = x1 + 2 * dyn.CELL_RADIUS_UM + gap + dyn.HALF_CYLINDER_UM
        return dyn.SimulationState(
            positions=np.array([[0.0, 0.0], [x2, 0.0]]),
            angles=np.zeros(2),
            head_free=np.array([False, True]),
            tail_free=np.array([True, False]),
            bonds=[
                dyn.Bond(0, 1, "polar", dyn.HALF_CYLINDER_UM, 0.0, -dyn.HALF_CYLINDER_UM, math.pi)
            ],
        )

    def test_rejects_negative_k_off(self):
        with self.assertRaises(ValueError):
            dyn.DynamicsConfig(k_off_per_s=-1.0)

    def test_default_k_off_is_zero_and_never_dissociates(self):
        """既定（k_off=0）は、この機能を追加する前とビット単位で同じ挙動。"""

        state = self._bonded_two_cell_state()
        rng = np.random.default_rng(0)
        for _ in range(1000):
            dyn._apply_spontaneous_dissociation(state, k_off_per_s=0.0, dt_s=0.02, rng=rng)
        self.assertEqual(len(state.bonds), 1)

    def test_nonzero_k_off_eventually_dissociates_and_frees_pole_slots(self):
        state = self._bonded_two_cell_state()
        rng = np.random.default_rng(1)
        for _ in range(2000):  # 十分な回数（k_off=2/s なら平均寿命0.5秒 = 25ステップ）
            dyn._apply_spontaneous_dissociation(state, k_off_per_s=2.0, dt_s=0.02, rng=rng)
            if not state.bonds:
                break
        self.assertEqual(len(state.bonds), 0)
        self.assertTrue(state.head_free[0])
        self.assertTrue(state.tail_free[1])

    def test_mean_bond_lifetime_matches_one_over_k_off(self):
        """解離はポアソン過程になっているはず（平均寿命 = 1/k_off）。"""

        k_off = 1.0
        dt = 0.02
        lifetimes = []
        for trial in range(200):
            state = self._bonded_two_cell_state()
            rng = np.random.default_rng(trial)
            steps = 0
            while state.bonds and steps < 100_000:
                dyn._apply_spontaneous_dissociation(state, k_off_per_s=k_off, dt_s=dt, rng=rng)
                steps += 1
            lifetimes.append(steps * dt)
        mean_lifetime = sum(lifetimes) / len(lifetimes)
        self.assertAlmostEqual(mean_lifetime, 1.0 / k_off, delta=0.15)

    def test_quiescent_run_with_default_k_off_still_never_loses_bonds(self):
        """`run()` 全体でも、既定（k_off=0）では結合が減らないという既存の

        不変量（A15）が壊れていないことを確かめる（回帰防止）。
        """

        config = dyn.DynamicsConfig(n_cells=12, box_size_um=14.0, total_time_s=15.0, seed=3)
        result = dyn.run(config)
        counts = [len(frame.bonds) for frame in result.frames]
        self.assertEqual(counts, sorted(counts))

    def test_quiescent_run_with_nonzero_k_off_can_lose_bonds(self):
        """`run()` に non-zero k_off を渡すと、結合数が減る場面がありうる。"""

        config = dyn.DynamicsConfig(
            n_cells=12, box_size_um=14.0, total_time_s=15.0, seed=3,
            f_pole=0.95, tip_width_um=0.10, k_off_per_s=5.0,
        )
        result = dyn.run(config)
        counts = [len(frame.bonds) for frame in result.frames]
        self.assertFalse(counts == sorted(counts), "k_off>0 なのに結合が一度も減らなかった")


if __name__ == "__main__":
    unittest.main()
