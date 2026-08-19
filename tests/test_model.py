"""モデルの振る舞いのテスト。

形態判定、スケーリング則、パラメータの分離（発現量は強さを、局在は形を決める）を
固定する。
"""

import unittest

from adhesion import literature as lit
from adhesion.model import (
    RUPTURE_FORCE_PN,
    UNIFORM_F_POLE,
    Morphology,
    morphology_threshold,
    predict,
)

N_MID = float(lit.ADHESINS_AT_MID_INDUCTION)


class MorphologyTests(unittest.TestCase):
    def test_tip_focused_bipolar_localization_gives_chains(self):
        result = predict(N_MID, 0.95, 0.5, tip_width_um=0.10)
        self.assertEqual(result.morphology, Morphology.CHAIN)

    def test_polar_enrichment_alone_is_not_enough(self):
        """極半球に広く集めるだけでは束のまま。先端に絞る必要がある。"""
        result = predict(N_MID, 0.95, 0.5, tip_width_um=None)
        self.assertEqual(result.morphology, Morphology.BUNDLE)

    def test_unipolar_localization_stops_at_dimers(self):
        """単極では1細胞1接続しか作れず、鎖に伸びない。"""
        result = predict(N_MID, 0.95, 1.0, tip_width_um=0.10)
        self.assertEqual(result.morphology, Morphology.DIMER)
        self.assertEqual(result.bonds_light_pole, 0.0)
        self.assertGreater(result.bonds_heavy_pole, 0.0)

    def test_unipolar_junction_is_stronger_than_bipolar(self):
        """単極は接合1つあたりでは強い。強い二量体か弱い鎖かのトレードオフ。"""
        bipolar = predict(N_MID, 0.95, 0.5, tip_width_um=0.10)
        unipolar = predict(N_MID, 0.95, 1.0, tip_width_um=0.10)
        self.assertGreater(unipolar.bonds_heavy_pole, bipolar.bonds_heavy_pole)

    def test_no_expression_gives_no_adhesion(self):
        self.assertEqual(predict(0.0, UNIFORM_F_POLE).morphology, Morphology.NO_ADHESION)

    def test_asymmetric_split_can_still_form_chains(self):
        """薄い側の極にも十分な分子が残っていれば鎖になる。"""
        result = predict(N_MID, 0.95, 0.75, tip_width_um=0.10)
        self.assertEqual(result.morphology, Morphology.CHAIN)
        self.assertGreater(result.bonds_heavy_pole, result.bonds_light_pole)


class ParameterSeparationTests(unittest.TestCase):
    """発現量は接合の強さを、局在は形を決める。両者は干渉しない。"""

    def test_morphology_threshold_does_not_depend_on_expression(self):
        threshold = morphology_threshold(0.20)
        self.assertIsNotNone(threshold)
        for n_total in (2000.0, 7300.0, 15300.0):
            with self.subTest(n_total=n_total):
                below = predict(n_total, threshold - 0.05, 0.5, 0.20)
                above = predict(n_total, threshold + 0.05, 0.5, 0.20)
                self.assertEqual(below.morphology, Morphology.BUNDLE)
                self.assertEqual(above.morphology, Morphology.CHAIN)

    def test_bond_count_is_linear_in_expression(self):
        single = predict(5000.0, 0.95, 0.5, 0.10)
        double = predict(10000.0, 0.95, 0.5, 0.10)
        self.assertAlmostEqual(
            double.bonds_heavy_pole / single.bonds_heavy_pole, 2.0, places=6
        )
        self.assertAlmostEqual(
            double.bonds_side_by_side / single.bonds_side_by_side, 2.0, places=6
        )

    def test_junction_force_is_bond_count_times_rupture_force(self):
        result = predict(N_MID, 0.95, 0.5, 0.10)
        self.assertAlmostEqual(
            result.junction_force_pn,
            result.bonds_per_junction * RUPTURE_FORCE_PN,
            places=9,
        )


class ThresholdTests(unittest.TestCase):
    def test_tighter_focusing_lowers_the_required_polar_fraction(self):
        widths = [0.30, 0.25, 0.20, 0.15, 0.10]
        thresholds = [morphology_threshold(w) for w in widths]
        self.assertTrue(all(t is not None for t in thresholds))
        self.assertEqual(thresholds, sorted(thresholds, reverse=True))

    def test_threshold_without_focusing_is_practically_unreachable(self):
        threshold = morphology_threshold(None)
        self.assertTrue(threshold is None or threshold > 0.9)

    def test_prediction_at_threshold_is_balanced(self):
        threshold = morphology_threshold(0.20)
        result = predict(N_MID, threshold, 0.5, 0.20)
        ratio = result.bonds_heavy_pole / result.bonds_side_by_side
        self.assertAlmostEqual(ratio, 1.0, delta=0.02)


class DensityTests(unittest.TestCase):
    def test_uniform_display_stays_far_below_membrane_capacity(self):
        """一様分布なら外膜のタンパク収容量（1e5 /um^2 のオーダー）に遠く届かない。"""
        result = predict(float(lit.ADHESINS_AT_MAX_INDUCTION), UNIFORM_F_POLE)
        self.assertLess(result.peak_density_per_um2, 1e4)

    def test_tight_focusing_approaches_membrane_capacity(self):
        """極先端に強く絞ると外膜の収容量に迫る。局在には物理的な上限がある。"""
        result = predict(N_MID, 0.95, 0.5, tip_width_um=0.05)
        self.assertGreater(result.peak_density_per_um2, 1e5)


if __name__ == "__main__":
    unittest.main()
