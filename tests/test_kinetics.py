"""荷重下での結合寿命のテスト。

値は anti-GFP / GFP ペアのもので、破断力を測ったペアとは異なる。
ここで固定するのは桁と単調性であって、本系の実測値ではない。
"""

import unittest

from adhesion import kinetics

#: 参考用。出典（Kirchhofer 2010）は未確認であることを承知で使う。
K = kinetics.UNVERIFIED_OFF_RATE_PER_S
from adhesion.model import predict
from adhesion.shear import chain_tension_pn


class BellModelTests(unittest.TestCase):
    def test_zero_load_lifetime_is_about_two_hours(self):
        self.assertAlmostEqual(kinetics.bond_lifetime_s(0.0, K) / 3600, 1.9, delta=0.2)

    def test_lifetime_decreases_with_load(self):
        lifetimes = [kinetics.bond_lifetime_s(f, K) for f in (0.0, 1.0, 5.0, 10.0)]
        self.assertEqual(lifetimes, sorted(lifetimes, reverse=True))

    def test_force_scale_is_a_few_piconewtons(self):
        """荷重がこのスケールを大きく下回るあいだ、寿命はほぼ変わらない。"""
        self.assertAlmostEqual(kinetics.FORCE_SCALE_PN, 2.6, delta=0.2)

    def test_negative_load_is_rejected(self):
        with self.assertRaises(ValueError):
            kinetics.bond_lifetime_s(-1.0, K)


class ShearDoesNotShortenBondsTests(unittest.TestCase):
    def test_shear_barely_changes_bond_lifetime(self):
        """震盪培養のせん断では結合寿命がほとんど変わらないこと。"""
        result = predict(7300, 0.95, 0.5, tip_width_um=0.10)
        load = chain_tension_pn(10, 300) / result.bonds_heavy_pole
        self.assertGreater(kinetics.lifetime_reduction(load), 0.95)


class DynamicForceSpectroscopyTests(unittest.TestCase):
    def test_rupture_force_grows_with_loading_rate(self):
        forces = [kinetics.rupture_force_pn(r, K) for r in (10, 100, 1000, 10000)]
        self.assertEqual(forces, sorted(forces))

    def test_reproduces_the_published_afm_range(self):
        """Klamecka らの AFM 荷重速度域で 28-45 pN 相当になること。"""
        self.assertAlmostEqual(kinetics.rupture_force_pn(20, K), 28, delta=3)
        self.assertAlmostEqual(kinetics.rupture_force_pn(1e4, K), 44, delta=3)

    def test_invalid_loading_rate_is_rejected(self):
        with self.assertRaises(ValueError):
            kinetics.rupture_force_pn(0.0, K)


if __name__ == "__main__":
    unittest.main()
