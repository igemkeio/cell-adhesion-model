"""全ての図・動画で共通の配色（`scripts/figstyle.py`）のテスト。

極（青）と側面（オレンジ）が、色覚の違いがあっても区別でき、文字用の濃い版が読めることを固定する。
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import figstyle  # noqa: E402

# Machado et al. (2009)、重症度 1.0。線形 RGB 上の色覚シミュレーション。
_CVD = {
    "normal": np.eye(3),
    "protan": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216],
                        [-0.003882, -0.048116, 1.051998]]),
    "deutan": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413],
                        [-0.011820, 0.042940, 0.968881]]),
    "tritan": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602],
                        [0.004733, 0.691367, 0.303900]]),
}


def _linear(hex_color: str) -> np.ndarray:
    rgb = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
    return np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)


def _oklab(linear_rgb: np.ndarray) -> np.ndarray:
    m1 = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                   [0.2119034982, 0.6806995451, 0.1073969566],
                   [0.0883024619, 0.2817188376, 0.6299787005]])
    m2 = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                   [1.9779984951, -2.4285922050, 0.4505937099],
                   [0.0259040371, 0.7827717662, -0.8086757660]])
    return m2 @ np.cbrt(m1 @ linear_rgb)


def delta_e(a: str, b: str, vision: str) -> float:
    """OKLab ΔE×100。色覚シミュレーション後の2色の距離。"""

    matrix = _CVD[vision]
    return 100.0 * float(np.linalg.norm(
        _oklab(np.clip(matrix @ _linear(a), 0, 1)) - _oklab(np.clip(matrix @ _linear(b), 0, 1))
    ))


def contrast(foreground: str, background: str = "#ffffff") -> float:
    def luminance(color: str) -> float:
        r, g, b = _linear(color)
        return 0.2126 * r + 0.7152 * g + 0.0722 * b

    high, low = sorted([luminance(foreground), luminance(background)], reverse=True)
    return (high + 0.05) / (low + 0.05)


class PaletteTests(unittest.TestCase):
    def test_polar_and_lateral_are_distinguishable_for_every_vision_type(self):
        for vision in _CVD:
            with self.subTest(vision=vision):
                self.assertGreaterEqual(delta_e(figstyle.POLAR, figstyle.LATERAL, vision), 8.0)

    def test_cell_fills_in_the_movies_are_distinguishable(self):
        """灰（単独）・青（鎖）・オレンジ（束）の3つが、どの色覚でも互いに区別できる。"""

        fills = [figstyle.MONOMER, figstyle.POLAR, figstyle.LATERAL]
        for vision in _CVD:
            for i, a in enumerate(fills):
                for b in fills[i + 1:]:
                    with self.subTest(vision=vision, a=a, b=b):
                        self.assertGreaterEqual(delta_e(a, b, vision), 8.0)

    def test_text_variants_are_readable_on_white(self):
        self.assertGreaterEqual(contrast(figstyle.POLAR_TEXT), 4.5)
        self.assertGreaterEqual(contrast(figstyle.LATERAL_TEXT), 4.5)
        self.assertGreaterEqual(contrast(figstyle.INK_2), 4.5)

    def test_text_variants_keep_the_hue_of_their_marks(self):
        """文字用の濃い版が、元の色よりも濃い（輝度が低い）こと。"""

        for mark, text in ((figstyle.POLAR, figstyle.POLAR_TEXT), (figstyle.LATERAL, figstyle.LATERAL_TEXT)):
            self.assertLess(sum(_linear(text)), sum(_linear(mark)))

    def test_the_two_roles_use_two_different_colors(self):
        self.assertNotEqual(figstyle.POLAR, figstyle.LATERAL)


class MovieColoringTests(unittest.TestCase):
    """動画の細胞の塗りは、その細胞がもつ結合の種類で決まる（かたまり全体ではなく細胞ごと）。"""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import simulate_chain_formation as movie
        from adhesion import dynamics as dyn

        cls.movie, cls.dyn = movie, dyn

    def _frame(self, bonds):
        n = 4
        return self.dyn.Frame(0.0, np.zeros((n, 2)), np.zeros(n), bonds)

    def _bond(self, i, j, kind):
        return self.dyn.Bond(i, j, kind, 0.0, 0.0, 0.0, 0.0)

    def test_cells_without_bonds_are_gray(self):
        colors = self.movie.cell_colors_for_frame(4, self._frame([]))
        self.assertEqual(colors, [figstyle.MONOMER] * 4)

    def test_cells_with_only_polar_bonds_are_blue(self):
        colors = self.movie.cell_colors_for_frame(4, self._frame([self._bond(0, 1, "polar")]))
        self.assertEqual(colors, [figstyle.POLAR, figstyle.POLAR, figstyle.MONOMER, figstyle.MONOMER])

    def test_one_lateral_bond_recolors_only_its_two_cells(self):
        """鎖の途中に側面接触が1本あっても、塗り替わるのはその2細胞だけ（鎖が束に見えてしまわない）。"""

        bonds = [self._bond(0, 1, "polar"), self._bond(1, 2, "polar"), self._bond(2, 3, "lateral")]
        colors = self.movie.cell_colors_for_frame(4, self._frame(bonds))
        self.assertEqual(colors, [figstyle.POLAR, figstyle.POLAR, figstyle.LATERAL, figstyle.LATERAL])

    def test_a_cell_with_both_kinds_is_lateral(self):
        bonds = [self._bond(0, 1, "polar"), self._bond(1, 2, "lateral")]
        colors = self.movie.cell_colors_for_frame(4, self._frame(bonds))
        self.assertEqual(colors[1], figstyle.LATERAL)

    def test_bond_lines_use_the_same_colors_as_the_cells(self):
        self.assertEqual(self.movie.POLAR_BOND_COLOR, figstyle.POLAR)
        self.assertEqual(self.movie.LATERAL_BOND_COLOR, figstyle.LATERAL)


class FocusViewTests(unittest.TestCase):
    """振盪培養の動画で、最初の 1/5 の時点で最も長く結合している群を追う視点。"""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        import types

        import simulate_chain_formation as movie
        from adhesion import dynamics as dyn

        cls.movie, cls.dyn, cls.types = movie, dyn, types

    def _result(self, bonds_by_frame, n_cells=6, box=16.0, positions_by_frame=None):
        """フレームごとの結合の組（i, j）の一覧から、最小限の結果オブジェクトを作る。"""

        frames = []
        for k, pairs in enumerate(bonds_by_frame):
            bonds = [self.dyn.Bond(i, j, "polar", 0.0, 0.0, 0.0, 0.0) for i, j in pairs]
            positions = np.zeros((n_cells, 2)) if positions_by_frame is None else positions_by_frame[k]
            frames.append(self.dyn.Frame(float(k), positions, np.zeros(n_cells), bonds))
        config = self.types.SimpleNamespace(n_cells=n_cells, box_size_um=box)
        return self.types.SimpleNamespace(frames=frames, config=config)

    def test_picks_the_group_with_the_oldest_bond(self):
        bonds = [[(0, 1)], [(0, 1), (2, 3), (3, 4)]] + [[(0, 1), (2, 3), (3, 4)]] * 9
        group = self.movie.find_focus_group(self._result(bonds), "oldest")
        self.assertEqual(group["members"], [0, 1])
        self.assertEqual(group["reference_index"], 2)
        self.assertAlmostEqual(group["age_s"], 2.0)

    def test_largest_rule_picks_the_biggest_group_instead(self):
        bonds = [[(0, 1)], [(0, 1), (2, 3), (3, 4)]] + [[(0, 1), (2, 3), (3, 4)]] * 9
        group = self.movie.find_focus_group(self._result(bonds), "largest")
        self.assertEqual(group["members"], [2, 3, 4])

    def test_a_bond_that_broke_and_re_formed_is_counted_from_the_re_formation(self):
        # 細胞 0-1 の結合は t=0 にでき、t=1 に切れ、t=2 につなぎ直された。細胞 2-3 は t=1 から続く。
        bonds = [[(0, 1)], [(2, 3)]] + [[(0, 1), (2, 3)]] * 9
        group = self.movie.find_focus_group(self._result(bonds), "oldest")
        self.assertEqual(group["members"], [2, 3])

    def test_returns_none_when_nothing_is_bonded_at_the_reference_time(self):
        bonds = [[]] * 11
        self.assertIsNone(self.movie.find_focus_group(self._result(bonds), "oldest"))

    def test_ties_go_to_the_larger_group(self):
        bonds = [[(0, 1), (2, 3), (3, 4)]] * 11
        group = self.movie.find_focus_group(self._result(bonds), "oldest")
        self.assertEqual(group["members"], [2, 3, 4])

    def test_camera_follows_a_group_that_straddles_the_periodic_boundary(self):
        """箱の端をまたぐ群の重心は、反対側（箱の中央）ではなく端にある。"""

        positions = np.zeros((6, 2))
        positions[0] = (0.5, 3.0)
        positions[1] = (15.5, 5.0)
        result = self._result([[(0, 1)]], positions_by_frame=[positions])
        cx, cy = self.movie.focus_camera(result, [0, 1])[0]
        wrapped = (cx + 8.0) % 16.0 - 8.0  # 箱の端からの距離（-8..8）
        self.assertAlmostEqual(wrapped, 0.0, delta=0.05)
        self.assertAlmostEqual(cy, 4.0)

    def test_tracked_cells_include_cells_that_joined_later(self):
        frame = self._result([[(0, 1), (1, 2)]]).frames[0]
        self.assertEqual(self.movie.tracked_cells(frame, [0, 1], 6), {0, 1, 2})
        self.assertEqual(self.movie.tracked_cells(frame, [4], 6), {4})


if __name__ == "__main__":
    unittest.main()
