"""主要結果の図を生成する。

    python3 scripts/make_figures.py

A: Costan et al. (2024) の報告値との照合
B: 局在に対する形態マップ（鎖と束の境界）
C: 接合強度と震盪培養のせん断荷重の比較
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import literature as lit
from adhesion.model import UNIFORM_F_POLE, predict
from adhesion.shear import chain_tension_pn
from figstyle import (
    INK, INK_2, LATERAL, LATERAL_TEXT, MUTED, POLAR, POLAR_MID, POLAR_TEXT, SURFACE,
    diverging_cmap, end_labels, style_axes,
)

# 日本語ラベル用のフォント。見つからなければ既定にフォールバックする。
for candidate in (
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
):
    try:
        font_manager.fontManager.addfont(candidate)
    except Exception:  # noqa: BLE001 - フォントが無い環境でも動かす
        pass
matplotlib.rcParams["font.family"] = [
    "Hiragino Sans GB",
    "Hiragino Sans",
    "Noto Sans CJK JP",
    "Arial Unicode MS",
    "DejaVu Sans",
]
matplotlib.rcParams["axes.unicode_minus"] = False

# 配色は全ての図・動画で共通（scripts/figstyle.py）: 極 = 青、側面 = オレンジ。
COLOR_POLE = POLAR
COLOR_SIDE = LATERAL
COLOR_REPORTED = MUTED  # 報告値

N_MID = float(lit.ADHESINS_AT_MID_INDUCTION)

# 形態マップは条件数が多いので、粗い求積で十分（境界の位置は変わらない）
MAP_POINTS = 6000


def panel_validation(axes) -> None:
    prediction = predict(N_MID, UNIFORM_F_POLE)
    model = [prediction.bonds_side_by_side, prediction.bonds_heavy_pole]
    reported = [lit.REPORTED_BONDS_PARALLEL.value, lit.REPORTED_BONDS_PERPENDICULAR.value]
    errors = [
        lit.REPORTED_BONDS_PARALLEL.uncertainty,
        lit.REPORTED_BONDS_PERPENDICULAR.uncertainty,
    ]

    positions = np.arange(2)
    width = 0.36
    # 色は接触の種類（横並び = 側面 = オレンジ、極 = 青）。報告値は灰。
    axes.bar(positions - width / 2, model, width, color=[COLOR_SIDE, COLOR_POLE])
    axes.bar(
        positions + width / 2, reported, width, yerr=errors, capsize=5,
        color=COLOR_REPORTED, error_kw=dict(ecolor=INK_2),
    )
    for position, value, text_color in zip(positions, model, (LATERAL_TEXT, POLAR_TEXT)):
        axes.text(position - width / 2, value + 4, f"このモデル\n{value:.0f}", ha="center",
                  fontsize=8.5, color=text_color)
    for position, value, error in zip(positions, reported, errors):
        axes.text(position + width / 2, value + error + 4, f"報告値\n{value:.0f}", ha="center",
                  fontsize=8.5, color=INK_2)

    axes.set_xticks(positions)
    axes.set_xticklabels(["横並び接触", "極接触"])
    axes.set_ylabel("2細胞間の結合数")
    axes.set_ylim(0, 200)
    axes.set_xlabel("報告値: Costan et al. (2024)", fontsize=9)
    axes.set_title("A  報告値との照合（N=7,300・一様分布）", loc="left", fontsize=11)
    style_axes(axes)
    axes.grid(False)


def panel_morphology_map(axes, figure) -> None:
    f_values = np.linspace(UNIFORM_F_POLE, 0.99, 60)
    w_values = np.linspace(0.06, 0.60, 55)

    grid = np.zeros((len(w_values), len(f_values)))
    for row, width in enumerate(w_values):
        for column, f_pole in enumerate(f_values):
            result = predict(N_MID, f_pole, 0.5, float(width), n_points=MAP_POINTS)
            grid[row, column] = np.log10(
                max(result.bonds_heavy_pole, 1e-9)
                / max(result.bonds_side_by_side, 1e-9)
            )

    mesh = axes.pcolormesh(
        f_values, w_values, grid, cmap=diverging_cmap(), vmin=-1.5, vmax=1.5, shading="auto"
    )
    axes.contour(f_values, w_values, grid, levels=[0], colors="k", linewidths=2)

    axes.plot(
        UNIFORM_F_POLE, 0.595, "o", ms=10, mfc=INK, mec="white", mew=1.8,
        zorder=5, clip_on=False,
    )
    axes.annotate(
        "実測の局在\n（一様・束）",
        xy=(UNIFORM_F_POLE, 0.595),
        xytext=(0.56, 0.53),
        fontsize=8.5,
        color=INK,
        weight="bold",
        arrowprops=dict(arrowstyle="-", color=INK, lw=1.2),
    )
    for f_pole, width, label in ((0.95, 0.20, "先端0.20µm"), (0.95, 0.10, "先端0.10µm")):
        axes.plot(f_pole, width, "o", ms=8, mfc="white", mec="k", mew=1.6, zorder=5)
        axes.annotate(
            label,
            xy=(f_pole, width),
            xytext=(-8, 0),
            textcoords="offset points",
            fontsize=8,
            ha="right",
            va="center",
        )

    axes.set_xlabel("極局在比  f_pole")
    axes.set_ylabel("先端集中度  w (µm)")
    axes.set_title("B  形態マップ（黒線 = 鎖と束の境界）", loc="left", fontsize=11)
    axes.text(0.62, 0.12, "鎖 chain", fontsize=11, color=POLAR_TEXT, weight="bold")
    axes.text(0.52, 0.44, "束 bundle", fontsize=11, color=LATERAL_TEXT, weight="bold")
    bar = figure.colorbar(mesh, ax=axes)
    bar.set_label("log₁₀(極接触 / 側面接触)", fontsize=9)


def panel_shear(axes) -> list[tuple[float, float, str]]:
    """せん断荷重の図。右端のラベルは、レイアウト確定後に `end_labels` で置く。"""

    cells = np.arange(2, 61)
    for shear_rate, style in ((100, ":"), (300, "-"), (1000, "--")):
        axes.plot(
            cells,
            [chain_tension_pn(int(n), shear_rate) for n in cells],
            style,
            color=INK_2,
            label=f"せん断荷重 γ̇={shear_rate} s⁻¹",
        )

    designs = [
        ("一様（束）", UNIFORM_F_POLE, None, COLOR_SIDE),
        ("極95% 先端0.20µm", 0.95, 0.20, POLAR_MID),
        ("極95% 先端0.10µm", 0.95, 0.10, COLOR_POLE),
    ]
    labels = []
    for label, f_pole, width, color in designs:
        result = predict(N_MID, f_pole, 0.5, width)
        axes.axhline(result.junction_force_pn, color=color, lw=2)
        labels.append((60, result.junction_force_pn, label))

    axes.set_yscale("log")
    axes.set_xlim(2, 60)
    axes.set_xlabel("鎖を構成する細胞数")
    axes.set_ylabel("力 (pN)")
    axes.set_title("C  接合強度 vs 震盪のせん断荷重", loc="left", fontsize=11)
    axes.legend(frameon=False, fontsize=8.5, loc="lower right")
    style_axes(axes)
    axes.grid(False)
    return labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "results"
        / "figures"
        / "main_results.png",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    panel_validation(axes[0])
    panel_morphology_map(axes[1], figure)
    labels = panel_shear(axes[2])
    figure.tight_layout()
    end_labels(axes[2], labels, min_gap_px=16.0)
    figure.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"図を {args.output.name} に保存した")


if __name__ == "__main__":
    main()
