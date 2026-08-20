"""用語集の解説図を生成する。

初見で分かりにくく、かつ結論を左右する用語を幾何学的に示す。

    python3 scripts/make_glossary_figure.py

A: 細胞の形状と測地距離
B: 隙間プロファイルと接触パッチの半径
C: 接触パッチと接触ストリップ（実寸比）
D: 極局在比と先端集中幅
E: 両極分配と、できる構造

図中の数値は末尾でモデルの計算値と突き合わせて検証する。
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Arc, Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion.geometry import (
    contact_patch_area,
    contact_strip_area,
    polar_area_fraction,
    total_surface_area,
)

for candidate in (
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
):
    try:
        font_manager.fontManager.addfont(candidate)
    except Exception:  # noqa: BLE001
        pass
matplotlib.rcParams["font.family"] = [
    "Hiragino Sans GB",
    "Hiragino Sans",
    "Noto Sans CJK JP",
    "Arial Unicode MS",
    "DejaVu Sans",
]
matplotlib.rcParams["axes.unicode_minus"] = False

LENGTH, RADIUS, BRIDGE = 1.80, 0.40, 0.008
HALF_CYLINDER = 0.5 * (LENGTH - 2 * RADIUS)
PATCH_RADIUS = math.sqrt(BRIDGE * RADIUS)
BRIDGE_NM = BRIDGE * 1000.0

COLOR_POLE = "#0E7C8A"
COLOR_SIDE = "#A6377C"
COLOR_ACCENT = "#A8722A"
COLOR_BODY = "#E4EAE6"
COLOR_PALE = "#C4CFCB"
COLOR_EDGE = "#6E7B77"
COLOR_INK = "#1B2422"
COLOR_MUTED = "#7B8783"


def capsule_outline(centre_x, centre_y, n_per_cap=120):
    """スフェロシリンダーを側面から見た輪郭の閉じた点列。

    全長は LENGTH、高さは 2 * RADIUS になる。丸めた端は円柱部の外側に付く。
    """
    right_angles = np.linspace(-math.pi / 2, math.pi / 2, n_per_cap)
    right = np.column_stack(
        (
            centre_x + HALF_CYLINDER + RADIUS * np.cos(right_angles),
            centre_y + RADIUS * np.sin(right_angles),
        )
    )
    left_angles = np.linspace(math.pi / 2, 3 * math.pi / 2, n_per_cap)
    left = np.column_stack(
        (
            centre_x - HALF_CYLINDER + RADIUS * np.cos(left_angles),
            centre_y + RADIUS * np.sin(left_angles),
        )
    )
    return np.vstack((right, left))


def draw_cell(axes, centre_x, centre_y, facecolor=COLOR_BODY, zorder=1):
    """スフェロシリンダーを側面から見た輪郭（全長 LENGTH、半径 RADIUS）。"""
    axes.add_patch(
        Polygon(
            capsule_outline(centre_x, centre_y),
            closed=True,
            linewidth=1.3,
            edgecolor=COLOR_EDGE,
            facecolor=facecolor,
            zorder=zorder,
        )
    )


def mark_pole(axes, x, y, colour, height=0.16):
    """極の位置に印をつける。"""
    axes.plot(
        [x, x], [y - height, y + height],
        color=colour, lw=5, solid_capstyle="round", zorder=6,
    )


def outline_samples(n_points=4000):
    """2次元輪郭上の点を、位置・領域・先端からの測地距離つきで返す。"""
    per_cap = per_side = n_points // 4

    x_side = np.linspace(-HALF_CYLINDER, HALF_CYLINDER, per_side)
    top = np.column_stack((x_side, np.full(per_side, RADIUS)))
    bottom = np.column_stack((x_side, np.full(per_side, -RADIUS)))

    theta = np.linspace(-math.pi / 2, math.pi / 2, per_cap)
    right = np.column_stack(
        (HALF_CYLINDER + RADIUS * np.cos(theta), RADIUS * np.sin(theta))
    )
    left = np.column_stack(
        (-HALF_CYLINDER - RADIUS * np.cos(theta), RADIUS * np.sin(theta))
    )

    points = np.vstack((top, bottom, right, left))
    regions = np.concatenate((np.zeros(2 * per_side), np.ones(per_cap), -np.ones(per_cap)))
    tip = RADIUS * np.abs(theta)
    tips = np.concatenate((np.zeros(2 * per_side), tip, tip))
    return points, regions, tips


def scatter_distribution(axes, centre_x, centre_y, f_pole, tip_width, n_dots=150, seed=3):
    """与えた局在でアドヘシンを輪郭上に散らす（模式図）。"""
    points, regions, tips = outline_samples()
    weights = np.zeros(len(points))
    for region, share in ((-1, f_pole / 2), (0, 1 - f_pole), (1, f_pole / 2)):
        mask = regions == region
        base = np.ones(mask.sum())
        if region != 0 and tip_width is not None:
            base = np.exp(-0.5 * (tips[mask] / tip_width) ** 2)
        if base.sum() > 0:
            weights[mask] = share * base / base.sum()

    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(points), size=n_dots, p=weights / weights.sum())
    dots = points[chosen] + rng.normal(scale=0.011, size=(n_dots, 2))
    dots += np.array([centre_x, centre_y])
    colours = np.where(regions[chosen] == 0, COLOR_SIDE, COLOR_POLE)
    axes.scatter(dots[:, 0], dots[:, 1], s=8, c=colours, zorder=5, linewidths=0)


def finish(axes, xlim, ylim, title):
    axes.set_xlim(*xlim)
    axes.set_ylim(*ylim)
    axes.set_aspect("equal")
    axes.axis("off")
    axes.set_title(title, loc="left", fontsize=12, weight="bold", pad=10)


# --------------------------------------------------------------------- パネル
def panel_anatomy(axes):
    draw_cell(axes, 0, 0)
    for sign in (-1, 1):
        axes.add_patch(
            Arc(
                (sign * HALF_CYLINDER, 0), 2 * RADIUS, 2 * RADIUS,
                theta1=-90 if sign > 0 else 90,
                theta2=90 if sign > 0 else 270,
                linewidth=5, edgecolor=COLOR_POLE, zorder=3,
            )
        )
    for sign in (-1, 1):
        axes.plot(
            [-HALF_CYLINDER, HALF_CYLINDER], [sign * RADIUS] * 2,
            color=COLOR_SIDE, lw=5, zorder=3,
        )

    # 寸法
    axes.annotate(
        "", xy=(-LENGTH / 2, -0.70), xytext=(LENGTH / 2, -0.70),
        arrowprops=dict(arrowstyle="<->", color=COLOR_INK, lw=1.1),
    )
    axes.text(0, -0.92, "細胞長 l = 1.80 µm", ha="center", fontsize=9)
    axes.annotate(
        "", xy=(-HALF_CYLINDER, 0), xytext=(-HALF_CYLINDER, RADIUS),
        arrowprops=dict(arrowstyle="<->", color=COLOR_INK, lw=1.1),
    )
    axes.text(-1.02, 0.20, "半径 r\n0.40 µm", fontsize=8.5, va="center", ha="right")

    # 領域
    axes.text(-0.15, 0.60, "円柱部（側面）　表面の 55.6 %",
              color=COLOR_SIDE, fontsize=9, weight="bold", ha="center")
    axes.text(1.62, -0.34, "極領域（半球）\n表面の 44.4 %",
              color=COLOR_POLE, fontsize=9, weight="bold", ha="center")

    # 測地距離
    theta = np.linspace(0, math.pi / 2.4, 60)
    arc_x = HALF_CYLINDER + (RADIUS + 0.05) * np.cos(theta)
    arc_y = (RADIUS + 0.05) * np.sin(theta)
    axes.plot(arc_x, arc_y, color=COLOR_ACCENT, lw=2.6, zorder=5)
    axes.plot([HALF_CYLINDER + RADIUS], [0], "o", ms=5, color=COLOR_ACCENT, zorder=6)
    axes.annotate(
        "測地距離 s\n先端から表面に沿って測る\n赤道まで最大 0.628 µm",
        xy=(arc_x[38], arc_y[38]), xytext=(0.78, 1.32),
        fontsize=8.5, color=COLOR_ACCENT, ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color=COLOR_ACCENT, lw=1.0),
    )

    finish(axes, (-2.30, 2.75), (-1.15, 1.72), "A　細胞の形状と測地距離")


def panel_gap(axes):
    """隙間プロファイル h = r²/R。縦方向を強調して描く。"""
    r_nm = np.linspace(-105, 105, 400)
    h_nm = r_nm**2 / (RADIUS * 1000)

    axes.axhline(0, color=COLOR_EDGE, lw=2)
    axes.plot(r_nm, h_nm, color=COLOR_EDGE, lw=2)
    axes.text(-102, 1.2, "下の細胞の表面", fontsize=8, color=COLOR_MUTED)
    axes.text(-102, 21.0, "上の細胞の表面", fontsize=8, color=COLOR_MUTED)

    inside = np.abs(r_nm) <= PATCH_RADIUS * 1000
    axes.fill_between(r_nm[inside], 0, h_nm[inside], color=COLOR_POLE, alpha=0.30, zorder=2)

    axes.axhline(BRIDGE_NM, color=COLOR_ACCENT, lw=1.8, ls="--")
    axes.text(104, BRIDGE_NM + 0.9, "架橋長 8 nm", ha="right",
              fontsize=9.5, color=COLOR_ACCENT, weight="bold")

    for sign in (-1, 1):
        axes.plot([sign * PATCH_RADIUS * 1000] * 2, [0, BRIDGE_NM],
                  color=COLOR_POLE, lw=1.2, ls=":")
    axes.annotate(
        "", xy=(-PATCH_RADIUS * 1000, -2.6), xytext=(PATCH_RADIUS * 1000, -2.6),
        arrowprops=dict(arrowstyle="<->", color=COLOR_POLE, lw=1.5),
    )
    axes.text(0, -5.4, "接触パッチの直径 113 nm", ha="center",
              fontsize=9.5, color=COLOR_POLE, weight="bold")
    axes.text(0, 12.5, "隙間  h = r² / r_cell\nh > 8 nm では架橋できない",
              ha="center", fontsize=9, color=COLOR_INK)

    axes.set_xlim(-105, 105)
    axes.set_ylim(-7.5, 25)
    axes.set_xlabel("接触点からの横方向距離 r (nm)", fontsize=9)
    axes.set_ylabel("隙間 h (nm)", fontsize=9)
    axes.set_title(
        "B　接触パッチの半径はどこから来るか\n"
        "　　　縁は r_c = √(架橋長 × 細胞半径) = 56.6 nm（縦方向を強調）",
        loc="left", fontsize=11, weight="bold", pad=8,
    )
    axes.spines[["top", "right"]].set_visible(False)
    axes.tick_params(labelsize=8.5)


def panel_patch_vs_strip(axes):
    """実寸比。パッチが線1本分にしか見えないことが要点。"""
    # 左: 極 - 極
    axes.text(-2.7, 1.42, "極 - 極（鎖の単位）", ha="center", fontsize=10, weight="bold")
    draw_cell(axes, -3.6, 0.30)
    draw_cell(axes, -1.8, 0.30)
    axes.plot([-2.7, -2.7], [0.30 - PATCH_RADIUS, 0.30 + PATCH_RADIUS],
              color=COLOR_POLE, lw=5, solid_capstyle="round", zorder=6)
    axes.annotate(
        "接触パッチ\n半径 56.6 nm・面積 0.0101 µm²\n表面の 0.22 %",
        xy=(-2.7, 0.30), xytext=(-2.7, -1.28),
        fontsize=9, color=COLOR_POLE, weight="bold", ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color=COLOR_POLE, lw=1.1),
    )

    # 右: 側面 - 側面
    axes.text(2.7, 1.42, "側面 - 側面（束の単位）", ha="center", fontsize=10, weight="bold")
    draw_cell(axes, 2.7, 0.30 + RADIUS)
    draw_cell(axes, 2.7, 0.30 - RADIUS)
    axes.plot([2.7 - HALF_CYLINDER, 2.7 + HALF_CYLINDER], [0.30, 0.30],
              color=COLOR_SIDE, lw=5, solid_capstyle="butt", zorder=6)
    axes.annotate(
        "接触ストリップ\n幅 113 nm × 長さ 1.0 µm・面積 0.1131 µm²\n表面の 2.50 %",
        xy=(2.7, 0.30), xytext=(2.7, -1.28),
        fontsize=9, color=COLOR_SIDE, weight="bold", ha="center", va="center",
        arrowprops=dict(arrowstyle="-", color=COLOR_SIDE, lw=1.1),
    )

    axes.text(
        0, 0.30, "面積比\n11.3 倍\n\n局在させなければ\n必ず束になる",
        fontsize=10, ha="center", va="center", weight="bold", color=COLOR_INK,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F2F4F1", edgecolor=COLOR_EDGE),
    )

    finish(axes, (-5.6, 5.6), (-2.05, 1.95),
           "C　接触パッチと接触ストリップ（細胞は実寸比）")


def panel_localization(axes):
    uniform = polar_area_fraction(LENGTH, RADIUS)
    cases = [
        (-3.4, uniform, None, f"f_pole = {uniform:.3f}\n先端集中なし",
         "一様 ＝ 局在なし\nこれが実測値", COLOR_SIDE),
        (0.0, 0.95, None, "f_pole = 0.95\n先端集中なし",
         "極に寄っているが\n半球全体に散る", COLOR_MUTED),
        (3.4, 0.95, 0.10, "f_pole = 0.95\nw = 0.10 µm",
         "先端に絞れている\nこれで鎖になる", COLOR_POLE),
    ]
    for x, f_pole, tip_width, label, comment, colour in cases:
        draw_cell(axes, x, 0.35)
        scatter_distribution(axes, x, 0.35, f_pole, tip_width)
        axes.text(x, 1.30, label, fontsize=9.5, ha="center", weight="bold")
        axes.text(x, -0.62, comment, fontsize=9, ha="center", va="top", color=colour)

    axes.text(
        0, -1.60,
        "注意：f_pole = 0 が「局在なし」ではない。極領域は表面の 44.4 % を占めるので、"
        "f_pole = 0.444 が一様分布にあたる。0.6 や 0.8 は「やや濃縮」にすぎない。",
        ha="center", va="center", fontsize=9, color=COLOR_ACCENT,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#F8F1E5", edgecolor=COLOR_ACCENT),
    )

    finish(axes, (-5.6, 5.6), (-2.25, 1.95),
           "D　極局在比 f_pole と先端集中幅 w　（● 極のアドヘシン　● 側面のアドヘシン）")


def panel_split(axes):
    # 対称二極 -> 鎖
    axes.text(-4.9, 0.72, "対称二極\npolar_split = 0.5",
              fontsize=9.5, ha="right", va="center", weight="bold")
    centres = [-2.7, -0.9, 0.9, 2.7]
    for x in centres:
        draw_cell(axes, x, 0.72)
    for x in (-3.6, -1.8, 0.0, 1.8, 3.6):
        mark_pole(axes, x, 0.72, COLOR_POLE)
    axes.text(4.3, 0.72, "鎖", fontsize=15, ha="left", va="center",
              weight="bold", color=COLOR_POLE)

    # 単極 -> 二量体
    axes.text(-4.9, -1.15, "単極\npolar_split = 1.0",
              fontsize=9.5, ha="right", va="center", weight="bold")
    for left_centre in (-3.1, 1.1):
        draw_cell(axes, left_centre, -1.15)
        draw_cell(axes, left_centre + LENGTH, -1.15)
        # 濃い極どうしが向かい合って二量体になる
        mark_pole(axes, left_centre + HALF_CYLINDER + RADIUS, -1.15, COLOR_POLE)
        # 外側の極は薄い
        mark_pole(axes, left_centre - HALF_CYLINDER - RADIUS, -1.15, COLOR_PALE)
        mark_pole(axes, left_centre + LENGTH + HALF_CYLINDER + RADIUS, -1.15, COLOR_PALE)
    axes.text(4.3, -1.15, "二量体\nどまり", fontsize=11, ha="left", va="center",
              weight="bold", color=COLOR_ACCENT)

    axes.text(
        0, -2.15,
        "濃い極（■ 濃色）と薄い側の極（■ 淡色）。鎖が伸びるには1細胞が2つの接続を要するので、"
        "律速は薄い側の極である。",
        ha="center", va="center", fontsize=9, color=COLOR_MUTED,
    )

    finish(axes, (-7.4, 5.6), (-2.65, 1.65),
           "E　両極分配 — 極に集めるだけでなく、両極に集める必要がある")


def verify(shown: dict[str, float]) -> bool:
    """図に書いた数値がモデルの計算値と一致することを確かめる。"""
    area = total_surface_area(LENGTH, RADIUS)
    computed = {
        "パッチ面積 (µm²)": contact_patch_area(RADIUS, BRIDGE),
        "ストリップ面積 (µm²)": contact_strip_area(LENGTH, RADIUS, BRIDGE),
        "パッチ半径 (nm)": PATCH_RADIUS * 1000,
        "面積比": contact_strip_area(LENGTH, RADIUS, BRIDGE)
        / contact_patch_area(RADIUS, BRIDGE),
        "一様の f_pole": polar_area_fraction(LENGTH, RADIUS),
        "パッチの表面比 (%)": contact_patch_area(RADIUS, BRIDGE) / area * 100,
        "ストリップの表面比 (%)": contact_strip_area(LENGTH, RADIUS, BRIDGE) / area * 100,
        "測地距離の最大 (µm)": math.pi * RADIUS / 2,
    }
    all_ok = True
    for name, value in computed.items():
        # 図中の表記は有効数字を丸めているので、その分を許容する
        ok = abs(value - shown[name]) / abs(shown[name]) < 0.025
        all_ok &= ok
        print(f"  {name:<22} 計算 {value:9.4f} / 図中 {shown[name]:9.4f}  "
              f"{'OK' if ok else '*** 不一致 ***'}")
    return all_ok


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "docs" / "figures" / "glossary_diagram.png",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    figure = plt.figure(figsize=(14.0, 16.6))
    grid = figure.add_gridspec(
        4, 2, height_ratios=[0.92, 1.00, 1.05, 1.00], hspace=0.34, wspace=0.16
    )
    panel_anatomy(figure.add_subplot(grid[0, 0]))
    panel_gap(figure.add_subplot(grid[0, 1]))
    panel_patch_vs_strip(figure.add_subplot(grid[1, :]))
    panel_localization(figure.add_subplot(grid[2, :]))
    panel_split(figure.add_subplot(grid[3, :]))

    figure.savefig(args.output, dpi=140, bbox_inches="tight", facecolor="white")
    print(f"図を {args.output.name} に保存した\n")

    ok = verify(
        {
            "パッチ面積 (µm²)": 0.0101,
            "ストリップ面積 (µm²)": 0.1131,
            "パッチ半径 (nm)": 56.6,
            "面積比": 11.3,
            "一様の f_pole": 0.444,
            "パッチの表面比 (%)": 0.22,
            "ストリップの表面比 (%)": 2.50,
            "測地距離の最大 (µm)": 0.628,
        }
    )
    if not ok:
        raise SystemExit("図中の数値がモデルの計算値と一致しません")


if __name__ == "__main__":
    main()
