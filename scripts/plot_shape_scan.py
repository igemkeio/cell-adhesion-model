"""`run_shape_scan.py` の結果を図と集計表にする。

    python3 scripts/plot_shape_scan.py

出力（`results/figures/`, `results/data/`）:
    shape_scan_static.png      静的モデルの結合数 P（極）・S（側面）と、その比
    shape_scan_bonds.png       時間発展: 20 秒時点の極結合数・側面結合数
    shape_scan_composition.png 時間発展: 結合の構成比（進行度を揃える）と凝集の速さ
    shape_scan_mechanism.png   時間発展: 形が効く経路の切り分けと、置き方への感度
    shape_scan_rotdiff.png     回転拡散の見積もり（A26）に対する円柱長の効果の感度
    shape_scan_summary.csv     動的な掃引の条件ごとの平均と標準誤差（図の元データ）

色の割り当ては全ての図・動画で共通（`scripts/figstyle.py`）: 極結合・鎖 = 青、側面結合・束 =
オレンジ。局在と一様、あるいは解析の変種のように「どの設計か」を表す系列は、色ではなく
線種・マーカー・直接ラベルで区別する（色だけに頼らない）。
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import model
from figstyle import (
    GRID, INK, INK_2, LATERAL, MUTED, POLAR, SURFACE, diverging_cmap, end_labels, style_axes,
)

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "results" / "data"
FIG_DIR = ROOT / "results" / "figures"
STATIC_CSV = DATA_DIR / "shape_scan_static.csv"
DYNAMIC_CSV = DATA_DIR / "shape_scan_dynamic.csv"
SUMMARY_CSV = DATA_DIR / "shape_scan_summary.csv"

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
    "Hiragino Sans GB", "Hiragino Sans", "Noto Sans CJK JP", "Arial Unicode MS", "DejaVu Sans",
]
matplotlib.rcParams["axes.unicode_minus"] = False

# 配色と体裁は全ての図・動画で共通（scripts/figstyle.py）。
DIVERGING = diverging_cmap()

DEFAULT_CYLINDER_UM = model.DEFAULT_SHAPE.cylinder_length_um
DEFAULT_RADIUS_UM = model.DEFAULT_SHAPE.radius_um


def panel_title(axes, text: str) -> None:
    axes.set_title(text, loc="left", fontsize=11, color=INK, fontweight="bold")


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


# ---------------------------------------------------------------------------
# 静的モデル
# ---------------------------------------------------------------------------
def static_figure() -> Path:
    rows = [{k: float(v) for k, v in r.items()} for r in read_csv(STATIC_CSV)]

    def at(cylinder: float, radius: float) -> dict:
        return next(
            r for r in rows
            if abs(r["cylinder_um"] - cylinder) < 1e-9 and abs(r["radius_um"] - radius) < 1e-9
        )

    cylinders = sorted({r["cylinder_um"] for r in rows})
    radii = sorted({r["radius_um"] for r in rows})

    figure, axes = plt.subplots(
        2, 2, figsize=(12.5, 8.6), facecolor=SURFACE, layout="constrained"
    )
    series = [
        # (キー, 線種, マーカー, ラベル, 色)
        ("localized_polar", "-", "o", "局在・極結合 P", POLAR),
        ("localized_side", "-", "o", "局在・側面結合 S", LATERAL),
        ("uniform_polar", "--", "s", "一様・極結合 P", POLAR),
        ("uniform_side", "--", "s", "一様・側面結合 S", LATERAL),
    ]

    def bond_panel(axes_, xs, picker, xlabel, title, default_x):
        labels = []
        for key, linestyle, marker, label, color in series:
            ys = [picker(x)[key] for x in xs]
            axes_.plot(xs, ys, linestyle, marker=marker, color=color, lw=2, ms=6,
                       mec=SURFACE, mew=1.5)
            labels.append((xs[-1], ys[-1], label))
        axes_.axvline(default_x, color=MUTED, lw=1, ls=":")
        axes_.text(default_x, 0.985, " 文献値", transform=axes_.get_xaxis_transform(),
                   fontsize=8, color=MUTED, va="top")
        axes_.set_yscale("log")
        axes_.set_xlabel(xlabel)
        axes_.set_ylabel("2細胞間の結合数（分子対）")
        axes_.set_xlim(xs[0], xs[-1] + 0.28 * (xs[-1] - xs[0]))
        panel_title(axes_, title)
        style_axes(axes_)
        end_labels(axes_, labels)

    bond_panel(
        axes[0][0], cylinders, lambda x: at(x, DEFAULT_RADIUS_UM),
        "円柱部分の長さ l_c (µm)   [半径 r = 0.40 µm 固定]",
        "A  円柱長を変える", DEFAULT_CYLINDER_UM,
    )
    bond_panel(
        axes[0][1], radii, lambda x: at(DEFAULT_CYLINDER_UM, x),
        "半球の半径 r (µm)   [円柱長 l_c = 1.0 µm 固定。曲率 = 1/r]",
        "B  半球の曲率を変える", DEFAULT_RADIUS_UM,
    )

    norm = TwoSlopeNorm(vmin=-2.0, vcenter=0.0, vmax=2.0)
    for axes_, scenario, title in (
        (axes[1][0], "localized", "C  P / S（局在: f_pole=0.95, 先端幅 0.10 µm）"),
        (axes[1][1], "uniform", "D  P / S（一様分布）"),
    ):
        grid = np.array([
            [np.log10(at(c, r)[f"{scenario}_polar"] / at(c, r)[f"{scenario}_side"]) for c in cylinders]
            for r in radii
        ])
        image = axes_.imshow(grid, origin="lower", cmap=DIVERGING, norm=norm, aspect="auto")
        axes_.set_xticks(range(len(cylinders)), [f"{c:g}" for c in cylinders])
        axes_.set_yticks(range(len(radii)), [f"{r:g}" for r in radii])
        axes_.set_xlabel("円柱部分の長さ l_c (µm)")
        axes_.set_ylabel("半球の半径 r (µm)")
        for (iy, ix), value in np.ndenumerate(grid):
            axes_.text(ix, iy, f"{10 ** value:.3g}", ha="center", va="center", fontsize=7.5,
                       color=INK if abs(value) < 1.2 else "#ffffff")
        # 文献値の形状を枠で示す
        ix = cylinders.index(DEFAULT_CYLINDER_UM)
        iy = radii.index(DEFAULT_RADIUS_UM)
        axes_.add_patch(plt.Rectangle((ix - 0.5, iy - 0.5), 1, 1, fill=False, ec=INK, lw=1.8))
        panel_title(axes_, title)
        for spine in axes_.spines.values():
            spine.set_visible(False)
        axes_.tick_params(colors=INK_2, labelsize=9, length=0)
        axes_.xaxis.label.set_color(INK_2)
        axes_.yaxis.label.set_color(INK_2)
        if scenario == "uniform":
            colorbar = figure.colorbar(
                image, ax=[axes[1][0], axes[1][1]], location="right", shrink=0.9, pad=0.02, aspect=28
            )
            colorbar.set_ticks(
                [-2, -1, 0, 1, 2],
                labels=["S が 100 倍", "S が 10 倍", "P = S", "P が 10 倍", "P が 100 倍"],
            )
            colorbar.ax.tick_params(labelsize=8, colors=INK_2)
            colorbar.outline.set_visible(False)

    figure.suptitle(
        "静的モデル: 細胞の形が結合数に与える影響（発現量 N = 7,300 / 細胞。枠は文献値の形状）",
        x=0.01, ha="left", fontsize=12.5, color=INK, fontweight="bold",
    )
    path = FIG_DIR / "shape_scan_static.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(figure)
    return path


# ---------------------------------------------------------------------------
# 動的な掃引の集計
# ---------------------------------------------------------------------------
N_CELLS = 24


def _number(value: str) -> float:
    return float(value) if value not in ("", None) else math.nan


def load_dynamic(path: Path) -> list[dict]:
    """1行1実行の CSV を読み、図に使う派生量を足す。"""

    rows = []
    for raw in read_csv(path):
        row = {k: _number(v) if k not in ("experiment", "scenario", "box_mode", "scaling",
                                          "first_bond_kind") else v for k, v in raw.items()}
        for t in (10, 20, 30):
            total = row[f"polar_t{t}"] + row[f"lateral_t{t}"]
            row[f"frac_t{t}"] = row[f"polar_t{t}"] / total if total > 0 else math.nan
        row["shape_key"] = (row["cylinder_um"], row["radius_um"])
        rows.append(row)
    return rows


def aggregate(rows: list[dict]) -> dict[tuple, dict]:
    """(実験, 設計, 円柱長, 半径) ごとに、各指標の平均・標準誤差・有効な実行数を返す。"""

    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["experiment"], row["scenario"], row["cylinder_um"], row["radius_um"])].append(row)

    metrics = [
        "polar_t10", "polar_t20", "polar_t30", "lateral_t10", "lateral_t20", "lateral_t30",
        "mono_t10", "mono_t20", "mono_t30", "largest_t20", "largest_t30",
        "frac_t10", "frac_t20", "frac_t30", "polar_frac_first10", "t_half_monomers_s",
        "dimer_t30", "chain_t30", "bundle_t30",
    ]
    result = {}
    for key, members in groups.items():
        entry = {
            "n_runs": len(members),
            "box_um": members[0]["box_um"],
            "area_fraction": members[0]["area_fraction"],
            "polar_bonds_P": members[0]["polar_bonds_P"],
            "side_bonds_S": members[0]["side_bonds_S"],
        }
        for metric in metrics:
            values = np.array([m[metric] for m in members], dtype=float)
            values = values[~np.isnan(values)]
            entry[f"{metric}_n"] = len(values)
            entry[f"{metric}_mean"] = float(values.mean()) if len(values) else math.nan
            entry[f"{metric}_sem"] = (
                float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else math.nan
            )
        result[key] = entry
    return result


def write_summary(table: dict[tuple, dict]) -> Path:
    keys = sorted(table)
    fields = ["experiment", "scenario", "cylinder_um", "radius_um"] + list(table[keys[0]])
    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for key in keys:
            writer.writerow(
                list(key) + [round(v, 4) if isinstance(v, float) else v for v in table[key].values()]
            )
    return SUMMARY_CSV


def scan_series(table, experiment, scenario, scan, metric):
    """1本の掃引（円柱長を振る / 半径を振る）の (x, 平均, 標準誤差) を x 昇順で返す。"""

    points = []
    for (exp, scen, cylinder, radius), entry in table.items():
        if exp != experiment or scen != scenario:
            continue
        if scan == "cylinder" and abs(radius - DEFAULT_RADIUS_UM) < 1e-9:
            points.append((cylinder, entry[f"{metric}_mean"], entry[f"{metric}_sem"]))
        if scan == "radius" and abs(cylinder - DEFAULT_CYLINDER_UM) < 1e-9:
            points.append((radius, entry[f"{metric}_mean"], entry[f"{metric}_sem"]))
    points.sort()
    if not points:
        return None
    xs, ys, es = (np.array(v, dtype=float) for v in zip(*points))
    return xs, ys, es


SCANS = {
    "cylinder": ("円柱部分の長さ l_c (µm)", DEFAULT_CYLINDER_UM, "円柱長を変える（r = 0.40 µm 固定）"),
    "radius": ("半球の半径 r (µm)  [曲率 = 1/r]", DEFAULT_RADIUS_UM, "半球の曲率を変える（l_c = 1.0 µm 固定）"),
}
SCENARIO_NAMES = {"localized": "局在（f_pole 0.95, 先端幅 0.10 µm）", "uniform": "一様分布"}


def draw_line(axes, series, *, color, linestyle="-", marker="o", label=None, band=True):
    """平均の折れ線と標準誤差の帯。系列が無ければ何もしない。"""

    if series is None:
        return None
    xs, ys, es = series
    axes.plot(xs, ys, linestyle, color=color, marker=marker, lw=2, ms=6, mec=SURFACE, mew=1.5, zorder=3)
    if band:
        axes.fill_between(xs, ys - es, ys + es, color=color, alpha=0.16, lw=0, zorder=2)
    return (xs[-1], ys[-1], label) if label else None


def mark_default(axes, x_default):
    axes.axvline(x_default, color=MUTED, lw=1, ls=":", zorder=1)
    axes.text(x_default, 0.985, " 文献値", transform=axes.get_xaxis_transform(),
              fontsize=8, color=MUTED, va="top")


def finish_axes(axes, xs_all, label_items, xlabel, ylabel, title):
    lo, hi = min(xs_all), max(xs_all)
    axes.set_xlim(lo - 0.03 * (hi - lo), hi + 0.30 * (hi - lo))
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    panel_title(axes, title)
    style_axes(axes)
    end_labels(axes, [item for item in label_items if item])


def bonds_figure(table, path: Path) -> Path:
    """20 秒時点での極結合（青）と側面結合（オレンジ）の本数。"""

    figure, axes = plt.subplots(2, 2, figsize=(12, 8), facecolor=SURFACE, layout="constrained",
                                sharey="row")
    letters = iter("ABCD")
    for row, scenario in enumerate(("localized", "uniform")):
        for col, scan in enumerate(("cylinder", "radius")):
            ax = axes[row][col]
            xlabel, default_x, scan_title = SCANS[scan]
            polar = scan_series(table, "main", scenario, scan, "polar_t20")
            lateral = scan_series(table, "main", scenario, scan, "lateral_t20")
            items = [
                draw_line(ax, polar, color=POLAR, label="極結合（鎖の背骨）"),
                draw_line(ax, lateral, color=LATERAL, marker="s", label="側面結合（束）"),
            ]
            mark_default(ax, default_x)
            xs_all = [] if polar is None else list(polar[0])
            finish_axes(ax, xs_all, items, xlabel, "20 秒時点の結合数（24 細胞あたり）",
                        f"{next(letters)}  {SCENARIO_NAMES[scenario].split('（')[0]} ─ {scan_title}")
        axes[row][0].set_ylim(bottom=0)
    figure.suptitle("時間発展: 形を変えたときに20秒までにできる結合の種類（平均 ± 標準誤差、面積占有率を一定に保つ）",
                    x=0.01, ha="left", fontsize=12, color=INK, fontweight="bold")
    figure.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(figure)
    return path


def composition_figure(table, path: Path) -> Path:
    """結合の構成比（10 本目の結合まで）と、単量体が半減するまでの時間。"""

    figure, axes = plt.subplots(2, 2, figsize=(12, 8), facecolor=SURFACE, layout="constrained",
                                sharey="row")
    styles = {
        "localized": dict(color=INK, linestyle="-", marker="o", label="局在"),
        "uniform": dict(color=MUTED, linestyle="--", marker="s", label="一様"),
    }
    metrics = [
        ("polar_frac_first10", "最初の10本のうち極結合の割合", (0, 1)),
        ("t_half_monomers_s", "単量体が半分になるまでの時間 (s)", None),
    ]
    letters = iter("ABCD")
    for row, (metric, ylabel, ylim) in enumerate(metrics):
        for col, scan in enumerate(("cylinder", "radius")):
            ax = axes[row][col]
            xlabel, default_x, scan_title = SCANS[scan]
            items, xs_all = [], []
            for scenario, style in styles.items():
                series = scan_series(table, "main", scenario, scan, metric)
                items.append(draw_line(ax, series, **style))
                if series is not None:
                    xs_all = list(series[0])
            if metric == "polar_frac_first10":
                ax.axhline(0.5, color=GRID, lw=1.2, zorder=1)
            mark_default(ax, default_x)
            if ylim:
                ax.set_ylim(*ylim)
            finish_axes(ax, xs_all, items, xlabel, ylabel, f"{next(letters)}  {scan_title}")
        if not ylim:
            axes[row][0].set_ylim(bottom=0)
    figure.suptitle("時間発展: 形が結合の構成比と凝集の速さに与える影響（平均 ± 標準誤差、面積占有率を一定に保つ）",
                    x=0.01, ha="left", fontsize=12, color=INK, fontweight="bold")
    figure.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(figure)
    return path


def mechanism_figure(table, path: Path) -> Path:
    """形が効く経路の切り分け（結合数を既定形状に固定）と、置き方への感度。"""

    figure, axes = plt.subplots(2, 2, figsize=(12, 8), facecolor=SURFACE, layout="constrained",
                                sharey=True)
    variants = [
        ("main", "全て反映", dict(color=INK, linestyle="-", marker="o")),
        ("decoupled", "結合数を既定形状に固定", dict(color=INK_2, linestyle="--", marker="D")),
        ("number_density", "箱を固定（細胞数密度一定）", dict(color=MUTED, linestyle="-.", marker="s")),
        ("per_area", "発現量を膜面積に比例", dict(color=MUTED, linestyle=":", marker="^")),
    ]
    letters = iter("ABCD")
    for row, scenario in enumerate(("localized", "uniform")):
        for col, scan in enumerate(("cylinder", "radius")):
            ax = axes[row][col]
            xlabel, default_x, scan_title = SCANS[scan]
            items, xs_all = [], []
            for experiment, label, style in variants:
                series = scan_series(table, experiment, scenario, scan, "polar_frac_first10")
                items.append(draw_line(ax, series, label=label, band=experiment in ("main", "decoupled"),
                                       **style))
                if series is not None and experiment == "main":
                    xs_all = list(series[0])
            ax.axhline(0.5, color=GRID, lw=1.2, zorder=1)
            mark_default(ax, default_x)
            ax.set_ylim(0, 1)
            finish_axes(ax, xs_all, items, xlabel, "最初の10本のうち極結合の割合",
                        f"{next(letters)}  {SCENARIO_NAMES[scenario].split('（')[0]} ─ {scan_title}")
    figure.suptitle("時間発展: 形の効果はどの経路を通るか（結合数 vs 衝突の幾何）と、置き方への感度",
                    x=0.01, ha="left", fontsize=12, color=INK, fontweight="bold")
    figure.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(figure)
    return path


def rotdiff_figure(table, path: Path) -> Path | None:
    """回転拡散係数の見積もり（A26）を変えても、円柱長への応答が保たれるか。"""

    variants = [
        ("main", "等価球近似（既定）", dict(color=INK, linestyle="-", marker="o")),
        ("rotdiff_0.5", "回転拡散 × 0.5", dict(color=INK_2, linestyle="--", marker="D")),
        ("rotdiff_0.25", "回転拡散 × 0.25", dict(color=MUTED, linestyle=":", marker="^")),
    ]
    if not any(scan_series(table, name, "localized", "cylinder", "polar_frac_first10") for name, *_ in variants[1:]):
        return None
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), facecolor=SURFACE, layout="constrained")
    shown = (0.1, 1.0, 3.0)
    for ax, metric, ylabel, ylim, title in (
        (axes[0], "polar_frac_first10", "最初の10本のうち極結合の割合", (0, 1), "A  結合の構成比"),
        (axes[1], "t_half_monomers_s", "単量体が半分になるまでの時間 (s)", None, "B  凝集の速さ"),
    ):
        items = []
        for name, label, style in variants:
            series = scan_series(table, name, "localized", "cylinder", metric)
            if series is None:
                continue
            xs, ys, es = (np.array([v[i] for i, x in enumerate(series[0]) if x in shown]) for v in series)
            ax.errorbar(xs, ys, yerr=es, capsize=3, lw=2, ms=6, mec=SURFACE, mew=1.5, **{
                "color": style["color"], "linestyle": style["linestyle"], "marker": style["marker"]})
            items.append((xs[-1], ys[-1], label))
        mark_default(ax, DEFAULT_CYLINDER_UM)
        if ylim:
            ax.axhline(0.5, color=GRID, lw=1.2, zorder=1)
            ax.set_ylim(*ylim)
        else:
            ax.set_ylim(bottom=0)
        ax.set_xlim(-0.05, 3.9)
        ax.set_xticks(shown)
        ax.set_xlabel("円柱部分の長さ l_c (µm)   [r = 0.40 µm、局在]")
        ax.set_ylabel(ylabel)
        panel_title(ax, title)
        style_axes(ax)
        end_labels(ax, items)
    figure.suptitle("回転拡散の見積もりを下げても、円柱が長いほど極結合が減る傾向は保たれるか（平均 ± 標準誤差）",
                    x=0.01, ha="left", fontsize=12, color=INK, fontweight="bold")
    figure.savefig(path, dpi=160, bbox_inches="tight", facecolor=SURFACE)
    plt.close(figure)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dynamic-csv", type=Path, default=DYNAMIC_CSV)
    parser.add_argument("--out-dir", type=Path, default=FIG_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def show(path: Path | None) -> None:
        if path is not None:
            print(f"図を保存した: {path}")

    show(static_figure())
    if not args.dynamic_csv.exists():
        print(f"{args.dynamic_csv} が無いので、動的な図は作らない（先に run_shape_scan.py dynamic）")
        return
    table = aggregate(load_dynamic(args.dynamic_csv))
    print(f"集計表: {write_summary(table)}")
    for make, name in (
        (bonds_figure, "shape_scan_bonds.png"),
        (composition_figure, "shape_scan_composition.png"),
        (mechanism_figure, "shape_scan_mechanism.png"),
        (rotdiff_figure, "shape_scan_rotdiff.png"),
    ):
        show(make(table, args.out_dir / name))


if __name__ == "__main__":
    main()
