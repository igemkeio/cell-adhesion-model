"""細胞の形を変えたときの鎖・束の形成を、並べて動画にする。

    python3 scripts/simulate_shape_comparison.py --scan cylinder   # 円柱長 l_c = 0.1 / 1.0 / 3.0 µm
    python3 scripts/simulate_shape_comparison.py --scan radius     # 半径 r = 0.25 / 0.40 / 0.60 µm
    python3 scripts/simulate_shape_comparison.py --shapes 1.0:0.4 2.0:0.3 --scenario uniform

半球の半径は円柱の半径と等しいので、半球の曲率 1/r は円柱の太さと切り離せない。
形は円柱長 l_c と半径 r（曲率 = 1/r）で指定する。

`scripts/simulate_chain_formation.py` と同じ描画（カプセル形の細胞。灰 = 結合なし、青 = 極結合
だけをもつ細胞、オレンジ = 側面結合をもつ細胞。● = Ag 極、○ = Nb 極）で、形ごとの
パネルを横に並べる。細胞の拡大率は
全パネルで揃えてあるので、形の違いがそのまま見える。下段は各形での接触の箇所数
（極接触 = 鎖の背骨、側面接触 = 束）の累積。

**この動画は1回の実行であり、統計ではない。** シードで結果はかなり変わる。
既定の `--seed representative` は、`run_shape_scan.py` の結果から、形ごとに「20 秒時点の
極結合の割合と単量体の半減時間が、その形の平均に最も近い実行」のシードを選ぶ（見栄えの
よい実行を選ぶ恣意性を避けるため）。同じ設定なら最初の 30 秒はスイープの実行と完全に
一致する。形の効果の統計は `results/figures/shape_scan_*.png` を見ること。

条件はスイープと同じ。箱の大きさは細胞が占める面積の割合を既定形状に揃え（`--box area`）、
時間刻みは全ての形で同じ（`DynamicsConfig` の既定、収束した値）。
"""

from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import run_shape_scan as scan
import simulate_chain_formation as base  # 配色・フォント設定・クラスタ着色を共有する
from adhesion import dynamics as dyn
from adhesion import model
from adhesion.shape import CellShape
from figstyle import INK

PRESETS = {
    "cylinder": [(0.1, 0.4), (1.0, 0.4), (3.0, 0.4)],
    "radius": [(1.0, 0.25), (1.0, 0.4), (1.0, 0.6)],
}
INCHES_PER_UM = 0.17


def capsule_outline(shape: CellShape, n_arc: int = 16) -> np.ndarray:
    """その形のカプセルの輪郭を、細胞中心基準で返す。"""

    return base._capsule_outline_local(shape.radius_um, shape.half_cylinder_um, n_arc)


def world_outline(local: np.ndarray, position: np.ndarray, theta: float) -> np.ndarray:
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    rotation = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
    return position + local @ rotation.T


def bond_segments(positions, angles, bonds, shape: CellShape):
    segments, colors = [], []
    for bond in bonds:
        p_i = dyn.attachment_point(positions[bond.i], angles[bond.i], bond.axial_i, bond.delta_i,
                                   shape.radius_um)
        p_j = dyn.attachment_point(positions[bond.j], angles[bond.j], bond.axial_j, bond.delta_j,
                                   shape.radius_um)
        segments.append(np.array([p_i, p_j]))
        colors.append(base.POLAR_BOND_COLOR if bond.kind == "polar" else base.LATERAL_BOND_COLOR)
    return segments, colors


def representative_seed(condition: scan.Condition) -> int | None:
    """スイープの実行のうち、その形の平均に最も近いもののシード。見つからなければ None。

    距離は、20 秒時点の極結合の割合と単量体の半減時間を、それぞれ標準偏差で割って足したもの。
    """

    if not scan.DYNAMIC_CSV.exists():
        return None
    rows = []
    for raw in scan.csv.DictReader(scan.DYNAMIC_CSV.open(encoding="utf-8")):
        if (raw["experiment"], raw["scenario"], raw["box_mode"], raw["scaling"], raw["decouple"]) != (
            "main", condition.scenario, condition.box, "per_cell", "0"
        ):
            continue
        if abs(float(raw["cylinder_um"]) - condition.cylinder_um) > 1e-9:
            continue
        if abs(float(raw["radius_um"]) - condition.radius_um) > 1e-9:
            continue
        total = float(raw["polar_t20"]) + float(raw["lateral_t20"])
        if total > 0:
            rows.append((int(raw["seed"]), float(raw["polar_t20"]) / total, float(raw["t_half_monomers_s"])))
    if len(rows) < 5:
        return None
    seeds, fractions, halves = (np.array(v, dtype=float) for v in zip(*rows))
    distance = np.zeros(len(rows))
    for values in (fractions, halves):
        spread = values.std() or 1.0
        distance += ((values - values.mean()) / spread) ** 2
    return int(seeds[int(np.argmin(distance))])


def movie_facts(results: list[dyn.SimulationResult], seeds: list[int], args, name: str) -> dict:
    """動画の条件と、最終フレームの状態。文書（docs/shape_movies.md）の元データになる。"""

    panels = []
    for result, seed in zip(results, seeds):
        config = result.config
        final = result.frames[-1]
        clusters = dyn.classify_clusters(config.n_cells, final.bonds, config.distinct_adhesins_per_pole)
        shape = config.shape
        panels.append({
            "cylinder_um": round(shape.cylinder_length_um, 4),
            "radius_um": shape.radius_um,
            "length_um": round(shape.length_um, 4),
            "curvature_per_um": round(shape.curvature_per_um, 4),
            "seed": seed,
            "box_um": config.box_size_um,
            "dt_s": config.dt_s,
            "polar_bonds_per_contact": round(result.polar_bond_count, 2),
            "side_bonds_per_contact": round(result.side_bond_count, 2),
            "final_polar_contacts": sum(b.kind == "polar" for b in final.bonds),
            "final_lateral_contacts": sum(b.kind == "lateral" for b in final.bonds),
            "final_monomers": sum(c.kind == "monomer" for c in clusters),
            "final_clusters": sum(c.kind != "monomer" for c in clusters),
            "final_largest_cluster": max(len(c.cells) for c in clusters),
        })
    return {
        "movie": name,
        "scenario": args.scenario,
        "seed_rule": str(args.seed),
        "box_mode": args.box,
        "n_cells": args.n_cells,
        "total_time_s": args.total_time_s,
        "fps": args.fps,
        "n_frames": len(results[0].frames),
        "panels": panels,
    }


def run_condition(job: tuple) -> dyn.SimulationResult:
    condition, seed, total_time_s, n_cells = job
    config = scan.build_config(condition, seed, total_time_s, None, n_cells)
    return dyn.run(config)


def build_animation(results: list[dyn.SimulationResult], fps: int) -> animation.FuncAnimation:
    n_panels = len(results)
    boxes = [r.config.box_size_um for r in results]
    n_frames = min(len(r.frames) for r in results)
    times = [f.time_s for f in results[0].frames[:n_frames]]
    shapes = [r.config.shape for r in results]
    n_cells = results[0].config.n_cells

    # 拡大率（1 µm あたりの長さ）を全パネルで厳密にそろえるため、軸は絶対寸法（インチ）で置く。
    # gridspec の比率配置だと、正方形に縮められたパネルの拡大率だけ小さくなる。
    left_margin, right_margin, gap = 0.75, 0.25, 0.85
    top_margin, middle, bottom_margin, count_height = 0.95, 0.95, 1.05, 1.7
    box_inches = [b * INCHES_PER_UM for b in boxes]
    width = max(6.0, left_margin + sum(box_inches) + gap * (n_panels - 1) + right_margin)
    height = top_margin + max(box_inches) + middle + count_height + bottom_margin
    figure = plt.figure(figsize=(width, height))
    lefts = [left_margin + sum(box_inches[:k]) + gap * k for k in range(n_panels)]
    default = model.DEFAULT_SHAPE

    # 接触の箇所数（極・側面）の時系列
    series = []
    for result in results:
        frames = result.frames[:n_frames]
        series.append((
            [sum(b.kind == "polar" for b in f.bonds) for f in frames],
            [sum(b.kind == "lateral" for b in f.bonds) for f in frames],
        ))
    y_max = max(max(p + l) for p, l in series) * 1.15 + 1

    panels = []
    for k, (result, shape) in enumerate(zip(results, shapes)):
        top_row_bottom = (bottom_margin + count_height + middle + max(box_inches) - box_inches[k]) / height
        ax_box = figure.add_axes(
            [lefts[k] / width, top_row_bottom, box_inches[k] / width, box_inches[k] / height]
        )
        ax_box.set_xlim(0, boxes[k])
        ax_box.set_ylim(0, boxes[k])
        ax_box.set_aspect("equal", adjustable="box")
        ax_box.set_xlabel("x (µm)", fontsize=8)
        ax_box.set_ylabel("y (µm)", fontsize=8)
        ax_box.tick_params(labelsize=7)
        is_default = (abs(shape.radius_um - default.radius_um) < 1e-9
                      and abs(shape.cylinder_length_um - default.cylinder_length_um) < 1e-9)
        mark = "（文献値）" if is_default else ""
        ax_box.set_title(
            f"円柱 l_c={shape.cylinder_length_um:.2g} µm, 半径 r={shape.radius_um:.2g} µm{mark}\n"
            f"1接触あたり 極 {result.polar_bond_count:.0f} 対 / 側面 {result.side_bond_count:.1f} 対",
            fontsize=9, loc="left",
        )
        local = capsule_outline(shape)
        patches = [
            Polygon(world_outline(local, result.frames[0].positions[i], result.frames[0].angles[i]),
                    closed=True, facecolor=base.MONOMER_COLOR, edgecolor="black", linewidth=0.4,
                    zorder=2)
            for i in range(n_cells)
        ]
        for patch in patches:
            ax_box.add_patch(patch)
        heads = ax_box.scatter([], [], s=10, color=INK, zorder=3)
        tails = ax_box.scatter([], [], s=10, facecolor="white", edgecolor=INK, linewidth=0.7, zorder=3)
        lines = LineCollection([], linewidths=1.6, zorder=1)
        ax_box.add_collection(lines)

        ax_count = figure.add_axes(
            [lefts[k] / width, bottom_margin / height, box_inches[k] / width, count_height / height]
        )
        status = ax_count.set_title("", loc="left", fontsize=8.5)
        ax_count.set_xlim(0, times[-1])
        ax_count.set_ylim(0, y_max)
        ax_count.set_xlabel("時間 (s)", fontsize=8)
        ax_count.set_ylabel("接触の箇所数", fontsize=8)
        ax_count.tick_params(labelsize=7)
        ax_count.spines[["top", "right"]].set_visible(False)
        (polar_line,) = ax_count.plot([], [], color=base.POLAR_BOND_COLOR, lw=2, label="極接触（鎖）")
        (lateral_line,) = ax_count.plot([], [], color=base.LATERAL_BOND_COLOR, lw=2, label="側面接触（束）")
        if k == 0:
            ax_count.legend(loc="upper left", fontsize=7.5, frameon=False)
        panels.append((patches, heads, tails, lines, status, polar_line, lateral_line, local))

    title = figure.suptitle("t = 0 s", x=left_margin / width, ha="left", fontsize=12, fontweight="bold")
    figure.legend(
        handles=base.legend_handles(), loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=5,
        fontsize=8.5, frameon=False, columnspacing=1.6,
    )

    def update(index: int):
        artists = []
        for k, result in enumerate(results):
            frame = result.frames[index]
            patches, heads, tails, lines, status, polar_line, lateral_line, local = panels[k]
            colors = base.cell_colors_for_frame(n_cells, frame)
            reach = shapes[k].half_cylinder_um + shapes[k].radius_um
            head_xy = np.zeros((n_cells, 2))
            tail_xy = np.zeros((n_cells, 2))
            for i in range(n_cells):
                patches[i].set_xy(world_outline(local, frame.positions[i], frame.angles[i]))
                patches[i].set_facecolor(colors[i])
                axis = np.array([np.cos(frame.angles[i]), np.sin(frame.angles[i])])
                head_xy[i] = frame.positions[i] + reach * axis
                tail_xy[i] = frame.positions[i] - reach * axis
            heads.set_offsets(head_xy)
            tails.set_offsets(tail_xy)
            segments, seg_colors = bond_segments(frame.positions, frame.angles, frame.bonds, shapes[k])
            lines.set_segments(segments)
            lines.set_color(seg_colors)

            n_polar, n_lateral = series[k][0][index], series[k][1][index]
            n_mono = sum(c.kind == "monomer" for c in dyn.classify_clusters(
                n_cells, frame.bonds, result.config.distinct_adhesins_per_pole))
            status.set_text(f"極接触 {n_polar}  側面接触 {n_lateral}  単独の細胞 {n_mono}")
            polar_line.set_data(times[: index + 1], series[k][0][: index + 1])
            lateral_line.set_data(times[: index + 1], series[k][1][: index + 1])
            artists += [*patches, heads, tails, lines, status, polar_line, lateral_line]
        title.set_text(f"t = {times[index]:.0f} s")
        return artists

    return animation.FuncAnimation(figure, update, frames=n_frames, interval=1000.0 / fps, blit=False)


def parse_shapes(values: list[str]) -> list[tuple[float, float]]:
    shapes = []
    for value in values:
        cylinder, radius = value.split(":")
        shapes.append((float(cylinder), float(radius)))
    return shapes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scan", choices=sorted(PRESETS), default=None,
                        help="既定の組: cylinder（円柱長を変える）/ radius（曲率を変える）")
    parser.add_argument("--shapes", nargs="+", metavar="L_C:R", default=None,
                        help="任意の形を 円柱長:半径（µm）で指定。例 1.0:0.4 3.0:0.4")
    parser.add_argument("--scenario", choices=("localized", "uniform"), default="localized")
    parser.add_argument("--box", choices=("area", "number"), default="area")
    parser.add_argument("--seed", default="representative",
                        help="整数、または representative（既定。スイープの結果から形ごとに選ぶ）")
    parser.add_argument("--n-cells", type=int, default=24)
    parser.add_argument("--total-time-s", type=float, default=60.0)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument("--dpi", type=int, default=80)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--output", type=Path, default=None,
                        help="既定では results/animations/shape_<scan>_<scenario>.gif")
    args = parser.parse_args()

    if args.shapes:
        shapes, name = parse_shapes(args.shapes), "custom"
    else:
        name = args.scan or "cylinder"
        shapes = PRESETS[name]
    conditions = [
        scan.Condition("movie", args.scenario, cylinder, radius, box=args.box)
        for cylinder, radius in shapes
    ]
    seeds = []
    for condition in conditions:
        if args.seed == "representative":
            seed = representative_seed(condition)
            if seed is None:
                print(f"  注意: l_c={condition.cylinder_um} r={condition.radius_um} は"
                      "スイープの結果が足りないので seed 0 を使う")
                seed = 0
        else:
            seed = int(args.seed)
        seeds.append(seed)
    jobs = [(c, seed, args.total_time_s, args.n_cells) for c, seed in zip(conditions, seeds)]

    print(f"{len(jobs)} 形状をシミュレーション中（{args.scenario}, seed {seeds}）...", flush=True)
    with Pool(min(args.jobs, len(jobs))) as pool:
        results = pool.map(run_condition, jobs)
    for (cylinder, radius), result in zip(shapes, results):
        final = result.frames[-1].bonds
        print(f"  l_c={cylinder} r={radius}: dt={result.config.dt_s:.5f} s, 箱 {result.config.box_size_um} µm, "
              f"P={result.polar_bond_count:.1f} S={result.side_bond_count:.1f}, "
              f"最終 極 {sum(b.kind == 'polar' for b in final)} / 側面 {sum(b.kind == 'lateral' for b in final)}")

    output = args.output or (
        Path(__file__).resolve().parent.parent / "results" / "animations"
        / f"shape_{name}_{args.scenario}.gif"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    anim = build_animation(results, args.fps)
    anim.save(output, writer=animation.PillowWriter(fps=args.fps), dpi=args.dpi)
    print(f"動画を {output} に保存した")

    facts_path = output.with_suffix(".json")
    facts_path.write_text(
        json.dumps(movie_facts(results, seeds, args, name), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(f"動画の条件と最終状態を {facts_path} に保存した（docs/shape_movies.md の元データ）")


if __name__ == "__main__":
    main()
