"""細胞集団が時間とともに衝突し、鎖・束を作っていく過程を動画にする。

    python3 scripts/simulate_chain_formation.py --scenario chain
    python3 scripts/simulate_chain_formation.py --scenario uniform
    python3 scripts/simulate_chain_formation.py --scenario chain --shear-rate-per-s 200 \
        --total-time-s 6 --n-cells 12 --box-size-um 16
    python3 scripts/simulate_chain_formation.py --scenario chain --well-depth-um 20000 \
        --total-time-s 3600 --n-cells 16 --box-size-um 20

`adhesion.dynamics` が返す時間発展を、2次元の細胞の形（カプセル）として描画し、
右側に極接触（鎖）・側面接触（束）の累積本数を時間とともに示す。沈降モードでは
さらに右側に「上清に残る細胞の割合」——実際のアッセイが測定する量——を表示する。

既定では静置（拡散のみ）条件を動画にする。`--shear-rate-per-s` を正の値にすると
震盪（平面クエット流）を課す。ただし Glass & Riedel-Kruse (2018) の実際の凝集
アッセイは震盪中には行われない（育成培養は震盪するが、凝集そのものは短時間の
ボルテックス混合のあと静置で進む）。**実際のアッセイに対応するのは
`--well-depth-um` の沈降モードである。** 震盪モードは「もし振盪し続けたら」
という仮想的な条件として残してある。

震盪モードは接触の持続時間を取りこぼさないよう dt を自動的に細かくするため
（`adhesion.dynamics._safe_dt_for_shear`）、同じ `--total-time-s` でも静置より
ずっと多くのステップ・計算時間を要する。まずは短い `--total-time-s` と
少ない `--n-cells` で試すことを勧める。

沈降モードでは、実際のアッセイの時間スケール（24時間 = 86400秒）を模すには
非常に多くのステップが要る。デモとして現実的な計算時間で済ませるには
`--total-time-s` を数十分〜数時間程度に短く取り、「その時間内で何が起こるか」
として見るとよい（大きな凝集塊ほど速く沈むという定性的な傾向は、時間を短く
取っても変わらない）。

母集団の規模・視野サイズ・総シミュレーション時間は、較正できる経時データが
文献にないため、「動画として見やすい」という基準で選んだ自由パラメータである
（`adhesion.dynamics.DynamicsConfig` の docstring を参照）。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import dynamics as dyn
from adhesion import model

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

CLUSTER_COLORS = [
    "#0E7C8A",
    "#A6377C",
    "#C98A2B",
    "#3F7D3F",
    "#5B5EA6",
    "#B4436C",
    "#4C8577",
    "#8A5A44",
]
MONOMER_COLOR = "#B9BDB6"
POLAR_BOND_COLOR = "#0E7C8A"
LATERAL_BOND_COLOR = "#A6377C"
HEAD_MARKER_COLOR = "#D94E3B"
TAIL_MARKER_COLOR = "#3B6FD9"
SUPERNATANT_COLOR = "#2E7D32"

SCENARIOS = {
    "chain": dict(f_pole=0.95, tip_width_um=0.10),
    "uniform": dict(f_pole=model.UNIFORM_F_POLE, tip_width_um=None),
}


def _capsule_outline_local(radius_um: float, half_length_um: float, n_arc: int = 16) -> np.ndarray:
    """カプセル（スフェロシリンダーの2次元断面）の輪郭を、細胞中心基準で返す。"""

    head_angles = np.linspace(-np.pi / 2, np.pi / 2, n_arc)
    head = np.column_stack(
        (half_length_um + radius_um * np.cos(head_angles), radius_um * np.sin(head_angles))
    )
    tail_angles = np.linspace(np.pi / 2, 3 * np.pi / 2, n_arc)
    tail = np.column_stack(
        (-half_length_um + radius_um * np.cos(tail_angles), radius_um * np.sin(tail_angles))
    )
    return np.vstack((head, tail))


_LOCAL_OUTLINE = _capsule_outline_local(dyn.CELL_RADIUS_UM, dyn.HALF_CYLINDER_UM)


def _capsule_world_outline(position: np.ndarray, theta: float) -> np.ndarray:
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    rotation = np.array([[cos_t, -sin_t], [sin_t, cos_t]])
    return position + _LOCAL_OUTLINE @ rotation.T


def _pole_marker_positions(position: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    axis = np.array([np.cos(theta), np.sin(theta)])
    reach = dyn.HALF_CYLINDER_UM + dyn.CELL_RADIUS_UM
    return position + reach * axis, position - reach * axis


def _unwrap_positions_for_periodic_x(
    n_cells: int, positions: np.ndarray, bonds: list[dyn.Bond], box_size_um: float
) -> np.ndarray:
    """周期境界（震盪モード）でも、結合でつながった細胞どうしが連続した
    x 座標になるように「アンラップ」する。

    そのままの座標で結合線を引くと、鎖が境界をまたいだ瞬間に線が箱の端から端
    まで伸びて見える（周期像の選び方が2細胞でバラバラになるため）。結合グラフを
    たどって、隣の細胞の x を「自分に最も近い像」へ揃えることでこれを防ぐ。
    結果として得られる座標は箱の外にはみ出しうるので、描画側で ±box_size_um
    ずらした「影」も一緒に描いて境界をまたぐ様子を自然に見せる。
    """

    unwrapped = positions.copy()
    adjacency: list[list[int]] = [[] for _ in range(n_cells)]
    for bond in bonds:
        adjacency[bond.i].append(bond.j)
        adjacency[bond.j].append(bond.i)

    visited = [False] * n_cells
    for start in range(n_cells):
        if visited[start]:
            continue
        visited[start] = True
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if visited[neighbor]:
                    continue
                dx = unwrapped[neighbor, 0] - unwrapped[current, 0]
                shift = round(dx / box_size_um) * box_size_um
                unwrapped[neighbor, 0] -= shift
                visited[neighbor] = True
                stack.append(neighbor)
    return unwrapped


def _cluster_colors_for_frame(n_cells: int, frame: dyn.Frame, directed: bool) -> list[str]:
    clusters = dyn.classify_clusters(n_cells, frame.bonds, directed)
    colors = [MONOMER_COLOR] * n_cells
    color_index = 0
    for cluster in clusters:
        if cluster.kind == "monomer":
            continue
        color = CLUSTER_COLORS[color_index % len(CLUSTER_COLORS)]
        color_index += 1
        for cell in cluster.cells:
            colors[cell] = color
    return colors


def _bond_line_segments(
    positions: np.ndarray, angles: np.ndarray, bonds: list[dyn.Bond]
) -> tuple[list[np.ndarray], list[str]]:
    segments = []
    colors = []
    for bond in bonds:
        p_i = dyn.attachment_point(positions[bond.i], angles[bond.i], bond.axial_i, bond.delta_i)
        p_j = dyn.attachment_point(positions[bond.j], angles[bond.j], bond.axial_j, bond.delta_j)
        segments.append(np.array([p_i, p_j]))
        colors.append(POLAR_BOND_COLOR if bond.kind == "polar" else LATERAL_BOND_COLOR)
    return segments, colors


def build_animation(result: dyn.SimulationResult, fps: int) -> animation.FuncAnimation:
    config = result.config
    n_cells = config.n_cells
    periodic_x = config.periodic_x
    settling = config.settling_enabled
    # 震盪モード（周期境界）では、境界をまたぐクラスタが自然に見えるよう、
    # ±box_size_um ずらした「影」も一緒に描く。静置モードでは影を描かず
    # （offsets=[0.0], n_copies=1）、これを追加する前とまったく同じ見た目になる。
    offsets = [-config.box_size_um, 0.0, config.box_size_um] if periodic_x else [0.0]
    n_copies = len(offsets)

    if settling:
        figure, (ax_main, ax_counts, ax_super) = plt.subplots(
            1, 3, figsize=(15.5, 5.0), gridspec_kw={"width_ratios": [1.2, 1.0, 1.0]}
        )
    else:
        figure, (ax_main, ax_counts) = plt.subplots(
            1, 2, figsize=(11.5, 5.4), gridspec_kw={"width_ratios": [1.4, 1.0]}
        )

    ax_main.set_xlim(0, config.box_size_um)
    ax_main.set_ylim(0, config.box_size_um)
    ax_main.set_aspect("equal")
    ax_main.set_xlabel("x (µm)")
    ax_main.set_ylabel("y (µm)")
    time_text = ax_main.set_title("t = 0 s")

    patches = [
        Polygon(_capsule_world_outline(result.frames[0].positions[k % n_cells], result.frames[0].angles[k % n_cells]),
                closed=True, facecolor=MONOMER_COLOR, edgecolor="black", linewidth=0.4, zorder=2)
        for k in range(n_cells * n_copies)
    ]
    for patch in patches:
        ax_main.add_patch(patch)

    head_scatter = ax_main.scatter(
        [], [], s=10, color=HEAD_MARKER_COLOR, zorder=3, label="Ag 極（頭）"
    )
    tail_scatter = ax_main.scatter(
        [], [], s=10, color=TAIL_MARKER_COLOR, zorder=3, label="Nb 極（尾）"
    )
    ax_main.legend(loc="upper right", fontsize=8, frameon=False)

    bond_lines = LineCollection([], linewidths=1.6, zorder=1)
    ax_main.add_collection(bond_lines)

    times = [frame.time_s for frame in result.frames]
    polar_counts = [sum(1 for b in frame.bonds if b.kind == "polar") for frame in result.frames]
    lateral_counts = [sum(1 for b in frame.bonds if b.kind == "lateral") for frame in result.frames]

    ax_counts.set_xlim(0, times[-1])
    ax_counts.set_ylim(0, max(max(polar_counts, default=1), max(lateral_counts, default=1)) + 2)
    ax_counts.set_xlabel("時間 (s)")
    ax_counts.set_ylabel("累積結合数")
    ax_counts.set_title("極接触（鎖）と側面接触（束）", fontsize=10)
    (polar_line,) = ax_counts.plot([], [], color=POLAR_BOND_COLOR, lw=2, label="極接触（鎖）")
    (lateral_line,) = ax_counts.plot([], [], color=LATERAL_BOND_COLOR, lw=2, label="側面接触（束）")
    progress_marker = ax_counts.axvline(0.0, color="black", lw=0.8, ls=":")
    ax_counts.legend(loc="upper left", fontsize=8, frameon=False)
    ax_counts.spines[["top", "right"]].set_visible(False)

    progress_marker_super = None
    if settling:
        supernatant = [
            dyn.supernatant_fraction(frame, config.well_depth_um, config.supernatant_fraction)
            for frame in result.frames
        ]
        ax_super.set_xlim(0, times[-1])
        ax_super.set_ylim(0, 1.05)
        ax_super.set_xlabel("時間 (s)")
        ax_super.set_ylabel("上清に残る細胞の割合")
        ax_super.set_title(
            f"上層 {config.supernatant_fraction:.0%} のOD相当\n"
            "（実際のアッセイの測定量）",
            fontsize=10,
        )
        (super_line,) = ax_super.plot([], [], color=SUPERNATANT_COLOR, lw=2)
        progress_marker_super = ax_super.axvline(0.0, color="black", lw=0.8, ls=":")
        ax_super.spines[["top", "right"]].set_visible(False)

    def update(frame_index: int):
        frame = result.frames[frame_index]
        colors = _cluster_colors_for_frame(n_cells, frame, config.distinct_adhesins_per_pole)
        if periodic_x:
            base_positions = _unwrap_positions_for_periodic_x(
                n_cells, frame.positions, frame.bonds, config.box_size_um
            )
        else:
            base_positions = frame.positions

        # 沈降モードでは、沈んだ細胞ほど淡く描き、実際の顕微鏡・OD測定で
        # 「上層から抜けていく」様子を視覚的に表す。
        if settling:
            alphas = 1.0 - 0.85 * np.clip(frame.depth_um / config.well_depth_um, 0.0, 1.0)
        else:
            alphas = np.ones(n_cells)

        heads = np.zeros((n_cells * n_copies, 2))
        tails = np.zeros((n_cells * n_copies, 2))
        for c, offset in enumerate(offsets):
            shift = np.array([offset, 0.0])
            for k in range(n_cells):
                position = base_positions[k] + shift
                idx = c * n_cells + k
                patches[idx].set_xy(_capsule_world_outline(position, frame.angles[k]))
                patches[idx].set_facecolor(colors[k])
                patches[idx].set_alpha(float(alphas[k]))
                heads[idx], tails[idx] = _pole_marker_positions(position, frame.angles[k])
        head_scatter.set_offsets(heads)
        tail_scatter.set_offsets(tails)

        all_segments: list[np.ndarray] = []
        all_colors: list[str] = []
        for offset in offsets:
            shift = np.array([offset, 0.0])
            segments, seg_colors = _bond_line_segments(
                base_positions + shift, frame.angles, frame.bonds
            )
            all_segments.extend(segments)
            all_colors.extend(seg_colors)
        bond_lines.set_segments(all_segments)
        bond_lines.set_color(all_colors)

        time_text.set_text(f"t = {frame.time_s:.0f} s")
        polar_line.set_data(times[: frame_index + 1], polar_counts[: frame_index + 1])
        lateral_line.set_data(times[: frame_index + 1], lateral_counts[: frame_index + 1])
        progress_marker.set_xdata([frame.time_s, frame.time_s])
        artists = [*patches, head_scatter, tail_scatter, bond_lines, polar_line, lateral_line]
        if settling:
            super_line.set_data(times[: frame_index + 1], supernatant[: frame_index + 1])
            progress_marker_super.set_xdata([frame.time_s, frame.time_s])
            artists.extend([super_line, progress_marker_super])
        return artists

    figure.tight_layout()
    return animation.FuncAnimation(
        figure, update, frames=len(result.frames), interval=1000.0 / fps, blit=False
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", choices=sorted(SCENARIOS), default="chain")
    parser.add_argument("--n-cells", type=int, default=24)
    parser.add_argument("--box-size-um", type=float, default=20.0)
    parser.add_argument("--total-time-s", type=float, default=120.0)
    parser.add_argument("--dt-s", type=float, default=0.02)
    parser.add_argument("--frame-stride", type=int, default=25)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--fps", type=int, default=15)
    parser.add_argument(
        "--shear-rate-per-s",
        type=float,
        default=0.0,
        help=(
            "震盪培養のせん断速度（1/s）。既定 0 は静置・拡散のみ。>0 にすると"
            "平面クエット流（流れ方向を周期境界）を課す。典型的なフラスコ震盪は"
            "100–1000 /s（shear.py 参照）。大きくするほど自動的に dt を細かく"
            "刻むため、計算時間が大きく伸びる点に注意。"
        ),
    )
    parser.add_argument(
        "--well-depth-um",
        type=float,
        default=None,
        help=(
            "沈降モード（A23）を有効にする液柱の深さ（µm）。Glass & "
            "Riedel-Kruse (2018) の実際の凝集アッセイ（震盪ではなく、短時間の"
            "ボルテックス混合 → 深底96穴プレートで静置 → 上層の光学密度で判定）"
            "に対応する。既定 None は水平面内の運動のみ（沈降を計算しない）。"
            "深底96穴プレートの物理的な深さは製品仕様でおよそ 30000–44000 µm。"
        ),
    )
    parser.add_argument(
        "--supernatant-fraction",
        type=float,
        default=0.25,
        help="上清として測定される上層の割合（既定 0.25、原著の上部 約25パーセント）。",
    )
    parser.add_argument(
        "--k-off-per-s",
        type=float,
        default=0.0,
        help=(
            "結合の自発的な解離速度定数（1/s）。既定 0 は Costan et al. (2024) の"
            "best estimate（k_off はゼロ）に合わせた不可逆凝集（A15）。ただしこの"
            "結論は狭い荷重速度範囲での検出限界に基づく主張であり、この結合"
            "（anti-EPEA/EPEA）は他の測定に使われている例外的に安定な"
            "anti-GFP/GFP ペアとは異なる、より低親和性のペアである（A24）。"
            "短時間のデモでは 0 のままでよいが、長時間の沈降モードでは"
            "1e-4〜1e-2 /s 程度の値を試すことを検討するとよい。"
        ),
    )
    parser.add_argument(
        "--swim-speed-um-s",
        type=float,
        default=0.0,
        help=(
            "遊泳速度（µm/s）。既定 0 は能動運動なし（受動拡散のみ）。正の値に"
            "すると、結合していない自由な細胞は run-and-tumble 運動をする"
            "（A25）。使用株 MG1655 は運動性を持つが、既定の速度・run/tumble"
            "時間・旋回角度は Berg & Brown (1972) の野生株 AW405 の値（全文"
            "照合済み）で、別の K12 系統・別の測定条件でのもの。目安は"
            "~14.2 µm/s。"
        ),
    )
    parser.add_argument(
        "--mean-run-time-s", type=float, default=0.86,
        help="run 相の平均継続時間（秒）。既定 0.86（Berg & Brown 1972、全文照合済み）。",
    )
    parser.add_argument(
        "--mean-tumble-time-s", type=float, default=0.14,
        help="tumble 相の平均継続時間（秒）。既定 0.14（同上）。",
    )
    parser.add_argument(
        "--tumble-turn-angle-deg", type=float, default=68.0,
        help="tumble で向きを変える角度の平均（度）。既定 68.0（同上）。",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="既定では results/animations/<scenario>.gif",
    )
    args = parser.parse_args()

    overrides = SCENARIOS[args.scenario]
    config = dyn.DynamicsConfig(
        n_cells=args.n_cells,
        box_size_um=args.box_size_um,
        total_time_s=args.total_time_s,
        dt_s=args.dt_s,
        frame_stride=args.frame_stride,
        seed=args.seed,
        shear_rate_per_s=args.shear_rate_per_s,
        well_depth_um=args.well_depth_um,
        supernatant_fraction=args.supernatant_fraction,
        k_off_per_s=args.k_off_per_s,
        swim_speed_um_s=args.swim_speed_um_s,
        mean_run_time_s=args.mean_run_time_s,
        mean_tumble_time_s=args.mean_tumble_time_s,
        tumble_turn_angle_deg=args.tumble_turn_angle_deg,
        **overrides,
    )

    print(f"シミュレーション実行中（{args.scenario}: f_pole={config.f_pole:.3f}, "
          f"tip_width_um={config.tip_width_um}）...")
    if config.periodic_x:
        effective_dt = dyn._safe_dt_for_shear(config)
        n_steps = int(round(config.total_time_s / effective_dt))
        print(
            f"震盪モード: shear_rate_per_s={config.shear_rate_per_s:.0f} のため "
            f"dt を {effective_dt:.2e} s に自動調整（{n_steps} ステップ）。"
            "流れ方向（x）は周期境界なので、境界を横切る細胞は反対側から"
            "現れて見える（既知の描画上の制約）。注意: Glass & Riedel-Kruse "
            "(2018) の実際の凝集アッセイは震盪中には行われない（--well-depth-um "
            "を参照）。"
        )
    if config.settling_enabled:
        print(
            f"沈降モード: well_depth_um={config.well_depth_um:.0f} "
            f"（実際の凝集アッセイに対応、震盪なし・静置）。"
            f"単量体の沈降速度={dyn._settling_velocity_um_s(1):.4f} µm/s、"
            f"20細胞の凝集塊={dyn._settling_velocity_um_s(20):.4f} µm/s。"
        )
    if config.motile:
        effective_dt_m = dyn._safe_dt_for_motility(config)
        print(
            f"運動性モード: swim_speed_um_s={config.swim_speed_um_s:.1f}（A25、"
            f"Berg & Brown 1972・野生株 AW405 の値、全文照合済み。ただし別の"
            f"K12系統・別の測定条件でのもの）。"
            f"dt を {effective_dt_m:.2e} s に自動調整。"
        )
    result = dyn.run(config)
    print(
        f"極接触の結合数定数={result.polar_bond_count:.1f}, "
        f"側面接触の結合数定数={result.side_bond_count:.1f}"
    )
    final_bonds = result.frames[-1].bonds
    print(
        f"最終フレーム: 結合 {len(final_bonds)} 本"
        f"（極 {sum(1 for b in final_bonds if b.kind == 'polar')},"
        f" 側面 {sum(1 for b in final_bonds if b.kind == 'lateral')}）"
    )

    output = args.output or (
        Path(__file__).resolve().parent.parent
        / "results"
        / "animations"
        / f"{args.scenario}.gif"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    anim = build_animation(result, fps=args.fps)
    anim.save(output, writer=animation.PillowWriter(fps=args.fps))
    print(f"動画を {output} に保存した")


if __name__ == "__main__":
    main()
