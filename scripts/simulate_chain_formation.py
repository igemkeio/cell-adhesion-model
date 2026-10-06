"""細胞集団が時間とともに衝突し、鎖・束を作っていく過程を動画にする。

    python3 scripts/simulate_chain_formation.py --scenario chain
    python3 scripts/simulate_chain_formation.py --scenario uniform
    python3 scripts/simulate_chain_formation.py --scenario chain --shear-rate-per-s 200 \
        --total-time-s 6 --n-cells 12 --box-size-um 16 --frame-interval-s 0.0125
    python3 scripts/simulate_chain_formation.py --scenario chain --well-depth-um 600 \
        --total-time-s 1400 --n-cells 18 --box-size-um 20 --frame-interval-s 5.8

`adhesion.dynamics` が返す時間発展を、2次元の細胞の形（カプセル）として描画する。
細胞の塗りは、その細胞がもつ結合の種類で決める（灰 = 結合なし、青 = 極結合だけをもつ、
オレンジ = 側面結合を1本でももつ）。● は Ag 極、○ は Nb 極。配色は全ての図・動画で共通
（`scripts/figstyle.py`）。右側に極接触（鎖）・側面接触（束）の累積結合数（分子対の数。接触の箇所数では
なく、静的モデルと同じ「結合数」の定義——1箇所の接触あたり `polar_bond_count`
または `side_bond_count` 個の分子対——で数える）を時間とともに示す。沈降モードでは
さらに右側に「上清に残る細胞の割合」——実際のアッセイが測定する量——を表示する。

既定では静置（拡散のみ）条件を動画にする。`--shear-rate-per-s` を正の値にすると
震盪（平面クエット流）を課す。ただし Glass & Riedel-Kruse (2018) の実際の凝集
アッセイは震盪中には行われない（育成培養は震盪するが、凝集そのものは短時間の
ボルテックス混合のあと静置で進む）。**実際のアッセイに対応するのは
`--well-depth-um` の沈降モードである。** 震盪モードは「もし振盪し続けたら」
という仮想的な条件として残してある。

震盪モードは接触の持続時間を取りこぼさないよう dt を自動的に細かくするため
（`adhesion.dynamics._safe_dt_for_shear`）、同じ `--total-time-s` でも静置より
多くのステップ・計算時間を要する。まずは短い `--total-time-s` と
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
import dataclasses
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Polygon, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import dynamics as dyn
from adhesion import model
from figstyle import INK, INK_2, LATERAL, MONOMER, POLAR

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

# 配色は全ての図・動画で共通（scripts/figstyle.py）。細胞の塗りは、その細胞がもつ結合の種類で
# 決める: 灰 = 結合なし（単独）、青 = 極結合だけをもつ、オレンジ = 側面結合を1本でももつ。
# 結合の種類も同じ色（極 = 青、側面 = オレンジ）。Ag 極・Nb 極の目印は色ではなく形で区別する
# （● = Ag 極、○ = Nb 極）。かたまり全体ではなく細胞ごとに塗るので、鎖の途中に側面接触が
# 1本あっても、塗り替わるのはその2細胞だけで、鎖がほとんど青のまま見える。
MONOMER_COLOR = MONOMER
CHAIN_COLOR = POLAR  # 極結合だけをもつ細胞
BUNDLE_COLOR = LATERAL  # 側面結合を1本でももつ細胞
POLAR_BOND_COLOR = POLAR
LATERAL_BOND_COLOR = LATERAL
SUPERNATANT_COLOR = INK

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


def cell_colors_for_frame(n_cells: int, frame: dyn.Frame) -> list[str]:
    """細胞の塗り。その細胞がもつ結合の種類で決める（結合なし = 灰、極だけ = 青、側面あり = オレンジ）。"""

    has_polar = [False] * n_cells
    has_lateral = [False] * n_cells
    for bond in frame.bonds:
        target = has_polar if bond.kind == "polar" else has_lateral
        target[bond.i] = target[bond.j] = True
    return [
        BUNDLE_COLOR if has_lateral[k] else (CHAIN_COLOR if has_polar[k] else MONOMER_COLOR)
        for k in range(n_cells)
    ]


def legend_handles() -> list:
    """動画の凡例。細胞の塗り（その細胞がもつ結合の種類）と、極の目印。"""

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    return [
        Patch(facecolor=MONOMER_COLOR, edgecolor="black", linewidth=0.4, label="結合なし（単独）"),
        Patch(facecolor=CHAIN_COLOR, edgecolor="black", linewidth=0.4, label="極結合だけをもつ細胞（鎖）"),
        Patch(facecolor=BUNDLE_COLOR, edgecolor="black", linewidth=0.4, label="側面結合をもつ細胞（束）"),
        Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor=INK, markeredgecolor=INK,
               label="Ag 極（頭）"),
        Line2D([], [], marker="o", ls="", markersize=5, markerfacecolor="white", markeredgecolor=INK,
               label="Nb 極（尾）"),
    ]


#: 追跡していない細胞を薄く描くときの不透明度。
FADED_ALPHA = 0.12
#: 追跡する群を決める時刻（総時間に対する割合）。
FOCUS_REFERENCE_FRACTION = 0.2


def find_focus_group(result: dyn.SimulationResult, rule: str = "oldest",
                     fraction: float = FOCUS_REFERENCE_FRACTION) -> dict | None:
    """追跡する細胞群を決める。

    総時間の最初の `fraction`（既定 1/5）の時点で存在する、2細胞以上のかたまりのうち、
    **最も長く結合しているもの**を選ぶ（rule="oldest"）。かたまりの結合期間は、その時点で
    まだ切れていない結合のうち最も古いものの経過時間。同点なら大きいかたまりを選ぶ。
    rule="largest" なら、その時点で最も大きいかたまりを選ぶ。
    結合した時刻は、フレームに結合が現れた最初の時刻（フレーム間隔の精度）。結合が一度切れて
    つなぎ直された場合は、つなぎ直した時刻から数える。
    該当するかたまりが無ければ None（追跡しない）。
    """

    frames = result.frames
    n_cells = result.config.n_cells
    reference = max(1, round((len(frames) - 1) * fraction))
    started: dict[tuple[int, int], float] = {}
    previous: set[tuple[int, int]] = set()
    for k in range(reference + 1):
        current = {(min(b.i, b.j), max(b.i, b.j)) for b in frames[k].bonds}
        for pair in current - previous:
            started[pair] = frames[k].time_s
        for pair in previous - current:
            started.pop(pair, None)
        previous = current
    reference_time = frames[reference].time_s

    best = None
    for component in dyn.connected_components(n_cells, frames[reference].bonds):
        if len(component) < 2:
            continue
        members = set(component)
        ages = [reference_time - started[pair] for pair in previous
                if pair[0] in members and pair[1] in members]
        if not ages:
            continue
        key = (max(ages), len(component)) if rule == "oldest" else (len(component), max(ages))
        if best is None or key > best[0]:
            best = (key, sorted(component), max(ages))
    if best is None:
        return None
    return {
        "members": best[1], "reference_index": reference, "reference_time_s": reference_time,
        "age_s": best[2], "rule": rule,
    }


def focus_camera(result: dyn.SimulationResult, members: list[int]) -> np.ndarray:
    """各フレームで画面の中心に置く点 (cx, cy)。追跡する細胞の重心。

    x は周期境界なので、円周上の平均（箱の長さを1周とみなした平均の向き）で求める。
    普通の平均だと、境界をまたいだ群の重心が箱の反対側に出てしまう。
    """

    box = result.config.box_size_um
    camera = np.zeros((len(result.frames), 2))
    for k, frame in enumerate(result.frames):
        angle = 2.0 * np.pi * frame.positions[members, 0] / box
        cx = (np.arctan2(np.sin(angle).mean(), np.cos(angle).mean()) / (2.0 * np.pi)) % 1.0 * box
        camera[k] = (cx, frame.positions[members, 1].mean())
    return camera


def tracked_cells(frame: dyn.Frame, members: list[int], n_cells: int) -> set[int]:
    """追跡する細胞。最初に決めた群の細胞がいま属しているかたまり（後から合流した細胞を含む）。"""

    wanted = set(members)
    tracked: set[int] = set()
    for component in dyn.connected_components(n_cells, frame.bonds):
        if wanted.intersection(component):
            tracked.update(component)
    return tracked


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


def build_animation(result: dyn.SimulationResult, fps: int, focus: dict | None = None,
                    focus_window_um: float = 12.0) -> animation.FuncAnimation:
    """動画を作る。

    focus（`find_focus_group` の結果）を渡すと、その群の重心が常に画面の中心に来る動く視点で描き、
    追跡している群以外の細胞は薄く描く。震盪培養では細胞が速く流れて固定した視点では見づらいため。
    座標は群の中心からの相対位置になり、x は周期境界を考慮して折り返す。壁（y = 0 と y = 箱の
    大きさ）は線で示す。
    """

    config = result.config
    n_cells = config.n_cells
    periodic_x = config.periodic_x
    settling = config.settling_enabled
    # 震盪モード（周期境界）で固定した視点のときは、境界をまたぐクラスタが自然に見えるよう、
    # ±box_size_um ずらした「影」も一緒に描く。静置モードや、群を追う視点では影を描かない
    # （offsets=[0.0], n_copies=1）。
    offsets = [-config.box_size_um, 0.0, config.box_size_um] if periodic_x and focus is None else [0.0]
    n_copies = len(offsets)
    camera = focus_camera(result, focus["members"]) if focus else None
    tracked_by_frame = (
        [tracked_cells(f, focus["members"], n_cells) for f in result.frames] if focus else None
    )
    box = config.box_size_um

    if settling:
        figure, (ax_main, ax_counts, ax_super) = plt.subplots(
            1, 3, figsize=(15.5, 5.0), gridspec_kw={"width_ratios": [1.2, 1.0, 1.0]}
        )
    else:
        figure, (ax_main, ax_counts) = plt.subplots(
            1, 2, figsize=(11.5, 5.4), gridspec_kw={"width_ratios": [1.4, 1.0]}
        )

    if focus:
        half = 0.5 * focus_window_um
        ax_main.set_xlim(-half, half)
        ax_main.set_ylim(-half, half)
        ax_main.set_xlabel("群の中心からの x (µm)")
        ax_main.set_ylabel("群の中心からの y (µm)")
        # 壁の外側（y < 0 と y > 箱の大きさ）。群が動くので、フレームごとに位置を更新する。
        wall_below = Rectangle((-box, -10 * box), 2 * box, 10 * box, facecolor="#efeeea",
                               edgecolor="none", hatch="///", zorder=0)
        wall_above = Rectangle((-box, 0), 2 * box, 10 * box, facecolor="#efeeea",
                               edgecolor="none", hatch="///", zorder=0)
        ax_main.add_patch(wall_below)
        ax_main.add_patch(wall_above)
        wall_lines = [ax_main.axhline(0.0, color=INK_2, lw=1.6, zorder=1) for _ in range(2)]
        ax_main.text(
            0.02, 0.02,
            f"追跡: t = {focus['reference_time_s']:.2f} s で最も長く結合していた群"
            if focus["rule"] == "oldest" else f"追跡: t = {focus['reference_time_s']:.2f} s で最大の群",
            transform=ax_main.transAxes, fontsize=8, color=INK_2, zorder=6,
        )
    else:
        ax_main.set_xlim(0, config.box_size_um)
        ax_main.set_ylim(0, config.box_size_um)
        ax_main.set_xlabel("x (µm)")
        ax_main.set_ylabel("y (µm)")
    ax_main.set_aspect("equal")
    time_text = ax_main.set_title("t = 0 s")

    patches = [
        Polygon(_capsule_world_outline(result.frames[0].positions[k % n_cells], result.frames[0].angles[k % n_cells]),
                closed=True, facecolor=MONOMER_COLOR, edgecolor="black", linewidth=0.4, zorder=2)
        for k in range(n_cells * n_copies)
    ]
    for patch in patches:
        ax_main.add_patch(patch)

    head_scatter = ax_main.scatter([], [], s=12, color=INK, zorder=3)
    tail_scatter = ax_main.scatter(
        [], [], s=12, facecolor="white", edgecolor=INK, linewidth=0.7, zorder=3
    )
    handles = legend_handles()
    if focus:
        from matplotlib.patches import Patch

        handles.append(Patch(facecolor=MONOMER_COLOR, edgecolor="black", linewidth=0.4,
                             alpha=FADED_ALPHA * 2, label="追跡していない細胞（薄く表示）"))
    ax_main.legend(
        handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
        fontsize=8, frameon=False, columnspacing=1.2, handletextpad=0.5,
    )

    bond_lines = LineCollection([], linewidths=1.6, zorder=1)
    ax_main.add_collection(bond_lines)

    times = [frame.time_s for frame in result.frames]
    # 「結合数」は静的モデルと同じ定義（接触1箇所あたりの分子対の数）で数える。
    # 接触の箇所数（＝くっついている細胞の対の数）をそのまま結合数として数える
    # と、局在によって1接触あたりの分子対数が大きく異なる（例: 先端集中で
    # 極接触1箇所が ~265 対）ことを反映できず、過小評価になる。
    polar_counts = [
        sum(1 for b in frame.bonds if b.kind == "polar") * result.polar_bond_count
        for frame in result.frames
    ]
    lateral_counts = [
        sum(1 for b in frame.bonds if b.kind == "lateral") * result.side_bond_count
        for frame in result.frames
    ]

    ax_counts.set_xlim(0, times[-1])
    ax_counts.set_ylim(0, max(max(polar_counts, default=1), max(lateral_counts, default=1)) * 1.1 + 2)
    ax_counts.set_xlabel("時間 (s)")
    ax_counts.set_ylabel("累積結合数（分子対）")
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

    def update_focused(frame_index: int):
        """群を追う視点。群の中心を原点に、x は周期境界で折り返して描く。"""

        frame = result.frames[frame_index]
        cx, cy = camera[frame_index]
        dx = frame.positions[:, 0] - cx
        shown = np.column_stack((dx - box * np.round(dx / box), frame.positions[:, 1] - cy))
        colors = cell_colors_for_frame(n_cells, frame)
        tracked = tracked_by_frame[frame_index]
        alphas = np.array([1.0 if k in tracked else FADED_ALPHA for k in range(n_cells)])

        heads = np.zeros((n_cells, 2))
        tails = np.zeros((n_cells, 2))
        for k in range(n_cells):
            patches[k].set_xy(_capsule_world_outline(shown[k], frame.angles[k]))
            patches[k].set_facecolor(colors[k])
            patches[k].set_alpha(float(alphas[k]))
            heads[k], tails[k] = _pole_marker_positions(shown[k], frame.angles[k])
        head_scatter.set_offsets(heads)
        head_scatter.set_facecolors([(*to_rgb(INK), a) for a in alphas])
        tail_scatter.set_offsets(tails)
        tail_scatter.set_facecolors([(1.0, 1.0, 1.0, a) for a in alphas])
        tail_scatter.set_edgecolors([(*to_rgb(INK), a) for a in alphas])
        bond_lines.set_segments([])  # 結合線は細胞の下でほぼ見えない。群を追う視点では描かない

        wall_below.set_y(-cy - 10 * box)
        wall_above.set_y(box - cy)
        wall_lines[0].set_ydata([-cy, -cy])
        wall_lines[1].set_ydata([box - cy, box - cy])

        time_text.set_text(f"t = {frame.time_s:.2f} s　追跡中の群 {len(tracked)} 細胞")
        polar_line.set_data(times[: frame_index + 1], polar_counts[: frame_index + 1])
        lateral_line.set_data(times[: frame_index + 1], lateral_counts[: frame_index + 1])
        progress_marker.set_xdata([frame.time_s, frame.time_s])
        return [*patches, head_scatter, tail_scatter, bond_lines, polar_line, lateral_line,
                wall_below, wall_above, *wall_lines]

    def update(frame_index: int):
        if focus:
            return update_focused(frame_index)
        frame = result.frames[frame_index]
        colors = cell_colors_for_frame(n_cells, frame)
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
    parser.add_argument(
        "--dt-s", type=float, default=None,
        help="時間刻み（s）。省略すると DynamicsConfig の既定（収束した値、0.0025 s）。",
    )
    parser.add_argument(
        "--frame-interval-s", type=float, default=0.5,
        help=(
            "動画の1コマに当たる模擬時間（s）。既定 0.5。震盪・運動性モードでは dt が自動で"
            "細かくなるので、見たい細かさに応じて小さくすること（震盪の例は README 参照）。"
        ),
    )
    parser.add_argument(
        "--focus", choices=("auto", "oldest", "largest", "off"), default="auto",
        help=(
            "群を追う視点にする。振盪培養は細胞が速く流れて固定した視点では見づらいので、"
            "総時間の最初の 1/5 の時点で『最も長く結合している群』（oldest）または『最大の群』"
            "（largest）を決め、その重心を画面の中心に置いて追い、他の細胞は薄く描く。"
            "既定の auto は、震盪培養（--shear-rate-per-s > 0）のときだけ oldest、それ以外は off。"
        ),
    )
    parser.add_argument(
        "--focus-window-um", type=float, default=12.0,
        help="群を追う視点での、画面の一辺（µm）。既定 12。",
    )
    parser.add_argument(
        "--frame-stride", type=int, default=None,
        help="フレームを記録するステップ間隔。指定すると --frame-interval-s より優先する。",
    )
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
            "遊泳速度（µm/s）。既定 0 は能動運動なし（受動拡散のみ）で、これが"
            "推奨。正の値にすると、結合していない自由な細胞は run-and-tumble "
            "運動をする（A25、任意の探索ツール）。実際のアッセイ濃度では拡散だけで"
            "凝集は数分〜数十分で進むため能動運動は必須ではなく、アッセイ条件"
            "（定常期・静置・アドヘシン負荷）で運動性が保たれるかも疑わしい。"
            "速度・run/tumble時間・旋回角度の既定値は Berg & Brown (1972) の野生株"
            "AW405（全文照合済み、別系統・別条件）。有効にする場合の目安は ~14.2。"
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
    dt_kwargs = {} if args.dt_s is None else {"dt_s": args.dt_s}
    config = dyn.DynamicsConfig(
        n_cells=args.n_cells,
        box_size_um=args.box_size_um,
        total_time_s=args.total_time_s,
        seed=args.seed,
        shear_rate_per_s=args.shear_rate_per_s,
        well_depth_um=args.well_depth_um,
        supernatant_fraction=args.supernatant_fraction,
        k_off_per_s=args.k_off_per_s,
        swim_speed_um_s=args.swim_speed_um_s,
        mean_run_time_s=args.mean_run_time_s,
        mean_tumble_time_s=args.mean_tumble_time_s,
        tumble_turn_angle_deg=args.tumble_turn_angle_deg,
        **dt_kwargs,
        **overrides,
    )
    # コマ間隔は模擬時間で指定する（実際に使う dt は震盪・運動性で自動的に細かくなる）。
    effective_dt = min(dyn._safe_dt_for_shear(config), dyn._safe_dt_for_motility(config))
    stride = args.frame_stride or max(1, round(args.frame_interval_s / effective_dt))
    config = dataclasses.replace(config, frame_stride=stride)
    print(f"時間刻み dt = {effective_dt:.2e} s、{stride} ステップ（{stride * effective_dt:.3g} s）ごとに1コマ")

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
    n_polar_contacts = sum(1 for b in final_bonds if b.kind == "polar")
    n_lateral_contacts = sum(1 for b in final_bonds if b.kind == "lateral")
    print(
        f"最終フレーム: 接触 {len(final_bonds)} 箇所"
        f"（極 {n_polar_contacts} 箇所, 側面 {n_lateral_contacts} 箇所）。"
        f"結合数（分子対）に換算すると"
        f" 極 {n_polar_contacts * result.polar_bond_count:.0f} 対,"
        f" 側面 {n_lateral_contacts * result.side_bond_count:.0f} 対"
        f"（1箇所あたり 極 {result.polar_bond_count:.1f} 対, 側面 {result.side_bond_count:.1f} 対）"
    )

    output = args.output or (
        Path(__file__).resolve().parent.parent
        / "results"
        / "animations"
        / f"{args.scenario}.gif"
    )
    output.parent.mkdir(parents=True, exist_ok=True)

    rule = args.focus
    if rule == "auto":
        rule = "oldest" if config.periodic_x else "off"
    focus = None
    if rule != "off":
        focus = find_focus_group(result, rule)
        if focus is None:
            print("追跡する群が無い（最初の 1/5 の時点で結合したかたまりが無い）ので、固定した視点で描く。")
        else:
            print(
                f"群を追う視点: t = {focus['reference_time_s']:.2f} s で"
                f"{'最も長く結合している' if rule == 'oldest' else '最大の'}群 "
                f"（細胞 {focus['members']}、結合から {focus['age_s']:.2f} s）の重心を画面の中心に置く。"
            )

    anim = build_animation(result, fps=args.fps, focus=focus, focus_window_um=args.focus_window_um)
    anim.save(output, writer=animation.PillowWriter(fps=args.fps))
    print(f"動画を {output} に保存した")


if __name__ == "__main__":
    main()
