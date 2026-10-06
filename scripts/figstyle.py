"""全ての図・動画で共通の配色と体裁。

1つの意味に1色を割り当て、どの図でも同じ意味には同じ色を使う。

    青      極（極接触・極結合・鎖・極のアドヘシン・極領域）
    オレンジ 側面（側面接触・側面結合・束・側面のアドヘシン・円柱側面）
    灰      単独の細胞、報告値、力・寸法・注意などの補助的な意味

青とオレンジは、4 種の色覚（通常・protan・deutan・tritan）のシミュレーション下でも
OKLab ΔE×100 が 24 以上離れる（基準は 8 以上）。以前の配色（ティール / マゼンタ）は
deutan で ΔE 2.6 しかなく区別しにくかったので置き換えた。第3の色相は増やさず、補助的な
意味は濃い灰と線種・ラベルで区別する（色だけに頼らない）。

細胞の塗りは、動画では「その細胞がもつ結合の種類」で決める（灰 = 結合なし、青 = 極結合だけを
もつ、オレンジ = 側面結合を1本でももつ）。かたまり全体ではなく細胞ごとなので、鎖の途中に側面接触が
あっても塗り替わるのはその細胞だけ。Ag 極・Nb 極の目印は色ではなく形で区別する
（● = Ag 極、○ = Nb 極）。
"""

from __future__ import annotations

import numpy as np

# --- 意味に割り当てた色 ----------------------------------------------------------
POLAR = "#2a78d6"  # 極
LATERAL = "#eb6834"  # 側面

#: 文字に使う濃い版（同じ色相。白地でコントラスト 4.5:1 以上）。塗りや線には POLAR / LATERAL。
POLAR_TEXT = "#1c5cab"
LATERAL_TEXT = "#b8430f"

#: 薄い塗り。「薄い側の極」など、同じ意味の弱い版。
POLAR_PALE = "#b9d3f3"
LATERAL_PALE = "#f8c9b5"
#: 青の中間の濃さ。同じ意味の中での強弱（先端の絞り方など）に使う。
POLAR_MID = "#6da7ec"

# --- 中立色 -----------------------------------------------------------------------
INK = "#0b0b0b"
INK_2 = "#52514e"  # 力・寸法・注意など、補助的な意味
MUTED = "#8a8985"  # 報告値、注記
GRID = "#e4e3df"
SURFACE = "#ffffff"
CELL_BODY = "#e8e7e3"  # 細胞の塗り（意味を持たないとき）
CELL_EDGE = "#6e6d69"
MONOMER = "#b9b8b3"  # 単独の細胞
NOTE_FILL = "#f3f2ef"

#: 発散型（P = S を中心に、極が勝つ側を青、側面が勝つ側をオレンジ。中点は中立灰）。
DIVERGING_STOPS = ["#9a3a12", LATERAL, LATERAL_PALE, "#f0efec", "#9ec5f4", POLAR, "#184f95"]


def diverging_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("polar_vs_lateral", DIVERGING_STOPS)


def style_axes(axes) -> None:
    """上と右の枠を消し、目盛りと格子を控えめにする。"""

    axes.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(MUTED)
    axes.tick_params(colors=INK_2, labelsize=9)
    axes.grid(axis="y", color=GRID, lw=0.8)
    axes.set_axisbelow(True)
    axes.xaxis.label.set_color(INK_2)
    axes.yaxis.label.set_color(INK_2)


def end_labels(axes, items: list[tuple[float, float, str]], min_gap_px: float = 13.0,
               color: str = INK_2, fontsize: float = 8.5) -> None:
    """系列の端に直接ラベルを置く。値が近い系列どうしのラベルが重ならないよう縦に散らす。

    items: (x, y, 文字列)。軸のスケールと範囲は先に確定させておくこと。
    """

    figure = axes.figure
    scale = 72.0 / figure.dpi
    displayed = [axes.transData.transform((x, y))[1] for x, y, _ in items]
    order = np.argsort(displayed)
    placed = [0.0] * len(items)
    previous = -np.inf
    for index in order:
        placed[index] = max(displayed[index], previous + min_gap_px)
        previous = placed[index]
    # 散らしたことで全体が上にずれるぶんを、中心が元の位置に戻るよう補正する
    shift = np.mean(np.array(displayed) - np.array(placed))
    for (x, y, text), shown, target in zip(items, displayed, placed):
        axes.annotate(
            text, (x, y), xytext=(7, (target + shift - shown) * scale), textcoords="offset points",
            va="center", fontsize=fontsize, color=color,
        )
