"""形状スイープの文書（`docs/shape_scan.md`, `docs/shape_movies.md`）をデータから生成する。

    python3 scripts/export_shape_docs.py

文章中の数値も表も、`results/data/shape_scan_*.csv`・`time_step_convergence.csv`・
`results/animations/shape_*.json` から埋め込む。手で書き写さないことで、データと文書の
食い違いをなくす。

文章には「頑健か」「向きが保たれるか」といった判断も含まれる。判断はデータから導けないので
本文に固定してあるが、データが判断と食い違ったときに気づけるよう、`check` で条件を確かめる。
食い違えばエラーで止まるので、データを見直して本文を直すこと。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import plot_shape_scan as ps  # noqa: E402
from adhesion import dynamics as dyn  # noqa: E402
from adhesion.shape import CellShape  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
ANIMATIONS = ROOT / "results" / "animations"
TIME_STEP_CSV = ROOT / "results" / "data" / "time_step_convergence.csv"

L_SHAPES = [(c, 0.4) for c in (0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0)]
R_SHAPES = [(1.0, r) for r in (0.25, 0.3, 0.4, 0.5, 0.6, 0.8)]
LITERATURE = (1.0, 0.4)


class ClaimContradicted(SystemExit):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise ClaimContradicted(f"本文の判断がデータと食い違う: {message}\n"
                                "データを見直して、scripts/export_shape_docs.py の本文を直すこと。")


class Data:
    def __init__(self, dynamic_csv: Path) -> None:
        self.t = ps.aggregate(ps.load_dynamic(dynamic_csv))
        self.static = [{k: float(v) for k, v in r.items()} for r in ps.read_csv(ps.STATIC_CSV)]

    def entry(self, exp: str, scen: str, c: float, r: float) -> dict | None:
        return self.t.get((exp, scen, c, r))

    def mean(self, exp, scen, c, r, metric) -> float:
        e = self.entry(exp, scen, c, r)
        return math.nan if e is None else e[f"{metric}_mean"]

    def sem(self, exp, scen, c, r, metric) -> float:
        e = self.entry(exp, scen, c, r)
        return math.nan if e is None else e[f"{metric}_sem"]

    def text(self, exp, scen, c, r, metric, digits=2) -> str:
        e = self.entry(exp, scen, c, r)
        if e is None or e[f"{metric}_n"] == 0:
            return "—"
        return f"{e[f'{metric}_mean']:.{digits}f}±{e[f'{metric}_sem']:.{digits}f}"

    def static_row(self, cylinder: float, radius: float) -> dict:
        return next(r for r in self.static
                    if abs(r["cylinder_um"] - cylinder) < 1e-9 and abs(r["radius_um"] - radius) < 1e-9)


def mark(c, r) -> str:
    return " **(文献値)**" if (c, r) == LITERATURE else ""


# ---------------------------------------------------------------------------
# 表
# ---------------------------------------------------------------------------
def composition_table(d: Data, shapes, scan: str) -> str:
    head = "| 円柱長 l_c (µm) | 全長 (µm) |" if scan == "L" else "| 半径 r (µm) | 曲率 1/r (1/µm) |"
    lines = [head + " P | S | 全て反映 | 結合数を既定形状に固定 | 箱を固定 | 発現量∝膜面積 | 20 s の極接触 | 20 s の側面接触 |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for c, r in shapes:
        e = d.entry("main", "localized", c, r)
        first = f"{c:g}" if scan == "L" else f"{r:g}"
        second = f"{c + 2 * r:.2f}" if scan == "L" else f"{1 / r:.2f}"
        lines.append(
            f"| {first}{mark(c, r)} | {second} | {e['polar_bonds_P']:.0f} | {e['side_bonds_S']:.1f} | "
            f"{d.text('main', 'localized', c, r, 'polar_frac_first10')} | "
            f"{d.text('decoupled', 'localized', c, r, 'polar_frac_first10')} | "
            f"{d.text('number_density', 'localized', c, r, 'polar_frac_first10')} | "
            f"{d.text('per_area', 'localized', c, r, 'polar_frac_first10')} | "
            f"{d.text('main', 'localized', c, r, 'polar_t20', 1)} | {d.text('main', 'localized', c, r, 'lateral_t20', 1)} |")
    return "\n".join(lines)


def uniform_table(d: Data, shapes, scan: str) -> str:
    head = "| 円柱長 l_c (µm) |" if scan == "L" else "| 半径 r (µm) |"
    lines = [head + " P | S | 最初の10本のうち極結合の割合 | 20 s の極接触 | 20 s の側面接触 |", "|---:|---:|---:|---:|---:|---:|"]
    for c, r in shapes:
        e = d.entry("main", "uniform", c, r)
        first = f"{c:g}" if scan == "L" else f"{r:g}"
        lines.append(f"| {first}{mark(c, r)} | {e['polar_bonds_P']:.1f} | {e['side_bonds_S']:.0f} | "
                     f"{d.text('main', 'uniform', c, r, 'polar_frac_first10')} | "
                     f"{d.text('main', 'uniform', c, r, 'polar_t20', 1)} | {d.text('main', 'uniform', c, r, 'lateral_t20', 1)} |")
    return "\n".join(lines)


def speed_table(d: Data, shapes, scan: str) -> str:
    head = "| 円柱長 l_c (µm) |" if scan == "L" else "| 半径 r (µm) |"
    lines = [head + " 面積占有率を固定（箱 µm） | 局在 | 一様 | 箱 20 µm 固定・局在 |", "|---:|---:|---:|---:|---:|"]
    for c, r in shapes:
        e = d.entry("main", "localized", c, r)
        first = f"{c:g}" if scan == "L" else f"{r:g}"
        lines.append(f"| {first}{mark(c, r)} | {e['box_um']:g} | {d.text('main', 'localized', c, r, 't_half_monomers_s', 1)} | "
                     f"{d.text('main', 'uniform', c, r, 't_half_monomers_s', 1)} | "
                     f"{d.text('number_density', 'localized', c, r, 't_half_monomers_s', 1)} |")
    return "\n".join(lines)


def rotdiff_table(d: Data) -> str:
    lines = ["| 円柱長 l_c (µm) | 等価球近似（既定） | 回転拡散 ×0.5 | 回転拡散 ×0.25 |", "|---:|---:|---:|---:|"]
    for c in (0.1, 1.0, 3.0):
        lines.append(f"| {c:g}{mark(c, 0.4)} | " + " | ".join(
            d.text(name, "localized", c, 0.4, "polar_frac_first10") for name in ("main", "rotdiff_0.5", "rotdiff_0.25")) + " |")
    return "\n".join(lines)


def structure_table(d: Data) -> str:
    lines = ["| 形状 | 純粋な鎖（3細胞以上） | 二量体 | 束（側面結合を含む） | 単独 | 最大クラスタ |", "|---|---:|---:|---:|---:|---:|"]
    for c, r in [(0.1, 0.4), (1.0, 0.4), (3.0, 0.4), (1.0, 0.25), (1.0, 0.8)]:
        e = d.entry("main", "localized", c, r)
        lines.append(f"| l_c={c:g}, r={r:g}{mark(c, r)} | {e['chain_t30_mean']:.1f} | {e['dimer_t30_mean']:.1f} | "
                     f"{e['bundle_t30_mean']:.1f} | {e['mono_t30_mean']:.1f} | {e['largest_t30_mean']:.1f} |")
    return "\n".join(lines)


def time_step_table() -> str:
    if not TIME_STEP_CSV.exists():
        print(f"警告: {TIME_STEP_CSV.name} が無い（先に scripts/check_time_step.py を実行すること）")
        return "（未実行: `python3 scripts/check_time_step.py`）"
    rows = list(csv.DictReader(TIME_STEP_CSV.open(encoding="utf-8")))
    cases = sorted({(float(r["cylinder_um"]), float(r["radius_um"])) for r in rows})
    n = len([r for r in rows if (r["cylinder_um"], r["radius_um"], r["dt_s"]) == (rows[0]["cylinder_um"], rows[0]["radius_um"], rows[0]["dt_s"])])
    lines = [f"| 形 (l_c, r) | dt (s) | δ | 20 s の側面結合 | 20 s の極結合の割合 |", "|---|---:|---:|---:|---:|"]
    for c, r in cases:
        dts = sorted({float(x["dt_s"]) for x in rows if (float(x["cylinder_um"]), float(x["radius_um"])) == (c, r)}, reverse=True)
        for dt in dts:
            group = [x for x in rows if (float(x["cylinder_um"]), float(x["radius_um"]), float(x["dt_s"])) == (c, r, dt)]
            lat = np.array([float(g["lateral"]) for g in group])
            frac = np.array([float(g["polar"]) / max(float(g["polar"]) + float(g["lateral"]), 1) for g in group])
            sem = lambda a: a.std(ddof=1) / math.sqrt(len(a))
            tag = "（既定）" if dt == dyn.DynamicsConfig().dt_s else ("（旧既定）" if dt == 0.02 else "")
            lines.append(f"| ({c:g}, {r:g}) | {dt:g}{tag} | {float(group[0]['delta']):.2f} | "
                         f"{lat.mean():.2f}±{sem(lat):.2f} | {frac.mean():.3f}±{sem(frac):.3f} |")
    lines.append(f"\n各行 {n} 回の平均 ± 標準誤差（局在、20 s 時点）。")
    return "\n".join(lines)


def time_step_commentary() -> str:
    """時間刻みの表の読み方。数値は CSV から計算する。"""

    if not TIME_STEP_CSV.exists():
        return ""
    rows = list(csv.DictReader(TIME_STEP_CSV.open(encoding="utf-8")))

    def stat(c, r, dt, key):
        values = np.array([float(x[key]) if key != "fraction" else
                           float(x["polar"]) / max(float(x["polar"]) + float(x["lateral"]), 1)
                           for x in rows if (float(x["cylinder_um"]), float(x["radius_um"]), float(x["dt_s"])) == (c, r, dt)])
        return values.mean(), values.std(ddof=1) / math.sqrt(len(values))

    fine = dyn.DynamicsConfig().dt_s
    lat_old, lat_old_se = stat(1.0, 0.4, 0.02, "lateral")
    lat_new, lat_new_se = stat(1.0, 0.4, fine, "lateral")
    mono_a, mono_a_se = stat(1.0, 0.25, fine, "monomers")
    mono_b, mono_b_se = stat(1.0, 0.25, fine / 2, "monomers")
    frac_a, _ = stat(1.0, 0.25, fine, "fraction")
    frac_b, _ = stat(1.0, 0.25, fine / 2, "fraction")
    check(lat_old > lat_new + lat_old_se, "旧既定の dt で側面結合が過大に数えられていない")
    return f"""読み方:

- 文献値の形では、旧既定（dt = 0.02 s、δ = 0.73）の側面結合は {lat_old:.2f}、既定（{fine:g} s）は {lat_new:.2f} で、
  旧既定は約 {100 * (lat_old / lat_new - 1):.0f}% 多い。dt を 0.005 s 以下にすると値が一致する。
- 細い形（l_c = 1.0, r = 0.25）では、{fine:g} s と {fine / 2:g} s で側面結合数と極結合の割合
  （{frac_a:.3f} と {frac_b:.3f}）は標準誤差内で一致する。**ただし単量体の数は {mono_a:.2f}±{mono_a_se:.2f} と
  {mono_b:.2f}±{mono_b_se:.2f} で、{fine:g} s でもまだ動く。** 結合の構成比は収束しているが、細い形の
  凝集の速さには dt の残差がありうる（速さの形依存性は元々結論にしていない）。
- 標準誤差は側面結合数で約 {lat_new_se:.1f}（5%）なので、5〜10% の差は読み取れない。"""


# ---------------------------------------------------------------------------
# 本文
# ---------------------------------------------------------------------------
def scan_doc(d: Data) -> str:
    dt = dyn.DynamicsConfig().dt_s
    f = lambda exp, c, r=0.4, scen="localized": d.mean(exp, scen, c, r, "polar_frac_first10")
    se = lambda exp, c, r=0.4, scen="localized": d.sem(exp, scen, c, r, "polar_frac_first10")

    # --- 判断の確認 ---------------------------------------------------------------------------
    f01, f10, f30 = f("main", 0.1), f("main", 1.0), f("main", 3.0)
    check(f01 - f10 > 2 * math.hypot(se("main", 0.1), se("main", 1.0)) and
          f10 - f30 > 2 * math.hypot(se("main", 1.0), se("main", 3.0)),
          f"円柱長 0.1 → 1.0 → 3.0 で極結合の割合が（標準誤差の2倍以上で）減らない: {f01:.2f}, {f10:.2f}, {f30:.2f}")
    dec_gap = max(abs(f("main", c) - f("decoupled", c)) for c, _ in L_SHAPES)
    check(dec_gap <= 0.12, f"結合数を固定した曲線が全て反映と一致しない（最大差 {dec_gap:.2f}）")
    for exp, label in (("number_density", "箱固定"), ("per_area", "発現量∝膜面積")):
        check(f(exp, 0.1) - f(exp, 3.0) >= 0.2, f"{label}で円柱長 0.1 → 3.0 の減少が 0.2 未満")
    for exp in ("main", "rotdiff_0.5", "rotdiff_0.25"):
        check(f(exp, 0.1) > f(exp, 1.0) > f(exp, 3.0), f"回転拡散の確認（{exp}）で円柱長に対して単調に減らない")
    uniform_max = max(d.mean("main", "uniform", c, r, "polar_frac_first10") for c, r in L_SHAPES + R_SHAPES)
    check(uniform_max <= 0.2, f"一様分布の極結合の割合が 0.2 を超える（{uniform_max:.2f}）")
    # 半径: 細い側（r = 0.25）では全ての置き方で極結合の割合が下がり、そのほとんどは結合数の変化による。
    # 太い側（r = 0.8）では、全て反映では平らで、置き方によって散る。
    variants = ("main", "decoupled", "number_density", "per_area")
    thin = {v: f(v, 1.0, 0.25) - f(v, 1.0, 0.4) for v in variants}
    fat = {v: f(v, 1.0, 0.8) - f(v, 1.0, 0.4) for v in variants}
    check(thin["main"] < -2 * math.hypot(se("main", 1.0, 0.25), se("main", 1.0, 0.4)),
          f"細い側（r = 0.25）で全て反映の極結合の割合が下がっていない: {thin['main']:+.2f}")
    check(all(thin[v] < 0 for v in variants), f"細い側で下がらない置き方がある: {thin}")
    check(thin["main"] - thin["decoupled"] <= -0.08,
          f"細い側の低下が結合数の変化によるとは言えない（全て反映 {thin['main']:+.2f}、結合数を固定 {thin['decoupled']:+.2f}）")
    check(abs(fat["main"]) <= 0.08, f"太い側（r = 0.8）で全て反映が平らでない: {fat['main']:+.2f}")
    check(fat["decoupled"] <= -0.1, f"結合数を固定すると太い側で下がる、という記述と食い違う: {fat['decoupled']:+.2f}")
    fat_values = [f(v, 1.0, 0.8) for v in variants]
    check(max(fat_values) - min(fat_values) >= 0.1, f"太い側で置き方による散らばりが 0.1 未満: {fat_values}")
    th = lambda exp, c, r=0.4: d.mean(exp, "localized", c, r, "t_half_monomers_s")
    sign_main = th("main", 3.0) - th("main", 0.1)
    sign_box = th("number_density", 3.0) - th("number_density", 0.1)
    check(sign_main * sign_box < 0, f"凝集の速さの円柱長依存の符号が箱の取り方で逆転しない: {sign_main:.1f}, {sign_box:.1f}")
    lit = d.entry("main", "localized", *LITERATURE)
    r_values = [f("main", 1.0, r) for _, r in R_SHAPES]
    check(lit["chain_t30_mean"] <= 5 and lit["bundle_t30_mean"] >= 14,
          f"文献値の形で純粋な鎖が少ない・束が多いという記述と食い違う（鎖 {lit['chain_t30_mean']:.1f}, 束 {lit['bundle_t30_mean']:.1f}）")

    # 最も円柱が長い形で、極結合と側面結合が同数に近づく（逆転はしない）
    p30 = d.mean("main", "localized", 3.0, 0.4, "polar_t20")
    l30 = d.mean("main", "localized", 3.0, 0.4, "lateral_t20")
    check(abs(f30 - 0.5) <= 0.1, f"円柱 3.0 µm で極結合の割合が 0.5 に近くない: {f30:.2f}")

    ps_ratio = lambda c, r: d.static_row(c, r)["localized_polar"] / d.static_row(c, r)["localized_side"]
    ratio_l = [ps_ratio(c, 0.4) for c in (0.1, 1.0, 4.0)]
    check(max(ratio_l) / min(ratio_l) < 1.05, "局在の P/S が円柱長で一定でない")
    uni_ratio = lambda c: d.static_row(c, 0.4)["uniform_polar"] / d.static_row(c, 0.4)["uniform_side"]
    pr = lambda r: d.static_row(1.0, r)["localized_polar"]
    sr = lambda r: d.static_row(1.0, r)["localized_side"]

    return f"""# 細胞の形（円柱の長さ・半球の曲率）を変えたときの結果

大腸菌の形を、時間発展モデル（`adhesion.dynamics`）のパラメータとして振った結果。
円柱部分の長さと半球の曲率を変えたとき、局在させたアドヘシンが作る構造（鎖か束か）と
凝集の速さがどう変わるかを調べた。結論の確かさは条件ごとに違うので、**頑健なものと
そうでないものを分けて**書く。動画の見方は [`shape_movies.md`](shape_movies.md)。

> このページの数値と表は `scripts/export_shape_docs.py` がデータから生成する。直接編集しない。

## 何を振ったか

細胞はスフェロシリンダー（円柱の両端に半球）で、独立な寸法は2つ。

| 量 | 意味 | 文献値 |
|---|---|---|
| 円柱長 l_c | 円柱部分の長さ。全長は l = l_c + 2r | 1.0 µm（全長 1.8 µm） |
| 半径 r | 円柱の半径であり、半球の曲率半径（曲率 = 1/r） | 0.40 µm（曲率 2.5 /µm） |

**半球の曲率は円柱の太さと切り離せない**（半球の半径は円柱の半径と等しい）。
そこで全長 l を固定せず、円柱長 l_c を固定して r だけを動かすことで、「円柱長」と
「曲率」を独立な軸にした。半球を円柱と違う曲率にする（楕円体の端など）形は扱っていない。

形は `adhesion.shape.CellShape` に集め、静的モデル（`model.predict(shape=...)`）と
時間発展（`DynamicsConfig(shape=...)`）が同じ寸法を使う。**形を指定しなければ、形状を
パラメータにする前と同じ計算になる**（旧既定の dt = 0.02 s を指定すれば、静置・一様・震盪・
沈降・運動性・解離の6つの運転モードで、変更前後の結果がビット単位で一致することを確認した）。

## 要点

**頑健な結論**（複数の対照条件で同じ向き）

1. **円柱が長いほど、局在させても極結合（鎖の背骨）の割合が下がり、側面結合（束）が増える。**
   最初の10本のうち極結合の割合は、l_c = 0.1 → 1.0 → 3.0 µm で **{f01:.2f} → {f10:.2f} → {f30:.2f}**。
   l_c = 3.0 µm（全長 3.8 µm）では極結合と側面結合がほぼ同数になる（20 s 時点の本数は
   極 {p30:.1f}・側面 {l30:.1f}。掃引の範囲では逆転しない）。
2. **この効果は結合数の変化ではなく、細胞の出会い方（衝突の幾何）による。** 局在させた
   場合の結合数 P・S は円柱長にほとんど依らない（P/S ≈ {np.mean(ratio_l):.0f} で一定）。しかも P・S を既定形状の値に
   固定して形だけを変えても、同じ曲線になる（全て反映との差は最大 {dec_gap:.2f}）。箱を固定した場合や、
   回転拡散を 0.5 倍・0.25 倍にした場合も、向きは変わらない。側面が長いほど出会いやすい標的になる、
   という解釈と整合するが、接触の種類ごとの頻度を直接測ったわけではない。
3. **一様分布は形に関わらず束になる**（極結合の割合は {uniform_max:.2f} 以下）。極結合が目立つのは
   ほぼ球（l_c = 0.1 µm）のときだけ。
4. **静的モデル**: 局在させた場合、極結合数 P は r にほぼ比例し（{pr(0.25):.0f} → {pr(0.4):.0f} → {pr(0.8):.0f}、r = 0.25 → 0.4 → 0.8）、
   側面結合数 S は r^(-1/2) に比例して減る（{sr(0.25):.1f} → {sr(0.4):.1f} → {sr(0.8):.1f}）。どちらも円柱長に依らない。
   よって P/S は r のおよそ 1.5 乗で増え（r = 0.25 → 0.8 で {pr(0.25) / sr(0.25):.0f} → {pr(0.8) / sr(0.8):.0f}）、l_c には依らない。
   一様分布では f_pole 自体が形で変わるため、P/S は l_c とともに急に下がる
   （l_c = 0.1 → 1.0 → 4.0 µm で {uni_ratio(0.1):.2f} → {uni_ratio(1.0):.3f} → {uni_ratio(4.0):.3f}）。

**頑健でない結論**（置き方で変わる。結論にしてはいけない）

5. **半球の曲率（半径）の効果は、細い側では一貫して極結合が減る向きだが、太い側は置き方で変わる。**
   - 細い側: r = 0.25 µm で、最初の10本のうち極結合の割合は {f("main", 1.0, 0.25):.2f}（文献値の半径 0.4 µm では {f("main", 1.0, 0.4):.2f}）。
     静的な P/S が r のおよそ 1.5 乗で下がる（{pr(0.4) / sr(0.4):.0f} → {pr(0.25) / sr(0.25):.0f}）ことと整合し、結合数を既定形状に固定すると
     {f("decoupled", 1.0, 0.25):.2f} まで戻るので、低下のほとんどは結合数の変化による。箱固定（{f("number_density", 1.0, 0.25):.2f}）、
     発現量∝膜面積（{f("per_area", 1.0, 0.25):.2f}）でも同じ向きに下がる。
   - 太い側: r = 0.8 µm で、全て反映では {f("main", 1.0, 0.8):.2f} と文献値の半径と同じ（r = 0.4–0.8 µm で {min(f("main", 1.0, r) for r in (0.4, 0.5, 0.6, 0.8)):.2f}–{max(f("main", 1.0, r) for r in (0.4, 0.5, 0.6, 0.8)):.2f} で平ら）。
     結合数を固定すると {f("decoupled", 1.0, 0.8):.2f} まで下がるので、結合数の増加（P/S の上昇）と、太いほど極結合が減る
     衝突の幾何が打ち消し合っている。箱固定では {f("number_density", 1.0, 0.8):.2f}、発現量∝膜面積では {f("per_area", 1.0, 0.8):.2f} と、
     置き方で {min(fat_values):.2f}–{max(fat_values):.2f} に散る。
   - したがって、半径を太くして鎖が増える根拠はない。標準誤差は {np.mean([se("main", 1.0, r) for _, r in R_SHAPES]):.2f} 程度。
6. **凝集の速さの形依存性は、箱の取り方で符号が逆転する。** 円柱長 0.1 → 3.0 µm で、
   面積占有率を固定すると単量体の半減時間は {th('main', 0.1):.1f} → {th('main', 3.0):.1f} s に{'伸び' if sign_main > 0 else '縮み'}、箱を 20 µm に固定すると
   {th('number_density', 0.1):.1f} → {th('number_density', 3.0):.1f} s に{'伸びる' if sign_box > 0 else '縮む'}。この2次元の箱は「動画として見やすい」基準で選んだ自由パラメータなので
   （`DynamicsConfig` の docstring）、形が速さに与える影響は言えない。

## 結果の表

最初の10本の結合のうち極結合が占める割合（平均 ± 標準誤差。局在 f_pole=0.95・先端幅 0.10 µm、
発現量 7,300 /細胞）。各条件 20 回（「箱を固定」「発現量∝膜面積」は 12 回）、30 s まで計算、
時間刻みは全て {dt:g} s。「全て反映」は、結合数も衝突の幾何も形に応じて変える通常の計算。

### 円柱長を変える（半径 0.40 µm 固定）

{composition_table(d, L_SHAPES, "L")}

### 半径（曲率）を変える（円柱長 1.0 µm 固定）

{composition_table(d, R_SHAPES, "R")}

### 一様分布（対照）

{uniform_table(d, L_SHAPES, "L")}

{uniform_table(d, R_SHAPES, "R")}

### 凝集の速さ（単量体が半分以下になるまでの s）

箱の取り方で傾向が逆転する（要点6）。

{speed_table(d, L_SHAPES, "L")}

{speed_table(d, R_SHAPES, "R")}

### 回転拡散の見積もりへの感度（A26）

等価球近似は細長い細胞の回転拡散を過大評価する。倍率を下げても、円柱が長いほど極結合が減る傾向は
保たれる。ただし文献値の形（l_c = 1.0 µm）での絶対値は {abs(f('main', 1.0) - f('rotdiff_0.5', 1.0)):.1f} 程度動く。

{rotdiff_table(d)}

### 局在させても純粋な鎖は少ない（形によらず）

30 s 時点で、24 細胞が入っているクラスタの種類（平均、20 回）。「束」は側面結合を1本でも
含むクラスタ（`classify_clusters` の定義）なので、鎖に側面結合が1本混ざると束に数える。

{structure_table(d)}

## 図と動画

| 図 | 内容 |
|---|---|
| `results/figures/shape_scan_static.png` | 静的モデルの結合数 P・S と P/S（円柱長×半径の格子） |
| `results/figures/shape_scan_bonds.png` | 20 s 時点の極結合・側面結合の本数 |
| `results/figures/shape_scan_composition.png` | 結合の構成比と凝集の速さ |
| `results/figures/shape_scan_mechanism.png` | 結合数を固定した対照・箱固定・発現量∝膜面積との比較 |
| `results/figures/shape_scan_rotdiff.png` | 回転拡散の見積もりへの感度 |

動画（4本）は [`shape_movies.md`](shape_movies.md)。集計表は `results/data/shape_scan_summary.csv`、
1回ごとの生データは `shape_scan_dynamic.csv`。

## 実験の設計と、その理由

- **時間刻み**: 全ての形・全ての実験で同じ dt = {dt:g} s（`DynamicsConfig` の既定）。旧既定の 0.02 s は
  側面結合を過大に数え、しかも偏りは形に依存する。1 ステップの先端の移動量（並進＋回転）を半径で
  割った δ（`dynamics.step_displacement_ratio`）で整理した（下の「時間刻みの収束」）。掃引した全ての
  形で δ ≤ 0.55。δ が {dyn.STEP_DISPLACEMENT_LIMIT} を超える形（細くて長い細胞など）には警告が出る。
- **混み具合の対照**: 円柱が長い・太い細胞は箱を占める面積が増え、形の効果と混み具合の効果が交絡する。
  主要な実験は、細胞が箱を占める面積の割合を文献値の形での値（7.8%）に揃えた。「箱を固定」実験は
  箱を 20 µm に固定（細胞数密度が一定）した場合。
- **局在と一様の対照**: 一様分布の f_pole は形で変わる（`model.uniform_f_pole(shape)`）。
  形ごとにその値を使った。
- **機構の切り分け**: 形が効く経路は「アドヘシンの接触幾何（結合数 P・S）」と「細胞どうしの
  衝突幾何」の2つ。`bond_counts_override` で P・S を既定形状の値に固定し、後者だけを残した。
- **発現量**: 既定は細胞あたりの分子数 n_total を固定（表面密度は面積に反比例して変わる）。
  「発現量∝膜面積」は n_total を表面積に比例させる別の置き方（A26）。
- **シード**: 20 回ずつ（感度の確認は 12 回）。標準誤差は 0.02–0.06 程度で、円柱長の効果
  （{f01:.2f} → {f30:.2f}）はこれより十分大きいが、半径の効果は誤差と同程度。

### 時間刻みの収束

{time_step_table()}

{time_step_commentary()}

## 限界

- **等価球近似（A16, A26）**: 拡散係数と沈降速度は同体積の球で見積もる。体積が同じなら形が違っても
  同じ値になり、細長い棒の異方性（回転が長さの3乗に反比例して遅い）を含まない。円柱を長くした側の
  動的な結果は特に怪しい。回転拡散への感度は上の表のとおり（向きは保たれ、絶対値は動く）。
- **2次元・24 細胞・30 s**: 平面拘束で、細胞数・箱・時間は較正できるデータがない自由パラメータ。
  絶対的な時間スケールは使えない。
- **形の範囲**: l_c は 0.1–3.0 µm、r は 0.25–0.8 µm。円柱長と半径を同時に変えた形（格子）は
  掃引していない。**細くて長い細胞（r = 0.25 µm かつ l_c ≥ 2 µm 程度）は、既定の dt でも粗く**
  （δ > {dyn.STEP_DISPLACEMENT_LIMIT}）、`DynamicsConfig` が警告する。半径 0.2 µm 以下も同様。
  l_c = 0 の球は、頭側・尾側の区別がつかなくなるので扱えない（`MIN_CYLINDER_LENGTH_UM`）。
- **半球を円柱と同じ曲率に限った**。端だけを尖らせる・平らにするといった形は別の形状族で、
  `CellShape` では表せない。

## 副次的に見つかったこと（形とは別に、既存モデルについて）

1. **旧既定の dt = 0.02 s は側面結合を過大に数えていた。** 既定を {dt:g} s に変更し、
   既存の動画（`chain.gif` ほか5本）と数値を作り直した。収束の確認は上の表。
   沈降の動画は、収束した dt では seed 0 の実行が 1500 µm の液柱で沈まなかったため、液柱を
   600 µm にした（README の沈降の節）。
2. **局在させても、30 s 時点で純粋な鎖に入っている細胞は 24 個中 {lit['chain_t30_mean']:.0f} 個ほど。** 上の表のとおり、
   文献値の形でも {lit['bundle_t30_mean']:.0f} 個前後が側面結合を含むクラスタ（束）に入る。「束」の定義が「側面結合を
   1本でも含む」ためで、鎖に側面結合が混ざったものも束に数える。README の
   「局在させると頭-尾の揃った有向な鎖が育つ」という記述と、この数字の関係は確認してほしい。
3. **配色を全ての図・動画で統一した**（極 = 青、側面 = オレンジ）。旧配色（極接触 = ティール、
   側面接触 = マゼンタ）は2型色覚（deutan）で区別しにくかった（OKLab ΔE×100 が 2.6。基準は 8 以上）。
   `scripts/figstyle.py` に集約し、`tests/test_figstyle.py` で色覚シミュレーション下の分離を固定した。

## 再現手順

```bash
python3 scripts/run_shape_scan.py static                       # 静的モデルの格子
python3 scripts/run_shape_scan.py dynamic --seeds 20 --jobs 8  # 時間発展（約 3 時間。再開可能）
python3 scripts/check_time_step.py --seeds 40 --jobs 8         # 時間刻みの収束
python3 scripts/plot_shape_scan.py                             # 図と集計表
python3 scripts/simulate_shape_comparison.py --scan cylinder   # 動画（--scan radius で曲率）
python3 scripts/export_shape_docs.py                           # このページと shape_movies.md
```

```python
from adhesion import dynamics as dyn
from adhesion.shape import CellShape

shape = CellShape.from_cylinder(cylinder_length_um=2.0, radius_um=0.3)   # 円柱長と半径で指定
dyn.run(dyn.DynamicsConfig(shape=shape, f_pole=0.95, tip_width_um=0.10))
```
"""


# ---------------------------------------------------------------------------
# 動画のページ
# ---------------------------------------------------------------------------
def _cp(p: dict) -> str:
    return f"極 {p['final_polar_contacts']} / 側面 {p['final_lateral_contacts']}"


#: 各動画で見えることの説明。数値は JSON から埋め込む。形容（数珠状、網目状など）は動画の最終フレームを
#: 見て書いたもので、動画を作り直したら見直すこと。
OBSERVATIONS = {
    "shape_cylinder_localized": lambda a: (
        f"見えること: 円柱が短い（{a[0]['cylinder_um']:g} µm、ほぼ球）細胞は、青い数珠状の鎖が長くつながり"
        f"（{_cp(a[0])}、最大のかたまりは {a[0]['final_largest_cluster']} 細胞）、オレンジはわずか。"
        f"文献値の形では、青い鎖の区間と、オレンジ（側面結合をもつ細胞）の区間が混ざったかたまりになる（{_cp(a[1])}）。"
        f"円柱が長い（{a[2]['cylinder_um']:g} µm）と、オレンジの細胞が増え、細胞が横に並んだ対や束が目立つ"
        f"（{_cp(a[2])}）。局在させていても、円柱が長いと側面接触が極接触と同じくらいできる。"),
    "shape_cylinder_uniform": lambda a: (
        f"見えること: どの形でも、ほとんどの細胞がオレンジ（側面結合をもつ）になり、細胞が横に並んだ束ができる"
        f"（{_cp(a[0])}、{_cp(a[1])}、{_cp(a[2])}）。青い細胞（極結合だけをもつ）は、あっても1〜2個。"
        "局在させた動画 1 と見比べると、青い鎖の区間があるかどうかが、分布の違いの見え方になる。"
        f"円柱が長いほど、側面の1接触あたりの分子対が増える（{a[0]['side_bonds_per_contact']:.0f} → "
        f"{a[1]['side_bonds_per_contact']:.0f} → {a[2]['side_bonds_per_contact']:.0f} 対）。"),
    "shape_radius_localized": lambda a: (
        f"見えること: 半径が小さい（細い、{a[0]['radius_um']:g} µm）細胞は、{a[0]['final_largest_cluster']} 細胞が"
        f"つながった1つの網目状のかたまりになり、側面接触が極接触より多い（{_cp(a[0])}）。"
        f"文献値の半径では、青い鎖とオレンジの束が混ざる（{_cp(a[1])}）。"
        f"半径が大きい（太い、{a[2]['radius_um']:g} µm）細胞は、{a[2]['final_clusters']} つのかたまり"
        f"（最大 {a[2]['final_largest_cluster']} 細胞）と単独の細胞 {a[2]['final_monomers']} 個に分かれたままで、"
        f"結合が増えるのがゆっくりに見える（{_cp(a[2])}）。ただし、この速さの違いは箱の取り方"
        "（細胞が箱を占める面積の割合を一定にした）に依存する見え方で、曲率そのものの効果とは言えない。"),
    "shape_radius_uniform": lambda a: (
        f"見えること: どの半径でも、オレンジの束になる（{_cp(a[0])}、{_cp(a[1])}、{_cp(a[2])}）。"
        f"細い（{a[0]['radius_um']:g} µm）細胞は {a[0]['final_largest_cluster']} 細胞が1つの大きな束に、"
        f"太い（{a[2]['radius_um']:g} µm）細胞は {a[2]['final_clusters']} つの小さなかたまりと単独の細胞"
        f" {a[2]['final_monomers']} 個に分かれて見える。"),
}

MOVIE_ORDER = [
    ("shape_cylinder_localized", "円柱長を変える（局在）"),
    ("shape_cylinder_uniform", "円柱長を変える（一様分布）"),
    ("shape_radius_localized", "半球の曲率（半径）を変える（局在）"),
    ("shape_radius_uniform", "半球の曲率（半径）を変える（一様分布）"),
]


def load_movie(name: str) -> dict:
    path = ANIMATIONS / f"{name}.json"
    if not path.exists():
        raise SystemExit(f"{path.name} が無い（先に scripts/simulate_shape_comparison.py で動画を作ること）")
    return json.loads(path.read_text(encoding="utf-8"))


def shape_label(panel: dict, scan: str) -> str:
    literature = (abs(panel["cylinder_um"] - 1.0) < 1e-9 and abs(panel["radius_um"] - 0.4) < 1e-9)
    base = f"円柱 {panel['cylinder_um']:g} µm" if scan == "cylinder" else f"半径 {panel['radius_um']:g} µm"
    return base + ("（文献値）" if literature else "")


def movie_section(number: int, name: str, title: str, facts: dict) -> str:
    scan = "cylinder" if "cylinder" in name else "radius"
    panels = facts["panels"]
    head = "| 形 | 全長 | " + ("" if scan == "cylinder" else "曲率 | ") + "極 / 側面 | 単独の細胞 | 1接触あたり 極 / 側面 | シード |"
    sep = "|---|---:|" + ("" if scan == "cylinder" else "---:|") + "---:|---:|---:|---:|"
    rows = []
    for p in panels:
        curvature = "" if scan == "cylinder" else f"{p['curvature_per_um']:.1f} /µm | "
        rows.append(
            f"| {shape_label(p, scan)} | {p['length_um']:.2f} µm | {curvature}"
            f"{p['final_polar_contacts']} / {p['final_lateral_contacts']} | {p['final_monomers']} | "
            f"{p['polar_bonds_per_contact']:.0f} / {p['side_bonds_per_contact']:.1f} 対 | {p['seed']} |")
    return f"""## {number}. {title}

![{title}](../results/animations/{name}.gif)

{facts['total_time_s']:g} 秒時点の接触数（極 / 側面）と、動画に使った実行のシード:

{head}
{sep}
""" + "\n".join(rows) + f"""

{OBSERVATIONS[name](panels)}
"""


def movies_doc() -> str:
    all_facts = {name: load_movie(name) for name, _ in MOVIE_ORDER}
    first = all_facts[MOVIE_ORDER[0][0]]
    n_frames, fps = first["n_frames"], first["fps"]
    total = first["total_time_s"]
    frame_interval = total / (n_frames - 1)
    from PIL import Image

    with Image.open(ANIMATIONS / f"{MOVIE_ORDER[0][0]}.gif") as gif:
        frame_ms = gif.info.get("duration", 1000 / fps)
        check(gif.n_frames == n_frames, "GIF のコマ数が JSON と違う")
    dt = first["panels"][0]["dt_s"]
    check(all(abs(p["dt_s"] - dt) < 1e-12 for f in all_facts.values() for p in f["panels"]),
          "動画によって dt が違う（全ての形で同じ dt のはず）")

    # 形の表（動画に出てくる形をまとめる）
    shapes = {}
    for facts in all_facts.values():
        for p in facts["panels"]:
            shapes[(p["cylinder_um"], p["radius_um"])] = p
    shape_rows = []
    for (c, r), p in sorted(shapes.items(), key=lambda kv: (kv[0][1] != 0.4, kv[0][0] if kv[0][1] == 0.4 else kv[0][1])):
        literature = (c, r) == LITERATURE
        shape_rows.append(f"| {'文献値' if literature else ('円柱 ' + format(c, 'g') + ' µm' if r == 0.4 else '半径 ' + format(r, 'g') + ' µm')} | "
                          f"{p['length_um']:.2f} µm | {c:g} µm | {r:g} µm | {p['curvature_per_um']:.1f} /µm | {p['box_um']:g} µm |")

    sections = "\n".join(movie_section(i + 1, name, title, all_facts[name]) for i, (name, title) in enumerate(MOVIE_ORDER))
    return f"""# 細胞の形を変えたときの時間発展動画

大腸菌の形（円柱の長さ、半球の曲率）を変えたとき、鎖や束ができる過程がどう変わるかを、
形の違う細胞集団を横に並べて同時に時間発展させた動画 4 本の見方と作り方。

このページは**動画だけ**を説明する。形ごとの統計と結論の確かさは
[`shape_scan.md`](shape_scan.md) にある。

> **動画は1回の実行で、統計ではない。** 同じ条件でもシードで結果はかなり変わる。
> 動画から言えるのは「こういう形のかたまりができうる」という絵であって、
> 「形を変えるとこうなる」という頻度や傾向ではない（[動画から言えること・言えないこと](#動画から言えること言えないこと)）。
>
> このページの表は `scripts/export_shape_docs.py` が、動画と同時に出力される JSON から生成する。

## 動画の一覧

| 動画 | 比べる形 | アドヘシンの分布 |
|---|---|---|
| [1. 円柱長・局在](#1-円柱長を変える局在) | 円柱長 l_c = 0.1 / 1.0（文献値）/ 3.0 µm、半径 0.40 µm 固定 | 局在（f_pole = 0.95、先端幅 0.10 µm） |
| [2. 円柱長・一様](#2-円柱長を変える一様分布) | 同上 | 一様（実測の分布、対照） |
| [3. 曲率・局在](#3-半球の曲率半径を変える局在) | 半径 r = 0.25 / 0.40（文献値）/ 0.60 µm、円柱長 1.0 µm 固定 | 局在 |
| [4. 曲率・一様](#4-半球の曲率半径を変える一様分布) | 同上 | 一様 |

どの動画も、**{first['n_cells']} 細胞を 0〜{total:g} 秒**時間発展させ、3 つの形を並べている。局在と一様は
同じ形・同じ条件で、アドヘシンの分布だけが違う。

## 画面の見方

![画面の構成](../results/animations/shape_cylinder_localized.gif)

**上段: 細胞が動く箱（2 次元）**

- 細胞はカプセル形（円柱の両端に半球）。**拡大率は 3 つのパネルで同じ**なので、細胞の大きさの違いが
  そのまま見える。箱の大きさはパネルごとに違う（下の「比べている形」）。
- **細胞の色は、その細胞がもつ結合の種類**を表す。**灰色**は結合なし（単独）、**青**は極結合だけを
  もつ細胞（鎖・二量体）、**オレンジ**は側面結合を1本でももつ細胞（束）。かたまり全体ではなく
  細胞ごとに塗るので、鎖の途中に側面接触が1本あっても、塗り替わるのはその2細胞だけ。
  この色（極 = 青、側面 = オレンジ）は、全ての図・動画で共通。（統計で使う「束」は、側面結合を
  1本でも含む*かたまり*という別の定義。）
- **● が Ag 極（頭）、○ が Nb 極（尾）**。極どうしで結合できるのは Ag 極と Nb 極だけなので、
  極接触でできる鎖は必ず頭-尾に揃う（● と ○ が向き合って接する）。
- 結合線は細胞の下に描かれ、架橋長が 8 nm と細胞に比べて極小のため、ほぼ見えない。
  **結合は「色のついたかたまり」と下段の曲線で読む。**
- パネルの見出しの2行目は、その形で**1回の接触あたりに結合する分子対の数**（極接触 P 対 /
  側面接触 S 対）。静的モデルが形から計算した値で、時間発展ではこの数で結合のしやすさが決まる。

**下段: 接触の箇所数の累積**

- 青が**極接触**（鎖の背骨）、オレンジが**側面接触**（束）。細胞の対がくっついた箇所の数で、
  分子対の数ではない。縦軸は同じ動画の3つのパネルで共通。
- グラフの上の文字は、その時刻の「極接触 n　側面接触 m　単独の細胞 k」。

**再生**: {frame_interval:g} 秒ごとの {n_frames} コマを 1 コマ {frame_ms:g} ms で再生するので、実時間の約 {frame_interval / (frame_ms / 1000):.0f} 倍速（1 周約 {n_frames * frame_ms / 1000:.0f} 秒）。

{sections}
## 比べている形

| 形 | 全長 l | 円柱長 l_c | 半径 r | 曲率 1/r | 箱の一辺 |
|---|---:|---:|---:|---:|---:|
""" + "\n".join(shape_rows) + f"""

- 半球の半径は円柱の半径と等しいので、**曲率 1/r は円柱の太さと切り離せない**。曲率だけを変えたい
  動画 3・4 では、全長ではなく円柱長 1.0 µm を固定して半径だけを動かしている（全長は 2r の分だけ変わる）。
- **箱の大きさ**は、細胞が箱（2 次元）を占める面積の割合を、文献値の形での値（7.8%）にそろえて決めた。
  細胞が大きいほど箱も大きい。そろえないと、円柱が長い・太い細胞ほど混み合って、形の効果と混み具合の
  効果が区別できなくなる。
- **時間刻み**は全ての動画・全ての形で同じ {dt:g} s（収束した値。`shape_scan.md` の「時間刻みの収束」）。

## 動画の条件

- {first['n_cells']} 細胞を、箱の中にランダムな位置・向きで置いて始める。細胞は 2 次元（平面）でブラウン運動する
  （並進・回転拡散は同体積の球で見積もる）。静置で、震盪・沈降・遊泳はなし。
- 結合は不可逆（k_off = 0）。局在は極に 95%、先端幅 0.10 µm。一様は、その形の表面積に応じた割合
  （形で変わる）。発現量は 1 細胞あたり 7,300 分子で、形を変えても変えない。
- 動画の 1 コマは {frame_interval:g} 秒。シミュレーションは {total:g} 秒まで。

## 代表的な実行の選び方

動画は 1 回の実行なので、見栄えのよい実行を選ぶと印象が偏る。そこで、形ごとに
[`run_shape_scan.py`](../scripts/run_shape_scan.py) が回した 20 回の実行から、**20 秒時点の極結合の割合と
単量体が半分になるまでの時間が、その形の平均に最も近い実行**を自動で選んだ（標準偏差で割った
距離の和が最小のもの。`--seed representative`）。上の表の「シード」がそれ。

同じ設定なので、動画の最初の 30 秒は、スイープの同じシードの実行と完全に一致する（4本の
動画に出てくる 10 通りの（分布・形・シード）全てで、10・20・30 秒の極接触・側面接触の数が
一致することを確認した）。30 秒以降は動画のために延長した部分。

## 動画から言えること・言えないこと

**言えること**

- 形ごとに、どんな種類のかたまり（数珠状の鎖、側面で並んだ束、網目状の大きなかたまり）が
  できうるかという絵。
- 局在と一様で、できるかたまりの種類が違うこと。

**言えないこと**

- 「円柱が長いと側面接触が増える」のような傾向。動画はパネルごとに別の 1 回の実行で、
  数本の接触の差はシードでいくらでも入れ替わる。傾向は 20 回ずつの統計で見ること
  （[`shape_scan.md`](shape_scan.md)）。
- 凝集の速さの比較。箱の取り方で逆転することが分かっている。
- 実際の凝集アッセイの時間スケール。細胞数・箱・時間は、較正できるデータがない自由パラメータ。

## 作り方

```bash
python3 scripts/run_shape_scan.py dynamic --seeds 20 --jobs 8   # 20 回ずつの実行（約 3 時間、再開可能）
python3 scripts/simulate_shape_comparison.py --scan cylinder    # 動画 1・2（--scenario uniform で一様）
python3 scripts/simulate_shape_comparison.py --scan radius      # 動画 3・4
python3 scripts/export_shape_docs.py                            # このページを再生成
```

任意の形を並べることもできる（円柱長:半径、µm）。

```bash
python3 scripts/simulate_shape_comparison.py --shapes 2.0:0.3 1.0:0.4 --scenario localized \\
    --total-time-s 60 --fps 15 --output my_shapes.gif
```

主なオプション: `--scenario localized|uniform`、`--box area|number`（箱の決め方）、
`--seed`（整数、または既定の `representative`）、`--total-time-s`、`--n-cells`。
実装は [`scripts/simulate_shape_comparison.py`](../scripts/simulate_shape_comparison.py)。
描画は [`simulate_chain_formation.py`](../scripts/simulate_chain_formation.py) と共通で、配色は
[`scripts/figstyle.py`](../scripts/figstyle.py)。
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dynamic-csv", type=Path, default=ps.DYNAMIC_CSV)
    parser.add_argument("--skip-movies", action="store_true", help="動画のページを作らない")
    args = parser.parse_args()

    data = Data(args.dynamic_csv)
    (DOCS / "shape_scan.md").write_text(scan_doc(data), encoding="utf-8")
    print("生成: docs/shape_scan.md")
    if not args.skip_movies:
        (DOCS / "shape_movies.md").write_text(movies_doc(), encoding="utf-8")
        print("生成: docs/shape_movies.md")


if __name__ == "__main__":
    main()
