"""細胞の形（円柱の長さ・半球の曲率）を振り、鎖・束の形成がどう変わるかを調べる。

    python3 scripts/run_shape_scan.py static
    python3 scripts/run_shape_scan.py dynamic --seeds 20 --jobs 10
    python3 scripts/plot_shape_scan.py

半球の半径は円柱の半径と等しいので、半球の曲率 1/r は円柱の太さと切り離せない。
そこで形を2つの独立な量で振る。

    円柱長 l_c   円柱部分の長さ。全長は l = l_c + 2r
    半径 r       円柱の半径であり、半球の曲率半径（曲率 1/r）

`static` は静的モデル（`adhesion.model`）が返す結合数 P（極）・S（側面）を解析的に
掃引する。決定的で速い。`dynamic` は時間発展（`adhesion.dynamics`）を多数のシードで
回し、結合の種類（極 = 鎖の背骨、側面 = 束）の構成比と凝集の速さを測る。

動的な掃引の設計上の選択
------------------------
- **時間刻み** 全ての形・全ての実験で同じ dt（`DynamicsConfig` の既定、0.0025 s）を使う。
  旧既定の 0.02 s は側面結合を過大に数え、しかも偏りは形に依存する（細い細胞ほど大きい）。
  0.0025 s は、掃引する全ての形で、結合の構成比が収束している（1 ステップの先端の移動量を
  半径で割った δ が最大 0.55。`dynamics.step_displacement_ratio`、`scripts/check_time_step.py`、
  `docs/shape_scan.md`）。細い形では単量体の数（凝集の速さ）に dt の残差がありうる。
- **混み具合の対照** 円柱が長い・太い細胞は箱を占める面積が増え、形の効果と混み具合の
  効果が交絡する。`--box area`（既定）は細胞が箱を占める面積の割合を既定形状の値に揃え、
  `--box number` は箱を 20 µm に固定する（細胞数密度を一定に保つ）。
- **局在と一様の対照** 一様分布の f_pole は形で変わる（`model.uniform_f_pole`）。
- **機構の切り分け** `decouple` は、形は変えるが結合数 P・S は既定形状の値に固定する。
  形が効く経路は「アドヘシンの接触幾何（結合数）」と「細胞どうしの衝突幾何」の2つなので、
  前者を止めると後者だけの効果が見える。
- **発現量の置き方** `per_area`（`per_area` 実験）は発現量を膜面積に比例させる。既定は
  n_total を固定（A26）。
- **回転拡散の感度** `rotdiff_0.5` / `rotdiff_0.25` は、等価球近似の回転拡散係数に倍率を
  掛ける（A26）。細長い棒の回転は等価球より遅いので、円柱長への応答が拡散係数の見積もりに
  依存するかを見る。
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import sys
import time
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import dynamics as dyn
from adhesion import literature as lit
from adhesion import model
from adhesion.shape import CellShape

DATA_DIR = Path(__file__).resolve().parent.parent / "results" / "data"
STATIC_CSV = DATA_DIR / "shape_scan_static.csv"
DYNAMIC_CSV = DATA_DIR / "shape_scan_dynamic.csv"

N_MID = float(lit.ADHESINS_AT_MID_INDUCTION)

#: 全ての形・実験で共通の時間刻み（収束した値。`DynamicsConfig` の既定）。
DEFAULT_DT_S = dyn.DynamicsConfig().dt_s

#: 局在させた場合（鎖を狙う設計）の条件。`scripts/simulate_chain_formation.py` の chain と同じ。
LOCALIZED = dict(f_pole=0.95, tip_width_um=0.10)

#: 形の掃引範囲。文献値は l = 1.80 ± 0.25 µm、r = 0.40 ± 0.10 µm（円柱長は 1.0 µm）。
#: 半径 0.25 未満は拡散が速く、dt の収束が悪化するので含めない。
CYLINDER_SCAN_UM = [0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0]
RADIUS_SCAN_UM = [0.25, 0.3, 0.4, 0.5, 0.6, 0.8]
DEFAULT_CYLINDER_UM = model.DEFAULT_SHAPE.cylinder_length_um
DEFAULT_RADIUS_UM = model.DEFAULT_SHAPE.radius_um

#: 結果を記録する時刻（s）。
CHECKPOINTS_S = (10.0, 20.0, 30.0)
FRAME_INTERVAL_S = 0.5

#: 回転拡散の感度確認（A26）。
ROT_DIFFUSION_SCALES = (0.5, 0.25)
ROT_DIFFUSION_CYLINDERS_UM = (0.1, 1.0, 3.0)

#: 実行順の周回ごとのシード数。
SEED_ROUND = 4


# ---------------------------------------------------------------------------
# 静的モデルの掃引
# ---------------------------------------------------------------------------
STATIC_CYLINDERS_UM = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
STATIC_RADII_UM = [0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]


def run_static() -> None:
    """結合数 P（極）・S（側面）を円柱長×半径の格子で求める。"""

    default_area = model.DEFAULT_SHAPE.surface_area_um2
    rows = []
    for cylinder, radius in itertools.product(STATIC_CYLINDERS_UM, STATIC_RADII_UM):
        shape = CellShape.from_cylinder(cylinder, radius)
        f_uniform = model.uniform_f_pole(shape)
        localized = model.predict(N_MID, shape=shape, **LOCALIZED)
        uniform = model.predict(N_MID, f_uniform, 0.5, None, shape=shape)
        area_scale = shape.surface_area_um2 / default_area
        localized_pa = model.predict(N_MID * area_scale, shape=shape, **LOCALIZED)
        uniform_pa = model.predict(N_MID * area_scale, f_uniform, 0.5, None, shape=shape)
        rows.append(
            {
                "cylinder_um": cylinder,
                "radius_um": radius,
                "length_um": round(shape.length_um, 4),
                "curvature_per_um": round(shape.curvature_per_um, 4),
                "aspect_ratio": round(shape.aspect_ratio, 4),
                "surface_area_um2": round(shape.surface_area_um2, 4),
                "volume_um3": round(shape.volume_um3, 4),
                "uniform_f_pole": round(f_uniform, 4),
                "localized_polar": round(localized.bonds_heavy_pole, 3),
                "localized_side": round(localized.bonds_side_by_side, 3),
                "uniform_polar": round(uniform.bonds_heavy_pole, 3),
                "uniform_side": round(uniform.bonds_side_by_side, 3),
                "localized_polar_per_area": round(localized_pa.bonds_heavy_pole, 3),
                "localized_side_per_area": round(localized_pa.bonds_side_by_side, 3),
                "uniform_polar_per_area": round(uniform_pa.bonds_heavy_pole, 3),
                "uniform_side_per_area": round(uniform_pa.bonds_side_by_side, 3),
            }
        )
    STATIC_CSV.parent.mkdir(parents=True, exist_ok=True)
    with STATIC_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} 形状を {STATIC_CSV.name} に書き出した")


# ---------------------------------------------------------------------------
# 動的な掃引
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Condition:
    """1つの実験条件（シード以外）。"""

    experiment: str
    scenario: str  # "localized" | "uniform"
    cylinder_um: float
    radius_um: float
    box: str = "area"  # "area" | "number"
    scaling: str = "per_cell"  # "per_cell" | "per_area"
    decouple: bool = False
    #: 回転拡散係数に掛ける倍率（A26 の感度確認）。実験名 `rotdiff_<倍率>` に対応する。
    rot_scale: float = 1.0

    @property
    def shape(self) -> CellShape:
        return CellShape.from_cylinder(self.cylinder_um, self.radius_um)

    def key(self, seed: int) -> tuple:
        return (
            self.experiment, self.scenario, self.cylinder_um, self.radius_um,
            self.box, self.scaling, self.decouple, seed,
        )


def _scenario_kwargs(scenario: str, shape: CellShape) -> dict:
    if scenario == "localized":
        return dict(LOCALIZED)
    if scenario == "uniform":
        return dict(f_pole=model.uniform_f_pole(shape), tip_width_um=None)
    raise ValueError(scenario)


def _box_size_um(shape: CellShape, n_cells: int, mode: str) -> float:
    """箱の一辺（µm）。"""

    base = dyn.DynamicsConfig().box_size_um
    if mode == "number":
        return base
    # 面積占有率を既定形状の値に揃える。
    reference_fraction = n_cells * model.DEFAULT_SHAPE.planar_footprint_um2 / base**2
    return round(math.sqrt(n_cells * shape.planar_footprint_um2 / reference_fraction), 2)


def _reference_counts(scenario: str) -> tuple[float, float]:
    """既定形状での (P, S)。decouple で固定する値。"""

    prediction = model.predict(N_MID, shape=model.DEFAULT_SHAPE, **_scenario_kwargs(
        scenario, model.DEFAULT_SHAPE
    ))
    return prediction.bonds_heavy_pole, prediction.bonds_side_by_side


def build_config(
    condition: Condition, seed: int, total_time_s: float, dt_s: float | None, n_cells: int
) -> dyn.DynamicsConfig:
    shape = condition.shape
    if dt_s is None:
        dt_s = DEFAULT_DT_S
    n_total = N_MID
    if condition.scaling == "per_area":
        n_total = N_MID * shape.surface_area_um2 / model.DEFAULT_SHAPE.surface_area_um2
    override = _reference_counts(condition.scenario) if condition.decouple else None
    return dyn.DynamicsConfig(
        shape=shape,
        n_cells=n_cells,
        box_size_um=_box_size_um(shape, n_cells, condition.box),
        n_total=n_total,
        dt_s=dt_s,
        total_time_s=total_time_s,
        frame_stride=max(1, round(FRAME_INTERVAL_S / dt_s)),
        seed=seed,
        bond_counts_override=override,
        diffusion_scale=(1.0, condition.rot_scale),
        **_scenario_kwargs(condition.scenario, shape),
    )


def _frame_stats(frame: dyn.Frame, n_cells: int, distinct: bool) -> dict:
    clusters = dyn.classify_clusters(n_cells, frame.bonds, distinct)
    sizes = [len(c.cells) for c in clusters]
    return {
        "mono": sum(1 for c in clusters if c.kind == "monomer"),
        "polar": sum(1 for b in frame.bonds if b.kind == "polar"),
        "lateral": sum(1 for b in frame.bonds if b.kind == "lateral"),
        "largest": max(sizes),
        # 極結合だけでできた2細胞（二量体）、3細胞以上（鎖）、側面結合を含む（束）
        "dimer": sum(len(c.cells) for c in clusters if c.kind == "chain" and len(c.cells) == 2),
        "chain": sum(len(c.cells) for c in clusters if c.kind == "chain" and len(c.cells) >= 3),
        "bundle": sum(len(c.cells) for c in clusters if c.kind == "bundle"),
    }


def run_one(task: tuple) -> dict:
    condition, seed, total_time_s, dt_s, n_cells = task
    config = build_config(condition, seed, total_time_s, dt_s, n_cells)
    started = time.time()
    result = dyn.run(config)

    stats = [(f.time_s, _frame_stats(f, n_cells, config.distinct_adhesins_per_pole))
             for f in result.frames]
    row = {
        "experiment": condition.experiment,
        "scenario": condition.scenario,
        "box_mode": condition.box,
        "scaling": condition.scaling,
        "decouple": int(condition.decouple),
        "cylinder_um": condition.cylinder_um,
        "radius_um": condition.radius_um,
        "length_um": round(condition.shape.length_um, 4),
        "seed": seed,
        "box_um": config.box_size_um,
        "area_fraction": round(n_cells * condition.shape.planar_footprint_um2 / config.box_size_um**2, 4),
        "n_total": round(config.n_total, 1),
        "f_pole": round(config.f_pole, 4),
        "polar_bonds_P": round(result.polar_bond_count, 3),
        "side_bonds_S": round(result.side_bond_count, 3),
        "dt_s": round(config.dt_s, 6),
        "total_time_s": total_time_s,
    }
    # 単量体が半分以下になる最初の時刻（打ち切りなら総時間）
    half = next((t for t, s in stats if s["mono"] <= n_cells / 2), total_time_s)
    row["t_half_monomers_s"] = round(half, 3)
    for checkpoint in CHECKPOINTS_S:
        _, s = min(stats, key=lambda item: abs(item[0] - checkpoint))
        for name, value in s.items():
            row[f"{name}_t{checkpoint:g}"] = value

    # 形成順の情報。bonds は成立した順に並ぶ（既定では解離しない）。
    final_bonds = result.frames[-1].bonds
    row["first_bond_kind"] = final_bonds[0].kind if final_bonds else ""
    first = final_bonds[:10]
    row["polar_frac_first10"] = (
        round(sum(b.kind == "polar" for b in first) / len(first), 4) if len(first) == 10 else ""
    )
    row["runtime_s"] = round(time.time() - started, 1)
    return row


def _line(experiment: str, scenario: str, cylinder: float, radius: float, **kwargs) -> Condition:
    return Condition(experiment, scenario, cylinder, radius, **kwargs)


def build_conditions(include_secondary: bool = True) -> list[tuple[Condition, str]]:
    """実験条件の一覧。(条件, 優先度) を返す。優先度 'primary' は多めのシードで回す。"""

    shapes = [(c, DEFAULT_RADIUS_UM) for c in CYLINDER_SCAN_UM]
    shapes += [(DEFAULT_CYLINDER_UM, r) for r in RADIUS_SCAN_UM if r != DEFAULT_RADIUS_UM]
    conditions: list[tuple[Condition, str]] = []
    for cylinder, radius in shapes:
        for scenario in ("localized", "uniform"):
            conditions.append((_line("main", scenario, cylinder, radius), "primary"))
            conditions.append(
                (_line("decoupled", scenario, cylinder, radius, decouple=True), "primary")
            )
    # A26 の感度確認: 細長い棒の回転拡散は等価球近似より遅い。円柱長への応答が
    # 回転拡散の見積もりに依存するかを、倍率 1（= main）・0.5・0.25 で比べる。
    for scale in ROT_DIFFUSION_SCALES:
        for cylinder in ROT_DIFFUSION_CYLINDERS_UM:
            conditions.append(
                (_line(f"rotdiff_{scale:g}", "localized", cylinder, DEFAULT_RADIUS_UM,
                       rot_scale=scale), "primary")
            )
    if include_secondary:
        for cylinder, radius in shapes:
            conditions.append(
                (_line("number_density", "localized", cylinder, radius, box="number"), "secondary")
            )
            conditions.append(
                (_line("per_area", "localized", cylinder, radius, scaling="per_area"), "secondary")
            )
    return conditions


def _existing_keys(path: Path) -> set[tuple]:
    if not path.exists():
        return set()
    keys = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            keys.add(
                (
                    row["experiment"], row["scenario"], float(row["cylinder_um"]),
                    float(row["radius_um"]), row["box_mode"], row["scaling"],
                    bool(int(row["decouple"])), int(row["seed"]),
                )
            )
    return keys


def run_dynamic(args: argparse.Namespace) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    done = _existing_keys(DYNAMIC_CSV)
    tasks = []
    for condition, priority in build_conditions(not args.primary_only):
        n_seeds = args.seeds if priority == "primary" else args.secondary_seeds
        for seed in range(n_seeds):
            if condition.key(seed) in done:
                continue
            tasks.append((condition, seed, args.total_time_s, args.dt_s, args.n_cells))
    # 重い（円柱が長く箱が大きい）ものから流して、最後に空きコアが出ないようにする。
    def cost(task) -> float:
        shape = task[0].shape
        dt = task[3] or DEFAULT_DT_S
        return (task[2] / dt) * (1.0 + 0.2 * shape.length_um)

    # シードを SEED_ROUND 個ずつの周回に分け、各周回の中で重いものから流す。途中で止めても
    # 全条件が同じシード数だけ揃っているので、そこまでの結果で比較できる。
    tasks.sort(key=lambda t: (t[1] // SEED_ROUND, -cost(t)))
    print(f"{len(tasks)} 回を実行する（既に完了: {len(done)} 回）。ジョブ数 {args.jobs}", flush=True)
    if not tasks:
        return

    started = time.time()
    new_file = not DYNAMIC_CSV.exists()
    with DYNAMIC_CSV.open("a", newline="", encoding="utf-8") as handle, Pool(args.jobs) as pool:
        writer = None
        for index, row in enumerate(pool.imap_unordered(run_one, tasks, chunksize=1), 1):
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                if new_file:
                    writer.writeheader()
            writer.writerow(row)
            handle.flush()
            if index % 20 == 0 or index == len(tasks):
                elapsed = time.time() - started
                eta = elapsed / index * (len(tasks) - index)
                print(f"  {index}/{len(tasks)}  経過 {elapsed/60:.1f} 分  残り約 {eta/60:.1f} 分", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("static", help="静的モデルの結合数を円柱長×半径の格子で求める")
    dynamic = sub.add_parser("dynamic", help="時間発展を多数のシードで回す")
    dynamic.add_argument("--seeds", type=int, default=20, help="主要な条件のシード数")
    dynamic.add_argument("--secondary-seeds", type=int, default=12, help="感度確認の条件のシード数")
    dynamic.add_argument("--primary-only", action="store_true", help="感度確認の条件を省く")
    dynamic.add_argument("--jobs", type=int, default=8)
    dynamic.add_argument("--total-time-s", type=float, default=30.0)
    dynamic.add_argument("--dt-s", type=float, default=None,
                         help="時間刻み（s）。省略すると DynamicsConfig の既定（収束した値）")
    dynamic.add_argument("--n-cells", type=int, default=24)
    args = parser.parse_args()

    if args.command == "static":
        run_static()
    else:
        run_dynamic(args)


if __name__ == "__main__":
    main()
