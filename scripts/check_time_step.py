"""時間刻み dt に対する結果の収束を確かめる。

    python3 scripts/check_time_step.py --seeds 40 --jobs 9

形ごとに、いくつかの dt で局在条件のシミュレーションを多数のシードで回し、20 秒時点の
単量体数・極結合数・側面結合数・極結合の割合を比べる。dt を小さくしても値が変わらなく
なれば収束している。結果は `results/data/time_step_convergence.csv` に書き、
`scripts/export_shape_docs.py` が `docs/shape_scan.md` の表に使う。

目安として、1 ステップの先端の移動量（並進＋回転）を半径で割った δ
（`adhesion.dynamics.step_displacement_ratio`）も併せて出す。実測では δ ≲ 0.55 で収束し、
δ ≳ 0.7 で側面結合を約 20% 過大に数える。
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import dynamics as dyn
from adhesion.shape import CellShape

OUTPUT = Path(__file__).resolve().parent.parent / "results" / "data" / "time_step_convergence.csv"
CHECK_TIME_S = 20.0
N_CELLS = 24

#: 形（円柱長:半径、µm）ごとに比べる dt。文献値の形は旧既定の 0.02 s まで、細い形・長い形は
#: 既定（0.0025 s）の前後を調べる。
DEFAULT_CASES = {
    (1.0, 0.4): [0.02, 0.01, 0.005, 0.0025],
    (1.0, 0.25): [0.005, 0.0025, 0.00125],
    (3.0, 0.4): [0.005, 0.0025, 0.00125],
}


def run_one(task: tuple) -> dict:
    cylinder, radius, dt, seed = task
    shape = CellShape.from_cylinder(cylinder, radius)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", dyn.CoarseTimeStepWarning)  # 粗い dt を調べるのが目的
        config = dyn.DynamicsConfig(
            shape=shape, n_cells=N_CELLS, seed=seed + 5000, total_time_s=CHECK_TIME_S,
            f_pole=0.95, tip_width_um=0.10, dt_s=dt, frame_stride=round(CHECK_TIME_S / dt),
        )
    final = dyn.run(config).frames[-1]
    clusters = dyn.connected_components(N_CELLS, final.bonds)
    return {
        "cylinder_um": cylinder, "radius_um": radius, "dt_s": dt, "seed": seed,
        "delta": round(dyn.step_displacement_ratio(shape, dt), 3),
        "monomers": sum(len(c) == 1 for c in clusters),
        "polar": sum(b.kind == "polar" for b in final.bonds),
        "lateral": sum(b.kind == "lateral" for b in final.bonds),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=int, default=40)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    tasks = [
        (c, r, dt, seed)
        for (c, r), dts in DEFAULT_CASES.items()
        for dt in sorted(dts)  # 細かい dt（重い）ほど後ろ。chunksize=1 で空きコアに回る
        for seed in range(args.seeds)
    ]
    tasks.sort(key=lambda t: t[2])
    print(f"{len(tasks)} 回を実行する", flush=True)
    started = time.time()
    with Pool(args.jobs) as pool:
        rows = pool.map(run_one, tasks, chunksize=1)
    print(f"経過 {(time.time() - started) / 60:.1f} 分")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"{args.output.name}\n形(l_c,r)   dt        δ     単量体   極結合   側面結合   極結合の割合（平均 ± 標準誤差、n={args.seeds}）")
    for (c, r), dts in DEFAULT_CASES.items():
        for dt in sorted(dts, reverse=True):
            group = [x for x in rows if (x["cylinder_um"], x["radius_um"], x["dt_s"]) == (c, r, dt)]
            def mean_sem(key):
                v = np.array([g[key] for g in group], float)
                return f"{v.mean():5.2f}±{v.std(ddof=1) / np.sqrt(len(v)):.2f}"
            fractions = np.array([g["polar"] / max(g["polar"] + g["lateral"], 1) for g in group])
            print(f"({c:3.1f},{r:4.2f}) {dt:8.5f} {group[0]['delta']:5.2f}  {mean_sem('monomers')}  "
                  f"{mean_sem('polar')}  {mean_sem('lateral')}  "
                  f"{fractions.mean():.3f}±{fractions.std(ddof=1) / np.sqrt(len(fractions)):.3f}")


if __name__ == "__main__":
    main()
