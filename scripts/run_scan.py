"""局在パターンを掃引し、どの局在がどの多細胞構造を作るかの対応表を出す。

発現量も併せて振るが、形態の閾値が発現量に依存しないことが結果から確認できる。

    python3 scripts/run_scan.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import literature as lit
from adhesion.model import UNIFORM_F_POLE, morphology_threshold, predict
from adhesion.shear import chain_tension_pn

TOTALS = [1000.0, float(lit.ADHESINS_AT_MID_INDUCTION), float(lit.ADHESINS_AT_MAX_INDUCTION)]
F_POLES = [UNIFORM_F_POLE, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]
TIP_WIDTHS = [None, 0.30, 0.25, 0.20, 0.15, 0.10]
SPLITS = [0.5, 0.75, 1.0]

#: 判定の基準にする震盪条件（散逸率 0.1 W/kg 相当）と鎖の長さ。
REFERENCE_SHEAR_PER_S = 300.0
REFERENCE_CHAIN_CELLS = 10


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "results"
        / "data"
        / "localization_scan.csv",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    load_pn = chain_tension_pn(REFERENCE_CHAIN_CELLS, REFERENCE_SHEAR_PER_S)
    rows = []
    for n_total in TOTALS:
        for f_pole in F_POLES:
            for tip_width in TIP_WIDTHS:
                for split in SPLITS:
                    result = predict(n_total, f_pole, split, tip_width)
                    rows.append(
                        {
                            "n_total": n_total,
                            "f_pole": round(f_pole, 4),
                            "polar_split": split,
                            "tip_width_um": "" if tip_width is None else tip_width,
                            "bonds_side_by_side": round(result.bonds_side_by_side, 2),
                            "bonds_heavy_pole": round(result.bonds_heavy_pole, 2),
                            "bonds_light_pole": round(result.bonds_light_pole, 2),
                            "peak_density_per_um2": round(result.peak_density_per_um2, 1),
                            "junction_force_pn": round(result.junction_force_pn, 1),
                            "morphology": result.morphology.name.lower(),
                            "morphology_ja": str(result.morphology),
                            "shear_load_pn": round(load_pn, 2),
                            "safety_factor": round(result.junction_force_pn / load_pn, 1),
                        }
                    )

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} 条件を {args.output.name} に書き出した\n")

    print("=" * 74)
    print("鎖ができるのに必要な f_pole（発現量に依存しない）")
    print("=" * 74)
    for tip_width in TIP_WIDTHS:
        label = "先端集中なし" if tip_width is None else f"先端 {tip_width:.2f} µm"
        threshold = morphology_threshold(tip_width)
        if threshold is None or threshold > 0.99:
            verdict = "事実上達成不能"
        else:
            verdict = f"f_pole > {threshold:.3f}"
        print(f"  {label:>14} : {verdict}")

    print()
    print("=" * 74)
    print(f"代表条件（N = {lit.ADHESINS_AT_MID_INDUCTION.value:,.0f}、100 ng/mL aTc）")
    print("=" * 74)
    print(f"{'局在':>30}{'横並び':>8}{'濃い極':>8}{'薄い極':>8}{'接合強度':>11}{'構造':>14}")
    print("-" * 74)
    highlights = [
        ("一様（実測値）", UNIFORM_F_POLE, 0.5, None),
        ("極80% 対称二極", 0.80, 0.5, None),
        ("極95% 対称二極", 0.95, 0.5, None),
        ("極95% 二極 + 先端0.20µm", 0.95, 0.5, 0.20),
        ("極95% 二極 + 先端0.10µm", 0.95, 0.5, 0.10),
        ("極95% 単極 + 先端0.10µm", 0.95, 1.0, 0.10),
    ]
    for label, f_pole, split, tip_width in highlights:
        result = predict(float(lit.ADHESINS_AT_MID_INDUCTION), f_pole, split, tip_width)
        print(
            f"{label:>30}{result.bonds_side_by_side:8.0f}{result.bonds_heavy_pole:8.0f}"
            f"{result.bonds_light_pole:8.0f}{result.junction_force_pn:9.0f}pN"
            f"{str(result.morphology):>14}"
        )

    print()
    print(
        f"参考: γ̇ = {REFERENCE_SHEAR_PER_S:.0f} /s で {REFERENCE_CHAIN_CELLS} 細胞の鎖に"
        f"かかる張力は {load_pn:.1f} pN。"
    )
    print(
        f"      単結合の破断力は {lit.RUPTURE_FORCE.formatted} なので、"
        "せん断は接合を壊さない。"
    )


if __name__ == "__main__":
    main()
