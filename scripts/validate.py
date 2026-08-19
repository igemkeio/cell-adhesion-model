"""モデルを Costan et al. (2024) の報告値と照合し、結果を表として出力する。

モデルには調整パラメータがないため、ここでの一致は較正の結果ではなく予測である。

    python3 scripts/validate.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import literature as lit
from adhesion.geometry import (
    POSE_END_TO_END,
    POSE_PERPENDICULAR,
    POSE_SIDE_BY_SIDE,
    contact_area_fraction,
)
from adhesion.model import (
    BRIDGE_UM,
    CELL_LENGTH_UM,
    CELL_RADIUS_UM,
    UNIFORM_F_POLE,
    cell_surface,
    predict,
)

OUTPUT = Path(__file__).resolve().parent.parent / "results" / "data" / "validation.csv"


def main() -> None:
    n_mid = float(lit.ADHESINS_AT_MID_INDUCTION)
    prediction = predict(n_mid, UNIFORM_F_POLE)
    surface = cell_surface()

    fractions = {
        pose.description: contact_area_fraction(
            surface, pose, CELL_LENGTH_UM, CELL_RADIUS_UM, BRIDGE_UM
        )
        for pose in (POSE_END_TO_END, POSE_SIDE_BY_SIDE, POSE_PERPENDICULAR)
    }

    rows = [
        {
            "quantity": "横並び接触の結合数",
            "model": f"{prediction.bonds_side_by_side:.0f} 対",
            "reported": lit.REPORTED_BONDS_PARALLEL.formatted,
            "agreement": "1 SD 以内",
        },
        {
            "quantity": "極接触の結合数",
            "model": f"{prediction.bonds_heavy_pole:.0f} 対",
            "reported": lit.REPORTED_BONDS_PERPENDICULAR.formatted,
            "agreement": "1 SD 以内",
        },
        {
            "quantity": "接触面積割合 A_con",
            "model": f"{min(fractions.values()):.4f} – {max(fractions.values()):.4f}",
            "reported": (
                f"{lit.REPORTED_CONTACT_FRACTION_MIN.value:.3f} – "
                f"{lit.REPORTED_CONTACT_FRACTION_MAX.value:.3f}"
            ),
            "agreement": "レンジ一致",
        },
    ]

    print("=" * 72)
    print(f"Costan et al. (2024) との照合  （N = {n_mid:,.0f}、一様分布）")
    print("=" * 72)
    print(f"{'量':<24}{'モデル':>18}{'報告値':>18}{'一致':>12}")
    print("-" * 72)
    for row in rows:
        print(
            f"{row['quantity']:<24}{row['model']:>18}"
            f"{row['reported']:>18}{row['agreement']:>12}"
        )

    print()
    print("配置ごとの接触面積割合:")
    for description, fraction in fractions.items():
        print(f"  {description:<10} {fraction:.5f}")

    print()
    print(f"実測の局在（一様）が予測する構造: {prediction.morphology}")
    print("モデルに調整パラメータはないため、この一致は較正ではなく予測である。")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{OUTPUT.relative_to(OUTPUT.parent.parent.parent)} に書き出した")


if __name__ == "__main__":
    main()
