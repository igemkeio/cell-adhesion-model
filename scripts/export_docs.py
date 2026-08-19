"""docs/ 以下のドキュメントをコードから生成する。

文献値の表と仮定の一覧は、`adhesion.literature` と `adhesion.assumptions` を
唯一の情報源として自動生成する。手で書き写さないことで、コードとドキュメントが
食い違う余地をなくしている。

    python3 scripts/export_docs.py
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from adhesion import literature as lit
from adhesion.assumptions import ALL_ASSUMPTIONS, Status, assumptions_markdown

DOCS = Path(__file__).resolve().parent.parent / "docs"

GENERATED_NOTICE = (
    "<!-- このファイルは scripts/export_docs.py が生成します。直接編集しないでください。 -->"
)


def write_parameters() -> Path:
    counts = {
        key: len(lit.values_from(key)) for key in lit.REFERENCES if lit.values_from(key)
    }
    summary = "、".join(
        f"{lit.REFERENCES[key].short} から {count} 件" for key, count in counts.items()
    )

    content = f"""{GENERATED_NOTICE}

# 使用した文献値

モデルが使う物理定数の一覧。すべて `adhesion/literature.py` で出典と一体で
定義されており、コード中に数値を直書きしている箇所はない。

内訳は{summary}。

{lit.parameters_markdown()}

## 値の性格について

上の表は2種類の値を含む。

- **モデルの入力**: 細胞寸法、架橋長、破断力、拡散係数、発現量など。
  これらはモデルに与える。
- **検証の対象**: 報告された結合数（`n_∥`, `n_⊥`）と接触面積割合（`A_con`）。
  これらはモデルに与えず、モデルが独立に再現できるかを確かめるために使う。
  検証結果は `tests/test_validation.py` に固定してある。

## 引用文献

{lit.bibliography_markdown()}

---
生成日: {date.today().isoformat()}
"""
    path = DOCS / "parameters.md"
    path.write_text(content, encoding="utf-8")
    return path


def write_assumptions() -> Path:
    tally = {status: 0 for status in Status}
    for item in ALL_ASSUMPTIONS:
        tally[item.status] += 1
    breakdown = "、".join(f"{status}が{count}件" for status, count in tally.items() if count)

    content = f"""{GENERATED_NOTICE}

# モデルが置いている仮定

このモデルは全部で {len(ALL_ASSUMPTIONS)} 個の仮定を置いている（{breakdown}）。
各仮定には、そう置いた根拠と、もし破れていたら結果がどう変わるかを併記してある。

仮定は `adhesion/assumptions.py` で定義され、各モジュールは自分が依存する仮定を
`uses(...)` で宣言している。存在しない ID を参照すると読み込み時に失敗するため、
仮定を消したり変えたりすれば参照側も必ず追随することになる。

{assumptions_markdown()}

---
生成日: {date.today().isoformat()}
"""
    path = DOCS / "assumptions.md"
    path.write_text(content, encoding="utf-8")
    return path


def main() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    for path in (write_parameters(), write_assumptions()):
        print(f"生成: {path.relative_to(path.parent.parent)}")


if __name__ == "__main__":
    main()
