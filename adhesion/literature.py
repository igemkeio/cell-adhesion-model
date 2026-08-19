"""文献から取った値を、出典付きの第一級オブジェクトとして定義する。

モデルが使う物理定数はすべてこのモジュールに集約し、値・単位・不確かさ・出典を
一体で持たせる。数値をコード中に直書きしないことで、「どの数字がどの文献から
来たのか」を後から必ず追跡できるようにしている。

`docs/parameters.md` はこのモジュールから自動生成される（`scripts/export_docs.py`）。
したがってコードとドキュメントが食い違うことはない。

文献値ではない量（モデル固有の仮定）は `adhesion.assumptions` に置く。
両者を混ぜないこと。
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Reference:
    """引用文献1件。"""

    key: str
    authors: str
    title: str
    journal: str
    year: int
    doi: str
    volume: str = ""
    pages: str = ""

    @property
    def short(self) -> str:
        """本文中で使う短い表記。例: 'Costan et al. (2024)'"""
        first = self.authors.split(",")[0].split(" and ")[0].strip()
        surname = first.split()[-1] if first else self.key
        multiple = (" and " in self.authors) or ("," in self.authors)
        return f"{surname} et al. ({self.year})" if multiple else f"{surname} ({self.year})"

    @property
    def url(self) -> str:
        return f"https://doi.org/{self.doi}"


REFERENCES: dict[str, Reference] = {
    "costan2024": Reference(
        key="costan2024",
        authors="Costan, Sara A. and Kim, Honesty and Others",
        title=(
            "Biophysical characterization of synthetic adhesins for predicting "
            "and tuning engineered living material properties"
        ),
        journal="Matter",
        year=2024,
        volume="7",
        pages="2836-2854",
        doi="10.1016/j.matt.2024.03.019",
    ),
    "glass2018": Reference(
        key="glass2018",
        authors="Glass, David S. and Riedel-Kruse, Ingmar H.",
        title=(
            "A Synthetic Bacterial Cell-Cell Adhesion Toolbox for Programming "
            "Multicellular Morphologies and Patterns"
        ),
        journal="Cell",
        year=2018,
        volume="174",
        pages="649-658",
        doi="10.1016/j.cell.2018.06.041",
    ),
    "rassam2015": Reference(
        key="rassam2015",
        authors="Rassam, Patrice and Copeland, Nicholas A. and Others",
        title="Supramolecular assemblies underpin turnover of outer membrane proteins in bacteria",
        journal="Nature",
        year=2015,
        volume="523",
        pages="333-336",
        doi="10.1038/nature14461",
    ),
    "tirado1979": Reference(
        key="tirado1979",
        authors="Tirado, Maria M. and Garcia de la Torre, Jose",
        title=(
            "Translational friction coefficients of rigid, symmetric top macromolecules. "
            "Application to circular cylinders"
        ),
        journal="The Journal of Chemical Physics",
        year=1979,
        volume="71",
        pages="2581-2587",
        doi="10.1063/1.438613",
    ),
}


@dataclass(frozen=True)
class LiteratureValue:
    """実測値1つ。値そのものではなく、出所ごと持ち回すための型。

    `float()` で数値になるので、計算式の中では通常の数値と同じように書ける。
    """

    symbol: str
    value: float
    unit: str
    quantity: str
    source: str
    uncertainty: float | None = None
    note: str = ""

    def __post_init__(self) -> None:
        if self.source not in REFERENCES:
            raise KeyError(f"未登録の出典です: {self.source}")

    def __float__(self) -> float:
        return float(self.value)

    @property
    def reference(self) -> Reference:
        return REFERENCES[self.source]

    @property
    def formatted(self) -> str:
        """例: '16.1 ± 0.4 pN'"""
        magnitude = abs(self.value)
        digits = 0 if magnitude >= 1000 else 2 if magnitude < 10 else 1
        text = f"{self.value:,.{digits}f}"
        if self.uncertainty is not None:
            text += f" ± {self.uncertainty:,.{digits}f}"
        return f"{text} {self.unit}".strip()


# --------------------------------------------------------------------------
# 細胞の形状
# --------------------------------------------------------------------------
CELL_LENGTH = LiteratureValue(
    symbol="l",
    value=1.80,
    uncertainty=0.25,
    unit="µm",
    quantity="細胞の全長",
    source="costan2024",
    note="スフェロシリンダー（円柱＋両端半球）として扱ったときの端から端まで。",
)

CELL_RADIUS = LiteratureValue(
    symbol="r",
    value=0.40,
    uncertainty=0.10,
    unit="µm",
    quantity="細胞の半径",
    source="costan2024",
)

# --------------------------------------------------------------------------
# アドヘシン分子
# --------------------------------------------------------------------------
BRIDGE_LENGTH = LiteratureValue(
    symbol="2M",
    value=8.0,
    unit="nm",
    quantity="アドヘシン対が架橋できる全長",
    source="costan2024",
    note=(
        "片方の細胞表面から、ナノボディと抗原の両方を経て、相手の細胞表面までの全長。"
        "2細胞をまたぐ総額なので、両細胞ぶんを足してはいけない。"
    ),
)

RUPTURE_FORCE = LiteratureValue(
    symbol="F_b",
    value=16.1,
    uncertainty=0.4,
    unit="pN",
    quantity="単一アドヘシン結合の破断力",
    source="costan2024",
    note="光ピンセットによる測定。荷重速度 5–100 pN/s の範囲で値は変わらなかった。",
)

MEMBRANE_DIFFUSION = LiteratureValue(
    symbol="D",
    value=0.36,
    uncertainty=0.06,
    unit="µm²/s",
    quantity="外膜内でのアドヘシンの拡散係数",
    source="costan2024",
)

MOBILE_FRACTION = LiteratureValue(
    symbol="M_f",
    value=0.37,
    uncertainty=0.04,
    unit="",
    quantity="拡散可能なアドヘシンの割合",
    source="costan2024",
    note="残り63%は不動。外膜タンパクの島構造と整合的。",
)

# --------------------------------------------------------------------------
# 発現量
# --------------------------------------------------------------------------
ADHESINS_AT_MID_INDUCTION = LiteratureValue(
    symbol="N_100",
    value=7300,
    uncertainty=1800,
    unit="分子/細胞",
    quantity="アドヘシン数（100 ng/mL aTc 誘導時）",
    source="costan2024",
)

ADHESINS_AT_MAX_INDUCTION = LiteratureValue(
    symbol="N_300",
    value=15300,
    uncertainty=4100,
    unit="分子/細胞",
    quantity="アドヘシン数（300 ng/mL aTc、飽和誘導時）",
    source="costan2024",
)

SURFACE_DENSITY_AT_MAX = LiteratureValue(
    symbol="σ_300",
    value=3400,
    uncertainty=1200,
    unit="分子/µm²",
    quantity="表面密度（飽和誘導時）",
    source="costan2024",
    note="外膜タンパク全種の総密度（10⁵ /µm² のオーダー）の数%にすぎない。",
)

PRODUCTION_RATE = LiteratureValue(
    symbol="a_max",
    value=9900,
    uncertainty=2100,
    unit="分子/細胞/h",
    quantity="最大産生速度",
    source="costan2024",
)

DECAY_RATE = LiteratureValue(
    symbol="b_eff",
    value=0.78,
    uncertainty=0.12,
    unit="1/h",
    quantity="アドヘシンの実効減衰速度",
    source="costan2024",
    note="大半は細胞分裂による希釈（0.82 ± 0.30 /h）。分解自体は 0.05 ± 0.01 /h と遅い。",
)

# --------------------------------------------------------------------------
# 検証に使う報告値（モデルの入力ではなく、答え合わせの対象）
# --------------------------------------------------------------------------
REPORTED_BONDS_PARALLEL = LiteratureValue(
    symbol="n_∥",
    value=130,
    uncertainty=30,
    unit="対",
    quantity="平行に並んだ2細胞間の結合数",
    source="costan2024",
    note="100 ng/mL 誘導時。モデルの検証対象。",
)

REPORTED_BONDS_PERPENDICULAR = LiteratureValue(
    symbol="n_⊥",
    value=7,
    uncertainty=2,
    unit="対",
    quantity="直交した2細胞間の結合数",
    source="costan2024",
    note="100 ng/mL 誘導時。接触面積が最小になる配置。モデルの検証対象。",
)

REPORTED_CONTACT_FRACTION_MIN = LiteratureValue(
    symbol="A_con,min",
    value=0.002,
    unit="",
    quantity="接触面積が細胞表面積に占める割合（下限）",
    source="costan2024",
)

REPORTED_CONTACT_FRACTION_MAX = LiteratureValue(
    symbol="A_con,max",
    value=0.036,
    unit="",
    quantity="接触面積が細胞表面積に占める割合（上限）",
    source="costan2024",
    note="配置に依存して 0.002–0.036 の範囲で変わると報告されている。",
)

# --------------------------------------------------------------------------
# 環境
# --------------------------------------------------------------------------
WATER_VISCOSITY = LiteratureValue(
    symbol="η",
    value=1.0,
    unit="mPa·s",
    quantity="培地の動粘性（水で近似）",
    source="tirado1979",
    note="Tirado らの抵抗係数式と組で使うため、便宜上ここに置く。",
)


ALL_VALUES: list[LiteratureValue] = [
    CELL_LENGTH,
    CELL_RADIUS,
    BRIDGE_LENGTH,
    RUPTURE_FORCE,
    MEMBRANE_DIFFUSION,
    MOBILE_FRACTION,
    ADHESINS_AT_MID_INDUCTION,
    ADHESINS_AT_MAX_INDUCTION,
    SURFACE_DENSITY_AT_MAX,
    PRODUCTION_RATE,
    DECAY_RATE,
    REPORTED_BONDS_PARALLEL,
    REPORTED_BONDS_PERPENDICULAR,
    REPORTED_CONTACT_FRACTION_MIN,
    REPORTED_CONTACT_FRACTION_MAX,
    WATER_VISCOSITY,
]


def values_from(source: str) -> list[LiteratureValue]:
    """特定の文献に由来する値だけを取り出す。"""
    if source not in REFERENCES:
        raise KeyError(f"未登録の出典です: {source}")
    return [value for value in ALL_VALUES if value.source == source]


def parameters_markdown() -> str:
    """文献値の対応表を Markdown で返す。docs/parameters.md の中身になる。"""

    lines = [
        "| 記号 | 量 | 値 | 出典 | 備考 |",
        "|---|---|---:|---|---|",
    ]
    for value in ALL_VALUES:
        reference = value.reference
        citation = f"[{reference.short}]({reference.url})"
        lines.append(
            f"| `{value.symbol}` | {value.quantity} | {value.formatted} "
            f"| {citation} | {value.note} |"
        )
    return "\n".join(lines)


def bibliography_markdown() -> str:
    """引用文献一覧を Markdown で返す。"""

    lines = []
    for reference in REFERENCES.values():
        detail = reference.journal
        if reference.volume:
            detail += f" **{reference.volume}**"
        if reference.pages:
            detail += f", {reference.pages}"
        lines.append(
            f"- {reference.authors} ({reference.year}). {reference.title}. "
            f"*{detail}*. <{reference.url}>"
        )
    return "\n".join(lines)
