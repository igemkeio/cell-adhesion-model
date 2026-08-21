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
    #: 全文を読んで値を照合したか。False の文献から値を取ってはならない。
    full_text_verified: bool = False
    #: 全文のどこまでを確認したかの記録。
    verification_note: str = ""

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
        authors=(
            "Costan, Stefana A. and Ryan, Paul M. and Kim, Honesty and "
            "Wolgemuth, Charles W. and Riedel-Kruse, Ingmar H."
        ),
        title=(
            "Biophysical characterization of synthetic adhesins for predicting "
            "and tuning engineered living material properties"
        ),
        journal="Matter",
        year=2024,
        volume="7",
        pages="2125-2143",
        doi="10.1016/j.matt.2024.03.019",
        full_text_verified=True,
        verification_note="本文20ページと Supplemental Information（Note 1-6）を全文照合。",
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
        full_text_verified=False,
        verification_note=(
            "未入手。背景として言及するのみで、値は一つも取っていない。"
            "値を取る前に全文を確認すること。"
        ),
    ),
    "klamecka2015": Reference(
        key="klamecka2015",
        authors=(
            "Klamecka, Kamila and Severin, Philip M. and Milles, Lukas F. and "
            "Gaub, Hermann E. and Leonhardt, Heinrich"
        ),
        title="Energy profile of nanobody-GFP complex under force",
        journal="Physical Biology",
        year=2015,
        volume="12",
        pages="056009",
        doi="10.1088/1478-3975/12/5/056009",
            full_text_verified=True,
        verification_note=(
            "本文9ページを全文照合。Δx と AFM 破断力は本論文自身の測定であることを"
            "確認した。一方 k_off と K_D は本論文の測定ではなく Kirchhofer et al. "
            "(2010) および Kubala et al. (2010) からの引用値であり、Bell–Evans 回帰に"
            "固定値として使われているだけなので、本論文の値としては採らない。"
        ),
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
        full_text_verified=True,
        verification_note="本文と Extended Data を全文照合。",
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
    note=(
        "光ピンセットによる測定（n = 14 細胞）。本文の値は 16.1 ± 0.4 pN、"
        "要約表では 16.2 と表記され系統誤差 1.8 pN が併記されている。本モデルは本文の値を採る。"
        "この値は anti-EPEA / EPEA ペアの測定である（EPEA は4アミノ酸のペプチドタグ）。"
        "原著の表が期待値として挙げる 28–45 pN は anti-GFP / GFP ペアの値なので、"
        "両者は別の分子ペアであり食い違いではない。原著自身も『パラメータは"
        "アドヘシンのペアごとに、とくに F_b は変わる』と註記している。"
        "すなわち破断力はペアの選択で変えられる設計変数である。"
    ),
)

MEMBRANE_DIFFUSION = LiteratureValue(
    symbol="D",
    value=0.36,
    uncertainty=0.06,
    unit="µm²/s",
    quantity="外膜内でのアドヘシンの拡散係数",
    source="costan2024",
    note=(
        "FRAP による測定。ただし内在性外膜タンパクについて Rassam et al. (2015) が"
        "報告した値（BtuB で 0.013、Cir で 0.019 µm²/s）より約 20 倍速い。"
        "可動割合が 0.37 にとどまることと合わせて、外膜内の輸送は速いとは言い切れない。"
    ),
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

COPY_NUMBER_CORRECTION = LiteratureValue(
    symbol="f_p",
    value=0.5,
    unit="",
    quantity="プラスミドコピー数の非対称に対する補正係数",
    source="costan2024",
    note=(
        "ナノボディ側は中コピー数プラスミド（20–30 コピー）、抗原側は低コピー数"
        "プラスミド（10–12 コピー）に載っている。抗原の数が約半分で結合数の律速に"
        "なるため、著者らは総分子数を 0.5 倍して結合ペア数を見積もっている。"
        "「幾何的に届く分子のうち結合する割合」ではないことに注意。"
    ),
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
    note="anti-GFP / GFP ペアでの測定。",
)

ADHESINS_AT_MAX_INDUCTION = LiteratureValue(
    symbol="N_300",
    value=15300,
    uncertainty=4100,
    unit="分子/細胞",
    quantity="アドヘシン数（300 ng/mL aTc、飽和誘導時）",
    source="costan2024",
    note="anti-GFP / GFP ペアでの測定（n = 20 細胞）。",
)

SURFACE_DENSITY_AT_MAX = LiteratureValue(
    symbol="σ_300",
    value=3400,
    uncertainty=1200,
    unit="分子/µm²",
    quantity="表面密度（飽和誘導時）",
    source="costan2024",
    note=(
        "参考までに、外膜タンパク全種の総密度は 10⁵ /µm² のオーダーとされるが、"
        "この比較値は一次文献を確認できていない（UNVERIFIED_ESTIMATES 参照）。"
    ),
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

BOND_FORMATION_HALF_TIME = LiteratureValue(
    symbol="t_0.5",
    value=0.1,
    unit="s",
    quantity="結合が形成される半減時間",
    source="costan2024",
    note=(
        "測定は装置の時間分解能（1 s）に阻まれ、上限 t_0.5 < 1 s しか決まっていない。"
        "原著は「典型的な on-rate は 100 ms のオーダーと予想される」と述べており、"
        "ここではその予想値を採る。桁の見積もりとして扱うこと。"
    ),
)

SPONTANEOUS_OFF_RATE = LiteratureValue(
    symbol="k_off",
    value=0.0,
    unit="1/s",
    quantity="荷重なしでの自発的な解離速度",
    source="costan2024",
    note=(
        "荷重速度に対する破断力の傾きが有意に 0 と異ならず、原著は「k_off はゼロと"
        "みなすのが最良の推定であり、アドヘシンは自発的には外れない」と結論している。"
        "したがって接着は事実上不可逆で、律速は形成の側にある。"
    ),
)

# --------------------------------------------------------------------------
# 外膜そのものの性質（別文献。モデルの制約を評価するために使う）
# --------------------------------------------------------------------------
OMP_ISLAND_DIAMETER = LiteratureValue(
    symbol="d_island",
    value=0.5,
    unit="µm",
    quantity="外膜タンパクが作る島の直径",
    source="rassam2015",
    note=(
        "外膜タンパクは数百分子からなる島に凝集し、島の内部に拡散が閉じ込められる。"
        "極半球の測地距離は 0.628 µm なので、島ひとつが極をほぼ覆う大きさにあたる。"
        "これより細かく局在を絞れるかは自明でない（仮定 A11）。"
    ),
)

NATIVE_OMP_DIFFUSION = LiteratureValue(
    symbol="D_native",
    value=0.013,
    unit="µm²/s",
    quantity="内在性外膜タンパク（BtuB）の拡散係数",
    source="rassam2015",
    note=(
        "単一分子追跡による。Cir では 0.019 µm²/s。コリシン標識した BtuB・Cir は"
        "3 分たっても FRAP による蛍光回復を示さず、事実上不動であった。"
        "Costan et al. が合成アドヘシンで測った 0.36 µm²/s とは約 20 倍の開きがある。"
    ),
)

BOND_BARRIER_POSITION = LiteratureValue(
    symbol="Δx",
    value=1.59,
    unit="nm",
    quantity="結合解離のエネルギー障壁までの距離",
    source="klamecka2015",
    note=(
        "単一末端でアンカーした野生型 GFP の値。eGFP では 1.36 nm、両末端"
        "アンカーでは 2.02 nm となり、引く向きによって約 50 % 変わる。"
        "細胞表面での引かれ方は AFM の配置と異なるため、この幅は無視できない。"
    ),
)

NANOBODY_RUPTURE_FORCE_WTGFP = LiteratureValue(
    symbol="F_AFM",
    value=45.0,
    unit="pN",
    quantity="anti-GFP ナノボディと野生型 GFP の破断力（AFM、範囲の上端）",
    source="klamecka2015",
    note=(
        "AFM による直接測定。野生型 GFP で 28–45 pN、eGFP で 41–56 pN。"
        "引く速度 300 nm/s–10 µm/s の範囲での最頻破断力。"
        "Costan et al. が期待値として挙げる 28–45 pN はこの野生型の値にあたる。"
        "本モデルは接着に anti-EPEA ペアの 16.1 pN を使っており、"
        "この値は比較の参照としてのみ保持する。"
    ),
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

ALL_VALUES: list[LiteratureValue] = [
    CELL_LENGTH,
    CELL_RADIUS,
    BRIDGE_LENGTH,
    RUPTURE_FORCE,
    MEMBRANE_DIFFUSION,
    MOBILE_FRACTION,
    COPY_NUMBER_CORRECTION,
    BOND_FORMATION_HALF_TIME,
    SPONTANEOUS_OFF_RATE,
    ADHESINS_AT_MID_INDUCTION,
    ADHESINS_AT_MAX_INDUCTION,
    SURFACE_DENSITY_AT_MAX,
    PRODUCTION_RATE,
    DECAY_RATE,
    REPORTED_BONDS_PARALLEL,
    REPORTED_BONDS_PERPENDICULAR,
    REPORTED_CONTACT_FRACTION_MIN,
    REPORTED_CONTACT_FRACTION_MAX,
    OMP_ISLAND_DIAMETER,
    NATIVE_OMP_DIFFUSION,
    BOND_BARRIER_POSITION,
    NANOBODY_RUPTURE_FORCE_WTGFP,
]


#: 出典を確認できていない見積もり。値として使う前に必ず一次文献を確認すること。
#: LiteratureValue にしていないのは、出典付きの実測値と混ざらないようにするため。
UNVERIFIED_ESTIMATES: dict[str, str] = {
    "外膜タンパクの総面密度": (
        "10^5 /µm^2 のオーダーと述べてきたが、一次文献を確認していない。"
        "極先端への濃縮が物理的に可能かを論じる際の比較対象に使っているため、"
        "仮定 A11 の議論の強さはこの値に依存する。要確認。"
    ),
}


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
