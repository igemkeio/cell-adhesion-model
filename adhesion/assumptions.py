"""このモデルが置いている仮定を、明示的に列挙する。

文献値（`adhesion.literature`）と違い、ここにあるのは「測られていないが、
計算を進めるために置いた」判断である。両者を混ぜないために別モジュールにしてある。

各仮定は次の3点をセットで持つ。

    rationale  なぜそう置いたか
    impact     もしこの仮定が破れていたら結果がどう変わるか
    status     文献に支持されているのか、単なる簡略化なのか、未検証なのか

`docs/assumptions.md` はこのモジュールから自動生成される。
モデル関数は `uses` でどの仮定に依存しているかを宣言する。
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Status(Enum):
    """仮定の性格。読む側が信頼度を判断できるようにする。"""

    SUPPORTED = "文献に支持されている"
    SIMPLIFICATION = "意図的な簡略化"
    UNTESTED = "未検証・要確認"

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True)
class Assumption:
    """モデルが置いている仮定1件。"""

    id: str
    statement: str
    rationale: str
    impact: str
    status: Status

    def __str__(self) -> str:
        return f"[{self.id}] {self.statement}"


A1_RIGID_SPHEROCYLINDER = Assumption(
    id="A1",
    statement="細胞を変形しない剛体スフェロシリンダー（円柱＋両端半球）として扱う。",
    rationale=(
        "大腸菌はペプチドグリカン層により形状を保つ。接着で生じる力は pN オーダーで、"
        "細胞壁を目に見えて変形させる大きさではない。"
    ),
    impact=(
        "接触部が潰れて平坦化する場合、接触面積は本モデルより大きくなり、"
        "結合数は過小評価になる。ただし平坦化は極接触・側面接触の双方を増やすので、"
        "鎖か束かという比の判定は大きくは動かない。"
    ),
    status=Status.SIMPLIFICATION,
)

A2_BRIDGE_IS_TOTAL_GAP = Assumption(
    id="A2",
    statement=(
        "2細胞の表面間距離が架橋長 2M = 8 nm 以内なら、そこにあるアドヘシンは"
        "結合しうるとみなす。"
    ),
    rationale=(
        "Costan et al. (2024) は 2M を「片方の細胞表面から相手の細胞表面までの全長」"
        "として定義している。したがってこれは2細胞をまたぐ総額であり、"
        "両細胞ぶんを足すのは誤りである。"
    ),
    impact=(
        "接触パッチ面積は 2M に比例するので、架橋長を取り違えると結合数が直接その比で"
        "スケールする。極接触は 2M に比例、側面接触は √(2M) に比例するため、"
        "架橋長を長くすると鎖が相対的に有利になる。"
    ),
    status=Status.SUPPORTED,
)

A3_FACING_NORMALS = Assumption(
    id="A3",
    statement=(
        "結合が成立するには、距離条件に加えて両表面の法線が向き合っている"
        "（なす角が 180° から 60° 以内）必要があるとする。"
    ),
    rationale=(
        "アドヘシンは膜から外向きに突き出るため、相手側を向いていない分子は"
        "距離が近くても架橋できない。"
    ),
    impact=(
        "この角度閾値は接触面積の縁をわずかに削るだけで、結合数への影響は数%に留まる。"
        "閾値を 45°–90° で振っても結論は変わらない。"
    ),
    status=Status.SIMPLIFICATION,
)

A4_FUNCTIONAL_FRACTION = Assumption(
    id="A4",
    statement=(
        "接触可能なアドヘシン数に補正係数 f_p = 0.5 を掛けて結合ペア数とする。"
    ),
    rationale=(
        "Costan et al. (2024) の構成では、ナノボディ側が中コピー数プラスミド"
        "（20–30 コピー）、抗原側が低コピー数プラスミド（10–12 コピー）に載っており、"
        "抗原の数が約半分で結合数の律速になる。著者らはこの非対称を 0.5 倍で補正している。"
        "すなわち f_p はプラスミド構成に由来する量であって、分子の機能率ではない。"
    ),
    impact=(
        "結合数に線形に効く。極接触と側面接触の双方に同じ係数がかかるため、鎖か束かの"
        "判定と形態の閾値 f_pole には影響せず、効くのは接合強度の絶対値だけである。"
        "ただしプラスミド構成が変われば f_p も変わる。両アドヘシンを同じコピー数で"
        "発現させれば f_p は 1.0 に近づき、結合数は倍になりうる。"
    ),
    status=Status.SUPPORTED,
)

A5_ONE_BOND_PER_MOLECULE = Assumption(
    id="A5",
    statement="1つのアドヘシン分子は同時に1本しか結合を作らない。",
    rationale="ナノボディ–抗原の結合は1対1であり、複数の相手を同時に掴むことはない。",
    impact=(
        "接触可能分子数を両細胞側で数え、少ない方を採用することで担保している。"
        "この制約を外すと結合数は上限寄りに過大評価される。"
    ),
    status=Status.SUPPORTED,
)

A6_GAUSSIAN_TIP_PROFILE = Assumption(
    id="A6",
    statement=(
        "極への局在は、極先端からの測地距離に対するガウス分布（幅 w）でモデル化する。"
    ),
    rationale=(
        "局在機構が未確定なため、集中の強さを1つの幅パラメータで表せる最も単純な形を選んだ。"
        "w を十分大きくすれば極半球への一様分布に連続的に移行する。"
    ),
    impact=(
        "分布の関数形が変われば接触パッチ内の分子数は変わる。ただし「先端に絞るほど"
        "パッチ内分子数が w⁻² で増える」という定性的な傾向は、先端に単峰性のピークを"
        "持つ分布であれば形によらず成立する。"
    ),
    status=Status.UNTESTED,
)

A7_QUASI_STATIC_CONTACT = Assumption(
    id="A7",
    statement=(
        "2細胞はすでに接触しているものとし、衝突するまでの待ち時間・結合の寿命・"
        "解離速度定数は扱わない。"
    ),
    rationale=(
        "本モデルの問いは「接触したときにどの配置が有利か」であり、"
        "時間発展を入れずに答えられる。速度論を入れるには k_off と 2D 親和性が必要だが、"
        "いずれも本系では測られていない。"
    ),
    impact=(
        "凝集にかかる時間や鎖長分布の時間発展は予測できない。"
        "予測できるのは、接触した2細胞がどの配置を取るかという相対的な有利不利のみ。"
    ),
    status=Status.SIMPLIFICATION,
)

A8_EQUAL_LOAD_SHARING = Assumption(
    id="A8",
    statement="接合部にかかる荷重は、そこにある結合に均等に分配されるとする。",
    rationale=(
        "剛体どうしが法線方向の引張だけを受ける場合、線形バネで結ばれた結合は"
        "同じ変位を受けるため、荷重は実際に均等に分配される。"
    ),
    impact=(
        "曲げモーメントが加わると縁の結合に応力が集中し、実際の最大荷重は"
        "本モデルの見積もりより大きくなる。ただし本系では1結合あたりの荷重が"
        "破断力より3桁小さいため、数倍の集中があっても結論は変わらない。"
    ),
    status=Status.SIMPLIFICATION,
)

A9_LAMINAR_SHEAR = Assumption(
    id="A9",
    statement=(
        "震盪培養中の細胞は、乱流ではなく局所的な層流せん断の中にいるとする。"
    ),
    rationale=(
        "散逸率 0.01–1 W/kg に対して Kolmogorov 長は 32–100 µm であり、"
        "細胞（1.8 µm）よりはるかに大きい。細胞は最小渦の内側にいる。"
    ),
    impact=(
        "この条件が崩れる（極端に激しい撹拌）と、細胞スケールの速度勾配が"
        "本モデルの見積もりより大きくなる。通常のフラスコ震盪では成立する。"
    ),
    status=Status.SUPPORTED,
)

A10_UNIFORM_IS_BASELINE = Assumption(
    id="A10",
    statement=(
        "局在を操作しない場合のアドヘシン分布は、細胞表面に一様であるとする。"
    ),
    rationale=(
        "Costan et al. (2024) は同一構成のアドヘシンについて分布を蛍光強度から直接測り、"
        "極と中央で有意差がないと報告している（t検定 p = 0.56）。"
    ),
    impact=(
        "これが本モデルの基準点であり、結論そのものを左右する。"
        "一様分布は側面接触を強く有利にするため、束を予測する。"
    ),
    status=Status.SUPPORTED,
)


ALL_ASSUMPTIONS: list[Assumption] = [
    A1_RIGID_SPHEROCYLINDER,
    A2_BRIDGE_IS_TOTAL_GAP,
    A3_FACING_NORMALS,
    A4_FUNCTIONAL_FRACTION,
    A5_ONE_BOND_PER_MOLECULE,
    A6_GAUSSIAN_TIP_PROFILE,
    A7_QUASI_STATIC_CONTACT,
    A8_EQUAL_LOAD_SHARING,
    A9_LAMINAR_SHEAR,
    A10_UNIFORM_IS_BASELINE,
]

BY_ID: dict[str, Assumption] = {item.id: item for item in ALL_ASSUMPTIONS}


def uses(*ids: str) -> list[Assumption]:
    """関数が依存する仮定を宣言するためのヘルパ。

    存在しない ID を書いた時点で失敗するので、仮定を消したり ID を変えたりすると
    参照側も必ず直すことになる。
    """

    return [BY_ID[identifier] for identifier in ids]


def assumptions_markdown() -> str:
    """仮定の一覧を Markdown で返す。docs/assumptions.md の中身になる。"""

    sections = [
        "| ID | 仮定 | 性格 |",
        "|---|---|---|",
    ]
    for item in ALL_ASSUMPTIONS:
        sections.append(f"| {item.id} | {item.statement} | {item.status} |")

    sections.append("")
    for item in ALL_ASSUMPTIONS:
        sections.extend(
            [
                f"### {item.id}　{item.statement}",
                "",
                f"**性格**　{item.status}",
                "",
                f"**根拠**　{item.rationale}",
                "",
                f"**破れた場合の影響**　{item.impact}",
                "",
            ]
        )
    return "\n".join(sections)
