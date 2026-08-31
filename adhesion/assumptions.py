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
        "本モデルの問いは「接触したときにどの配置が有利か」であり、時間発展を"
        "入れずに答えられる。本系のアドヘシンペアについて解離速度定数も2D親和性も"
        "測られていないため、速度論を入れても値を埋められない。"
    ),
    impact=(
        "凝集にかかる時間や鎖長分布の時間発展は予測できない。"
        "予測できるのは、接触した2細胞がどの配置を取るかという相対的な有利不利のみ。"
        "ただし別のナノボディ–抗原ペア（anti-GFP / GFP）については解離速度定数と"
        "障壁位置が測られており、`adhesion.kinetics` で桁の目安を出せる。"
        "それによれば無荷重での結合寿命は約2時間で、震盪培養のせん断はこれをほとんど"
        "変えない。時間を入れても力学が律速にならないという結論は変わらない見込みである。"
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
        "なお Rassam et al. (2015) が示した外膜の空間構造は、古い外膜タンパクが"
        "古い極へ押しやられるというもので、細胞あたり片側の極に偏る。"
        "したがって自然に生じる偏りは対称二極ではなく単極であり、鎖には使えない。"
        "ただしこれは古いタンパクについての話で、恒常発現の定常状態では"
        "上記の実測どおり一様になる。"
    ),
    status=Status.SUPPORTED,
)


A11_ISLAND_CONFINEMENT = Assumption(
    id="A11",
    statement=(
        "アドヘシンは極先端まわりに任意の幅で分布させられるとし、外膜タンパクが作る"
        "島構造による下限は課していない。"
    ),
    rationale=(
        "局在機構が未確定なため、まず幾何がどこまで要求するかを見るのが本モデルの"
        "目的である。そのため到達可能性の制約を先に課さず、必要条件を先に求めている。"
    ),
    impact=(
        "Rassam et al. (2015) によれば外膜タンパクは直径およそ 0.5 µm の島に凝集し、"
        "島の内部に拡散が閉じ込められる。極半球の測地距離は 0.628 µm なので、"
        "島ひとつが極をほぼ覆う大きさにあたる。本モデルが鎖の条件として要求する"
        "先端集中幅 w ≤ 0.20 µm は島より細かく、島が局在の単位であるならば"
        "到達できない可能性がある。すなわち本モデルの「鎖ができる条件」は"
        "必要条件であって、達成可能性を保証しない。"
        "ただし閉じ込めの大きさはタンパクによって変わる。Spector et al. (2010) は"
        "OmpF が直径 100 nm の領域に閉じ込められていると報告しており、"
        "0.5 µm が下限とは限らない。この懸念は当初考えたほど強くない。"
    ),
    status=Status.UNTESTED,
)

A12_STATIC_DISTRIBUTION = Assumption(
    id="A12",
    statement=(
        "接触している間もアドヘシンの分布は変化しないとし、拡散による接触部への"
        "集積は考えない。"
    ),
    rationale=(
        "本モデルは時間を持たないため、接触の継続時間を定義できない（A7）。"
    ),
    impact=(
        "実際には可動なアドヘシンが接触部へ拡散して集まりうるので、結合数は"
        "本モデルより増える可能性がある。Costan et al. 自身もこれを予測誤差の"
        "要因として挙げている。ただし外膜内の拡散は遅く（内在性 OMP で"
        "0.013–0.019 µm²/s、可動割合も 0.37）、島に閉じ込められている場合は"
        "この集積が起きるかどうか自体が不確かである。"
    ),
    status=Status.SIMPLIFICATION,
)


A13_BRIDGE_MUCH_SMALLER_THAN_CAP = Assumption(
    id="A13",
    statement=(
        "架橋長は極半球の大きさに比べて十分小さいとし、極接触と側面接触を"
        "別々の事象として扱う。"
    ),
    rationale=(
        "実測の架橋長 2M = 8 nm では接触パッチの半径が 57 nm で、極半球の測地距離"
        "0.628 µm の 1 割にも満たない。両者は明確に分かれる。"
    ),
    impact=(
        "架橋長を仮想的に大きくすると接触パッチが極半球に匹敵する大きさになり、"
        "極に置いた分子が横並び接触でも拾われるようになる。2M = 128 nm では"
        "全分子を極に置いても極接触と側面接触がほぼ同数になり、判定が数値的な"
        "同着になって意味を失う。リンカー長の設計を検討する場合、本モデルは"
        "2M が 50 nm 程度までしか信頼できない。それ以上を扱うなら、極と側面を"
        "面積で二分する枠組み自体を作り直す必要がある。"
    ),
    status=Status.SIMPLIFICATION,
)


A14_TETHERING_WORKS_FOR_BETA_BARRELS = Assumption(
    id="A14",
    statement=(
        "ペリプラズムでの係留により、βバレル外膜タンパクを極へ集められるとする。"
    ),
    rationale=(
        "Lim & Bernhardt (2019) は内膜に固定した融合体がペリプラズムの相手を極へ"
        "引き寄せること（POLAR）を実証しており、CoilE/CoilR の組でも成立している。"
    ),
    impact=(
        "ただし同論文が外膜側で試したのはリポタンパク（LpoP）であって、βバレル"
        "外膜タンパクではない（本文に β-barrel の語が一度も現れない）。Intimin は"
        "βバレルで外膜に埋め込まれ、島構造に閉じ込められている（Rassam 2015）。"
        "リポタンパクを引き寄せられることは、βバレルを引き寄せられることを"
        "意味しない。この仮定が破れれば局在そのものが成立せず、本モデルの入力が"
        "得られない。資料の実験⑥（係留 vs 全長発現）がこれを直接確かめる実験にあたる。"
    ),
    status=Status.UNTESTED,
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
    A11_ISLAND_CONFINEMENT,
    A12_STATIC_DISTRIBUTION,
    A13_BRIDGE_MUCH_SMALLER_THAN_CAP,
    A14_TETHERING_WORKS_FOR_BETA_BARRELS,
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
