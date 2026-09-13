"""大腸菌の細胞接着モデル。

極局在させたアドヘシンで大腸菌を接着させる設計について、できあがる多細胞構造が
鎖になるか束になるかを予測する。物理定数はすべて Costan et al., Matter (2024) の
実測値に固定してある。

構成
----
`literature`   文献値。値・単位・不確かさ・出典を一体で持つ
`assumptions`  モデルが置いている仮定。根拠と、破れた場合の影響つき
`geometry`     接触幾何。表面の求積と、架橋しうる分子の数え上げ
`model`        発現量と局在から接着と形態を予測する（時間を持たない）
`dynamics`     多数の細胞が時間とともに出会い、鎖・束へ育っていく過程の時間発展
`shear`        震盪培養での衝突頻度とせん断荷重
`kinetics`     荷重下での結合寿命（別ペアの値を使う補助的な見積もり）

使い方
------
>>> from adhesion.model import predict, UNIFORM_F_POLE
>>> predict(7300, UNIFORM_F_POLE).morphology          # 実測の局在
<Morphology.BUNDLE: '束'>
>>> predict(7300, 0.95, tip_width_um=0.10).morphology  # 極先端に集中させた場合
<Morphology.CHAIN: '鎖'>
"""

from . import assumptions, dynamics, geometry, kinetics, literature, model, shear

__all__ = [
    "assumptions",
    "dynamics",
    "geometry",
    "kinetics",
    "literature",
    "model",
    "shear",
]

__version__ = "1.0.0"
