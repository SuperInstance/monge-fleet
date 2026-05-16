"""
monge-fleet — Monge's theorem for fleet mathematics

Mathematical foundation:
- Monge's theorem: 3 circles → 3 external homothetic centers collinear
- Lassak generalization: n+1 bodies in ℝ^n → (n-1)-flat
- Menelaus ratios: O(V³) β₁ computation via homothetic centers
- Area-based Lyapunov function for consensus convergence
- Wasserstein geometry for optimal transport bounds
"""

from .homothetic import HomotheticCenter, ExternalCenter, external_center, internal_center
from .radical_axis import RadicalAxis, RadicalCochain
from .cohomology import MongeCohomology, menelaus_beta1
from .consensus import MongeConsensus, homonomy_area, homothetic_center
from .pythagorean48 import P48Monge, verify_monge_p48
from .geometry import lift_to_sphere, sandwich_planes, monge_flat

__all__ = [
    "HomotheticCenter", "ExternalCenter", "external_center", "internal_center",
    "RadicalAxis", "RadicalCochain", 
    "MongeCohomology", "menelaus_beta1",
    "MongeConsensus", "homonomy_area", "homothetic_center",
    "P48Monge", "verify_monge_p48",
    "lift_to_sphere", "sandwich_planes", "monge_flat",
]
