"""
monge-fleet — Monge's theorem for fleet mathematics

Mathematical foundation:
- Monge's theorem: 3 circles → 3 external homothetic centers collinear
- Lassak generalization: n+1 bodies in ℝ^n → (n-1)-flat
- Menelaus ratios: O(V³) β₁ computation via homothetic centers
- Area-based Lyapunov function for consensus convergence
- Wasserstein geometry for optimal transport bounds
- Zero holonomy ↔ Monge collinearity ↔ Byzantine detection
- Pythagorean48 zero drift via Monge consistency

Modules:
    geometry.py          — 3D lifting, sandwich planes, Monge flat
    homothetic.py        — External/internal homothetic centers
    radical_axis.py      — Radical axis as 1-cocycle
    cohomology.py        — H¹ computation via Menelaus
    consensus.py         — MongeConsensus, Zero Holonomy Consensus
    pythagorean48.py     — Original P48 module (top-level)
    pythagorean48/       — Deep dive: directions, collinearity, drift analysis
    transport/           — Wasserstein geometry, optimal transport
    bridge/              — Monge-CT integration
    proofs/              — Formal proofs (convergence, Byzantine, zero holonomy)
"""

from .homothetic import HomotheticCenter, ExternalCenter, external_center, internal_center
from .radical_axis import RadicalAxis, RadicalCochain
from .cohomology import MongeCohomology, menelaus_beta1
from .consensus import MongeConsensus, homonomy_area, homothetic_center
from .pythagorean48 import P48Monge, verify_monge_p48
from .geometry import lift_to_sphere, sandwich_planes, monge_flat

# Re-export transport modules for convenience
from .transport import (
    TransportPlan, wasserstein_2, compute_lambda, wasserstein_decay_bound,
    EntropicTransport, sinkhorn_knopp,
    FleetTransport, fleet_consensus_w2,
    c_transform, kantorovich_dual,
)

__all__ = [
    # Core geometry
    "HomotheticCenter", "ExternalCenter", "external_center", "internal_center",
    "RadicalAxis", "RadicalCochain",
    "MongeCohomology", "menelaus_beta1",
    "MongeConsensus", "homonomy_area", "homothetic_center",
    "P48Monge", "verify_monge_p48",
    "lift_to_sphere", "sandwich_planes", "monge_flat",
    # Transport
    "TransportPlan", "wasserstein_2", "compute_lambda", "wasserstein_decay_bound",
    "EntropicTransport", "sinkhorn_knopp",
    "FleetTransport", "fleet_consensus_w2",
    "c_transform", "kantorovich_dual",
]
