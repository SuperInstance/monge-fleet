"""
Optimal transport layer for fleet consensus.

The Monge problem: find the minimum-cost transport plan between distributions.
In fleet context: models are distributions, coupling is transport plan,
conservation law is the budget constraint.

Key result: Wasserstein distance gives convergence bound for consensus.

Modules:
    monge_map       — Discrete Monge map and basic transport plan computation
    consensus_transport — ConsensusTransport and fleet_consensus
    wasserstein     — Formal W₂ decay proof, TransportPlan with Sinkhorn
    entropy         — Entropic regularization, Sinkhorn-Knopp algorithm
    geometry        — Dual potentials, c-transform, Kantorovich duality
    fleet           — Fleet consensus as optimal transport (orchestrator)
"""

from .monge_map import MongeMap, transport_plan, wasserstein_bound
from .consensus_transport import ConsensusTransport, fleet_consensus
from .wasserstein import (
    TransportPlan, wasserstein_2, wasserstein_to_consensus,
    compute_lambda, wasserstein_decay_bound, convergence_time_ms,
    verify_w2_decay_bound, wasserstein_barycenter,
)
from .entropy import (
    EntropicTransport, sinkhorn_knopp, epsilon_schedule,
)
from .geometry import (
    c_transform, c_conjugate, kantorovich_dual,
    monge_map_from_potential, c_concave_check,
    transport_plan_from_potentials,
)
from .fleet import (
    FleetTransport, fleet_consensus_w2, monge_optimal_plan,
)

__all__ = [
    # From monge_map
    "MongeMap", "transport_plan", "wasserstein_bound",
    # From consensus_transport
    "ConsensusTransport", "fleet_consensus",
    # From wasserstein
    "TransportPlan", "wasserstein_2", "wasserstein_to_consensus",
    "compute_lambda", "wasserstein_decay_bound", "convergence_time_ms",
    "verify_w2_decay_bound", "wasserstein_barycenter",
    # From entropy
    "EntropicTransport", "sinkhorn_knopp", "epsilon_schedule",
    # From geometry
    "c_transform", "c_conjugate", "kantorovich_dual",
    "monge_map_from_potential", "c_concave_check",
    "transport_plan_from_potentials",
    # From fleet
    "FleetTransport", "fleet_consensus_w2", "monge_optimal_plan",
]
