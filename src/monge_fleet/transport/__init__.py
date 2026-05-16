"""
Optimal transport layer for fleet consensus.

The Monge problem: find the minimum-cost transport plan between distributions.
In fleet context: models are distributions, coupling is transport plan,
conservation law is the budget constraint.

Key result: Wasserstein distance gives convergence bound for consensus.
"""

from .monge_map import MongeMap, transport_plan, wasserstein_bound
from .consensus_transport import ConsensusTransport, fleet_consensus

__all__ = [
    "MongeMap", "transport_plan", "wasserstein_bound",
    "ConsensusTransport", "fleet_consensus",
]
