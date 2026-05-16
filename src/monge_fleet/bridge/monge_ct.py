"""
Monge-CT Bridge — Monge cohomology + constraint theory integration.

This bridge connects Monge's theorem to the existing jc1-ct-bridge constraint theory:
- H1 cohomology (Laman rigidity) + Monge collinearity → enhanced emergence detection
- Zero Holonomy Consensus + Monge area → tighter convergence bounds
- Pythagorean48 + Monge line → gauge-invariant direction encoding

The Monge layer adds a PROJECTIVE INVARIANT to the constraint graph.
Deviation from Monge collinearity = emergence signal.

Usage:
    from monge_fleet.bridge import MongeCTBridge, monge_emergence_detection
    
    bridge = MongeCTBridge(n_agents=5, dim=2)
    bridge.set_agent(0, [0, 0], radius=1.0)
    # ... set all agents ...
    
    signal = bridge.emergence_signal()
    convergence_time = bridge.convergence_time_ms(f=0)
"""

import numpy as np
import math
from typing import Optional

from ..homothetic import ExternalCenter, homothetic_center
from ..cohomology import MongeCohomology
from ..consensus import MongeConsensus
from ..transport.monge_map import consensus_convergence_rate
from ..pythagorean48 import P48Monge, verify_monge_p48


class MongeCTBridge:
    """
    Bridge between Monge cohomology and constraint theory.
    
    Integrates:
    - H1 cohomology via Monge-Menelaus (beta1 computation)
    - Zero Holonomy Consensus via Monge-Lyapunov (convergence bounds)
    - Pythagorean48 via Monge line (direction encoding)
    - Optimal transport via Wasserstein geometry (consensus dynamics)
    """
    
    def __init__(self, n_agents: int, dim: int = 2):
        """
        Args:
            n_agents: Number of agents in the fleet
            dim: Dimension of agent state space (default 2)
        """
        self.n = n_agents
        self.dim = dim
        self.positions = [np.zeros(dim) for _ in range(n_agents)]
        self.radii = [1.0] * n_agents
        
        # Sub-systems
        self._consensus: Optional[MongeConsensus] = None
        self._cohomology: Optional[MongeCohomology] = None
    
    def set_agent(self, i: int, position: list[float], radius: float = 1.0):
        """Set agent i's position and trust radius."""
        self.positions[i] = np.asarray(position, dtype=float)
        self.radii[i] = max(radius, 1e-10)
    
    def _ensure_consensus(self):
        """Lazy initialization of MongeConsensus."""
        if self._consensus is None:
            self._consensus = MongeConsensus(self.n, self.dim)
            for i in range(self.n):
                self._consensus.set_agent(i, self.positions[i], self.radii[i])
    
    def _ensure_cohomology(self):
        """Lazy initialization of MongeCohomology."""
        if self._cohomology is None:
            # Build edges from complete graph (all pairs)
            edges = [(i, j) for i in range(self.n) for j in range(i + 1, self.n)]
            # Circles centered at positions with radii as "constraint size"
            circles = [(self.positions[i], self.radii[i]) for i in range(self.n)]
            self._cohomology = MongeCohomology(self.n, edges, circles)
    
    def update_agent(self, i: int, position: list[float], radius: Optional[float] = None):
        """Update agent i's state (and recompute derived state)."""
        self.positions[i] = np.asarray(position, dtype=float)
        if radius is not None:
            self.radii[i] = max(radius, 1e-10)
        
        # Invalidate derived state
        if self._consensus is not None:
            self._consensus.update(i, self.positions[i], radius)
        self._cohomology = None  # Will be rebuilt on next access
    
    def emergence_signal(self) -> float:
        """
        Compute the emergence signal via Monge deviation.
        
        Emergence is measured as deviation from Monge collinearity.
        For all triples (i, j, k), the homothetic centers S_ij, S_jk, S_ki
        should be collinear. The maximum area of these triangles is the
        emergence signal.
        
        Returns:
            float: Emergence signal (0 = no emergence, > 0 = emergence detected)
        """
        self._ensure_consensus()
        return self._consensus.max_area()
    
    def beta1(self) -> int:
        """
        First Betti number β₁ via Monge-Menelaus.
        
        For minimally rigid graphs: β₁ = 0
        For over-constrained: β₁ > 0 (count of independent cycles)
        
        The Monge enhancement adds Menelaus violation counting.
        """
        self._ensure_cohomology()
        return self._cohomology.beta1()
    
    def convergence_time_ms(self, f: int = 0, ε: float = 1e-6) -> float:
        """
        Predicted convergence time with optional Byzantine tolerance.
        
        Args:
            f: Number of Byzantine (faulty) agents
            ε: Desired tolerance (Wasserstein distance)
        
        Returns:
            float: Predicted convergence time in milliseconds
        """
        self._ensure_consensus()
        
        if f == 0:
            return self._consensus.convergence_time(ε)
        return self._consensus.byzantine_convergence_time(f, ε)
    
    def lambda_bound(self) -> float:
        """
        Convergence rate bound λ from trust radii.
        
        λ = max_{i,j,k} (r_i + r_j) / (r_i + r_j + r_k)
        
        This is the theoretical maximum contraction ratio per step.
        """
        self._ensure_consensus()
        return self._consensus.lambda_bound()
    
    def all_triple_areas(self) -> dict:
        """Get holonomy area for all triples (for debugging/visualization)."""
        self._ensure_consensus()
        return self._consensus.all_triple_areas()
    
    def p48_monke_check(self) -> dict:
        """
        Verify Monge consistency of Pythagorean48 directions.
        
        Returns:
            dict with verification results (zero_drift_verified, etc.)
        """
        return verify_monge_p48()


def monge_emergence_detection(positions: list[list[float]], 
                              radii: list[float]) -> dict:
    """
    Standalone function for Monge-based emergence detection.
    
    Args:
        positions: Agent positions [[x,y], ...]
        radii: Trust radii per agent
    
    Returns:
        dict with emergence metrics
    """
    n = len(positions)
    bridge = MongeCTBridge(n, dim=len(positions[0]))
    
    for i, (pos, r) in enumerate(zip(positions, radii)):
        bridge.set_agent(i, pos, r)
    
    return {
        'emergence_signal': bridge.emergence_signal(),
        'beta1': bridge.beta1(),
        'lambda_bound': bridge.lambda_bound(),
        'convergence_time_ms': bridge.convergence_time_ms(),
        'max_triple_area': bridge.emergence_signal(),
    }


def monge_consensus(positions: list[list[float]], 
                    radii: list[float],
                    f: int = 0,
                    tol: float = 1e-6,
                    max_steps: int = 100) -> dict:
    """
    Run Monge-enhanced consensus.
    
    Args:
        positions: Initial positions
        radii: Trust radii
        f: Number of Byzantine agents
        tol: Convergence tolerance
        max_steps: Maximum iterations
    
    Returns:
        dict with consensus results
    """
    n = len(positions)
    bridge = MongeCTBridge(n, dim=len(positions[0]))
    
    for i, (pos, r) in enumerate(zip(positions, radii)):
        bridge.set_agent(i, pos, r)
    
    # Simulate convergence
    for step in range(max_steps):
        signal = bridge.emergence_signal()
        if signal < tol:
            return {
                'converged': True,
                'steps': step,
                'final_signal': signal,
                'convergence_time_ms': bridge.convergence_time_ms(f, tol),
                'lambda_bound': bridge.lambda_bound(),
            }
        
        # Move all agents toward centroid
        centroid = np.mean([np.asarray(p) for p in positions], axis=0)
        alpha = 0.3
        new_positions = []
        for i, pos in enumerate(positions):
            new_pos = np.asarray(pos) + alpha * (centroid - np.asarray(pos))
            new_positions.append(new_pos.tolist())
            bridge.update_agent(i, new_pos.tolist())
        positions = new_positions
    
    return {
        'converged': False,
        'steps': max_steps,
        'final_signal': bridge.emergence_signal(),
        'convergence_time_ms': bridge.convergence_time_ms(f, tol),
        'lambda_bound': bridge.lambda_bound(),
    }


if __name__ == "__main__":
    print("=== Monge-CT Bridge Demo ===\n")
    
    # 5 agents in pentagon with varying radii
    n = 5
    positions = [[math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)] for i in range(n)]
    radii = [1.0, 1.2, 0.8, 1.5, 1.0]
    
    bridge = MongeCTBridge(n, dim=2)
    for i, (pos, r) in enumerate(zip(positions, radii)):
        bridge.set_agent(i, pos, r)
    
    print(f"Emergence signal: {bridge.emergence_signal():.6f}")
    print(f"β₁: {bridge.beta1()}")
    print(f"λ bound: {bridge.lambda_bound():.4f}")
    print(f"Convergence time (f=0): {bridge.convergence_time_ms(f=0):.1f}ms")
    print(f"Convergence time (f=1): {bridge.convergence_time_ms(f=1):.1f}ms")
    
    # P48 check
    print("\nPythagorean48 Monge verification:")
    p48 = bridge.p48_monke_check()
    print(f"  Triples: {p48['n_triples_checked']:,}")
    print(f"  Violations: {p48['violations']}")
    print(f"  Zero drift: {p48['zero_drift_verified']}")
    
    # Standalone function test
    print("\nStandalone emergence detection:")
    result = monge_emergence_detection(positions, radii)
    print(f"  emergence_signal: {result['emergence_signal']:.6f}")
    print(f"  convergence_time_ms: {result['convergence_time_ms']:.1f}ms")
