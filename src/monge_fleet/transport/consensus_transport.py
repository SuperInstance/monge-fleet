"""
ConsensusTransport — Zero Holonomy Consensus via Optimal Transport.

Key insight: Consensus is optimal transport.
- Agents = mass points with different costs (radii)
- Consensus target = transport to centroid
- Convergence = transport plan contraction

The Monge line is the axis of the transport plan.
The holonomy area is the "wasted" transport (not along the Monge line).
Zero holonomy = all transport along the Monge line = perfect consensus.
"""

import numpy as np
import math
from typing import Optional

Vector2 = np.ndarray


class ConsensusTransport:
    """
    Fleet consensus via optimal transport.
    
    Each consensus step is a Monge map from current positions to
    updated positions (toward consensus). The transport cost is the
    "energy" spent on consensus.
    
    The convergence is bounded by the Wasserstein distance decay:
      W₂(μ_t, μ_∞) ≤ W₂(μ_0, μ_∞) · λ^t
    where λ comes from the Monge collinearity condition.
    """
    
    def __init__(self, n_agents: int, dim: int = 2):
        self.n = n_agents
        self.dim = dim
        self.positions = [np.zeros(dim) for _ in range(n_agents)]
        self.radii = [1.0] * n_agents
        self.history: list[dict] = []
    
    def set_agent(self, i: int, position: Vector2, radius: float = 1.0):
        self.positions[i] = np.asarray(position, dtype=float)
        self.radii[i] = max(radius, 1e-10)
    
    def centroid(self) -> np.ndarray:
        return np.mean(self.positions, axis=0)
    
    def wasserstein_from_consensus(self) -> float:
        """W₂ distance from current positions to consensus centroid."""
        cent = self.centroid()
        dist_sq = 0.0
        for i in range(self.n):
            diff = self.positions[i] - cent
            dist_sq += np.dot(diff, diff) / (self.radii[i] ** 2)
        return float(np.sqrt(dist_sq))
    
    def step(self, alpha: float = 0.3) -> dict:
        """
        Take one consensus step.
        
        All agents move toward the centroid by fraction alpha.
        Returns metrics about the transport step.
        """
        cent = self.centroid()
        
        # Compute transport
        old_positions = self.positions.copy()
        new_positions = []
        transport_costs = []
        
        for i in range(self.n):
            diff = cent - self.positions[i]
            new_pos = self.positions[i] + alpha * diff
            new_positions.append(new_pos)
            
            diff2 = new_pos - self.positions[i]
            cost = np.dot(diff2, diff2) / (self.radii[i] ** 2)
            transport_costs.append(cost)
        
        self.positions = new_positions
        
        # Compute metrics
        metrics = {
            'step': len(self.history),
            'transport_cost': sum(transport_costs),
            'wasserstein': self.wasserstein_from_consensus(),
            'centroid': cent.copy(),
        }
        self.history.append(metrics)
        
        return metrics
    
    def converge(self, tol: float = 1e-6, max_steps: int = 100, alpha: float = 0.3) -> dict:
        """
        Run consensus until convergence.
        
        Returns:
            dict with convergence stats
        """
        for step in range(max_steps):
            w = self.wasserstein_from_consensus()
            if w < tol:
                return {
                    'converged': True,
                    'steps': step,
                    'final_wasserstein': w,
                    'total_transport_cost': sum(h['transport_cost'] for h in self.history),
                }
            
            self.step(alpha=alpha)
        
        return {
            'converged': False,
            'steps': max_steps,
            'final_wasserstein': self.wasserstein_from_consensus(),
            'total_transport_cost': sum(h['transport_cost'] for h in self.history),
        }
    
    def lambda_bound(self) -> float:
        """
        Compute convergence rate bound λ from radii.
        
        λ = max_{i,j,k} (r_i + r_j) / (r_i + r_j + r_k)
        
        This is the theoretical bound on convergence rate.
        Actual convergence may be faster.
        """
        λ = 0.0
        n = self.n
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    r_i, r_j, r_k = self.radii[i], self.radii[j], self.radii[k]
                    bound = (r_i + r_j) / (r_i + r_j + r_k)
                    λ = max(λ, bound)
        return λ if λ > 0 else 0.5
    
    def convergence_time_ms(self, ε: float = 1e-6) -> float:
        """
        Predicted convergence time in milliseconds.
        
        T ≤ 38ms · log(1/ε) / log(1/λ) · (n / (n - 2f))
        
        For f=0 (no Byzantine agents):
        T = 38 · log(1/ε) / log(1/λ)
        """
        λ = self.lambda_bound()
        
        if λ >= 1.0:
            return float('inf')
        
        iterations = math.log(1.0 / ε) / math.log(1.0 / λ)
        base_time_ms = 38.0  # Network round-trip
        
        return base_time_ms * iterations
    
    def byzantine_convergence_time_ms(self, f: int, ε: float = 1e-6) -> float:
        """
        Convergence time with Byzantine agents.
        
        Byzantine tolerance requires n > 3f.
        Safety factor = n / (n - 2f)
        """
        n = self.n
        if n <= 2 * f:
            return float('inf')
        
        base = self.convergence_time_ms(ε)
        return base * (n / (n - 2 * f))


def fleet_consensus(positions: list[Vector2], radii: list[float],
                    α: float = 0.3, tol: float = 1e-6, max_steps: int = 100) -> dict:
    """
    Run fleet consensus via optimal transport.
    
    Args:
        positions: Initial agent positions
        radii: Trust radii per agent
        α: Step size (fraction toward centroid)
        tol: Convergence tolerance (Wasserstein distance)
        max_steps: Maximum iterations
    
    Returns:
        dict with convergence metrics
    """
    ct = ConsensusTransport(len(positions))
    for i, (pos, r) in enumerate(zip(positions, radii)):
        ct.set_agent(i, pos, r)
    
    result = ct.converge(tol=tol, max_steps=max_steps, alpha=α)
    result['final_positions'] = ct.positions.copy()
    result['final_centroid'] = ct.centroid().copy()
    result['lambda_bound'] = ct.lambda_bound()
    result['predicted_time_ms'] = ct.convergence_time_ms(tol)
    
    return result


if __name__ == "__main__":
    print("=== ConsensusTransport Demo ===\n")
    
    # 5 agents in pentagon with varying radii
    n = 5
    positions = [np.array([math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)]) for i in range(n)]
    radii = [1.0, 1.2, 0.8, 1.5, 1.0]
    
    result = fleet_consensus(positions, radii, α=0.3)
    
    print(f"Converged: {result['converged']}")
    print(f"Steps: {result['steps']}")
    print(f"Final Wasserstein: {result['final_wasserstein']:.6f}")
    print(f"λ bound: {result['lambda_bound']:.4f}")
    print(f"Predicted time: {result['predicted_time_ms']:.1f}ms")
    print(f"Total transport cost: {result['total_transport_cost']:.4f}")
    
    # Show convergence curve
    ct = ConsensusTransport(n)
    for i, (pos, r) in enumerate(zip(positions, radii)):
        ct.set_agent(i, pos, r)
    
    print("\nConvergence curve:")
    for step in range(12):
        m = ct.step(alpha=0.3)
        print(f"  step {step:2d}: W₂={m['wasserstein']:.6f}, cost={m['transport_cost']:.4f}")
