"""
Fleet consensus as optimal transport — the Monge-Wasserstein layer.

Key insight: Fleet consensus IS optimal transport.

Each consensus step moves agents from their current positions to new
positions (closer to the centroid). This is exactly a Monge map —
a transport plan that minimizes cost while respecting the marginal
constraints (each agent has equal "mass" in the decision process).

The Wasserstein geometry provides:
1. A metric on the space of fleet configurations: W₂(μ, ν)
2. A convergence bound: W₂(μ_t, μ*) ≤ λ^t · W₂(μ_0, μ*)
3. A barycenter: the consensus limit μ* = argmin_μ Σ_i W₂²(μ_i, μ)

The Monge geometry provides:
1. The homothetic centers as transport fixed points
2. The Monge line as the optimal transport axis
3. The holonomy area as the "wasted" transport cost

For a fleet of n agents with trust radii r_i:
- Mass distribution: μ = (1/n) Σ_i δ_{x_i}
- Transport cost: c(x, y) = ||x - y||² / r²
- Consensus step: T(x_i) = x_i + α · (centroid - x_i)
- Convergence rate: λ = max (r_i + r_j)/(r_i + r_j + r_k)
"""

import numpy as np
import math
from typing import Optional

from ..consensus import homothetic_center
from .wasserstein import (
    TransportPlan, wasserstein_2, wasserstein_to_consensus,
    compute_lambda, wasserstein_decay_bound, convergence_time_ms
)
from .entropy import EntropicTransport, epsilon_schedule
from .geometry import c_transform, monge_map_from_potential

Vector = np.ndarray


class FleetTransport:
    """
    Fleet consensus as optimal transport.

    Orchestrates:
    - Wasserstein distance computation
    - Entropic transport between configurations
    - Monge map extraction from potentials
    - Convergence monitoring via W₂ decay
    - Byzantine fault tolerance via safety factor
    """

    def __init__(self, n_agents: int, dim: int = 2):
        self.n = n_agents
        self.dim = dim
        self.positions = [np.zeros(dim) for _ in range(n_agents)]
        self.radii = [1.0] * n_agents
        self.history: list[dict] = []

    def set_agent(self, i: int, position: Vector, radius: float = 1.0):
        """Set agent i's state."""
        self.positions[i] = np.asarray(position, dtype=float)
        self.radii[i] = max(radius, 1e-10)

    def centroid(self) -> np.ndarray:
        """Current centroid of all agents."""
        return np.mean(self.positions, axis=0)

    def w2_to_consensus(self) -> float:
        """W₂ distance from current configuration to consensus."""
        return wasserstein_to_consensus(self.positions, self.centroid(), self.radii)

    def lambda_bound(self) -> float:
        """Convergence contraction factor λ from trust radii."""
        return compute_lambda(self.radii)

    def convergence_time(self, ε: float = 1e-6, f: int = 0) -> float:
        """Predicted convergence time in ms."""
        return convergence_time_ms(self.radii, self.n, ε, f=f)

    def step(self, alpha: float = 0.3) -> dict:
        """
        Take one consensus step.

        All agents move toward the centroid by fraction alpha.
        The transport plan is constructed and analyzed.

        Args:
            alpha: Step size (fraction toward centroid)

        Returns:
            Metrics dict for this step
        """
        cent = self.centroid()
        old_positions = [p.copy() for p in self.positions]

        # Move toward centroid
        new_positions = []
        for i in range(self.n):
            diff = cent - self.positions[i]
            new_pos = self.positions[i] + alpha * diff
            new_positions.append(new_pos)

        self.positions = new_positions

        # Compute transport plan between old and new
        plan = TransportPlan(
            old_positions, new_positions,
            cost_fn=lambda x, y: float(np.dot(x - y, x - y))
        )
        W2 = plan.wasserstein_distance()

        # Compute entropic transport cost
        et = EntropicTransport(
            old_positions, new_positions,
            cost_fn=lambda x, y: float(np.dot(x - y, x - y))
        )
        pi, info = et.solve(reg=0.05)

        # Compute c-transform of the potential
        f, g = et.dual_potentials()

        # W₂ to final consensus
        W2_cons = self.w2_to_consensus()

        metrics = {
            'step': len(self.history),
            'alpha': alpha,
            'W2_step': W2,
            'W2_to_consensus': W2_cons,
            'transport_cost': info['total_cost'],
            'entropy': info['entropy'],
            'potential_f': f.copy(),
            'potential_g': g.copy(),
            'centroid': cent.copy(),
        }
        self.history.append(metrics)
        return metrics

    def converge(self, tol: float = 1e-6, max_steps: int = 100,
                 alpha: float = 0.3, anneal: bool = False) -> dict:
        """
        Run consensus until W₂ convergence.

        With annealing, the step size decreases over time for smoother
        convergence near the limit.

        Args:
            tol: Convergence tolerance on W₂ to consensus
            max_steps: Maximum iterations
            alpha: Step size (if not annealing)
            anneal: If True, decrease alpha over time

        Returns:
            Convergence summary
        """
        for step in range(max_steps):
            current_alpha = alpha
            if anneal:
                # Decay alpha to avoid oscillations near equilibrium
                current_alpha = alpha * (0.9 ** step)

            m = self.step(alpha=current_alpha)

            if m['W2_to_consensus'] < tol:
                return {
                    'converged': True,
                    'steps': step + 1,
                    'final_W2': m['W2_to_consensus'],
                    'final_cost': m['transport_cost'],
                    'total_transport_cost': sum(
                        h['transport_cost'] for h in self.history
                    ),
                    'λ': self.lambda_bound(),
                    'predicted_time_ms': self.convergence_time(tol),
                }

        return {
            'converged': False,
            'steps': max_steps,
            'final_W2': self.w2_to_consensus(),
            'total_transport_cost': sum(
                h['transport_cost'] for h in self.history
            ),
            'λ': self.lambda_bound(),
        }

    def byzantine_convergence(self, f: int, tol: float = 1e-6) -> dict:
        """
        Convergence estimate with Byzantine agents.

        Byzantine tolerance requires n > 3f.
        The convergence time scales as n/(n - 2f).

        Args:
            f: Number of Byzantine (faulty) agents
            tol: Convergence tolerance

        Returns:
            Convergence estimate
        """
        if self.n <= 2 * f:
            return {
                'possible': False,
                'reason': f'n={self.n} ≤ 2f={2*f}',
                'predicted_time_ms': float('inf'),
            }

        base_time = self.convergence_time(tol, f=0)
        safety = self.n / (self.n - 2 * f)
        byz_time = base_time * safety

        return {
            'possible': True,
            'n': self.n,
            'f': f,
            'n_over_3f': self.n > 3 * f,  # Byzantine agreement threshold
            'base_time_ms': base_time,
            'safety_factor': safety,
            'byz_time_ms': byz_time,
        }

    def gradient_flow(self, dt: float = 0.1, steps: int = 50) -> list[Vector]:
        """
        Wasserstein gradient flow for fleet consensus.

        The gradient flow of the Wasserstein distance to consensus:
            ∂_t μ = -∇_{W₂} F(μ)
        where F(μ) = ½·W₂²(μ, μ*) is the "energy" of being away from consensus.

        In discrete time:
            x_i(t + dt) = x_i(t) - dt · ∇_{x_i} F(μ)
                        = x_i(t) + dt · (x* - x_i(t)) / r_i²

        Returns:
            Final positions after gradient flow
        """
        positions = [p.copy() for p in self.positions]
        x_star = self.centroid()

        for step in range(steps):
            new_positions = []
            for i in range(self.n):
                # Gradient step in Wasserstein space
                grad = (positions[i] - x_star) / (self.radii[i] ** 2)
                new_pos = positions[i] - dt * grad
                new_positions.append(new_pos)
            positions = new_positions

            # Update centroid after each step
            x_star = np.mean(positions, axis=0)

        return positions

    def transport_efficiency(self) -> float:
        """
        Efficiency of transport relative to ideal Monge map.

        Efficiency = W₂(ideal) / W₂(actual)

        For the fleet consensus, the ideal Monge map moves each agent
        directly to the centroid (cost = Σ ||x_i - x*||² / r_i²).
        The actual step may have small deviations.

        Returns 1.0 for perfect efficiency.
        """
        if len(self.history) < 1:
            return 1.0

        last = self.history[-1]
        cent = last['centroid']

        ideal_cost = 0.0
        for i in range(self.n):
            diff = self.positions[i] - cent
            ideal_cost += np.dot(diff, diff) / (self.radii[i] ** 2)

        actual_cost = last['transport_cost']
        if actual_cost < 1e-15:
            return 1.0

        return min(1.0, float(ideal_cost / actual_cost))


def fleet_consensus_w2(positions: list[list[float]],
                       radii: list[float],
                       alpha: float = 0.3,
                       tol: float = 1e-6,
                       max_steps: int = 100,
                       return_history: bool = False) -> dict:
    """
    Run fleet consensus with Wasserstein monitoring.

    Args:
        positions: Initial agent positions
        radii: Trust radii
        alpha: Step size
        tol: Convergence tolerance
        max_steps: Maximum iterations
        return_history: If True, include full convergence history

    Returns:
        Dict with convergence metrics
    """
    n = len(positions)
    ft = FleetTransport(n, dim=len(positions[0]))

    for i, (pos, r) in enumerate(zip(positions, radii)):
        ft.set_agent(i, pos, r)

    result = ft.converge(tol=tol, max_steps=max_steps, alpha=alpha)
    result['final_positions'] = [p.copy() for p in ft.positions]
    result['λ'] = ft.lambda_bound()

    if return_history:
        result['history'] = ft.history

    return result


def monge_optimal_plan(positions: list[Vector],
                        targets: list[Vector],
                        radii: Optional[list[float]] = None) -> dict:
    """
    Compute the Monge optimal transport plan for fleet movement.

    Args:
        positions: Current agent positions
        targets: Target positions (e.g., after consensus step)
        radii: Trust radii (default all 1.0)

    Returns:
        Dict with coupling matrix, W₂ distance, costs
    """
    n = len(positions)
    if radii is None:
        radii = [1.0] * n

    # Mahalanobis cost: ||x - y||² / r²
    def cost_fn(x, y, r=1.0):
        return float(np.dot(x - y, x - y)) / (r ** 2)

    plan = TransportPlan(positions, targets)

    # Compute exact OT
    pi = plan.exact_ot()
    W2 = plan.wasserstein_distance()

    # Per-agent costs
    costs = []
    for i in range(n):
        c = cost_fn(positions[i], targets[i], radii[i])
        costs.append(c)

    return {
        'coupling_matrix': pi,
        'W2_distance': W2,
        'per_agent_costs': costs,
        'total_cost': sum(costs),
    }


# ─── Tests ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import math

    print("=== Fleet Consensus as Optimal Transport ===\n")

    # Test 1: Basic fleet transport
    print("Test 1: Basic fleet transport")
    n = 5
    pos = [np.array([math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)]) for i in range(n)]
    radii = [1.0, 1.2, 0.8, 1.5, 1.0]

    ft = FleetTransport(n)
    for i, (p, r) in enumerate(zip(pos, radii)):
        ft.set_agent(i, p, r)

    print(f"  Initial W₂ to consensus: {ft.w2_to_consensus():.6f}")
    print(f"  λ bound: {ft.lambda_bound():.4f}")
    print(f"  Convergence time (f=0): {ft.convergence_time():.1f}ms")

    # Test 2: Convergence check
    print("\nTest 2: Convergence check")
    result = ft.converge(tol=1e-6, max_steps=50)
    print(f"  Converged: {result['converged']}")
    print(f"  Steps needed: {result['steps']}")
    print(f"  Final W₂: {result['final_W2']:.6f}")
    print(f"  Total transport cost: {result['total_transport_cost']:.4f}")
    print(f"  Predicted time: {result.get('predicted_time_ms', 0):.1f}ms")

    # Test 3: Byzantine tolerance
    print("\nTest 3: Byzantine convergence")
    ft2 = FleetTransport(7)  # 7 agents
    for i in range(7):
        angle = 2 * math.pi * i / 7
        ft2.set_agent(i, [math.cos(angle), math.sin(angle)], 1.0)

    for f in [0, 1, 2]:
        bc = ft2.byzantine_convergence(f)
        print(f"  f={f}: possible={bc['possible']}, time={bc.get('byz_time_ms', float('inf')):.1f}ms")

    # Test 4: Gradient flow
    print("\nTest 4: Wasserstein gradient flow")
    ft3 = FleetTransport(n)
    for i, (p, r) in enumerate(zip(pos, radii)):
        ft3.set_agent(i, p, r)
    final = ft3.gradient_flow(dt=0.1, steps=30)
    print(f"  Final positions (first 3): {final[:3]}")
    print(f"  Final centroid: {np.mean(final, axis=0)}")

    # Test 5: Transport efficiency
    print("\nTest 5: Transport efficiency")
    ft4 = FleetTransport(3)
    ft4.set_agent(0, [0.0, 0.0], 1.0)
    ft4.set_agent(1, [4.0, 0.0], 1.0)
    ft4.set_agent(2, [2.0, 3.0], 1.0)
    ft4.step()
    print(f"  Step efficiency: {ft4.transport_efficiency():.4f}")
