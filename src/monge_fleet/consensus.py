"""
Monge Consensus — Zero Holonomy Consensus with Monge collinearity invariants.

Key insight: Monge's theorem says for any 3 agents with trust radii,
the 3 pairwise homothetic centers are COLLINEAR.
Deviation from collinearity = holonomy = consensus failure.

The Lyapunov function is the AREA of the triangle formed by
the three homothetic centers S_ij, S_jk, S_ki:
  H_ijk = area(S_ij, S_jk, S_ki)

Zero holonomy means H_ijk = 0 for all triples → perfect collinearity.
In practice, H_ijk > 0 indicates the consensus hasn't converged.

Convergence rate bound:
  area(t+1) <= λ · area(t)  where λ depends on trust radius ratios
  
With Byzantine fault tolerance:
  T_converge <= 38ms · log(1/ε) / log(1/λ) · (n / (n - 2f))
"""

import numpy as np
from typing import Optional

Vector2 = np.ndarray


def homonomy_area(c1: Vector2, r1: float, c2: Vector2, r2: float,
                  c3: Vector2, r3: float) -> float:
    """
    Compute the holonomy measure for a triple of agents.
    
    For agents with circles (Oᵢ, rᵢ), compute the external homothetic
    centers S_ij, S_jk, S_ki and return the area of their triangle.
    
    Args:
        c1, c2, c3: Agent center positions
        r1, r2, r3: Agent trust radii
    
    Returns:
        Area of triangle S_ij S_jk S_ki (zero = perfect consensus)
    """
    def ext_center(ai, ra, aj, rb):
        if abs(rb - ra) < 1e-12:
            return (ai + aj) / 2.0
        return (rb * ai - ra * aj) / (rb - ra)
    
    S12 = ext_center(c1, r1, c2, r2)
    S23 = ext_center(c2, r2, c3, r3)
    S31 = ext_center(c3, r3, c1, r1)
    
    # Cross product gives 2x signed area
    v1 = S23 - S12
    v2 = S31 - S12
    area = abs(np.cross(v1, v2))
    
    return float(area)


def homothetic_center(c1: Vector2, r1: float, c2: Vector2, r2: float) -> np.ndarray:
    """External homothetic center of two circles."""
    if abs(r2 - r1) < 1e-12:
        return (np.asarray(c1) + np.asarray(c2)) / 2.0
    return (r2 * np.asarray(c1) - r1 * np.asarray(c2)) / (r2 - r1)


class MongeConsensus:
    """
    Zero Holonomy Consensus with Monge collinearity monitoring.
    
    Each agent i has:
      - position c_i (in ℝ² or ℝⁿ)
      - trust radius r_i (uncertainty scale)
    
    The homothetic center S_ij is the fixed point of the consensus
    between agents i and j, accounting for their different scales.
    
    Monge invariant: For any triple (i, j, k), S_ij, S_jk, S_ki are collinear.
    Deviation from collinearity (triangle area) is the holonomy measure.
    """
    
    def __init__(self, n_agents: int, dim: int = 2):
        """
        Args:
            n_agents: Number of agents in the consensus group
            dim: Dimension of the agent state space
        """
        self.n = n_agents
        self.dim = dim
        
        # State: agent positions and trust radii
        self.positions = [np.zeros(dim) for _ in range(n_agents)]
        self.radii = [1.0] * n_agents
        
        # Homothetic centers for each pair
        self._S: dict[tuple[int, int], np.ndarray] = {}
        
        # History for convergence monitoring
        self.area_history: list[float] = []
        
        self._update_centers()
    
    def set_agent(self, i: int, position: Vector2, radius: float = 1.0):
        """Set the position and trust radius for agent i."""
        self.positions[i] = np.asarray(position, dtype=float)
        self.radii[i] = max(radius, 1e-10)
    
    def _update_centers(self):
        """Recompute all pairwise homothetic centers."""
        for i in range(self.n):
            for j in range(i + 1, self.n):
                self._S[(i, j)] = homothetic_center(
                    self.positions[i], self.radii[i],
                    self.positions[j], self.radii[j]
                )
                self._S[(j, i)] = self._S[(i, j)]
    
    def S(self, i: int, j: int) -> np.ndarray:
        """Get homothetic center between agents i and j."""
        if (i, j) in self._S:
            return self._S[(i, j)]
        return self._S.get((j, i))
    
    def triple_area(self, i: int, j: int, k: int) -> float:
        """Area of triangle formed by S_ij, S_jk, S_ki (holonomy measure)."""
        Sij = self.S(i, j)
        Sjk = self.S(j, k)
        Ski = self.S(k, i)
        
        if Sij is None or Sjk is None or Ski is None:
            return float('nan')
        
        v1 = Sjk - Sij
        v2 = Ski - Sij
        return float(abs(np.cross(v1, v2)))
    
    def all_triple_areas(self) -> dict[tuple[int, int, int], float]:
        """Compute holonomy area for all triples."""
        areas = {}
        for i in range(self.n):
            for j in range(i + 1, self.n):
                for k in range(j + 1, self.n):
                    areas[(i, j, k)] = self.triple_area(i, j, k)
        return areas
    
    def max_area(self) -> float:
        """Maximum holonomy area across all triples (convergence metric)."""
        areas = self.all_triple_areas()
        return max(areas.values()) if areas else 0.0
    
    def convergence_ratio(self, prev_area: float, curr_area: float) -> float:
        """Ratio of consecutive areas — λ in the convergence bound."""
        if prev_area < 1e-12:
            return 1.0 if curr_area < 1e-12 else float('inf')
        return curr_area / prev_area
    
    def lambda_bound(self) -> float:
        """
        Compute the convergence rate bound λ from trust radii.
        
        λ = max_{i,j,k} (r_i + r_j) / (r_i + r_j + r_k)
        
        This comes from the Monge update rule.
        """
        λ = 0.0
        for i in range(self.n):
            for j in range(i + 1, self.n):
                for k in range(j + 1, self.n):
                    r_i, r_j, r_k = self.radii[i], self.radii[j], self.radii[k]
                    bound = (r_i + r_j) / (r_i + r_j + r_k)
                    λ = max(λ, bound)
        return λ if λ > 0 else 0.5  # Default to 0.5 if no triples
    
    def convergence_time(self, ε: float = 1e-6) -> float:
        """
        Predicted convergence time in milliseconds.
        
        T <= 38ms · log(1/ε) / log(1/λ)
        
        Args:
            ε: Desired tolerance (max remaining holonomy area)
        
        Returns:
            Predicted convergence time in ms
        """
        λ = self.lambda_bound()
        
        if λ >= 1.0:
            return float('inf')  # Won't converge
        
        if λ < 1e-10:
            return 38.0  # Very fast convergence
        
        base_time_ms = 38.0  # Network round-trip base
        rate = np.log(1.0 / λ)
        iterations = np.log(1.0 / ε)
        
        return base_time_ms * iterations / rate
    
    def byzantine_convergence_time(self, f: int, ε: float = 1e-6) -> float:
        """
        Convergence time with Byzantine fault tolerance.
        
        T <= T_base · (n / (n - 2f))
        
        Args:
            f: Number of Byzantine (faulty) agents
            ε: Desired tolerance
        """
        n = self.n
        if n <= 2 * f:
            return float('inf')  # Can't tolerate this many faults
        
        base = self.convergence_time(ε)
        safety_factor = n / (n - 2 * f)
        
        return base * safety_factor
    
    def update(self, i: int, new_position: Vector2, new_radius: Optional[float] = None):
        """
        Update agent i's position and/or radius.
        Recomputes all homothetic centers.
        """
        self.positions[i] = np.asarray(new_position, dtype=float)
        if new_radius is not None:
            self.radii[i] = max(new_radius, 1e-10)
        
        self._update_centers()
        
        # Record area for convergence tracking
        area = self.max_area()
        self.area_history.append(area)
    
    def monge_line(self, i: int, j: int, k: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Get the Monge line (collinearity axis) for a triple.
        
        Returns two points on the Monge line (for visualization).
        """
        Sij = self.S(i, j)
        Sjk = self.S(j, k)
        
        if Sij is None or Sjk is None:
            return np.zeros(2), np.zeros(2)
        
        # Direction along the Monge line
        direction = Sjk - Sij
        norm = np.linalg.norm(direction)
        
        if norm < 1e-10:
            return Sij, Sij + np.array([1.0, 0.0])  # Degenerate
        
        direction = direction / norm
        # Extend the line in both directions
        extension = 10.0  # Arbitrary extension distance
        return Sij - extension * direction, Sjk + extension * direction


# --- Emergency verification ---
if __name__ == "__main__":
    import math
    
    print("=== Monge Consensus Verification ===\n")
    
    # Test 1: Three agents with arbitrary radii
    mc = MongeConsensus(3)
    mc.set_agent(0, [0.0, 0.0], 1.0)
    mc.set_agent(1, [4.0, 0.0], 2.0)
    mc.set_agent(2, [1.5, 3.0], 1.5)
    
    print("Test 1: Three agents")
    print(f"  Triple areas: {mc.all_triple_areas()}")
    print(f"  Max area (holonomy): {mc.max_area():.6f}")
    print(f"  Lambda bound: {mc.lambda_bound():.4f}")
    print(f"  Convergence time: {mc.convergence_time():.1f}ms")
    
    # Test 2: Convergence simulation — move agents toward consensus
    print("\nTest 2: Convergence simulation")
    
    for step in range(10):
        # Move each agent slightly toward the centroid
        centroid = np.mean(mc.positions, axis=0)
        for i in range(3):
            mc.positions[i] += 0.1 * (centroid - mc.positions[i])
        mc._update_centers()
        area = mc.max_area()
        print(f"  Step {step}: area = {area:.6f}")
    
    print(f"\n  Final convergence time prediction: {mc.convergence_time():.1f}ms")
    print(f"  With Byzantine f=1: {mc.byzantine_convergence_time(1):.1f}ms")
    
    # Test 3: 5 agents — all triples
    print("\nTest 3: 5 agents, all triples")
    mc5 = MongeConsensus(5)
    np.random.seed(42)
    for i in range(5):
        angle = 2 * math.pi * i / 5
        mc5.set_agent(i, [math.cos(angle), math.sin(angle)], 1.0 + i * 0.2)
    
    areas = mc5.all_triple_areas()
    print(f"  {len(areas)} triples")
    print(f"  Max area: {mc5.max_area():.6f}")
    print(f"  Lambda: {mc5.lambda_bound():.4f}")
    print(f"  Convergence: {mc5.convergence_time():.1f}ms")
    
    # Identify the most problematic triple
    worst = max(areas, key=areas.get)
    print(f"  Worst triple: {worst} with area {areas[worst]:.6f}")
