"""
Monge Cohomology — H¹ computation via homothetic centers and Menelaus ratios.

Key insight: The first cohomology group H¹ of a rigidity matroid can be
computed geometrically via homothetic center collinearity.

For a Laman graph (E = 2V - 3, minimally rigid in ℝ²):
  β₁ = |E| - (2|V| - 3) = 0 for minimally rigid graphs

For a graph with E > 2V - 3 (over-constrained, emergence signal):
  β₁ = E - (2V - 3) > 0 counts the number of linearly independent cycles

Monge enhancement:
  Each triple (i, j, k) defines a Monge line L_ijk
  The collection of all Monge lines spans a vector space
  dim(span{L_ijk}) relates to β₁

The Menelaus theorem gives O(V³) β₁ computation:
  For triangle (i, j, k) with points on sides at ratios α, β, γ:
    (AF/FB)·(BD/DC)·(CE/EA) = -1  (Menelaus)
  
  For homothetic centers S_ij on side (i, k) of triangle i-j-k:
    We can use Menelaus ratios to check collinearity of S_ij, S_jk, S_ki
"""

import numpy as np
from itertools import combinations

Vector2 = np.ndarray


def menelaus_ratio(a: Vector2, b: Vector2, p: Vector2) -> float:
    """
    Signed ratio of AP/PB where P is on line AB.
    
    Returns λ = |AP|/|AB| where sign indicates which side of A.
    """
    A = np.asarray(a, dtype=float)
    B = np.asarray(b, dtype=float)
    P = np.asarray(p, dtype=float)
    
    AB = B - A
    AP = P - A
    
    denom = np.dot(AB, AB)
    if denom < 1e-10:
        return 0.5  # degenerate
    
    return np.dot(AP, AB) / denom


def menelaus_check(p1: Vector2, q1: Vector2, p2: Vector2, q2: Vector2, 
                   p3: Vector2, q3: Vector2, r1: Vector2, r2: Vector2, r3: Vector2,
                   tol: float = 1e-9) -> bool:
    """
    Menelaus theorem: points R1 on P1Q1, R2 on P2Q2, R3 on P3Q3 are collinear
    iff (P1R1/R1Q1)·(Q1R2/R2P2)·(P2R3/R3Q3) = -1.
    
    Args:
        p1, q1: Endpoints of side 1, r1: point on p1-q1 line
        p2, q2: Endpoints of side 2, r2: point on p2-q2 line
        p3, q3: Endpoints of side 3, r3: point on p3-q3 line
    Returns:
        True if Menelaus collinearity holds (product ≈ -1)
    """
    λ1 = menelaus_ratio(p1, q1, r1)
    λ2 = menelaus_ratio(p2, q2, r2)
    λ3 = menelaus_ratio(p3, q3, r3)
    
    # Menelaus: (λ1/(1-λ1)) * ((1-λ2)/λ2) * (λ3/(1-λ3)) = -1
    # Handle edge cases (point at vertex)
    if abs(λ1) < tol or abs(1-λ1) < tol:
        λ1 = 1e-10
    if abs(λ2) < tol or abs(1-λ2) < tol:
        λ2 = 1e-10
    if abs(λ3) < tol or abs(1-λ3) < tol:
        λ3 = 1e-10
    
    product = (λ1 / (1 - λ1)) * ((1 - λ2) / λ2) * (λ3 / (1 - λ3))
    
    return abs(product + 1) < tol


def menelaus_beta1(V: int, edges: list[tuple[int, int]], 
                   circles: list[tuple[Vector2, float]]) -> int:
    """
    Compute β₁ (first Betti number) via Menelaus theorem on homothetic centers.
    
    For each triangle (i, j, k) in the graph:
      - The homothetic center S_ij lies on the line O_i-O_j
      - Check if S_ij, S_jk, S_ki satisfy Menelaus collinearity
      - Each violated triple contributes to β₁
    
    Args:
        V: Number of vertices
        edges: List of (i, j) edge pairs
        circles: List of (center, radius) for each vertex
    
    Returns:
        β₁ = E - (2V - 3) for Laman graphs (should be 0 for minimally rigid)
    """
    E = len(edges)
    expected_E = 2 * V - 3
    
    if E < expected_E:
        # Under-constrained: β₁ = 0 (no rigidity)
        return 0
    
    # Count triples that FAIL Menelaus collinearity
    violations = 0
    total_triples = 0
    
    # Build vertex positions and radii
    centers = [np.asarray(c, dtype=float) for c, r in circles]
    radii = [r for c, r in circles]
    
    # Find all triangles in the graph
    vertex_list = list(range(V))
    
    for combo in combinations(vertex_list, 3):
        i, j, k = combo
        
        # Check if edges (i,j), (j,k), (k,i) exist in the graph
        edge_set = set(map(tuple, edges)) | set(map(lambda e: (e[1], e[0]), edges))
        
        if (i, j) not in edge_set or (j, k) not in edge_set or (k, i) not in edge_set:
            continue  # Not a triangle in this graph
        
        total_triples += 1
        
        # Compute homothetic centers S_ij, S_jk, S_ki
        # Using external center formula: S = (r_j*c_i - r_i*c_j)/(r_j - r_i)
        def ext_center(a, ra, b, rb):
            if abs(rb - ra) < 1e-10:
                return (a + b) / 2.0
            return (rb * a - ra * b) / (rb - ra)
        
        Sij = ext_center(centers[i], radii[i], centers[j], radii[j])
        Sjk = ext_center(centers[j], radii[j], centers[k], radii[k])
        Ski = ext_center(centers[k], radii[k], centers[i], radii[i])
        
        # Check Menelaus: S_ij on side (i,j)?, S_jk on (j,k)?, S_ki on (k,i)?
        # For collinearity: (i,j) side has S_ij, (j,k) side has S_jk, (k,i) side has S_ki
        collinear = menelaus_check(
            centers[i], centers[j],   # side (i,j)
            centers[j], centers[k],   # side (j,k)  
            centers[k], centers[i],   # side (k,i)
            Sij, Sjk, Ski
        )
        
        if not collinear:
            violations += 1
    
    # β₁ = E - (2V - 3), enhanced by Menelaus violations
    # Each Menelaus violation indicates a cycle not rigidly constrained
    base_beta = E - expected_E
    menelaus_contribution = violations
    
    return base_beta + menelaus_contribution


class MongeCohomology:
    """
    Compute cohomology invariants for a fleet constraint graph using Monge geometry.
    
    Key method: beta1() — the first Betti number of the rigidity matroid
    Key method: monge_rank() — the dimension of the Monge line span
    """
    
    def __init__(self, V: int, edges: list[tuple[int, int]], 
                 circles: list[tuple[Vector2, float]]):
        """
        Args:
            V: Number of vertices
            edges: Graph edges
            circles: (center, radius) for each vertex — constraint "size"
        """
        self.V = V
        self.edges = edges
        self.circles = [(np.asarray(c, dtype=float), r) for c, r in circles]
        self.centers = [c for c, r in self.circles]
        self.radii = [r for c, r in self.circles]
        
        # Precompute all homothetic centers
        self._S: dict[tuple[int, int], np.ndarray] = {}
        for i in range(V):
            for j in range(i + 1, V):
                c1, r1 = self.circles[i]
                c2, r2 = self.circles[j]
                if abs(r2 - r1) < 1e-10:
                    S = (c1 + c2) / 2.0
                else:
                    S = (r2 * c1 - r1 * c2) / (r2 - r1)
                self._S[(i, j)] = S
                self._S[(j, i)] = S
    
    def S(self, i: int, j: int) -> np.ndarray:
        """Get homothetic center for vertices i and j."""
        if (i, j) in self._S: return self._S[(i, j)]; return self._S.get((j, i))
    
    def beta1(self) -> int:
        """
        First Betti number β₁ = E - (2V - 3) for Laman rigidity.
        
        In the Monge formulation, β₁ also counts Menelaus violations
        — triples where the three homothetic centers fail collinearity.
        """
        E = len(self.edges)
        base = E - (2 * self.V - 3)
        
        # Count Menelaus violations (Monge deviations)
        violations = 0
        edge_set = set(map(tuple, self.edges)) | set(map(lambda e: (e[1], e[0]), self.edges))
        
        for combo in combinations(range(self.V), 3):
            i, j, k = combo
            if (i, j) not in edge_set or (j, k) not in edge_set or (k, i) not in edge_set:
                continue
            
            Sij = self.S(i, j)
            Sjk = self.S(j, k)
            Ski = self.S(k, i)
            
            if Sij is None or Sjk is None or Ski is None:
                continue
            
            # Cross product area (2x triangle area)
            v1 = Sjk - Sij
            v2 = Ski - Sij
            area = abs(np.cross(v1, v2))
            
            if area > 1e-9:
                violations += 1
        
        return base + violations
    
    def monge_rank(self) -> int:
        """
        Dimension of the span of all Monge lines.
        
        For n vertices in ℝ², Monge lines live in the projective plane P².
        The rank is the dimension of the span of all line directions.
        
        If all Monge lines for all triples pass through a common point → rank 1
        If all Monge lines are parallel → rank 1
        If Monge lines fill the plane → rank 3 (in projective sense)
        """
        # Collect all Monge line directions
        edge_set = set(map(tuple, self.edges)) | set(map(lambda e: (e[1], e[0]), self.edges))
        directions = []
        
        for combo in combinations(range(self.V), 3):
            i, j, k = combo
            if (i, j) not in edge_set or (j, k) not in edge_set or (k, i) not in edge_set:
                continue
            
            Sij = self.S(i, j)
            Sjk = self.S(j, k)
            Ski = self.S(k, i)
            
            if Sij is None or Sjk is None or Ski is None:
                continue
            
            # Direction of Monge line = direction of Sjk - Sij (or any two points)
            d1 = Sjk - Sij
            d2 = Ski - Sij
            if np.linalg.norm(d1) > 1e-10:
                directions.append(d1 / np.linalg.norm(d1))
            if np.linalg.norm(d2) > 1e-10:
                directions.append(d2 / np.linalg.norm(d2))
        
        if not directions:
            return 0
        
        # Compute rank via QR decomposition
        D = np.array(directions)
        # Project to 2D (each direction is a 2D unit vector)
        # For rank computation, use SVD
        U, s, Vt = np.linalg.svd(D[:len(directions)//3] if len(directions) > 3 else D)
        
        # Count singular values above threshold
        rank = sum(1 for sv in s if sv > 1e-6)
        return rank
    
    def emergence_signal(self) -> float:
        """
        Measure of emergence: how much β₁ exceeds the minimally rigid case.
        
        Returns 0 for perfectly rigid (β₁ = 0), scales with over-constraint.
        """
        return float(max(0, self.beta1()))


if __name__ == "__main__":
    # Test: Laman graph (minimally rigid)
    # Triangle with 3 vertices, 3 edges
    V = 3
    edges = [(0, 1), (1, 2), (2, 0)]
    circles = [
        (np.array([0.0, 0.0]), 1.0),
        (np.array([3.0, 0.0]), 2.0),
        (np.array([1.0, 2.5]), 1.5),
    ]
    
    mc = MongeCohomology(V, edges, circles)
    print(f"Laman triangle: V={V}, E={len(edges)}, β₁={mc.beta1()}, monge_rank={mc.monge_rank()}")
    print(f"  Emergence signal: {mc.emergence_signal()}")
    
    # Test: Add an extra edge (over-constrained)
    edges_over = [(0, 1), (1, 2), (2, 0), (0, 2)]  # 4 edges, E > 2V-3
    mc_over = MongeCohomology(V, edges_over, circles)
    print(f"\nOver-constrained: V={V}, E={len(edges_over)}, β₁={mc_over.beta1()}")
    print(f"  Emergence signal: {mc_over.emergence_signal()}")
    
    # Test: 4-vertex Laman graph
    V4 = 4
    E4 = 2 * V4 - 3  # = 5
    edges4 = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]  # K4 minus one edge
    circles4 = [
        (np.array([0.0, 0.0]), 1.0),
        (np.array([4.0, 0.0]), 1.5),
        (np.array([4.0, 3.0]), 1.0),
        (np.array([0.0, 3.0]), 2.0),
    ]
    mc4 = MongeCohomology(V4, edges4, circles4)
    print(f"\n4-vertex Laman: V={V4}, E={E4}, β₁={mc4.beta1()}, monge_rank={mc4.monge_rank()}")
