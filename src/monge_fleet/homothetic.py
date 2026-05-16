"""
Homothetic Centers — external and internal similarity centers of two circles.

For two circles (O₁, r₁) and (O₂, r₂):
  - External homothetic center Sₑ: lies on O₁O₂ line, satisfies OS₁/r₁ = OS₂/r₂
    Sₑ = (r₂·O₁ - r₁·O₂)/(r₂ - r₁)
  - Internal homothetic center Sᵢ: lies between O₁ and O₂
    Sᵢ = (r₂·O₁ + r₁·O₂)/(r₁ + r₂)

Monge's theorem: for any 3 circles, the 3 external homothetic centers are collinear.
This is a PROJECTIVE invariant — holds for any configuration of circles.
"""

import numpy as np
from numpy.typing import ArrayLike

Vector2 = tuple[float, float] | np.ndarray


def external_center(c1: Vector2, r1: float, c2: Vector2, r2: float) -> np.ndarray:
    """
    External homothetic center of two circles.
    
    The external center lies on the line connecting the two centers,
    on the side of the larger circle.
    
    Args:
        c1: Center of circle 1 (x, y)
        r1: Radius of circle 1
        c2: Center of circle 2 (x, y)
        r2: Radius of circle 2
    
    Returns:
        External homothetic center coordinates
    
    Formula:
        Sₑ = (r₂·O₁ - r₁·O₂)/(r₂ - r₁)
    
    Note:
        When r₁ ≈ r₂, use limiting position at midpoint.
    """
    O1 = np.asarray(c1, dtype=float)
    O2 = np.asarray(c2, dtype=float)
    
    # Handle near-equal radii (degenerate case → use midpoint)
    if abs(r2 - r1) < 1e-10:
        return (O1 + O2) / 2.0
    
    Se = (r2 * O1 - r1 * O2) / (r2 - r1)
    return Se


def internal_center(c1: Vector2, r1: float, c2: Vector2, r2: float) -> np.ndarray:
    """
    Internal homothetic center of two circles.
    
    The internal center always lies between the two centers.
    
    Args:
        c1: Center of circle 1 (x, y)
        r1: Radius of circle 1
        c2: Center of circle 2 (x, y)
        r2: Radius of circle 2
    
    Returns:
        Internal homothetic center coordinates
    
    Formula:
        Sᵢ = (r₂·O₁ + r₁·O₂)/(r₁ + r₂)
    """
    O1 = np.asarray(c1, dtype=float)
    O2 = np.asarray(c2, dtype=float)
    Si = (r2 * O1 + r1 * O2) / (r1 + r2)
    return Si


class HomotheticCenter:
    """
    Represents the homothetic (similarity) center structure of two circles.
    
    A pair of circles has both an external and internal homothetic center.
    The external center is used in Monge's theorem.
    """
    
    def __init__(self, c1: Vector2, r1: float, c2: Vector2, r2: float):
        self.c1 = np.asarray(c1, dtype=float)
        self.c2 = np.asarray(c2, dtype=float)
        self.r1 = r1
        self.r2 = r2
        self.external = external_center(c1, r1, c2, r2)
        self.internal = internal_center(c1, r1, c2, r2)
    
    def monge_ratio(self) -> float:
        """
        Returns the ratio that places the external center on O1-O2 line.
        Used in Menelaus theorem for collinearity checking.
        
        Returns:
            λ = r₂/(r₂ - r₁) for external center placement
        """
        if abs(self.r2 - self.r1) < 1e-10:
            return 0.5  # midpoint
        return self.r2 / (self.r2 - self.r1)
    
    def __repr__(self) -> str:
        return f"HomotheticCenter(external={self.external}, internal={self.internal})"


class ExternalCenter:
    """External homothetic center only — used in Monge's theorem."""
    
    def __init__(self, circles: list[tuple[Vector2, float]]):
        """
        Args:
            circles: List of (center, radius) tuples for n circles
        """
        self.circles = [(np.asarray(c, dtype=float), r) for c, r in circles]
        self.n = len(circles)
        # Compute all pairwise external centers
        self._centers: dict[tuple[int, int], np.ndarray] = {}
        for i in range(self.n):
            for j in range(i + 1, self.n):
                c1, r1 = self.circles[i]
                c2, r2 = self.circles[j]
                self._centers[(i, j)] = external_center(c1, r1, c2, r2)
    
    def get(self, i: int, j: int) -> np.ndarray:
        """Get external center for circles i and j."""
        if (i, j) in self._centers:
            return self._centers[(i, j)]
        return self._centers.get((j, i), None)
    
    def monge_line_deviation(self, i: int, j: int, k: int) -> float:
        """
        Measure deviation from Monge collinearity for triple (i, j, k).
        
        Returns the signed area of triangle formed by the three external centers.
        Zero area = perfect Monge collinearity.
        
        Args:
            i, j, k: Indices of three circles
        Returns:
            Area of triangle Sₑᵢⱼ Sₑⱼₖ Sₑₖᵢ (2x the signed area)
        """
        Sij = self.get(i, j)
        Sjk = self.get(j, k)
        Ski = self.get(k, i)
        
        if Sij is None or Sjk is None or Ski is None:
            return float('nan')
        
        # Cross product for signed area
        v1 = Sjk - Sij
        v2 = Ski - Sij
        return float(np.abs(np.cross(v1, v2)))
    
    def all_triples_collinear(self, tol: float = 1e-9) -> bool:
        """
        Check if all triples of circles satisfy Monge collinearity.
        
        This is the Monge consistency condition for a set of circles.
        
        Args:
            tol: Numerical tolerance for area check
        Returns:
            True if all triples are collinear (Monge holds)
        """
        for i in range(self.n):
            for j in range(i + 1, self.n):
                for k in range(j + 1, self.n):
                    if self.monge_line_deviation(i, j, k) > tol:
                        return False
        return True
    
    def menelaus_ratio(self, i: int, j: int, k: int, point: Vector2) -> float:
        """
        Compute Menelaus ratio for point on line through circles i and j.
        
        Menelaus theorem: for collinear points A, B, C on sides of triangle,
        (AF/FB)·(BD/DC)·(CE/EA) = -1
        
        Here used to check collinearity of homothetic centers.
        
        Args:
            i, j: Circle indices defining the line
            k: Third circle index (forming triangle with i and j)
            point: Test point
        Returns:
            Menelaus ratio (should be ±1 for collinear)
        """
        Sij = self.get(i, j)
        if Sij is None:
            return float('nan')
        
        P = np.asarray(point, dtype=float)
        c_i, r_i = self.circles[i]
        c_j, r_j = self.circles[j]
        
        # Signed ratio along O_i-O_j line
        # Using directed segment ratios
        dot_ij = np.dot(P - Sij, c_j - c_i)
        len_ij_sq = np.dot(c_j - c_i, c_j - c_i)
        
        return dot_ij / len_ij_sq if len_ij_sq > 1e-10 else 0.0


# --- Verification / Tests ---
if __name__ == "__main__":
    import math
    
    # Test 1: Three circles — verify Monge collinearity
    circles = [
        (np.array([0.0, 0.0]), 1.0),
        (np.array([3.0, 0.0]), 2.0),
        (np.array([1.0, 2.0]), 1.5),
    ]
    
    ec = ExternalCenter(circles)
    print("Test 1: Monge collinearity for 3 arbitrary circles")
    for i, j, k in [(0, 1, 2)]:
        dev = ec.monge_line_deviation(i, j, k)
        print(f"  Triple ({i},{j},{k}) deviation area: {dev:.2e}")
    
    print(f"  All triples collinear: {ec.all_triples_collinear()}")
    
    # Test 2: 48 Pythagorean directions — verify Monge holds
    from .pythagorean48 import generate_p48_directions
    
    p48 = generate_p48_directions()
    p48_circles = [(np.array([math.cos(a), math.sin(a)]), 1.0) for a in p48]
    
    ec48 = ExternalCenter(p48_circles)
    
    # Sample check: first 10 triples
    print("\nTest 2: Pythagorean48 Monge check (sample triples)")
    violations = 0
    for i in range(min(10, len(p48))):
        for j in range(i + 1, min(i + 5, len(p48))):
            for k in range(j + 1, min(j + 3, len(p48))):
                dev = ec48.monge_line_deviation(i, j, k)
                if dev > 1e-9:
                    violations += 1
                    print(f"  VIOLATION: ({i},{j},{k}) area={dev:.2e}")
    
    if violations == 0:
        print("  No violations in sample — Monge holds for P48")
    
    # Test 3: Collinearity deviation as emergence detector
    print("\nTest 3: Emergence detection via Monge deviation")
    # Add noise to one circle to break rigidity
    noisy_circles = circles.copy()
    noisy_circles[1] = (np.array([3.0 + 0.3, 0.0 + 0.1]), 2.0 + 0.2)
    
    ec_noisy = ExternalCenter(noisy_circles)
    dev = ec_noisy.monge_line_deviation(0, 1, 2)
    print(f"  Clean deviation: {ec.monge_line_deviation(0,1,2):.2e}")
    print(f"  Noisy deviation:  {dev:.6f}")
    print(f"  Emergence signal: {dev > 1e-6}")
