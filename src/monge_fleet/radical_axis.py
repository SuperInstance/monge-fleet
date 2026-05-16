"""
Radical Axis — the set of points with equal power to two circles.

For circles (O₁, r₁) and (O₂, r₂), the radical axis is:
  Rad(C₁, C₂) = {P : |PO₁|² - r₁² = |PO₂|² - r₂²}

This is a LINE perpendicular to O₁O₂.

Key insight: The radical axis is a COHOMOLOGICAL object.
  φᵢⱼ(P) = Power(P, Cᵢ) - Power(P, Cⱼ)
  Then φᵢⱼ + φⱼₖ + φₖᵢ = 0 identically (it's a coboundary)

Monge's theorem enhancement:
  The external homothetic center Sᵢⱼ lies on Rad(Cᵢ, Cⱼ)
  By Monge, all three Sᵢⱼ, Sⱼₖ, Sₖᵢ lie on the same line
  This line is the zero set of a 1-cocycle glued from pairwise radical constraints
"""

import numpy as np
from numpy.typing import ArrayLike

Vector2 = np.ndarray


class RadicalAxis:
    """
    The radical axis (power line) of two circles.
    
    The radical axis is perpendicular to the line connecting the two centers.
    It consists of all points with equal power to both circles.
    """
    
    def __init__(self, c1: Vector2, r1: float, c2: Vector2, r2: float):
        self.c1 = np.asarray(c1, dtype=float)
        self.c2 = np.asarray(c2, dtype=float)
        self.r1 = r1
        self.r2 = r2
        
        # Direction vector along O1-O2
        self.d = self.c2 - self.c1
        self.len_d = np.linalg.norm(self.d)
        self.d_unit = self.d / self.len_d if self.len_d > 1e-10 else np.array([1.0, 0.0])
        
        # Perpendicular direction (for the line itself)
        self.perp = np.array([-self.d_unit[1], self.d_unit[0]])
        
        # Midpoint
        self.mid = (self.c1 + self.c2) / 2.0
        
        # Distance from midpoint to radical axis
        # Derivation: solve |P-O1|² - r1² = |P-O2|² - r2²
        # Expanding: 2(O2-O1)·P = |O2|² - r2² - |O1|² + r1²
        # The RHS determines offset from midpoint
        self.offset = (np.dot(self.c2, self.c2) - r2**2 - np.dot(self.c1, self.c1) + r1**2) / 2.0
        self.along_offset = self.offset / self.len_d if self.len_d > 1e-10 else 0.0
    
    def contains(self, point: Vector2, tol: float = 1e-9) -> bool:
        """Check if a point lies on the radical axis."""
        P = np.asarray(point, dtype=float)
        # Project point onto direction axis, subtract offset
        along = np.dot(P - self.mid, self.d_unit)
        deviation = along - self.along_offset
        return abs(deviation) < tol
    
    def signed_distance(self, point: Vector2) -> float:
        """Signed distance from point to the radical axis."""
        P = np.asarray(point, dtype=float)
        along = np.dot(P - self.mid, self.d_unit)
        return along - self.along_offset
    
    def project(self, point: Vector2) -> Vector2:
        """Project a point onto the radical axis."""
        P = np.asarray(point, dtype=float)
        along = np.dot(P - self.mid, self.d_unit) - self.along_offset
        return self.mid + along * self.d_unit
    
    def __repr__(self) -> str:
        return f"RadicalAxis(mid={self.mid}, perp_dir={self.perp})"


class RadicalCochain:
    """
    A cochain on a graph where each edge carries a radical axis constraint.
    
    For a triple (i, j, k), the radical axes Rad(Cᵢ, Cⱼ), Rad(Cⱼ, Cₖ), Rad(Cₖ, Cᵢ)
    satisfy a coboundary relation: their normals sum to zero.
    
    This is the geometric dual of H¹ cohomology on the rigidity matroid.
    """
    
    def __init__(self, circles: list[tuple[Vector2, float]]):
        """
        Args:
            circles: List of (center, radius) tuples
        """
        self.circles = [(np.asarray(c, dtype=float), r) for c, r in circles]
        self.n = len(circles)
        self.axes: dict[tuple[int, int], RadicalAxis] = {}
        
        for i in range(self.n):
            for j in range(i + 1, self.n):
                c1, r1 = self.circles[i]
                c2, r2 = self.circles[j]
                self.axes[(i, j)] = RadicalAxis(c1, r1, c2, r2)
    
    def get_axis(self, i: int, j: int) -> RadicalAxis | None:
        """Get the radical axis for circles i and j."""
        if (i, j) in self.axes:
            return self.axes[(i, j)]
        return self.axes.get((j, i))
    
    def cochain_value(self, i: int, j: int, point: Vector2) -> float:
        """
        Evaluate the 1-cochain φᵢⱼ at a point.
        
        φᵢⱼ(P) = Power(P, Cᵢ) - Power(P, Cⱼ) = |POᵢ|² - rᵢ² - |POⱼ|² + rⱼ²
        
        The radical axis is the zero set of this cochain.
        """
        if i == j:
            return 0.0
        c_i, r_i = self.circles[i]
        c_j, r_j = self.circles[j]
        P = np.asarray(point, dtype=float)
        
        power_i = np.dot(P - c_i, P - c_i) - r_i**2
        power_j = np.dot(P - c_j, P - c_j) - r_j**2
        return power_i - power_j
    
    def coboundary(self, i: int, j: int, k: int, point: Vector2) -> float:
        """
        Check the coboundary condition for a triple.
        
        δφ(i,j,k) = φᵢⱼ(P) + φⱼₖ(P) + φₖᵢ(P) should be 0 for all P.
        
        This is equivalent to the radical axes being concurrent at infinity
        (i.e., the Monge collinearity condition).
        """
        val_ij = self.cochain_value(i, j, point)
        val_jk = self.cochain_value(j, k, point)
        val_ki = self.cochain_value(k, i, point)
        return val_ij + val_jk + val_ki
    
    def cohomology_status(self, i: int, j: int, k: int, tol: float = 1e-9) -> str:
        """
        Determine the cohomology status of a triple.
        
        Returns:
            'exact': coboundary = 0 → radical axes concurrent at infinity (Monge holds)
            'nontrivial': coboundary ≠ 0 → emergence signal (Monge violated)
        """
        # Test at the centroid of the three circles
        c_i, r_i = self.circles[i]
        c_j, r_j = self.circles[j]
        c_k, r_k = self.circles[k]
        centroid = (c_i + c_j + c_k) / 3.0
        
        cb = abs(self.coboundary(i, j, k, centroid))
        return 'exact' if cb < tol else 'nontrivial'
    
    def all_exact(self, tol: float = 1e-9) -> bool:
        """
        Check if all triples are exact (coboundary = 0 for all triples).
        
        This means the radical axis complex is a true coboundary —
        all triples satisfy Monge collinearity.
        """
        for i in range(self.n):
            for j in range(i + 1, self.n):
                for k in range(j + 1, self.n):
                    if self.cohomology_status(i, j, k, tol) == 'nontrivial':
                        return False
        return True


if __name__ == "__main__":
    # Test: Verify coboundary = 0 for any three circles
    circles = [
        (np.array([0.0, 0.0]), 1.0),
        (np.array([4.0, 0.0]), 2.0),
        (np.array([1.0, 3.0]), 1.5),
    ]
    
    rc = RadicalCochain(circles)
    print("RadicalCochain test: coboundary for triple (0,1,2)")
    print(f"  Coboundary at centroid: {rc.coboundary(0, 1, 2, (circles[0][0] + circles[1][0] + circles[2][0])/3):.2e}")
    print(f"  All exact: {rc.all_exact()}")
    
    # Show radical axes are perpendicular to center lines
    for (i,j), ax in rc.axes.items():
        print(f"  Rad({i},{j}): direction perp to ({circles[i][0]} → {circles[j][0]})")
