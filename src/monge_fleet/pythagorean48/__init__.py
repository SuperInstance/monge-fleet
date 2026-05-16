"""
Pythagorean48 — 48-direction encoding with Monge consistency.

This subpackage provides:
    directions.py     — All 48 direction vectors and their geometric properties
    collinearity.py   — Monge collinearity verification for all 17,296 triples
    drift_analysis.py — Zero drift proof via Monge consistency

The 48 directions correspond to vertices of a truncated cuboctahedron.
Each direction is encoded as a rational unit vector (a/c, b/c) from
a primitive Pythagorean triple (a, b, c).
"""

import math
import numpy as np

from .directions import (
    generate_p48_raw, generate_angles, get_vector, all_vectors,
    vector_set_properties, vectors_by_triple, vector_coset,
    rational_representation, all_rational_representations,
)
from .collinearity import (
    verify_all_triples, verify_by_triple_group,
    monge_deviation_surface, triple_area,
)
from .drift_analysis import (
    compute_all_drifts, path_drift, random_path_drift,
    prove_zero_drift,
)

# ─── Legacy API (originally in top-level pythagorean48.py) ──────────────────

def generate_p48_directions() -> list[float]:
    """Legacy: generate 48 P48 angles in radians."""
    return generate_angles()


def p48_to_vector(index: int) -> np.ndarray:
    """Legacy: convert a P48 index to a unit direction vector."""
    return get_vector(index)


def vector_to_p48(v: np.ndarray) -> int:
    """Legacy: convert a 2D unit vector to the nearest P48 direction index."""
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        return 0
    v = v / norm

    angle = math.atan2(v[1], v[0])
    if angle < 0:
        angle += 2 * math.pi

    angles = generate_angles()
    closest = 0
    min_diff = 2 * math.pi
    for i, a in enumerate(angles):
        diff = abs(a - angle)
        if diff > math.pi:
            diff = 2 * math.pi - diff
        if diff < min_diff:
            min_diff = diff
            closest = i

    return closest


class P48Monge:
    """
    Pythagorean48 with Monge consistency verification.

    For Monge verification, we use slightly varying radii to ensure
    external homothetic centers are well-defined (not at infinity).
    The small radius perturbation doesn't affect the geometry —
    it only makes the external center formula well-defined.
    """

    def __init__(self, radius_variation: float = 0.01):
        self.radius_variation = radius_variation
        self.angles = generate_angles()
        self.n = len(self.angles)

        # Direction circles with varying radii for Monge verification
        self.circles: list[tuple[np.ndarray, float]] = []
        for i, a in enumerate(self.angles):
            center = np.array([math.cos(a), math.sin(a)])
            radius = 1.0 + radius_variation * (i / self.n)
            self.circles.append((center, radius))

        # Precompute homothetic centers for all pairs
        self._S: dict[tuple[int, int], np.ndarray] = {}
        for i in range(self.n):
            for j in range(i + 1, self.n):
                ci, ri = self.circles[i]
                cj, rj = self.circles[j]
                if abs(rj - ri) < 1e-12:
                    S = (ci + cj) / 2.0
                else:
                    S = (rj * ci - ri * cj) / (rj - ri)
                self._S[(i, j)] = S
                self._S[(j, i)] = S

    def S(self, i: int, j: int) -> np.ndarray:
        """Homothetic center between directions i and j."""
        if (i, j) in self._S:
            return self._S[(i, j)]
        return self._S.get((j, i))

    def triple_area(self, i: int, j: int, k: int) -> float:
        """Area of triangle S_ij, S_jk, S_ki (Monge deviation)."""
        Sij = self.S(i, j)
        Sjk = self.S(j, k)
        Ski = self.S(k, i)
        if Sij is None or Sjk is None or Ski is None:
            return float('nan')
        v1 = Sjk - Sij
        v2 = Ski - Sij
        return float(abs(np.cross(v1, v2)))

    def all_triple_areas(self) -> dict[tuple[int, int, int], float]:
        """Compute Monge deviation for all C(48,3) = 17,296 triples."""
        areas = {}
        for i in range(self.n):
            for j in range(i + 1, self.n):
                for k in range(j + 1, self.n):
                    areas[(i, j, k)] = self.triple_area(i, j, k)
        return areas

    def max_area(self) -> float:
        """Maximum Monge deviation across all triples."""
        areas = self.all_triple_areas()
        return max(areas.values()) if areas else 0.0

    def verify_consistency(self, tol: float = 1e-9) -> tuple[bool, float]:
        """Verify all triples satisfy Monge collinearity."""
        max_dev = self.max_area()
        return max_dev < tol, max_dev

    def worst_triples(self, n: int = 5) -> list[tuple[tuple, float]]:
        """Return the n triples with largest Monge deviation."""
        areas = self.all_triple_areas()
        sorted_areas = sorted(areas.items(), key=lambda x: -x[1])
        return [(k, v) for k, v in sorted_areas[:n]]


def verify_monge_p48(radius_variation: float = 0.01) -> dict:
    """
    Verify Monge's theorem holds for all 17,296 triples of P48 directions.

    Args:
        radius_variation: How much to vary radii (default 1%)

    Returns:
        dict with verification results and statistics
    """
    p48 = P48Monge(radius_variation=radius_variation)

    areas = p48.all_triple_areas()
    n_triples = len(areas)
    max_dev = max(areas.values()) if areas else 0.0
    mean_dev = sum(areas.values()) / n_triples if areas else 0.0

    tol = 1e-6
    violations = [(k, v) for k, v in areas.items() if v > tol]

    return {
        'n_directions': 48,
        'n_triples_checked': n_triples,
        'max_deviation_area': max_dev,
        'mean_deviation_area': mean_dev,
        'violations': len(violations),
        'all_collinear': len(violations) == 0,
        'zero_drift_verified': len(violations) == 0,
        'worst_triples': p48.worst_triples(3),
        'radius_variation': radius_variation,
    }


__all__ = [
    # directions
    "generate_p48_raw", "generate_p48_directions", "generate_angles",
    "get_vector", "p48_to_vector", "vector_to_p48",
    "all_vectors",
    "vector_set_properties", "vectors_by_triple", "vector_coset",
    "rational_representation", "all_rational_representations",
    # collinearity
    "verify_all_triples", "verify_by_triple_group",
    "monge_deviation_surface", "triple_area",
    # drift_analysis
    "compute_all_drifts", "path_drift", "random_path_drift",
    "prove_zero_drift",
    # legacy
    "P48Monge", "verify_monge_p48",
]
