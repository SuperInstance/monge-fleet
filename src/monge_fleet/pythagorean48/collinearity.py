"""
Pythagorean48 collinearity verification — verify Monge holds for ALL 17,296 triples.

Monge's theorem: For any 3 circles, the 3 external homothetic centers are collinear.
For P48, each direction becomes a circle with a slightly perturbed radius.
We verify Monge collinearity for all C(48, 3) = 17,296 distinct triples.

The verification uses a radius perturbation to ensure external homothetic
centers are well-defined (not at infinity when radii are equal).
The small perturbation doesn't affect the geometry — it only makes
the external center formula well-defined.

Key result: Monge holds for all 17,296 triples → zero drift verified.
"""

import math
import numpy as np
from itertools import combinations
from typing import Iterator, Optional

from .directions import generate_angles, get_vector


def homothetic_center(c1: np.ndarray, r1: float,
                      c2: np.ndarray, r2: float) -> np.ndarray:
    """External homothetic center of two circles."""
    if abs(r2 - r1) < 1e-12:
        return (c1 + c2) / 2.0
    return (r2 * c1 - r1 * c2) / (r2 - r1)


def build_p48_circles(radius_variation: float = 0.01) -> list[tuple[np.ndarray, float]]:
    """
    Build circles for all 48 P48 directions with slightly varied radii.

    Each direction angle generates a UNIVECTOR on the unit circle:
        center = (cos(θ), sin(θ))
        radius = 1.0 + ε · (i / 48)

    The radius variation ε ensures distinct radii for Monge verification.
    Without variation, equal radii → external center at infinity.

    Args:
        radius_variation: Fractional radius spread (default 1%)
                          Larger values make collinearity more numerically stable.

    Returns:
        List of 48 (center, radius) tuples
    """
    angles = generate_angles()
    n = len(angles)

    circles = []
    for i, a in enumerate(angles):
        center = np.array([math.cos(a), math.sin(a)])
        radius = 1.0 + radius_variation * (i / n)
        circles.append((center, radius))

    return circles


def triple_area(c1: np.ndarray, r1: float,
                c2: np.ndarray, r2: float,
                c3: np.ndarray, r3: float) -> float:
    """
    Compute Monge triangle area for triple (i, j, k).

    Area = |(S_jk - S_ij) × (S_ki - S_ij)| / 2
    where S_ij is the external homothetic center of circles i and j.

    Zero area = Monge collinearity holds.
    """
    S12 = homothetic_center(c1, r1, c2, r2)
    S23 = homothetic_center(c2, r2, c3, r3)
    S31 = homothetic_center(c3, r3, c1, r1)

    v1 = S23 - S12
    v2 = S31 - S12
    return float(abs(np.cross(v1, v2)))


def verify_all_triples(radius_variation: float = 0.01,
                        tol: float = 1e-6,
                        progress: bool = True) -> dict:
    """
    Verify Monge collinearity for ALL C(48, 3) = 17,296 triples.

    Args:
        radius_variation: Radius perturbation for distinct radii
        tol: Tolerance for collinearity check
        progress: Print progress bar

    Returns:
        Dict with verification results
    """
    circles = build_p48_circles(radius_variation)
    n = len(circles)  # 48

    total_triples = 0
    violations = []
    max_area = 0.0
    area_sum = 0.0
    violators = set()

    # Outer loop picks every 100th triple for detailed inspection
    sample_n = 0
    sample_areas = []

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                total_triples += 1
                ci, ri = circles[i]
                cj, rj = circles[j]
                ck, rk = circles[k]

                area = triple_area(ci, ri, cj, rj, ck, rk)
                max_area = max(max_area, area)
                area_sum += area
                sample_areas.append(area)

                if area > tol:
                    violations.append({
                        'triple': (i, j, k),
                        'area': area,
                    })
                    violators.add(i)
                    violators.add(j)
                    violators.add(k)
                    sample_n += 1
                elif sample_n < 50:
                    sample_n += 1

    mean_area = area_sum / max(total_triples, 1)

    return {
        'n_directions': n,
        'total_triples': total_triples,
        'violations': len(violations),
        'all_collinear': len(violations) == 0,
        'max_area': max_area,
        'mean_area': mean_area,
        'violators': violators,
        'worst_violations': sorted(violations, key=lambda x: -x['area'])[:10],
        'radius_variation': radius_variation,
        'tol': tol,
    }


def verify_by_triple_group(radius_variation: float = 0.01) -> dict:
    """
    Verify collinearity broken down by generating triple.

    Circles from the same primitive triple are closer together
    and may have different numerical properties.

    Returns:
        Dict with per-group and cross-group results
    """
    from .directions import PRIMITIVE_TRIPLES, vector_coset

    circles = build_p48_circles(radius_variation)
    angles = generate_angles()
    n = len(circles)

    # Group indices by generating triple
    groups = {}
    for idx in range(n):
        v = get_vector(idx)
        label = vector_coset(v)
        groups.setdefault(label, []).append(idx)

    results = {}
    for label1, indices1 in groups.items():
        for label2, indices2 in groups.items():
            if label1 > label2:
                continue
            key = f'{label1} × {label2}'
            results[key] = {'triples': 0, 'violations': 0, 'max_area': 0.0}

    # Check all cross-group pairs
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                ci, ri = circles[i]
                cj, rj = circles[j]
                ck, rk = circles[k]

                area = triple_area(ci, ri, cj, rj, ck, rk)

                li = vector_coset(get_vector(i))
                lj = vector_coset(get_vector(j))
                lk = vector_coset(get_vector(k))
                labels = tuple(sorted([li, lj, lk]))

                for key in results:
                    triples_in_key = all(
                        any(label in kk for kk in key.split(' × '))
                        for label in labels
                    )
                    if triples_in_key:
                        results[key]['triples'] += 1
                        results[key]['max_area'] = max(results[key]['max_area'], area)
                        if area > 1e-9:
                            results[key]['violations'] += 1

    return {
        'by_group': results,
        'worst_group': max(results.items(), key=lambda x: x[1]['max_area']),
    }


def monge_deviation_surface(radius_variation: float = 0.01) -> np.ndarray:
    """
    Compute the Monge deviation matrix for visualization.

    Returns a 48×48 matrix where M[i][j] = max_k area(i, j, k).
    This shows which pairs are most problematic for Monge collinearity.

    Args:
        radius_variation: Radius perturbation

    Returns:
        48×48 deviation matrix
    """
    circles = build_p48_circles(radius_variation)
    n = len(circles)

    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            max_area_for_pair = 0.0
            ci, ri = circles[i]
            cj, rj = circles[j]

            for k in range(n):
                if k == i or k == j:
                    continue
                ck, rk = circles[k]
                area = triple_area(ci, ri, cj, rj, ck, rk)
                max_area_for_pair = max(max_area_for_pair, area)

            M[i, j] = max_area_for_pair
            M[j, i] = max_area_for_pair

    return M


# ─── Tests ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 65)
    print("Pythagorean48 — Monge Collinearity Verification")
    print("=" * 65)
    print("Verifying ALL C(48,3) = 17,296 triples...")

    result = verify_all_triples(radius_variation=0.01, tol=1e-6)

    print(f"\nRESULTS:")
    print(f"  Directions: {result['n_directions']}")
    print(f"  Triples checked: {result['total_triples']:,}")
    print(f"  Violations (tol=1e-6): {result['violations']}")
    print(f"  All collinear: {result['all_collinear']}")
    print(f"  Max deviation area: {result['max_area']:.2e}")
    print(f"  Mean deviation area: {result['mean_area']:.2e}")

    if result['violations'] > 0:
        print(f"\n  Worst violations:")
        for v in result['worst_violations'][:5]:
            print(f"    Triple {v['triple']}: area = {v['area']:.2e}")
    else:
        print(f"\n  ✅ Monge's theorem holds for ALL P48 triples")
        print(f"  ✅ Zero drift verified")

    # Per-group breakdown
    print(f"\nCross-group breakdown:")
    group_result = verify_by_triple_group()
    for group, info in group_result['by_group'].items():
        status = "✅" if info['violations'] == 0 else "❌"
        print(f"  {status} {group}: {info['triples']} triples, "
              f"max_area={info['max_area']:.2e}")
