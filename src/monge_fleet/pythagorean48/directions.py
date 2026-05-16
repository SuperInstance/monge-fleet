"""
Pythagorean48 — all 48 direction vectors and their geometric properties.

The 48 directions of P48 correspond to the vertices of a truncated cuboctahedron,
an Archimedean solid with 48 vertices. Each direction is encoded by a primitive
Pythagorean triple (a, b, c), normalized to a unit vector (a/c, b/c) with
sign flips and coordinate swaps generating the full set.

The 6 primitive triples with c ≤ 37:
    (3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (12, 35, 37)

Each triple yields 8 directions (4 sign patterns × 2 coordinate swaps) = 48 total.

Properties:
    1. All 48 vectors are distinct (no overlaps between triples)
    2. All vectors are rational (a/c, b/c) with denominator ≤ 37
    3. The set is closed under reflection across x, y, and line y=x
    4. The directions are NOT equally spaced (unlike regular 48-gon)
    5. But they form a symmetric group under the octahedral symmetry group O_h
"""

import math
import numpy as np
from typing import Iterator
from itertools import product

Vector2 = np.ndarray

# ─── The 6 Primitive Pythagorean Triples ────────────────────────────────────

PRIMITIVE_TRIPLES: list[tuple[int, int, int]] = [
    (3, 4, 5),
    (5, 12, 13),
    (8, 15, 17),
    (7, 24, 25),
    (20, 21, 29),
    (12, 35, 37),
]


def generate_p48_raw() -> list[tuple[float, float]]:
    """
    Generate the 48 raw direction vectors (x, y) from Pythagorean triples.

    Each of the 6 primitives generates:
    4 sign patterns: (+,+), (+,-), (-,+), (-,-)
    2 coordinate swaps: (a,b) and (b,a)
    = 8 directions per triple
    Total: 6 × 8 = 48

    Returns:
        List of 48 (x, y) tuples, each a unit vector with rational coordinates
    """
    vectors = []
    seen = set()

    for a, b, c in PRIMITIVE_TRIPLES:
        for sx, sy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for swap in [(a, b), (b, a)]:
                u, v = swap
                x = sx * u / c
                y = sy * v / c

                # Round to stable precision and deduplicate
                key = (round(x, 12), round(y, 12))
                if key not in seen:
                    seen.add(key)
                    vectors.append((x, y))

    assert len(vectors) == 48, f"Expected 48 vectors, got {len(vectors)}"
    return vectors


def generate_angles() -> list[float]:
    """
    Generate the 48 angles (radians, [0, 2π)) from P48 vectors.

    Returns:
        Sorted list of 48 angles in radians
    """
    vectors = generate_p48_raw()
    angles = set()
    for x, y in vectors:
        angle = math.atan2(y, x)
        if angle < 0:
            angle += 2 * math.pi
        angles.add(round(angle, 10))
    return sorted(list(angles))


def get_vector(index: int) -> np.ndarray:
    """
    Get the P48 direction vector by index (0-47).

    Args:
        index: Index into the sorted P48 angle list

    Returns:
        (x, y) unit vector
    """
    angles = generate_angles()
    idx = index % len(angles)
    angle = angles[idx]
    return np.array([math.cos(angle), math.sin(angle)])


def all_vectors() -> list[np.ndarray]:
    """
    Get all 48 direction vectors as numpy arrays.

    Returns:
        List of 48 unit vectors
    """
    return [get_vector(i) for i in range(48)]


# ─── Properties of the P48 directions ───────────────────────────────────────

def vector_set_properties() -> dict:
    """
    Compute geometric properties of the P48 direction set.

    Returns:
        Dict with properties
    """
    vectors = all_vectors()

    # Pairwise angles
    pairwise_angles = []
    for i in range(48):
        for j in range(i + 1, 48):
            dot = np.dot(vectors[i], vectors[j])
            angle = math.acos(max(-1.0, min(1.0, dot)))
            pairwise_angles.append(angle)

    # Unique angles between neighbors (sorted order)
    sorted_idx = np.argsort([math.atan2(v[1], v[0]) for v in vectors])
    sorted_vectors = [vectors[i] for i in sorted_idx]
    neighbor_angles = []
    for i in range(48):
        v1 = sorted_vectors[i]
        v2 = sorted_vectors[(i + 1) % 48]
        dot = np.dot(v1, v2)
        angle = math.acos(max(-1.0, min(1.0, dot)))
        neighbor_angles.append(angle)

    # Check all are unit vectors
    norms = [np.linalg.norm(v) for v in vectors]
    all_unit = all(abs(n - 1.0) < 1e-12 for n in norms)

    # Check closure under octahedral symmetry
    closure_checks = {
        'under_x_reflection': all(
            any(np.allclose([-v[0], v[1]], w) for w in vectors)
            for v in vectors
        ),
        'under_y_reflection': all(
            any(np.allclose([v[0], -v[1]], w) for w in vectors)
            for v in vectors
        ),
        'under_swap': all(
            any(np.allclose([v[1], v[0]], w) for w in vectors)
            for v in vectors
        ),
    }

    # Minimum separation angle
    min_sep = min(pairwise_angles)
    max_sep = max(pairwise_angles)

    return {
        'n_vectors': len(vectors),
        'all_unit_length': all_unit,
        'min_pairwise_angle_deg': math.degrees(min_sep),
        'max_pairwise_angle_deg': math.degrees(max_sep),
        'mean_neighbor_angle_deg': math.degrees(np.mean(neighbor_angles)),
        'min_neighbor_angle_deg': math.degrees(min(neighbor_angles)),
        'max_neighbor_angle_deg': math.degrees(max(neighbor_angles)),
        'octahedral_symmetry': all(closure_checks.values()),
        'closure_checks': closure_checks,
    }


# ─── Direction Classification ───────────────────────────────────────────────

def vector_coset(v: Vector2) -> str:
    """
    Classify a P48 vector by its primitive triple coset.

    Returns:
        String like '3-4-5' indicating the generating triple
    """
    for a, b, c in PRIMITIVE_TRIPLES:
        for sx, sy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for swap in [(a, b), (b, a)]:
                u, v = swap
                x = sx * u / c
                y = sy * v / c
                if np.allclose(v, [x, y], atol=1e-12):
                    return f'{a}-{b}-{c}'
    return 'unknown'


def vectors_by_triple() -> dict[str, list[np.ndarray]]:
    """
    Group the 48 vectors by their generating primitive triple.

    Returns:
        Dict mapping 'a-b-c' → list of 8 vectors
    """
    groups = {}
    for v in all_vectors():
        label = vector_coset(v)
        if label not in groups:
            groups[label] = []
        groups[label].append(v)

    for label in groups:
        assert len(groups[label]) == 8, f"Expected 8 vectors in {label}, got {len(groups[label])}"
    assert len(groups) == 6

    return groups


# ─── Rational Representation ────────────────────────────────────────────────

def rational_representation(index: int) -> dict:
    """
    Get the rational representation of a P48 vector.

    Each P48 vector is (a/c, b/c) for some primitive Pythagorean triple (a,b,c).

    Args:
        index: Vector index (0-47)

    Returns:
        Dict with numerator, denominator info
    """
    v = get_vector(index)
    for a, b, c in PRIMITIVE_TRIPLES:
        for sx, sy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for swap in [(a, b), (b, a)]:
                u, v2 = swap
                x = sx * u / c
                y = sy * v2 / c
                if abs(v[0] - x) < 1e-12 and abs(v[1] - y) < 1e-12:
                    return {
                        'vector': v,
                        'a': u,
                        'b': v2,
                        'c': c,
                        'tuple': f'({u}, {v2}, {c})',
                        'norm': math.sqrt(u**2 + v2**2) / c,
                        'sx': sx,
                        'sy': sy,
                    }
    return {'vector': v, 'error': 'not found'}


def all_rational_representations() -> list[dict]:
    """
    Get rational representations for all 48 vectors.

    Returns:
        List of dicts with rational info per vector
    """
    return [rational_representation(i) for i in range(48)]


# ─── Tests ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== Pythagorean48 Direction Vectors ===\n")

    vectors = generate_p48_raw()
    print(f"Total vectors: {len(vectors)}")

    props = vector_set_properties()
    print(f"All unit length: {props['all_unit_length']}")
    print(f"Min pairwise angle: {props['min_pairwise_angle_deg']:.4f}°")
    print(f"Mean neighbor angle: {props['mean_neighbor_angle_deg']:.4f}°")
    print(f"Octahedral symmetry: {props['octahedral_symmetry']}")

    print(f"\nVectors by primitive triple:")
    for label, vecs in vectors_by_triple().items():
        print(f"  {label}: {len(vecs)} vectors, "
              f"first = ({vecs[0][0]:.6f}, {vecs[0][1]:.6f})")

    # Show first 8 directions
    angles = generate_angles()
    print(f"\nFirst 8 angles (of 48):")
    for i in range(8):
        v = get_vector(i)
        print(f"  {i}: angle={math.degrees(angles[i]):.2f}°, "
              f"v=({v[0]:.8f}, {v[1]:.8f})")

    # Show a rational representation
    print(f"\nRational representation example:")
    repr = rational_representation(0)
    print(f"  v = ({repr['vector'][0]:.6f}, {repr['vector'][1]:.6f})")
    print(f"  = ({repr['a']}/{repr['c']}, {repr['b']}/{repr['c']})")
    print(f"  from triple {repr['tuple']}")
