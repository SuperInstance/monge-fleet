"""
Geometric utilities for Monge's theorem and fleet mathematics.

Includes:
- lift_to_sphere: lift circles to spheres for the 3D Monge proof
- sandwich_planes: two planes that enclose all spheres
- monge_flat: compute the (n-1)-dimensional flat for n+1 bodies (Lassak generalization)
"""

import numpy as np
from typing import Tuple

Vector3 = np.ndarray


def lift_to_sphere(c2: np.ndarray, r: float) -> Tuple[np.ndarray, float]:
    """
    Lift a 2D circle (c, r) to a 3D sphere.
    
    The sphere has the same center and radius as the circle,
    but lives in ℝ³. For Monge's theorem proof.
    """
    return np.asarray(c2, dtype=float), float(r)


def sandwich_planes(spheres: list[tuple[Vector3, float]]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find two parallel planes that sandwich all spheres.
    
    For n spheres in ℝ³, there exist two parallel planes P₁ and P₂
    such that all spheres are between them.
    The planes are perpendicular to the line connecting the centers
    of the two extremal spheres.
    
    Returns:
        (plane1_point, plane1_normal), (plane2_point, plane2_normal)
        Both planes share the same normal.
    """
    if not spheres:
        raise ValueError("No spheres provided")
    
    centers = [np.asarray(c, dtype=float) for c, r in spheres]
    
    # Find the two spheres with maximum separation along each axis
    # The sandwich planes are perpendicular to the axis of maximum spread
    all_points = np.array(centers)
    
    # For Monge's proof, we want planes perpendicular to the line
    # connecting the two furthest-apart centers
    max_dist = 0.0
    p1, p2 = centers[0], centers[1]
    
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            dist = np.linalg.norm(centers[i] - centers[j])
            if dist > max_dist:
                max_dist = dist
                p1, p2 = centers[i], centers[j]
    
    # Direction from p1 to p2
    direction = p2 - p1
    norm_dir = np.linalg.norm(direction)
    
    if norm_dir < 1e-10:
        # All centers coincide — use x-axis as fallback
        normal = np.array([1.0, 0.0, 0.0])
        # Planes at ± max radius from centroid
        centroid = np.mean(centers, axis=0)
        max_r = max(r for c, r in spheres)
        return (centroid - max_r * normal, normal), (centroid + max_r * normal, normal)
    
    normal = direction / norm_dir
    
    # Compute the two extreme positions along this normal
    # Using support function of the sphere union
    support_min = float('inf')
    support_max = float('-inf')
    
    for c, r in spheres:
        proj = np.dot(c, normal)
        if proj - r < support_min:
            support_min = proj - r
        if proj + r > support_max:
            support_max = proj + r
    
    # Two parallel planes at these extremes
    plane1_point = support_min * normal
    plane2_point = support_max * normal
    
    return (plane1_point, normal), (plane2_point, normal)


def monge_flat(spheres: list[tuple[Vector3, float]], dim: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    General Monge flat computation (Lassak generalization).
    
    For n+1 homothetic bodies in ℝ^n, the homothetic centers lie in
    an (n-1)-dimensional affine subspace (the Monge flat).
    
    Args:
        spheres: List of (center, radius) for n+1 spheres
        dim: Dimension of the space (default 2 for plane, 3 for space)
    
    Returns:
        (point_on_flat, direction_basis) — parametric description of the flat
    """
    n = len(spheres)
    if n < 2:
        raise ValueError("Need at least 2 spheres")
    
    if dim == 2:
        # Trivial case: one point (the Monge line in ℝ²)
        # Compute all pairwise homothetic centers
        centers = [np.asarray(c, dtype=float) for c, r in spheres]
        radii = [r for c, r in spheres]
        
        S_list = []
        for i in range(n):
            for j in range(i + 1, n):
                if abs(radii[j] - radii[i]) < 1e-12:
                    S_list.append((centers[i] + centers[j]) / 2.0)
                else:
                    S = (radii[j] * centers[i] - radii[i] * centers[j]) / (radii[j] - radii[i])
                    S_list.append(S)
        
        if not S_list:
            # Degenerate — return arbitrary point and zero direction
            return np.zeros(2), np.zeros((0, 2))
        
        # The Monge "flat" is a line (1D). Fit it through the homothetic centers.
        # Use SVD to find the best-fit line
        S_arr = np.array(S_list)
        centroid = np.mean(S_arr, axis=0)
        
        # Centered points
        centered = S_arr - centroid
        
        if len(centered) < 2:
            return centroid, np.zeros((0, 2))
        
        # SVD — first singular vector gives the line direction
        U, s, Vt = np.linalg.svd(centered, full_matrices=False)
        
        if len(s) == 0:
            return centroid, np.zeros((0, 2))
        
        direction = Vt[0] if len(Vt) > 0 else np.array([1.0, 0.0])
        return centroid, direction.reshape(1, -1)
    
    elif dim == 3:
        # For 4 spheres in ℝ³, the Monge flat is a plane
        if n != 4:
            raise ValueError("Need exactly 4 spheres for 3D Monge plane")
        
        # Compute all 6 pairwise external homothetic centers
        centers = [np.asarray(c, dtype=float) for c, r in spheres]
        radii = [r for c, r in spheres]
        
        S_list = []
        for i in range(4):
            for j in range(i + 1, 4):
                if abs(radii[j] - radii[i]) < 1e-12:
                    S_list.append((centers[i] + centers[j]) / 2.0)
                else:
                    S = (radii[j] * centers[i] - radii[i] * centers[j]) / (radii[j] - radii[i])
                    S_list.append(S)
        
        if len(S_list) < 4:
            raise ValueError("Need at least 4 homothetic centers")
        
        # Fit a plane to the 6 points using SVD
        S_arr = np.array(S_list)
        centroid = np.mean(S_arr, axis=0)
        centered = S_arr - centroid
        
        U, s, Vt = np.linalg.svd(centered, full_matrices=False)
        
        # Normal to the plane = last right singular vector (smallest singular value)
        plane_normal = Vt[-1] if len(Vt) > 0 else np.array([0.0, 0.0, 1.0])
        
        return centroid, plane_normal.reshape(1, -1)
    
    else:
        raise ValueError(f"Unsupported dimension: {dim}")


def monge_projection_matrix(n: int) -> np.ndarray:
    """
    The n×n projection matrix onto the Monge hyperplane.
    
    For a vector x ∈ ℝⁿ, the Monge projection P·x gives the coordinates
    in the (n-1)-dimensional Monge flat (constraint subspace).
    
    For the Lassak theorem with n+1 bodies in ℝⁿ:
    the Monge flat is an (n-1)-dimensional affine subspace.
    The projection matrix has rank n-1.
    
    Returns:
        n×n projection matrix P with P² = P (idempotent)
    """
    I = np.eye(n)
    # The Monge projection removes one degree of freedom (the "height")
    # For n dimensions, we project onto the subspace where coordinates sum to zero
    # This is the standard (n-1)-dim hyperplane
    J = np.ones((n, n))
    P = I - J / n
    return P


# --- Tests ---
if __name__ == "__main__":
    print("=== Geometry utilities test ===\n")
    
    # Test 1: Lift circles to spheres
    print("Test 1: Lifting circles to spheres")
    c1, r1 = np.array([1.0, 0.0]), 1.0
    c2, r2 = np.array([4.0, 0.0]), 2.0
    c3, r3 = np.array([0.0, 3.0]), 1.5
    
    s1 = lift_to_sphere(c1, r1)
    s2 = lift_to_sphere(c2, r2)
    s3 = lift_to_sphere(c3, r3)
    print(f"  Circle (1,0) r=1 → Sphere center={s1[0]}, r={s1[1]}")
    
    # Test 2: Sandwich planes for three spheres
    print("\nTest 2: Sandwich planes")
    planes = sandwich_planes([s1, s2, s3])
    print(f"  Plane 1: point={planes[0][0]}, normal={planes[0][1]}")
    print(f"  Plane 2: point={planes[1][0]}, normal={planes[1][1]}")
    
    # Test 3: Monge flat for three circles (line in ℝ²)
    print("\nTest 3: Monge flat (line) for 3 circles")
    circles_2d = [(c1, r1), (c2, r2), (c3, r3)]
    point, direction = monge_flat(circles_2d, dim=2)
    print(f"  Line: point={point}, direction={direction}")
    
    # Test 4: Monge flat for 4 spheres (plane in ℝ³)
    print("\nTest 4: Monge flat (plane) for 4 spheres")
    spheres_3d = [
        (np.array([0.0, 0.0, 0.0]), 1.0),
        (np.array([4.0, 0.0, 0.0]), 1.5),
        (np.array([0.0, 4.0, 0.0]), 1.0),
        (np.array([0.0, 0.0, 4.0]), 2.0),
    ]
    try:
        p3d, n3d = monge_flat(spheres_3d, dim=3)
        print(f"  Plane: point={p3d}, normal={n3d}")
    except Exception as e:
        print(f"  Error: {e}")
    
    # Test 5: Projection matrix
    print("\nTest 5: Monge projection matrices")
    for n in [2, 3, 4, 5]:
        P = monge_projection_matrix(n)
        print(f"  dim={n}: P²-P=0? {np.allclose(P @ P, P)}, rank(P)={np.linalg.matrix_rank(P)}")
