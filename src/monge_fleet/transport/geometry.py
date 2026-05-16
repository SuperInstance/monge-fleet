"""
Transport geometry — dual potentials, c-transform, and Kantorovich duality.

The Kantorovich dual formulation of optimal transport:
    W₂²(μ, ν) = max_{f, g} ∫ f dμ + ∫ g dν
    subject to f(x) + g(y) ≤ c(x, y) for all x, y

For c(x, y) = ||x - y||², the c-transform is:
    g(y) = inf_x (c(x, y) - f(x))

The Stone duality theorem connects OT to Monge geometry:
    The c-concave potential f is the support function of the Monge map.
    For the Monge problem, ∇f(x) = 2(x - T(x)) where T is the transport.

In the fleet context:
    - The dual potentials encode the "value" of being at a position
    - The c-transform constraint says: moving costs can't exceed shared info
    - The optimal potentials give the Monge map (consensus update)
    - The Monge line (collinearity) emerges from the duality gap = 0 condition
"""

import numpy as np
from typing import Optional, Callable, Tuple
from scipy.optimize import minimize_scalar

Vector = np.ndarray


def c_transform(f: np.ndarray,
                source: list[Vector],
                target: list[Vector]) -> np.ndarray:
    """
    Compute the c-transform of potential f.

    c*convex conjugate: g(y_j) = min_i (c(x_i, y_j) - f_i)

    For c(x, y) = ||x - y||²/2:
        g(y) = inf_x (||x - y||²/2 - f(x))

    This is the Legendre-Fenchel transform of f under the cost c.

    Args:
        f: Potential values at source points (n,)
        source: Source positions [x_1, ..., x_n]
        target: Target positions [y_1, ..., y_m]

    Returns:
        c-transform g at target points (m,)
    """
    n = len(source)
    m = len(target)
    g = np.full(m, -np.inf)

    for j in range(m):
        yj = target[j]
        values = np.array([
            float(np.dot(xi - yj, xi - yj)) - f[i]
            for i, xi in enumerate(source)
        ])
        g[j] = np.min(values)

    return g


def c_conjugate(f: np.ndarray,
                points: list[Vector],
                cost_fn: Optional[Callable[[Vector, Vector], float]] = None) -> Callable:
    """
    Create a callable c-conjugate function.

    Returns a function g(y) = inf_x (c(x, y) - f(x))
    that can be evaluated at arbitrary points.

    Args:
        f: Potential values at source points
        points: Source points
        cost_fn: Cost function (default squared Euclidean / 2)

    Returns:
        Callable g: ℝ^d → ℝ
    """
    cost_fn = cost_fn or (lambda x, y: 0.5 * float(np.dot(x - y, x - y)))

    def g(y: Vector) -> float:
        y = np.asarray(y, dtype=float)
        values = np.array([
            cost_fn(xi, y) - f[i]
            for i, xi in enumerate(points)
        ])
        return float(np.min(values))

    return g


def kantorovich_dual(positions1: list[Vector],
                     positions2: list[Vector],
                     max_iter: int = 100,
                     lr: float = 0.1) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Solve the Kantorovich dual problem numerically.

    maximize_{f, g} ⟨f, a⟩ + ⟨g, b⟩
    subject to f_i + g_j ≤ C_{ij}

    This is the dual of the optimal transport problem.
    At optimality, the duality gap is zero.

    Args:
        positions1: Source positions
        positions2: Target positions
        max_iter: Maximum iterations
        lr: Learning rate for gradient ascent

    Returns:
        (f, g, dual_value) optimal dual potentials and dual value
    """
    n = len(positions1)
    m = len(positions2)

    # Cost matrix
    C = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            diff = np.asarray(positions1[i]) - np.asarray(positions2[j])
            C[i, j] = float(np.dot(diff, diff))

    # Initialize potentials
    f = np.zeros(n)
    g = np.zeros(m)

    for iteration in range(max_iter):
        # Constraint satisfaction: f_i + g_j ≤ C_{ij}
        # Violation: when f_i + g_j > C_{ij}, we need to decrease
        violation = f[:, np.newaxis] + g[np.newaxis, :] - C
        violation_mask = violation > 0

        # Gradient: d(dual)/df_i = a_i - Σ_j π_{ij}
        # where π is the induced plan
        if np.any(violation_mask):
            # Penalty approach: push potentials down where violated
            grad_f = np.zeros(n)
            grad_g = np.zeros(m)

            for i in range(n):
                grad_f[i] = 1.0 / n - np.mean(violation_mask[i, :])

            for j in range(m):
                grad_g[j] = 1.0 / m - np.mean(violation_mask[:, j])

            f += lr * grad_f
            g += lr * grad_g

            # Project back to constraint satisfaction
            sat = f[:, np.newaxis] + g[np.newaxis, :] - C
            # f = f - max(0, sat) / 2 (heuristic adjustment)
            overshoot = np.maximum(sat, 0)
            f -= 0.5 * lr * np.mean(overshoot, axis=1)
            g -= 0.5 * lr * np.mean(overshoot, axis=0)

    dual_value = np.mean(f) + np.mean(g)
    return f, g, float(dual_value)


def monge_map_from_potential(f: np.ndarray,
                              source: list[Vector],
                              target: list[Vector]) -> list[Vector]:
    """
    Extract the Monge map from a dual potential f.

    The Monge map is: T(x_i) = argmin_y (c(x_i, y) - f_i)

    For squared Euclidean cost:
        T(x_i) = x_i - ∇f_i / 2  (when f is c-concave)

    Args:
        f: Dual potential at source points
        source: Source positions
        target: Target positions

    Returns:
        Monge map images of source points
    """
    n = len(source)
    m = len(target)
    mapped = []

    for i in range(n):
        # Find best target
        best_j = 0
        best_cost = float('inf')
        for j in range(m):
            diff = source[i] - target[j]
            cost = float(np.dot(diff, diff)) - f[i]
            if cost < best_cost:
                best_cost = cost
                best_j = j
        mapped.append(target[best_j].copy())

    return mapped


def c_concave_check(f: np.ndarray,
                     points: list[Vector],
                     tol: float = 1e-9) -> bool:
    """
    Check if f is c-concave.

    A function f is c-concave if (f^c)^c = f,
    where f^c(y) = inf_x (c(x, y) - f(x)).

    For the squared distance cost, c-concave = concave up to quadratic:
        f(x) = ||x||²/2 - g(x) where g is convex.

    Args:
        f: Potential values at points
        points: Grid points where f is defined
        tol: Numerical tolerance

    Returns:
        True if f is c-concave
    """
    # Compute c-transform twice
    g = c_transform(f, points, points)
    f_double = c_transform(g, points, points)

    return np.allclose(f, f_double, atol=tol)


def kantorovich_rkhs(f: np.ndarray,
                     points: list[Vector],
                     kernel_scale: float = 1.0) -> Callable:
    """
    Extend the dual potential f to a function on ℝ^d.

    Uses a kernel approximation:
        f_hat(x) = Σ_i k(x, x_i) · α_i

    where k is a Gaussian kernel and α are chosen to match f at sample points.

    Args:
        f: Potential values at sample points
        points: Sample point locations
        kernel_scale: Gaussian kernel bandwidth

    Returns:
        Callable function f_hat: ℝ^d → ℝ
    """
    n = len(points)
    points_arr = np.array([np.asarray(p, dtype=float) for p in points])

    # Build kernel matrix K_{ij} = k(x_i, x_j)
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            diff = points_arr[i] - points_arr[j]
            K[i, j] = np.exp(-np.dot(diff, diff) / (2 * kernel_scale ** 2))

    # Solve K · α = f
    alpha = np.linalg.solve(K + 1e-10 * np.eye(n), f)

    def f_hat(x: Vector) -> float:
        x = np.asarray(x, dtype=float)
        vals = np.array([
            np.exp(-np.dot(x - points_arr[i], x - points_arr[i])
                   / (2 * kernel_scale ** 2))
            for i in range(n)
        ])
        return float(np.dot(vals, alpha))

    return f_hat


def transport_plan_from_potentials(f: np.ndarray,
                                    g: np.ndarray,
                                    C: np.ndarray,
                                    reg: float = 0.0) -> np.ndarray:
    """
    Build coupling matrix from dual potentials.

    At optimality:
        π_{ij} = a_i · b_j · exp((f_i + g_j - C_{ij}) / ε)

    For ε = 0 (exact OT):
        π_{ij} > 0 only if f_i + g_j = C_{ij}

    Args:
        f: Source dual potentials (n,)
        g: Target dual potentials (m,)
        C: Cost matrix (n, m)
        reg: Entropic regularization (default 0 = exact)

    Returns:
        Coupling matrix π (n, m)
    """
    n, m = C.shape
    a = np.ones(n) / n
    b = np.ones(m) / m

    if reg > 0:
        # Entropic case
        pi = a[:, np.newaxis] * b[np.newaxis, :] * \
             np.exp((f[:, np.newaxis] + g[np.newaxis, :] - C) / reg)
        return pi
    else:
        # Exact case: support on C - f - g = 0
        slack = C - f[:, np.newaxis] - g[np.newaxis, :]
        pi = np.zeros((n, m))
        for i in range(n):
            # Find active constraints (within tolerance)
            active = np.where(np.abs(slack[i, :]) < 1e-9)[0]
            if len(active) > 0:
                pi[i, active] = a[i] / len(active)
        return pi


# ─── Monge-Ampère equation (fleet consensus geometry) ───────────────────────

def monge_ampere_operator(positions: list[Vector],
                           potential: np.ndarray,
                           epsilon: float = 1e-6) -> float:
    """
    Discrete Monge-Ampère operator.

    The continuous Monge-Ampère equation:
        det(∇²φ) = ρ(∇φ) / σ

    describes the optimal transport map ∇φ between densities ρ and σ.
    In the fleet context, φ is the Kantorovich potential, and
    the map ∇φ gives the consensus update.

    The discrete operator approximates:
        (Δ_h φ)(x_i) ≈ (2/h²) · (φ_{i+1} + φ_{i-1} - 2φ_i)

    Args:
        positions: Agent positions
        potential: Potential values at positions
        epsilon: Smoothing parameter

    Returns:
        Monge-Ampère residual (0 = perfect transport map)
    """
    n = len(positions)
    pos_arr = np.array([np.asarray(p, dtype=float) for p in positions])

    # Compute all pairwise distances
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            diff = pos_arr[i] - pos_arr[j]
            D[i, j] = np.sqrt(np.dot(diff, diff) + epsilon)

    # Discrete Laplacian via weighted graph
    laplacian = np.zeros(n)
    for i in range(n):
        weights = np.exp(-D[i, :] / epsilon)
        weights[i] = 0
        weights /= np.sum(weights) + epsilon
        laplacian[i] = np.sum(weights * (potential[i] - potential))

    return float(np.sqrt(np.mean(laplacian ** 2)))


# ─── Tests ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import math

    print("=== Transport Geometry — Dual Potentials and c-Transform ===\n")

    # Test 1: c-transform
    print("Test 1: c-transform")
    source = [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([0.5, 1.0])]
    target = [np.array([0.1, 0.1]), np.array([0.9, 0.0]), np.array([0.4, 0.9])]

    f_test = np.array([0.0, 0.1, 0.2])
    g_test = c_transform(f_test, source, target)
    print(f"  f: {f_test}")
    print(f"  c-transform g: {g_test}")

    # Test 2: c-concave check
    print("\nTest 2: c-concave check")
    # A linear function should be c-concave for squared cost
    f_lin = np.array([0.0, 0.5, 1.0])
    is_cc = c_concave_check(f_lin, source)
    print(f"  Linear f is c-concave: {is_cc}")

    # Test 3: Monge map from potential
    print("\nTest 3: Monge map from potential")
    mapped = monge_map_from_potential(f_test, source, target)
    for i, m in enumerate(mapped):
        print(f"  {source[i]} → {m}")

    # Test 4: Kantorovich dual
    print("\nTest 4: Kantorovich dual")
    f_dual, g_dual, value = kantorovich_dual(source, target)
    print(f"  f: {f_dual}")
    print(f"  g: {g_dual}")
    print(f"  Dual value: {value:.6f}")

    # Test 5: RKHS extension
    print("\nTest 5: RKHS extension of potential")
    f_hat = kantorovich_rkhs(f_test, source)
    for i, s in enumerate(source):
        print(f"  f({s}) = {f_test[i]:.4f}, f_hat({s}) = {f_hat(s):.4f}")
