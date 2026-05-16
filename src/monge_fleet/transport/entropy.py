"""
Entropic regularization for optimal transport in fleet consensus.

The entropic optimal transport problem:
    min_π ⟨π, C⟩ + ε · H(π)
    subject to π · 1 = a, π^T · 1 = b, π ≥ 0

where H(π) = Σ_{ij} π_{ij} · log π_{ij} is the negative entropy.

The entropic regularization makes the transport plan smoother and
enables the Sinkhorn algorithm (a simple matrix scaling iteration).

For fleet consensus, ε controls the "spread" of the coupling:
- ε → 0: exact Monge map (hard assignment, potentially unstable)
- ε → ∞: independent coupling π = a · b^T (no transport information)
- Moderate ε: smooth transport plan with well-posed gradient

Key result: The entropic transport plan is a Gibbs kernel:
    π_{ij} = u_i · exp(-C_{ij}/ε) · v_j
where u, v are found by Sinkhorn scaling.

The connection to Monge:
    The entropic regularized transport plan converges to the Monge map
    as ε → 0. The Monge line (collinearity) gives the geometric structure
    of the limiting transport plan.
"""

import numpy as np
from typing import Optional, Callable, Tuple

Vector = np.ndarray


def sinkhorn_knopp(a: np.ndarray,
                   b: np.ndarray,
                   C: np.ndarray,
                   reg: float = 0.01,
                   max_iter: int = 1000,
                   tol: float = 1e-9,
                   log: bool = False) -> Tuple[np.ndarray, dict]:
    """
    Entropic Optimal Transport via Sinkhorn-Knopp algorithm.

    Solves:
        min_π ⟨π, C⟩ + ε · KL(π || a·b^T)
    subject to π·1 = a, π^T·1 = b, π ≥ 0

    The solution is: π_{ij} = u_i · K_{ij} · v_j
    where K_{ij} = exp(-C_{ij}/ε) and u, v are the Sinkhorn scaling vectors.

    The Sinkhorn-Knopp iterations:
        u^{(k+1)} = a / (K · v^{(k)})
        v^{(k+1)} = b / (K^T · u^{(k+1)})

    Args:
        a: Source marginal (n,)
        b: Target marginal (m,)
        C: Cost matrix (n, m)
        reg: Entropic regularization coefficient ε
        max_iter: Maximum iterations
        tol: Convergence tolerance on marginals
        log: If True, compute dual potentials (log-domain stabilization)

    Returns:
        (π, info) where π is the coupling matrix and info contains metadata
    """
    n, m = C.shape

    if log:
        return _sinkhorn_log(a, b, C, reg, max_iter, tol)

    # Gibbs kernel K
    K = np.exp(-C / reg)

    # Initialize
    u = np.ones(n) / n
    v = np.ones(m) / m

    converged = False
    iterations = 0

    for iteration in range(max_iter):
        u_prev = u.copy()
        v_prev = v.copy()

        # Update
        u = a / (K @ v + 1e-15)
        v = b / (K.T @ u + 1e-15)

        # Check convergence
        marginal_u = u * (K @ v)
        marginal_v = v * (K.T @ u)

        err_u = np.max(np.abs(marginal_u - a))
        err_v = np.max(np.abs(marginal_v - b))

        if err_u < tol and err_v < tol:
            converged = True
            iterations = iteration
            break

    # Build coupling matrix
    pi = u[:, np.newaxis] * K * v[np.newaxis, :]

    # Compute dual potentials f_i = ε · log u_i, g_j = ε · log v_j
    f = reg * np.log(np.maximum(u, 1e-15))
    g = reg * np.log(np.maximum(v, 1e-15))

    info = {
        'converged': converged,
        'iterations': iterations,
        'reg': reg,
        'f_potential': f,
        'g_potential': g,
        'total_cost': float(np.sum(pi * C)),
        'entropy': float(np.sum(pi * np.log(np.maximum(pi, 1e-15)))),
    }

    return pi, info


def _sinkhorn_log(a: np.ndarray,
                  b: np.ndarray,
                  C: np.ndarray,
                  reg: float = 0.01,
                  max_iter: int = 1000,
                  tol: float = 1e-9) -> Tuple[np.ndarray, dict]:
    """
    Log-domain Sinkhorn for numerical stability with very small ε.

    Works with f, g dual potentials instead of u, v:
        f^{(k+1)} = -reg · log(a / (exp(f^{(k)}/reg) · K))
        g^{(k+1)} = -reg · log(b / (exp(g^{(k+1)}/reg)^T · K))

    More stable for small ε but slower (requires exp/log each iteration).
    """
    n, m = C.shape
    K = np.exp(-C / reg)

    f = np.zeros(n)
    g = np.zeros(m)

    converged = False
    iterations = 0

    for iteration in range(max_iter):
        f_prev = f.copy()
        g_prev = g.copy()

        # Log-domain update
        # u = exp(f/reg), so K^T · u = K^T · exp(f/reg)
        # g = -reg · log(b / (K^T · exp(f/reg)))
        KTu = K.T @ np.exp(f / reg)
        g = -reg * np.log(np.maximum(b, 1e-15) / np.maximum(KTu, 1e-15))

        Kv = K @ np.exp(g / reg)
        f = -reg * np.log(np.maximum(a, 1e-15) / np.maximum(Kv, 1e-15))

        # Check convergence on potentials
        df = np.max(np.abs(f - f_prev))
        dg = np.max(np.abs(g - g_prev))

        if df < tol and dg < tol:
            converged = True
            iterations = iteration
            break

    # Build coupling from potentials
    u = np.exp(f / reg)
    v = np.exp(g / reg)
    pi = u[:, np.newaxis] * K * v[np.newaxis, :]

    info = {
        'converged': converged,
        'iterations': iterations,
        'reg': reg,
        'f_potential': f,
        'g_potential': g,
        'total_cost': float(np.sum(pi * C)),
        'entropy': float(np.sum(pi * np.log(np.maximum(pi, 1e-15)))),
    }

    return pi, info


class EntropicTransport:
    """
    Entropy-regularized optimal transport for fleet consensus.

    The entropic regularization parameter ε controls the tradeoff:
    - Small ε: nearly-exact transport (close to Monge map)
      → Better for precise cost computation
      → Slower Sinkhorn convergence
    - Large ε: smoother transport (more diffusive)
      → Faster convergence
      → Less accurate for sharp distributions

    In the fleet context:
    - Agents with very different trust radii → need small ε
      (sharp coupling required)
    - Agents with similar trust radii → can use larger ε
      (coupling is naturally smoother)
    """

    def __init__(self,
                 source: list[Vector],
                 target: list[Vector],
                 source_masses: Optional[list[float]] = None,
                 target_masses: Optional[list[float]] = None,
                 cost_fn: Optional[Callable[[Vector, Vector], float]] = None):
        """
        Args:
            source: Source positions
            target: Target positions
            source_masses: Source marginals (default uniform)
            target_masses: Target marginals (default uniform)
            cost_fn: Cost function (default squared Euclidean)
        """
        self.n = len(source)
        self.m = len(target)
        self.source = [np.asarray(s, dtype=float) for s in source]
        self.target = [np.asarray(t, dtype=float) for t in target]
        self.source_masses = np.ones(self.n) / self.n if source_masses is None else np.array(source_masses)
        self.target_masses = np.ones(self.m) / self.m if target_masses is None else np.array(target_masses)

        self._cost_fn = cost_fn or (lambda x, y: float(np.dot(x - y, x - y)))

        # Build cost matrix
        self.C = np.zeros((self.n, self.m))
        for i in range(self.n):
            for j in range(self.m):
                self.C[i, j] = self._cost_fn(self.source[i], self.target[j])

        # Cache
        self._last_pi: Optional[np.ndarray] = None
        self._last_reg: Optional[float] = None
        self._last_info: Optional[dict] = None

    def solve(self, reg: float = 0.01, max_iter: int = 1000, tol: float = 1e-9) -> Tuple[np.ndarray, dict]:
        """
        Solve entropic optimal transport.

        Args:
            reg: Entropic regularization coefficient
            max_iter: Maximum Sinkhorn iterations
            tol: Marginal convergence tolerance

        Returns:
            (π, info) coupling matrix and metadata
        """
        pi, info = sinkhorn_knopp(
            self.source_masses,
            self.target_masses,
            self.C,
            reg=reg,
            max_iter=max_iter,
            tol=tol,
        )
        self._last_pi = pi
        self._last_reg = reg
        self._last_info = info
        return pi, info

    def dual_potentials(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the dual optimal potentials (f, g).

        The dual problem is:
            max_{f, g} ⟨f, a⟩ + ⟨g, b⟩ - ε · ⟨exp((f⊕g - C)/ε), 1⟩

        The optimal f, g satisfy:
            π_{ij} = exp((f_i + g_j - C_{ij}) / ε)

        Returns:
            (f, g) dual potential vectors
        """
        if self._last_info is None:
            self.solve()

        assert self._last_info is not None
        return self._last_info['f_potential'], self._last_info['g_potential']

    def primal_value(self, reg: Optional[float] = None) -> float:
        """
        Compute the primal objective value:
            ⟨π, C⟩ + ε · Σ_{ij} π_{ij} · log π_{ij}

        At optimality, this equals the dual objective value.

        Args:
            reg: Regularization (uses last solve's reg if None)

        Returns:
            Primal value
        """
        if self._last_pi is None:
            self.solve(reg=reg or 0.01)

        assert self._last_pi is not None
        pi = self._last_pi
        r = self._last_reg if reg is None else reg

        cost = float(np.sum(pi * self.C))
        entropy_parts = pi * np.log(np.maximum(pi, 1e-15))
        entropy = float(np.sum(entropy_parts))

        return cost + r * entropy

    def entropic_wasserstein(self, p: int = 2) -> float:
        """
        Compute entropic Wasserstein distance.

        The entropic W₂ distance is:
            W₂²(μ, ν) ≈ Σ_{ij} π_{ij} · ||x_i - y_j||²

        where π is the entropic coupling.

        Args:
            p: Order (1 or 2)

        Returns:
            Entropic p-Wasserstein distance
        """
        if self._last_pi is None:
            self.solve()

        assert self._last_pi is not None
        pi = self._last_pi

        if p == 2:
            W_sq = float(np.sum(pi * self.C))
            return float(np.sqrt(W_sq))
        elif p == 1:
            # Use sqrt(C) for L1 = expected Euclidean distance
            C_l1 = np.sqrt(np.maximum(self.C, 1e-15))
            W1 = float(np.sum(self._last_pi * C_l1))
            return W1
        else:
            raise ValueError(f"Unsupported p={p}")

    def gibbs_kernel(self, reg: float = 0.01) -> np.ndarray:
        """
        Compute the Gibbs kernel K_{ij} = exp(-C_{ij} / ε).

        This is the kernel used by Sinkhorn scaling.
        """
        return np.exp(-self.C / reg)

    def entropic_map(self) -> list[np.ndarray]:
        """
        Compute the entropic Monge map (barycentric projection).

        The entropic map is the barycentric projection of the coupling:
            T(x_i) = (Σ_j π_{ij} · y_j) / (Σ_j π_{ij})

        This is a smoothed approximation to the exact Monge map.

        Returns:
            Target positions under the entropic map
        """
        if self._last_pi is None:
            self.solve()

        assert self._last_pi is not None
        pi = self._last_pi

        mapped = []
        for i in range(self.n):
            total = np.sum(pi[i, :])
            if total > 1e-15:
                weighted = sum(pi[i, j] * self.target[j] for j in range(self.m))
                mapped.append(weighted / total)
            else:
                mapped.append(self.source[i].copy())

        return mapped

    def entropic_planar(self, n_grid: int = 50) -> np.ndarray:
        """
        Compute the entropic transport plan as a 2D density for visualization.

        Only works for 1D source/target spaces embedded in the first coordinate.

        Args:
            n_grid: Number of grid points per dimension

        Returns:
            n_grid × n_grid density matrix
        """
        if self._last_pi is None:
            self.solve()

        assert self._last_pi is not None
        plan = self._last_pi

        # Extract 1D coordinates
        src_1d = np.array([s[0] for s in self.source])
        tgt_1d = np.array([t[0] for t in self.target])

        # Create grid
        lo = min(src_1d.min(), tgt_1d.min())
        hi = max(src_1d.max(), tgt_1d.max())
        grid = np.linspace(lo, hi, n_grid)

        density = np.zeros((n_grid, n_grid))
        for i in range(self.n):
            for j in range(self.m):
                ix = np.argmin(np.abs(grid - src_1d[i]))
                jy = np.argmin(np.abs(grid - tgt_1d[j]))
                density[ix, jy] += plan[i, j]

        return density

    def entropy_rate(self) -> float:
        """
        Entropy rate of the transport plan.

        H(π) = -Σ_{ij} π_{ij} · log π_{ij}

        This measures the "fuzziness" of the coupling.
        Maximum at independence: H(π) = H(a) + H(b)
        Minimum at deterministic coupling: H(π) = 0
        """
        if self._last_pi is None:
            self.solve()

        assert self._last_pi is not None
        pi = self._last_pi
        return -float(np.sum(pi * np.log(np.maximum(pi, 1e-15))))


def epsilon_schedule(step: int, initial: float = 0.1, min_val: float = 1e-4, decay: float = 0.9) -> float:
    """
    Annealing schedule for entropic regularization.

    Start with large ε (fast, smooth) and decrease to small ε (accurate).

    Args:
        step: Current iteration
        initial: Initial ε
        min_val: Minimum ε
        decay: Per-step decay factor

    Returns:
        ε at given step
    """
    return max(initial * (decay ** step), min_val)


# ─── Tests ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import math

    print("=== Entropic Transport for Fleet Consensus ===\n")

    # Test 1: Entropic transport between two agent configurations
    print("Test 1: Entropic transport")
    source = [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([0.5, 1.0])]
    target = [np.array([0.1, 0.1]), np.array([0.9, 0.0]), np.array([0.4, 0.9])]

    et = EntropicTransport(source, target)
    pi, info = et.solve(reg=0.05)
    print(f"  Converged: {info['converged']}")
    print(f"  Iterations: {info['iterations']}")
    print(f"  Total cost: {info['total_cost']:.6f}")
    print(f"  Entropy: {info['entropy']:.6f}")
    print(f"  Row marginals: {np.sum(pi, axis=1)}")
    print(f"  Col marginals: {np.sum(pi, axis=0)}")

    # Test 2: Entropic map (smooth Monge map)
    print("\nTest 2: Entropic Monge map")
    mapped = et.entropic_map()
    for i, m in enumerate(mapped):
        print(f"  Agent {i}: {source[i]} → {m}")

    # Test 3: Dual potentials
    print("\nTest 3: Dual potentials")
    f, g = et.dual_potentials()
    print(f"  f: {f}")
    print(f"  g: {g}")

    # Test 4: ε-schedule for annealing
    print("\nTest 4: Epsilon annealing")
    for step in range(10):
        eps = epsilon_schedule(step)
        print(f"  step {step}: ε = {eps:.6f}")

    # Test 5: Log-domain Sinkhorn for stability
    print("\nTest 5: Log-domain Sinkhorn (small ε)")
    pi_log, info_log = sinkhorn_knopp(
        np.ones(3)/3, np.ones(3)/3, et.C,
        reg=1e-3, max_iter=2000, tol=1e-9, log=True
    )
    print(f"  Converged: {info_log['converged']}")
    print(f"  Iterations: {info_log['iterations']}")
    print(f"  Total cost: {info_log['total_cost']:.6f}")
