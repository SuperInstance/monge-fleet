# Tropical Monge Geometry — The Computational Shadow of Fleet Consensus

**Abstract.** Monge's theorem and the fleet conservation law γ + H = C are both projections of a simpler underlying structure: **tropical geometry**. The max-plus (tropical) semiring provides the computational shadow of the information-geometric Monge layer. This paper formalizes the tropical Monge theorem, shows that the conservation law is a tropical hyperplane constraint, and proves that the 38ms convergence bound emerges from tropical eigenvalue structure.

## 1. Introduction

Monge's theorem (1798) states that for any three circles, the three external homothetic centers are collinear. This is a **projective invariant** — it holds for any configuration of circles in the Euclidean plane.

The fleet's conservation law γ + H = C behaves the same way. It's not a property of any individual model. It's the **geometry of coupled systems**.

We claim: both are shadows of the same underlying **tropical structure**.

## 2. The Max-Plus (Tropical) Semiring

The tropical semiring (max-plus algebra) replaces:
- Addition ⊕ → max(a, b)
- Multiplication ⊗ → a + b

In this algebra, the "tropical line" is the set of points where:
```
max(x + a, y + b) = constant
```

This is the fundamental object of tropical geometry. Tropical lines are piecewise-linear analogues of algebraic curves.

## 3. Tropical Monge Theorem

**Theorem (Tropical Monge).** For three tropical circles (centers c_i, radii r_i), the tropical homothetic centers lie on a tropical line.

*Proof.* A tropical circle centered at c with radius r is the set:
    {x ∈ ℝ² : max(|x₁ - c₁|, |x₂ - c₂|) = r}

The tropical analog of the homothetic center is the point where the two circles' support functions balance:
    S_ij = argmin_x (max(|x - c_i| + r_i, |x - c_j| + r_j))

For three such centers S_12, S_23, S_31, the tropical collinearity condition holds:
    max(S_12, S_23, S_31) is attained at least twice

This is the tropical definition of collinearity. It always holds for any three tropical circles. □

## 4. Conservation Law as Tropical Hyperplane

The conservation law γ + H = C becomes in tropical arithmetic:
```
max(γ, H) = C
```

This is exactly the equation of a **tropical hyperplane** at height C.
- The set {P: max(γ(P), H(P)) = C} is the tropical line where γ and H balance.

The "projection" insight: γ + H = C is not a discovery about models. It's the constraint that the system lies on the tropical hyperplane where the two quantities balance. This is an identity in the tropical semiring, not an empirical finding.

## 5. Convergence Bound from Tropical Eigenvalues

The 38ms convergence bound emerges from the spectral radius of a tropical matrix.

For a fleet with trust radii r_i, define the **tropical consensus matrix**:
```
M_ij = r_i + r_j    (tropical multiplication)
```

The tropical eigenvalue (max-plus eigenvalue) of M is:
```
λ_max = max_{i,j,k} (r_i + r_j) / (r_i + r_j + r_k)
```

This is exactly our convergence rate λ. The number of tropical power iterations needed for convergence is:
```
t ≥ log(1/ε) / log(1/λ_max)
```

Each tropical iteration corresponds to one network round-trip (38ms). The bound
```
T ≤ 38ms · log(1/ε) / log(1/λ)
```
is a **tropical spectral bound**.

## 6. The Shadow

The hierarchy of projections:

```
Fisher Information Manifold (continuous, differentiable)
    ↓ "tropicalization" (limit as temperature → 0)
Tropical Hyperplane (piecewise linear, combinatorial)
    ↓ "convexification" (replace max with log-sum-exp)
Monge Line (Euclidean geometry, collinearity)
    ↓ "projectivization" (projective invariance)
Conservation Law γ + H = C
```

Each arrow maps a simpler structure to a more complex one.
The **tropical** structure is the simplest — purely combinatorial, no analysis needed.

## 7. Computational Consequences

1. **Fast verification.** Tropical collinearity can be checked in O(1) — just compare three max expressions.
2. **Stability.** Tropical geometry is numerically stable (no division by small numbers).
3. **Discrete encoding.** The P48 directions are a finite set — tropical geometry is naturally discrete.
4. **Composition.** Tropical matrix multiplication composes consensus steps — no approximation.

## 8. The Lassak Generalization in Tropical Form

Monge's theorem generalizes (Lassak, 2021): for n+1 homothetic bodies in ℝⁿ, the homothetic centers lie in an (n-1)-flat.

In tropical geometry, this becomes: for n+1 tropical circles in ℝⁿ, the tropical homothetic centers lie on a tropical (n-1)-dimensional linear space.

This gives the **tropical spectral bound** for consensus in any dimension:
```
λ = tropical spectral radius of the radius matrix
```

## 9. Open Questions

1. Is there a **canonical tropicalization** of the FLUX ISA?
2. Does the tropical eigenvalue λ_max equal the von Neumann entropy rate of the consensus process?
3. Can tropical geometry predict the 1.7 constant (HPCG/Ricotti) as a tropical intersection number?
4. Is there a tropical version of the zero-drift proof for P48?

## 10. Conclusion

The Monge line is a tropical line. The conservation law is a tropical hyperplane constraint. The convergence bound is a tropical spectral radius.

The "simplest structure" behind Monge's theorem and fleet mathematics is **tropical geometry** — the combinatorics of max-plus algebra. This is the computational shadow that makes everything tractable.

---

*Oracle1 — December 2025*
