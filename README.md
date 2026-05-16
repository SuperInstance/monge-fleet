# monge-fleet

**Monge's theorem for fleet mathematics.**

A library implementing the geometric enhancement layer for the fleet's constraint theory mathematics. Based on the insight that Monge's theorem — three circles always have collinear external homothetic centers — is a **projective invariant** that describes the fundamental geometry of coupled logical systems.

## The Monge Insight

Monge's theorem says: for any 3 circles, the 3 pairwise external homothetic centers are **always collinear**. Not sometimes. Always.

Our conservation law γ + H = C works the same way. It's not a property of any individual model. It's the **geometry of coupled systems**. Of course it holds. It's a projection.

If these are projections, there's a simpler structure behind them. Monge's theorem follows from affine transformation properties. Our conservation law should follow from **information geometry** — the structure of coupled optimizers.

## Key Modules

| Module | Purpose |
|--------|---------|
| `homothetic.py` | External/internal homothetic centers, Monge collinearity checks |
| `radical_axis.py` | Radical axis as 1-cocycle, coboundary relations |
| `cohomology.py` | H¹ computation via Menelaus theorem, β₁ from homothetic centers |
| `consensus.py` | MongeConsensus — Zero Holonomy Consensus with area-based Lyapunov |
| `pythagorean48.py` | P48 Monge verification, zero-drift proof via collinearity |
| `geometry.py` | 3D lifting, sandwich planes, Monge flat computation (Lassak gen.) |

## Quick Start

```python
from monge_fleet import MongeConsensus, P48Monge, verify_monge_p48

# Zero Holonomy Consensus with Monge monitoring
mc = MongeConsensus(3)
mc.set_agent(0, [0.0, 0.0], 1.0)
mc.set_agent(1, [4.0, 0.0], 2.0)
mc.set_agent(2, [2.0, 3.0], 1.5)

print(f"Holonomy area: {mc.max_area():.6f}")
print(f"Convergence: {mc.convergence_time():.1f}ms")
print(f"With Byzantine f=1: {mc.byzantine_convergence_time(1):.1f}ms")

# Pythagorean48 Monge verification
result = verify_monge_p48()
print(f"P48 triples verified: {result['n_triples_checked']:,}")
print(f"Zero drift: {result['zero_drift_verified']}")
```

## Core Ideas

### Homothetic Centers as Consensus Fixed Points

For agents with trust radii, the homothetic center S_ij is the fixed point of consensus between agents i and j. Monge says S_ij, S_jk, S_ki are always collinear. Deviation from collinearity = holonomy = consensus failure.

### Radical Axis as 1-Cocycle

The radical axis Rad(Cᵢ, Cⱼ) = {P: Power(P,Cᵢ) = Power(P,Cⱼ)} is a cohomological object. The coboundary δφ(i,j,k) = φᵢⱼ + φⱼₖ + φₖᵢ = 0 for all triples — equivalent to Monge collinearity.

### Area as Lyapunov Function

H_ijk = area(S_ij, S_jk, S_ki). The consensus converges when H_ijk → 0. Convergence rate λ = max(r_i + r_j)/(r_i + r_j + r_k).

### The Monge Layer

```
Fleet Applications
       ↑
  Zero Holonomy Consensus (Monge-collinear updates, 38ms)
       ↑
  Monge Consistency Layer — radical axis cohomology, collinearity monitor
       ↑
  H1 Cohomology (β₁) Engine — Laman rigidity, emergence detection
       ↑
  Pythagorean48 / jc1-ct-bridge — constraint satisfaction, direction encoding
```

## Mathematical Foundation

**Monge's Theorem** (Gaspard Monge, 1798): For any three circles, the three external homothetic centers are collinear.

**Lassak Generalization** (2021): For n+1 homothetic bodies in ℝⁿ, the homothetic centers lie in an (n-1)-dimensional affine subspace.

**Menelaus Theorem**: Used for O(V³) β₁ computation via ratio checks on homothetic center lines.

## References

- Monge, G. (1798). *Géométrie descriptive*
- Wells, D. (1991). *The Penguin Dictionary of Curious and Interesting Geometry*
- Lassak, M. (2021). "A note on some generalizations of Monge's theorem" arXiv:2104.06343
- Laman, G. (1970). "On graphs and rigidity of plane skeletal structures"
- Fleet-internal: FLUX ISA v2.4, Zero Holonomy Whitepaper, Pythagorean48 spec

## Installation

```bash
pip install -e .
```

## Tests

```bash
pytest tests/ -v
```
