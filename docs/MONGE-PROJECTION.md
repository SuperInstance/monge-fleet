# The Monge Projection — What Lies Behind the Conservation Law

## The Core Insight

Casey said: "We keep finding these patterns and treating them as discoveries. 
But they're not discoveries. They're PROJECTIONS."

Monge's theorem is a projection of affine geometry.
The conservation law γ + H = C is a projection of information geometry.

What is being projected?

## The Simpler Structure

### Option 1: Information Geometry (Fisher Information)

Coupled models minimize KL divergence:
  KL(p∥q) = ∫ p(x) log(p(x)/q(x)) dx

The "similarity center" between two models p and q is the distribution
that lies at the geodesic midpoint in the Fisher information metric.

For three models p, q, r, the geodesic triangle has a collinearity property
analogous to Monge's theorem in the Fisher information manifold.

### Option 2: Tropical Geometry (Max-Plus Algebra)

The "shadow" of optimization is tropical arithmetic:
  a ⊕ b = max(a, b)
  a ⊗ b = a + b

The Monge line is the tropical hyperplane — the set of points where
max operations balance.

### Option 3: Optimal Transport (Wasserstein Geometry)

The Monge map: push-forward that minimizes transport cost.
For three measures μ, ν, ρ, the Monge maps between them satisfy
a collinearity condition in the Wasserstein space.

### Option 4: Categorical Limits (Colimits)

The "similarity center" is the pushout of two objects in the category
of models. Monge's theorem = colimit of 3 objects is a line in the 
projective completion of the category.

## Which is the Simplest?

The simplest structure is likely **tropical geometry**:

1. Models = points in ℝⁿ (parameter space)
2. Coupling = tropical sum (max)
3. Similarity center = tropical midpoint
4. Collinearity = tropical line equation

The Monge line is the set of points satisfying:
  max(x·a + r₁, x·b + r₂) = constant  (for circles centered at a, b with radii r₁, r₂)

For three circles, the three centers where these maxima balance are always collinear.
This is the tropical version of Monge's theorem.

## The Conservation Law as Tropical Projection

γ + H = C

In tropical arithmetic:
  γ ⊕ H = C
  
This is the condition that γ and H balance to give C.
The "projection" is the tropical line where this balance holds.

The eigenvalue concentration isn't a surprise — it's what happens
when you constrain a system to the tropical hyperplane.

## Next Steps

1. Formalize the tropical Monge theorem for fleet constraints
2. Show that γ + H = C is equivalent to tropical collinearity
3. Prove convergence via tropical geometry (max-plus semiring)
4. Connect to Lassak generalization via tropicalization
