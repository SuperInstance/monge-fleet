#!/usr/bin/env python3
"""
Demo: Monge's theorem for fleet mathematics.

Run: python3 demo.py

This demo shows:
1. Monge's theorem for three arbitrary circles
2. Zero Holonomy Consensus with area-based Lyapunov
3. Pythagorean48 zero drift verification
4. Emergence detection via Monge deviation
"""

import math
import numpy as np

from monge_fleet import (
    ExternalCenter, RadicalCochain, MongeCohomology,
    MongeConsensus, P48Monge, verify_monge_p48,
    lift_to_sphere, sandwich_planes, monge_flat
)
from monge_fleet.transport import ConsensusTransport, fleet_consensus
from monge_fleet.bridge import MongeCTBridge


def demo_mongeTheorem():
    print("=" * 60)
    print("Demo 1: Monge's Theorem for Three Arbitrary Circles")
    print("=" * 60)
    
    circles = [
        (np.array([0.0, 0.0]), 1.0),
        (np.array([4.0, 0.0]), 2.5),
        (np.array([1.5, 3.0]), 1.8),
    ]
    
    ec = ExternalCenter(circles)
    
    # Get the three external homothetic centers
    S01 = ec.get(0, 1)
    S12 = ec.get(1, 2)
    S20 = ec.get(2, 0)
    
    print(f"\nCircle 0: center={circles[0][0]}, radius={circles[0][1]}")
    print(f"Circle 1: center={circles[1][0]}, radius={circles[1][1]}")
    print(f"Circle 2: center={circles[2][0]}, radius={circles[2][1]}")
    print(f"\nExternal homothetic centers:")
    print(f"  S01 = {S01}")
    print(f"  S12 = {S12}")
    print(f"  S20 = {S20}")
    
    # Check collinearity
    print(f"\nMonge verification: all triples collinear = {ec.all_triples_collinear()}")
    dev = ec.monge_line_deviation(0, 1, 2)
    print(f"Collinearity deviation (should be ~0): {dev:.2e}")
    
    # 3D lifting proof
    print("\n--- 3D Lifting Proof ---")
    spheres = [lift_to_sphere(c, r) for c, r in circles]
    (p1, n1), (p2, n2) = sandwich_planes(spheres)
    print(f"Sandwich planes: normal={n1}")
    print(f"  Plane 1 at: {p1}")
    print(f"  Plane 2 at: {p2}")
    print("All cone apexes lie on the intersection line of the two planes.")


def demo_consensus():
    print("\n" + "=" * 60)
    print("Demo 2: Zero Holonomy Consensus")
    print("=" * 60)
    
    mc = MongeConsensus(3)
    mc.set_agent(0, [0.0, 0.0], 1.0)
    mc.set_agent(1, [4.0, 0.0], 2.0)
    mc.set_agent(2, [2.0, 3.0], 1.5)
    
    print("\n3 agents with different trust radii:")
    print(f"  Agent 0: pos=[0,0], radius=1.0")
    print(f"  Agent 1: pos=[4,0], radius=2.0")
    print(f"  Agent 2: pos=[2,3], radius=1.5")
    
    print(f"\nInitial state:")
    print(f"  Holonomy area: {mc.max_area():.6f}")
    print(f"  λ bound: {mc.lambda_bound():.4f}")
    print(f"  Convergence time: {mc.convergence_time():.1f}ms")
    print(f"  Byzantine (f=1): {mc.byzantine_convergence_time(1):.1f}ms")
    
    # Simulate convergence
    print("\nConvergence (20 steps, alpha=0.3):")
    for step in range(20):
        cent = np.mean(mc.positions, axis=0)
        for i in range(3):
            mc.positions[i] = mc.positions[i] + 0.3 * (cent - mc.positions[i])
        mc._update_centers()
        
        if step < 5 or step == 19:
            print(f"  Step {step:2d}: area={mc.max_area():.2e}, centroid={cent[:2]}")
    
    print(f"\nFinal: convergence time = {mc.convergence_time():.1f}ms")


def demo_p48():
    print("\n" + "=" * 60)
    print("Demo 3: Pythagorean48 Zero Drift Verification")
    print("=" * 60)
    
    result = verify_monge_p48()
    
    print(f"\n48 directions from 6 primitive triples (c ≤ 37)")
    print(f"Triples checked: {result['n_triples_checked']:,}")
    print(f"Violations (tol=1e-6): {result['violations']}")
    print(f"Max deviation area: {result['max_deviation_area']:.2e}")
    print(f"Zero drift verified: {result['zero_drift_verified']}")
    
    if result['zero_drift_verified']:
        print("\n✅ All 17,296 triples satisfy Monge — zero drift proven!")
        print("   The P48 encoding induces no accumulated error.")


def demo_transport():
    print("\n" + "=" * 60)
    print("Demo 4: Fleet Consensus via Optimal Transport")
    print("=" * 60)
    
    n = 5
    positions = [[math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)] for i in range(n)]
    radii = [1.0, 1.2, 0.8, 1.5, 1.0]
    
    result = fleet_consensus(positions, radii, α=0.3, tol=1e-6)
    
    print(f"\n5 agents in pentagon, varying radii")
    print(f"λ bound: {result['lambda_bound']:.4f}")
    print(f"Predicted convergence: {result['predicted_time_ms']:.1f}ms")
    print(f"Actual steps: {result['steps']}")
    print(f"Converged: {result['converged']}")
    print(f"Final W₂: {result['final_wasserstein']:.6f}")


def demo_bridge():
    print("\n" + "=" * 60)
    print("Demo 5: Monge-CT Bridge")
    print("=" * 60)
    
    n = 5
    positions = [[math.cos(2*math.pi*i/n), math.sin(2*math.pi*i/n)] for i in range(n)]
    radii = [1.0, 1.2, 0.8, 1.5, 1.0]
    
    bridge = MongeCTBridge(n, dim=2)
    for i, (pos, r) in enumerate(zip(positions, radii)):
        bridge.set_agent(i, pos, r)
    
    print(f"\nFleet with {n} agents:")
    print(f"  Emergence signal: {bridge.emergence_signal():.6f}")
    print(f"  β₁ (Monge-Menelaus): {bridge.beta1()}")
    print(f"  λ bound: {bridge.lambda_bound():.4f}")
    print(f"  T(f=0): {bridge.convergence_time_ms(f=0):.1f}ms")
    print(f"  T(f=1): {bridge.convergence_time_ms(f=1):.1f}ms")


if __name__ == "__main__":
    print("MONGE-FLEET DEMO")
    print("Monge's theorem for fleet mathematics")
    print("=" * 60)
    
    demo_mongeTheorem()
    demo_consensus()
    demo_p48()
    demo_transport()
    demo_bridge()
    
    print("\n" + "=" * 60)
    print("All demos complete!")
    print("=" * 60)
    print("""
For more:
  pip install monge-fleet
  python -c "import monge_fleet; help(monge_fleet)"

GitHub: github.com/SuperInstance/monge-fleet
PyPI: pypi.org/project/monge-fleet
""")
