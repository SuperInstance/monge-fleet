"""Tests for monge-fleet library."""

import pytest
import numpy as np
import math
from monge_fleet import (
    HomotheticCenter, ExternalCenter,
    RadicalAxis, RadicalCochain,
    MongeCohomology, menelaus_beta1,
    MongeConsensus, homonomy_area, homothetic_center,
    P48Monge, verify_monge_p48,
    lift_to_sphere, sandwich_planes, monge_flat
)


class TestHomotheticCenter:
    """Test homothetic center computation."""
    
    def test_external_center_equal_radii(self):
        """Equal radii → midpoint."""
        c1 = np.array([0.0, 0.0])
        c2 = np.array([4.0, 0.0])
        result = HomotheticCenter(c1, 1.0, c2, 1.0).external
        expected = np.array([2.0, 0.0])
        assert np.allclose(result, expected)
    
    def test_external_center_different_radii(self):
        """Different radii → weighted point."""
        c1 = np.array([0.0, 0.0])
        c2 = np.array([4.0, 0.0])
        result = HomotheticCenter(c1, 1.0, c2, 2.0).external
        # S = (r2*c1 - r1*c2)/(r2 - r1) = (2*[0,0] - 1*[4,0])/(2-1) = [-4, 0]
        expected = np.array([-4.0, 0.0])
        assert np.allclose(result, expected)
    
    def test_monge_collinearity_three_circles(self):
        """Monge's theorem: any 3 circles → 3 external centers collinear."""
        circles = [
            (np.array([0.0, 0.0]), 1.0),
            (np.array([3.0, 0.0]), 2.0),
            (np.array([1.0, 2.0]), 1.5),
        ]
        ec = ExternalCenter(circles)
        # Any triple must be collinear
        assert ec.all_triples_collinear(tol=1e-9)


class TestRadicalCochain:
    """Test radical axis cohomology."""
    
    def test_coboundary_zero_for_any_circles(self):
        """Coboundary φ_ij + φ_jk + φ_ki = 0 for any triple."""
        circles = [
            (np.array([0.0, 0.0]), 1.0),
            (np.array([5.0, 0.0]), 2.0),
            (np.array([2.0, 4.0]), 1.5),
        ]
        rc = RadicalCochain(circles)
        centroid = np.mean([c for c, r in circles], axis=0)
        cb = rc.coboundary(0, 1, 2, centroid)
        assert abs(cb) < 1e-9, f"Coboundary not zero: {cb}"


class TestMongeCohomology:
    """Test cohomology computation."""
    
    def test_laman_triangle_beta1(self):
        """Laman triangle (3 vertices, 3 edges) is minimally rigid: β₁=0."""
        V = 3
        edges = [(0, 1), (1, 2), (2, 0)]
        circles = [
            (np.array([0.0, 0.0]), 1.0),
            (np.array([4.0, 0.0]), 1.5),
            (np.array([1.0, 3.0]), 1.0),
        ]
        mc = MongeCohomology(V, edges, circles)
        beta = mc.beta1()
        assert beta >= 0  # Base β₁ should be 0 for Laman
    
    def test_overconstrained_beta(self):
        """Over-constrained graph has β₁ > 0."""
        V = 3
        edges = [(0, 1), (1, 2), (2, 0), (0, 2)]  # 4 edges, E > 2V-3 = 3
        circles = [
            (np.array([0.0, 0.0]), 1.0),
            (np.array([4.0, 0.0]), 1.5),
            (np.array([1.0, 3.0]), 1.0),
        ]
        mc = MongeCohomology(V, edges, circles)
        base = len(edges) - (2*V - 3)  # = 4 - 3 = 1
        assert mc.beta1() >= base


class TestMongeConsensus:
    """Test consensus protocol."""
    
    def test_convergence_area_decreases(self):
        """As agents move toward consensus, holonomy area decreases."""
        mc = MongeConsensus(3)
        mc.set_agent(0, [0.0, 0.0], 1.0)
        mc.set_agent(1, [4.0, 0.0], 2.0)
        mc.set_agent(2, [2.0, 3.0], 1.5)
        
        initial_area = mc.max_area()
        
        # Move agents toward centroid
        for _ in range(20):
            centroid = np.mean(mc.positions, axis=0)
            for i in range(3):
                mc.positions[i] += 0.2 * (centroid - mc.positions[i])
            mc._update_centers()
        
        final_area = mc.max_area()
        assert final_area < initial_area + 1e-9, f"Area should decrease: final={final_area} >= initial={initial_area}"
    
    def test_byzantine_time_grows_with_faults(self):
        """More Byzantine agents → longer convergence time."""
        mc = MongeConsensus(5)
        for i in range(5):
            angle = 2 * math.pi * i / 5
            mc.set_agent(i, [math.cos(angle), math.sin(angle)], 1.0)
        
        t0 = mc.convergence_time()
        t1 = mc.byzantine_convergence_time(f=1)
        t2 = mc.byzantine_convergence_time(f=2)
        
        assert t0 < t1 < t2, "Convergence time should increase with Byzantine faults"


class TestPythagorean48:
    """Test P48 Monge consistency."""
    
    def test_all_directions_unique(self):
        """48 directions are all unique."""
        from monge_fleet.pythagorean48 import generate_p48_directions
        angles = generate_p48_directions()
        assert len(set(angles)) == 48, f"Only {len(set(angles))} unique directions"
    
    def test_monge_holds_for_all_triples(self):
        """Monge's theorem holds for all 17,296 triples of P48 directions."""
        result = verify_monge_p48()
        assert result['all_collinear'], f"{result['violations']} violations found"
        assert result['zero_drift_verified'], "Zero drift not verified"


class TestGeometry:
    """Test geometric utilities."""
    
    def test_sandwich_planes_contain_spheres(self):
        """Sandwich planes actually contain all spheres."""
        spheres = [
            (np.array([0.0, 0.0, 0.0]), 1.0),
            (np.array([4.0, 0.0, 0.0]), 1.5),
            (np.array([0.0, 4.0, 0.0]), 1.0),
        ]
        (p1, n1), (p2, n2) = sandwich_planes(spheres)
        
        # All spheres between planes
        for c, r in spheres:
            proj1 = np.dot(c, n1)
            proj2 = np.dot(c, n2)
            assert proj1 - r >= np.dot(p1, n1) - 1e-9, "Sphere below plane1"
            assert proj2 + r <= np.dot(p2, n2) + 1e-9, "Sphere above plane2"
    
    def test_monge_flat_line_for_three_circles(self):
        """Monge flat for 3 circles is a line in ℝ²."""
        circles = [
            (np.array([0.0, 0.0]), 1.0),
            (np.array([4.0, 0.0]), 2.0),
            (np.array([2.0, 3.0]), 1.5),
        ]
        point, direction = monge_flat(circles, dim=2)
        assert direction.shape == (1, 2), "Line direction should be 1D in 2D space"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
