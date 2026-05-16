"""
Monge-constraint theory bridge.

Connects Monge cohomology to the existing jc1-ct-bridge constraint theory.
The Monge layer sits above CT and provides the geometric invariants
that enhance the emergence detection and consensus protocols.
"""

from .monge_ct import MongeCTBridge, monge_emergence_detection, monge_consensus

__all__ = ["MongeCTBridge", "monge_emergence_detection", "monge_consensus"]
