# SPDX-License-Identifier: LGPL-3.0-or-later
"""rzfield: 2D axisymmetric (r-z) nonlinear magnetostatic solver for SoLID solenoid field maps.

Modules: amdeck (Poisson .am reader), bh (B-H curves), mesh (gmsh P2 meshing), fem (P2 assembly + Newton),
field (solution I/O and B(r, z) evaluation), gemcmap (GEMC map reader / writer), geometry (polygon helpers),
driver (deck -> solution pipeline used by the command-line tools).  Only mesh needs gmsh.
"""
__version__ = "0.1.0"
