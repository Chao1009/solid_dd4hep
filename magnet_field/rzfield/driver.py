# SPDX-License-Identifier: LGPL-3.0-or-later
"""High-level helpers shared by the command-line tools: deck -> mesh -> nonlinear solve -> RZSolution."""
from __future__ import annotations

import time

import numpy as np

from .amdeck import read_deck
from .bh import material_from_spec
from .field import RZSolution


def join_negative_values(argv, options):
    """argparse treats '-500:700:1' as an option: rewrite '--opt -500:700:1' to '--opt=-500:700:1'."""
    out = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in options and i + 1 < len(argv) and argv[i + 1].startswith("-") and argv[i + 1][1:2].isdigit():
            out.append(a + "=" + argv[i + 1])
            i += 2
            continue
        out.append(a)
        i += 1
    return out


def parse_range(text: str):
    """'min:max:step' (cm) -> (min, max, step)."""
    parts = [float(x) for x in text.split(":")]
    if len(parts) != 3 or parts[2] <= 0 or parts[1] <= parts[0]:
        raise ValueError("range must be min:max:step with max > min and step > 0, got %r" % text)
    return parts[0], parts[1], parts[2]


def zone_materials(deck, zones, overrides: dict | None = None, interp: str = "poisson"):
    """Material law per mesh zone: the deck's (Poisson rules) unless overridden per material number."""
    overrides = overrides or {}
    cache = {}
    out = []
    for z in zones:
        mat = int(z["mat"])
        if mat not in cache:
            if mat >= 2 and mat in overrides:
                cache[mat] = material_from_spec(overrides[mat], interp=interp)
            elif mat >= 2 and "all" in overrides:
                cache[mat] = material_from_spec(overrides["all"], interp=interp)
            else:
                cache[mat] = deck.material(mat, interp=interp)
        out.append(cache[mat])
    return out


def solve_deck(deck_path: str, mode: str = "poisson", mesh_scale: float = 1.0, open_factor: float = 10.0,
               overrides: dict | None = None, interp: str = "poisson", refine_box=None, sizes: dict | None = None,
               rtol: float = 1e-9, atol: float = 1e-9, maxit: int = 60, threads: int = 1, linear_solver: str = "lu",
               log=print):
    from .fem import AxiMagnetostatics
    from .mesh import build_mesh

    t0 = time.time()
    deck = read_deck(deck_path)
    log(deck.summary())
    mesh = build_mesh(deck, mode=mode, scale=mesh_scale, open_factor=open_factor, refine_box=refine_box,
                      sizes=sizes, threads=threads)
    log("mesh (%s mode): %d quadratic triangles, %d nodes, %.1f s%s" % (
        mode, mesh.info["n_elements"], mesh.info["n_nodes"], mesh.info["mesh_seconds"],
        (", outer radius %.0f cm" % mesh.info["R_out_cm"]) if mode == "open" else ""))
    mats = zone_materials(deck, mesh.zones, overrides, interp)
    for z, m in zip(mesh.zones, mats):
        z["material"] = m.describe()
    names = sorted({m.describe() for z, m in zip(mesh.zones, mats) if z["mat"] >= 2})
    for n in names:
        log("iron material: %s" % n)
    fem = AxiMagnetostatics(mesh, mats, log=log)
    A = fem.solve(rtol=rtol, atol=atol, maxit=maxit, linear_solver=linear_solver)
    energy = fem.energy()
    zb = fem.zone_b_stats()
    for z, s in zip(mesh.zones, zb):
        if s is not None:
            z["B_mean_T"] = s["mean_T"]
            z["B_max_T"] = s["max_T"]
    info = dict(mesh.info)
    info.update({"deck": deck_path, "title": deck.title, "interp": interp, "overrides": overrides or {},
                 "newton": fem.history, "converged": fem.converged, "rtol": rtol, "atol": atol,
                 "linear_solver": linear_solver,
                 "energy_J": energy, "timing": fem.timing, "total_seconds": time.time() - t0,
                 "total_current_At": float(sum(r.cur for r in deck.regions[1:]))})
    log("stored energy %.4f MJ, Newton %s in %d iterations, total %.1f s" % (
        energy / 1e6, "converged" if fem.converged else "NOT converged", len(fem.history), info["total_seconds"]))
    sol = RZSolution(mesh.nodes, mesh.tri, A, zone=mesh.zone, zones=mesh.zones, info=info)
    return sol, deck


def standard_points(sol, points_cm):
    out = []
    for r, z in points_cm:
        Br, Bz = sol.B_cm_gauss(np.array([r]), np.array([z]))
        out.append((r, z, float(Br[0]), float(Bz[0])))
    return out
