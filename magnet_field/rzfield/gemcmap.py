# SPDX-License-Identifier: LGPL-3.0-or-later
"""GEMC 2D (cylindrical-z) field maps: rows 'R Z Br Bz' in cm and gauss, R-major with Z fastest, after an
<mfield> header laid out exactly like solid_gemc/field/solenoid_CLEOv10.dat (read by SoLID_FieldMapB)."""
from __future__ import annotations

import re

import numpy as np

_HEADER = (
    "<mfield>\n"
    '       <description name="{name}" factory="ASCII" comment="{comment}"/>\n'
    '       <symmetry type="cylindrical-z" format="map" integration="ClassicalRK4" minStep="0.01*mm"/>\n'
    "       <map>\n"
    "               <coordinate>\n"
    '                       <first  name="transverse" npoints="{nr}"   min="{rmin}"  max="{rmax}"  units="cm"/>\n'
    '                       <second name="longitudinal"   npoints="{nz}"  min="{zmin}" max="{zmax}"  units="cm"/>\n'
    "               </coordinate>\n"
    '               <field unit="gauss"/>\n'
    '               <interpolation type="linear"/>\n'
    '               <shift z="0" units="cm"/>\n'
    "       </map>\n"
    "</mfield>\n"
)


def grid_axis(vmin: float, vmax: float, step: float) -> np.ndarray:
    n = int(round((vmax - vmin) / step)) + 1
    if n < 2 or abs(vmin + (n - 1) * step - vmax) > 1e-6 * max(1.0, abs(step)):
        raise ValueError("range %g..%g is not a multiple of the step %g" % (vmin, vmax, step))
    return np.round(vmin + step * np.arange(n), 9) + 0.0


def _num(v: float) -> str:
    return "%d" % int(round(v)) if abs(v - round(v)) < 1e-9 else "%g" % v


def _attr(text: str) -> str:
    """text usable inside a double-quoted header attribute"""
    return text.replace('"', "'").replace("<", "(").replace(">", ")").replace("&", "+")


def header(name: str, r, z, comment: str | None = None) -> str:
    return _HEADER.format(name=_attr(name), comment=_attr(comment if comment is not None else name),
                          nr=len(r), nz=len(z), rmin=_num(r[0]), rmax=_num(r[-1]), zmin=_num(z[0]),
                          zmax=_num(z[-1]))


def write_map(path: str, name: str, r_cm, z_cm, Br_G, Bz_G, comment: str | None = None):
    """Write a GEMC map; Br_G, Bz_G have shape (len(r_cm), len(z_cm))."""
    r_cm = np.asarray(r_cm, dtype=float)
    z_cm = np.asarray(z_cm, dtype=float)
    Br_G = np.asarray(Br_G, dtype=float)
    Bz_G = np.asarray(Bz_G, dtype=float)
    if Br_G.shape != (len(r_cm), len(z_cm)) or Bz_G.shape != Br_G.shape:
        raise ValueError("field arrays must have shape (nR, nZ)")
    if not (np.all(np.isfinite(Br_G)) and np.all(np.isfinite(Bz_G))):
        raise ValueError("map contains non-finite values (points outside the mesh?)")
    R, Z = np.meshgrid(r_cm, z_cm, indexing="ij")
    # avoid '-0' in the output
    cols = [R.ravel() + 0.0, Z.ravel() + 0.0, Br_G.ravel() + 0.0, Bz_G.ravel() + 0.0]
    with open(path, "w") as f:
        f.write(header(name, r_cm, z_cm, comment))
        lines = ["%#.6g\t%#.6g\t%.6E\t%.6E\n" % row for row in zip(*[c.tolist() for c in cols])]
        f.writelines(lines)


def read_map(path: str) -> dict:
    """Read a GEMC 2D map (with <mfield> header) -> dict(name, r, z, Br, Bz) with fields shaped (nR, nZ)."""
    head = []
    nhead = 0
    with open(path) as f:
        for line in f:
            if not head and not line.strip():
                nhead += 1          # blank lines before the header
                continue
            if not head and not line.lstrip().startswith("<mfield>"):
                break               # no header
            head.append(line)
            nhead += 1
            if "</mfield>" in line:
                break
    text = "".join(head)
    if "<mfield>" not in text:
        raise ValueError("%s: no <mfield> header" % path)

    def axis(which):
        m = re.search(r"<%s\b([^>]*)>" % which, text)
        if not m:
            raise ValueError("%s: no <%s> axis in the header" % (path, which))
        at = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        return int(at["npoints"]), float(at["min"]), float(at["max"]), at.get("units", "cm")

    nr, rmin, rmax, ur = axis("first")
    nz, zmin, zmax, uz = axis("second")
    if ur != "cm" or uz != "cm":
        raise ValueError("only cm maps are supported")
    m = re.search(r'<field\s+unit="([^"]*)"', text)
    unit = m.group(1).lower() if m else "gauss"
    fac = {"gauss": 1.0, "kilogauss": 1000.0, "tesla": 1.0e4, "t": 1.0e4, "kg": 1000.0, "g": 1.0}[unit]
    m = re.search(r'<description\s+name="([^"]*)"', text)
    name = m.group(1) if m else ""
    data = np.loadtxt(path, skiprows=nhead, dtype=float)
    if data.shape != (nr * nz, 4):
        raise ValueError("%s: %s data rows, header says %d x %d" % (path, data.shape, nr, nz))
    r = data[::nz, 0]
    z = data[:nz, 1]
    if not (np.allclose(r, np.linspace(rmin, rmax, nr)) and np.allclose(z, np.linspace(zmin, zmax, nz))):
        raise ValueError("%s: rows are not R-major / Z-fastest on the header grid" % path)
    return {"name": name, "r": r, "z": z, "Br": data[:, 2].reshape(nr, nz) * fac,
            "Bz": data[:, 3].reshape(nr, nz) * fac}


def map_from_solution(sol, r_cm, z_cm, air_side: bool = False, info: dict | None = None):
    """Evaluate an RZSolution on the grid: (Br, Bz) in gauss, shape (len(r), len(z)).

    Default: the field at each node (air-side value for nodes exactly on an iron face, see RZSolution.locate).
    A cell cut by an iron face then interpolates between iron and air values: within one cell of the face the map
    is off by up to the jump of the tangential field (kG next to saturated iron), on both sides, as in POISSON/SF7
    maps. air_side=True replaces the values at the iron nodes of the cells that contain air by an extrapolation of
    the air field (air_side_values), so that the bilinear interpolation is accurate in air up to the iron faces; in
    turn the map is wrong in the iron within about two cells of its faces. info (a dict) receives statistics."""
    R, Z = np.meshgrid(np.asarray(r_cm, float), np.asarray(z_cm, float), indexing="ij")
    Br, Bz, in_air = sol.B_cm_gauss(R, Z, return_air=True)
    if info is not None:
        info.update(n_nodes=int(R.size), n_iron_nodes=int((~in_air).sum()), air_side=bool(air_side))
    if air_side:
        Br, Bz = air_side_values(sol, np.asarray(r_cm, float), np.asarray(z_cm, float), Br, Bz, in_air, info=info)
    return Br, Bz


def air_side_values(sol, r, z, Br, Bz, in_air, sub: int = 4, lam: float = 0.01, info: dict | None = None):
    """Map values (gauss) with the iron nodes next to air replaced by an extrapolation of the air field.

    Every grid cell with an iron corner is sampled on a sub x sub grid of interior points. The unknowns are the
    values at the iron corners of the cells with at least one air sample; they minimise, in the least-squares
    sense over all air samples, the difference between the bilinear interpolation of the map (as SoLID_FieldMapB
    computes it) and the air field of the solution. A weak regularisation (weight lam per unknown, against about
    1 per air sample) pulls each unknown towards the air field at the nearest air sample, which fixes the nodes
    whose cells hold almost no air. Nodes in air keep their values. Needs a uniform grid."""
    from scipy.sparse import coo_matrix, identity
    from scipy.sparse.linalg import splu
    from scipy.spatial import cKDTree

    nr, nz = len(r), len(z)
    hr, hz = r[1] - r[0], z[1] - z[0]
    if not (np.allclose(np.diff(r), hr) and np.allclose(np.diff(z), hz)):
        raise ValueError("air_side needs a uniform grid")
    if info is not None:
        info.update(n_air_side_nodes=0, sub=sub, lam=lam)
    iron = ~np.asarray(in_air, bool)
    out_r, out_z = np.array(Br, float), np.array(Bz, float)
    cell = iron[:-1, :-1] | iron[1:, :-1] | iron[:-1, 1:] | iron[1:, 1:]
    ci, cj = np.nonzero(cell)
    if ci.size == 0:
        return out_r, out_z
    t = (np.arange(sub) + 0.5) / sub
    ta, tb = (x.ravel() for x in np.meshgrid(t, t, indexing="ij"))
    ns = ta.size
    sr = (r[ci][:, None] + hr * ta[None, :]).ravel()
    sz = (z[cj][:, None] + hz * tb[None, :]).ravel()
    sbr, sbz, sair = sol.B_cm_gauss(sr, sz, return_air=True)
    keep = sair & np.isfinite(sbr) & np.isfinite(sbz)
    c = np.repeat(np.arange(ci.size), ns)[keep]
    a, b = np.tile(ta, ci.size)[keep], np.tile(tb, ci.size)[keep]
    i0, j0 = ci[c], cj[c]
    corners = ((i0, j0, (1 - a) * (1 - b)), (i0, j0 + 1, (1 - a) * b), (i0 + 1, j0, a * (1 - b)),
               (i0 + 1, j0 + 1, a * b))
    unknown = np.zeros((nr, nz), dtype=bool)
    for ii, jj, _ in corners:
        unknown[ii, jj] |= iron[ii, jj]
    nu = int(unknown.sum())
    if nu == 0:
        return out_r, out_z
    uid = np.full((nr, nz), -1, dtype=np.int64)
    uid[unknown] = np.arange(nu)
    nsamp = c.size
    rows, cols, vals = [], [], []
    rhs = np.column_stack([sbr[keep], sbz[keep]])
    for ii, jj, w in corners:
        u = uid[ii, jj]
        isu = u >= 0
        rows.append(np.nonzero(isu)[0])
        cols.append(u[isu])
        vals.append(w[isu])
        rhs[~isu, 0] -= w[~isu] * out_r[ii[~isu], jj[~isu]]      # air corners keep their own values
        rhs[~isu, 1] -= w[~isu] * out_z[ii[~isu], jj[~isu]]
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nsamp, nu)).tocsr()
    # regularisation target: the air field at the nearest air sample (np.nonzero is in the order of uid)
    ui, uj = np.nonzero(unknown)
    _, near = cKDTree(np.column_stack([sr[keep], sz[keep]])).query(np.column_stack([r[ui], z[uj]]))
    target = np.column_stack([sbr[keep][near], sbz[keep][near]])
    N = (A.T @ A + lam * identity(nu, format="csr")).tocsc()
    v = splu(N).solve(A.T @ rhs + lam * target)
    v[np.abs(r[ui]) < 1e-9, 0] = 0.0                               # Br = 0 on the axis
    out_r[ui, uj] = v[:, 0]
    out_z[ui, uj] = v[:, 1]
    if info is not None:
        info.update(n_air_side_nodes=nu, n_cells_sampled=int(ci.size), n_air_samples=int(nsamp))
    return out_r, out_z
