# SPDX-License-Identifier: LGPL-3.0-or-later
"""Quadratic-triangle finite elements for axisymmetric magnetostatics with nonlinear iron.

Unknown: the azimuthal vector potential A = A_phi(r, z) at the P2 nodes.  With B = curl(A e_phi):
    B_r = -dA/dz,   B_z = dA/dr + A/r,
and the weak form of curl(nu curl A) = J_phi (per radian, dV = r dr dz):
    R_i(A) = int nu(|B|) B . c_i r dr dz - int J N_i r dr dz = 0,   c_i = (-dN_i/dz, dN_i/dr + N_i/r).
A = 0 on the axis (symmetry) and on the outer boundary.  Because A vanishes along the whole axis edge of an
element, A/r is a polynomial there and B_z(0, z) = 2 dA/dr is finite; B_r = 0 on the axis exactly.
Newton:  K dA = -R  with  K_ij = int [nu c_i.c_j + (dnu/d|B|)/|B| (B.c_i)(B.c_j)] r dr dz  (symmetric positive
definite for a monotonic B-H curve), followed by a line search on the energy (a root of dA . R(A + a dA) in a).
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

from .bh import MU0

# Dunavant degree-5 rule, 7 points (barycentric l0, l1, l2; weights sum to 1)
_A1, _B1, _W1 = 0.059715871789770, 0.470142064105115, 0.132394152788506
_A2, _B2, _W2 = 0.797426985353087, 0.101286507323456, 0.125939180544827
QUAD7 = np.array([
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.225],
    [_A1, _B1, _B1, _W1], [_B1, _A1, _B1, _W1], [_B1, _B1, _A1, _W1],
    [_A2, _B2, _B2, _W2], [_B2, _A2, _B2, _W2], [_B2, _B2, _A2, _W2],
])


def p2_basis(l1, l2):
    """P2 basis and reference derivatives at reference coordinates (xi, eta) = (l1, l2).

    Returns N (..., 6), dN/dxi (..., 6), dN/deta (..., 6) for the node order v0, v1, v2, m01, m12, m20."""
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    l0 = 1.0 - l1 - l2
    N = np.stack([l0 * (2 * l0 - 1), l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), 4 * l0 * l1, 4 * l1 * l2, 4 * l2 * l0], -1)
    z = np.zeros_like(l1)
    dxi = np.stack([-(4 * l0 - 1), 4 * l1 - 1, z, 4 * (l0 - l1), 4 * l2, -4 * l2], -1)
    deta = np.stack([-(4 * l0 - 1), z, 4 * l2 - 1, -4 * l1, 4 * l1, 4 * (l0 - l2)], -1)
    return N, dxi, deta


def element_geometry(nodes, tri):
    """Affine map data: vertices x (ne,3), y (ne,3), det (ne,), inverse-Jacobian terms."""
    x = nodes[tri[:, :3], 0]
    y = nodes[tri[:, :3], 1]
    j00 = x[:, 1] - x[:, 0]
    j01 = x[:, 2] - x[:, 0]
    j10 = y[:, 1] - y[:, 0]
    j11 = y[:, 2] - y[:, 0]
    det = j00 * j11 - j01 * j10
    # d(xi)/dr, d(xi)/dz, d(eta)/dr, d(eta)/dz
    inv = np.stack([j11 / det, -j01 / det, -j10 / det, j00 / det], -1)
    return x, y, det, inv


class AxiMagnetostatics:
    """Assembles and solves the nonlinear r-z problem on a rzfield.mesh.Mesh (coordinates in metres)."""

    def __init__(self, mesh, zone_materials: list, quad=QUAD7, log=print):
        self.mesh = mesh
        self.log = log or (lambda *a, **k: None)
        t0 = time.time()
        nodes, tri = mesh.nodes, mesh.tri
        ne = len(tri)
        nq = len(quad)
        self.mats = zone_materials
        x, y, det, inv = element_geometry(nodes, tri)
        if np.any(det <= 0):
            raise RuntimeError("mesh has %d degenerate or inverted elements" % int(np.sum(det <= 0)))
        N, dxi, deta = p2_basis(quad[:, 1], quad[:, 2])          # (nq, 6)
        lam = quad[:, :3]
        rq = x @ lam.T                                            # (ne, nq) radius of the quadrature points
        if np.any(rq <= 0):
            raise RuntimeError("quadrature point on or left of the axis")
        dNdr = inv[:, None, None, 0] * dxi[None] + inv[:, None, None, 2] * deta[None]   # (ne, nq, 6)
        dNdz = inv[:, None, None, 1] * dxi[None] + inv[:, None, None, 3] * deta[None]
        self.cr = -dNdz
        self.cz = dNdr + N[None] / rq[:, :, None]
        del dNdr, dNdz
        self.wr = 0.5 * det[:, None] * quad[None, :, 3] * rq      # (ne, nq): weight * |J| * r
        self.Nq = N
        # load vector J N_i
        J = np.array([z.get("J", 0.0) for z in mesh.zones])[mesh.zone]
        Fe = (self.wr * J[:, None]) @ N                           # (ne, 6)
        # degrees of freedom
        nn = mesh.n_nodes
        free = np.ones(nn, bool)
        free[mesh.dirichlet] = False
        self.free = free
        dof = np.full(nn, -1, dtype=np.int64)
        dof[free] = np.arange(int(free.sum()))
        self.dof = dof
        self.nf = int(free.sum())
        de = dof[tri]                                             # (ne, 6)
        self.F = np.bincount(de[de >= 0], weights=Fe[de >= 0], minlength=self.nf)
        # sparsity of the reduced matrix
        rows = np.repeat(de[:, :, None], 6, axis=2)
        cols = np.repeat(de[:, None, :], 6, axis=1)
        valid = (rows >= 0) & (cols >= 0)
        keys = rows[valid] * self.nf + cols[valid]
        uk, invv = np.unique(keys, return_inverse=True)
        del keys, rows, cols
        self.nnz = len(uk)
        self.indptr = np.searchsorted(uk // self.nf, np.arange(self.nf + 1)).astype(np.int64)
        self.indices = (uk % self.nf).astype(np.int32)
        del uk
        pos = np.full((ne, 36), -1, dtype=np.int64)
        pos[valid.reshape(ne, 36)] = invv.ravel()
        self.pos = pos
        # element groups by material: linear ones are assembled once
        self.groups = []                                          # (element indices, material)
        lin_elems = []
        nu_lin = np.zeros(ne)
        for zi, m in enumerate(zone_materials):
            el = np.where(mesh.zone == zi)[0]
            if len(el) == 0:
                continue
            if m.nonlinear:
                self.groups.append((el, m))
            else:
                lin_elems.append(el)
                nu_lin[el] = 1.0 / (MU0 * m.mu_r)
        self.lin = np.concatenate(lin_elems) if lin_elems else np.zeros(0, dtype=np.int64)
        self.nu_lin = nu_lin
        self.data_lin = self._assemble_data(self.lin, self._ke_linear(self.lin, nu_lin[self.lin]))
        self.nonlin = (np.concatenate([g[0] for g in self.groups]) if self.groups
                       else np.zeros(0, dtype=np.int64))
        self.timing = {"setup": time.time() - t0}
        self.log("FEM setup: %d elements, %d nodes, %d free DOFs, %d nonzeros, %d nonlinear elements (%.1f s)"
                 % (ne, nn, self.nf, self.nnz, len(self.nonlin), self.timing["setup"]))

    # ------------------------------------------------------------------------------------------ element kernels
    def _ke_linear(self, el, nu):
        w = self.wr[el] * nu[:, None]                             # (n, nq)
        cr, cz = self.cr[el], self.cz[el]
        return np.einsum("eq,eqi,eqj->eij", w, cr, cr, optimize=True) + \
            np.einsum("eq,eqi,eqj->eij", w, cz, cz, optimize=True)

    def _assemble_data(self, el, ke):
        p = self.pos[el].ravel()
        v = ke.reshape(-1)
        ok = p >= 0
        return np.bincount(p[ok], weights=v[ok], minlength=self.nnz)

    def _b_at_quad(self, Afull, el):
        Ae = Afull[self.mesh.tri[el]]                              # (n, 6)
        Br = np.einsum("eqi,ei->eq", self.cr[el], Ae)
        Bz = np.einsum("eqi,ei->eq", self.cz[el], Ae)
        return Br, Bz

    def full(self, Afree):
        A = np.zeros(self.mesh.n_nodes)
        A[self.free] = Afree
        return A

    def residual(self, Afree):
        """R = K_lin A + sum_nonlinear int nu B.c - F (reduced)."""
        K_lin = self._k_lin()
        R = K_lin @ Afree - self.F
        if len(self.nonlin):
            A = self.full(Afree)
            for el, m in self.groups:
                Br, Bz = self._b_at_quad(A, el)
                nu, _ = m.nu(Br * Br + Bz * Bz)
                w = self.wr[el] * nu
                re = np.einsum("eq,eqi->ei", w * Br, self.cr[el]) + np.einsum("eq,eqi->ei", w * Bz, self.cz[el])
                de = self.dof[self.mesh.tri[el]]
                ok = de >= 0
                R += np.bincount(de[ok], weights=re[ok], minlength=self.nf)
        return R

    def _k_lin(self):
        if not hasattr(self, "_Klin"):
            self._Klin = sp.csr_matrix((self.data_lin, self.indices, self.indptr), shape=(self.nf, self.nf))
        return self._Klin

    def jacobian(self, Afree):
        data = self.data_lin.copy()
        if len(self.nonlin):
            A = self.full(Afree)
            for el, m in self.groups:
                Br, Bz = self._b_at_quad(A, el)
                nu, dnub = m.nu(Br * Br + Bz * Bz)
                cr, cz = self.cr[el], self.cz[el]
                w = self.wr[el] * nu
                ke = np.einsum("eq,eqi,eqj->eij", w, cr, cr, optimize=True) + \
                    np.einsum("eq,eqi,eqj->eij", w, cz, cz, optimize=True)
                bc = Br[:, :, None] * cr + Bz[:, :, None] * cz        # (n, nq, 6) = B . c_i
                ke += np.einsum("eq,eqi,eqj->eij", self.wr[el] * dnub, bc, bc, optimize=True)
                data += self._assemble_data(el, ke)
        return sp.csr_matrix((data, self.indices, self.indptr), shape=(self.nf, self.nf))

    # ---------------------------------------------------------------------------------------------- solution
    def solve(self, rtol: float = 1e-9, atol: float = 1e-9, maxit: int = 60, pcg_maxit: int = 40,
              linear_solver: str = "lu"):
        """Newton iteration from A = 0.  Converged when ||R|| <= rtol ||F|| and max|dA| <= atol max|A|.

        Each Newton system K dA = -R is solved by conjugate gradients preconditioned with the sparse LU of an
        earlier Jacobian (inexact Newton, forcing term min(1e-2, 0.1 ||R||/||F||)); the Jacobian is factorized
        again (SuperLU, MMD on A^T + A) after a large Newton step (max|dA|/max|A| > 0.05) and whenever CG does not
        converge within pcg_maxit iterations.
        linear_solver = "amg": smoothed-aggregation AMG (pyamg) + CG instead, rebuilt every Newton step, with the
        same forcing term floored at 1e-8; slower for the default meshes but needs far less memory (refined
        meshes with > 1e6 nodes).
        """
        if linear_solver not in ("lu", "amg"):
            raise ValueError("linear_solver must be 'lu' or 'amg'")
        t0 = time.time()
        Af = np.zeros(self.nf)
        Fn = np.linalg.norm(self.F)
        hist = []
        R = self.residual(Af)
        converged = False
        t_fact = 0.0
        lu = None
        nfact = 0
        for it in range(1, maxit + 1):
            K = self.jacobian(Af)
            eta = max(1e-13, min(1e-2, 0.1 * float(np.linalg.norm(R)) / Fn))
            d, nit, ok = None, 0, False
            if linear_solver == "amg":
                t1 = time.time()
                d, nit = _amg_cg(K, -R, max(eta, 1e-8), log=self.log)
                t_fact += time.time() - t1
                ok = True
            elif lu is not None and (not hist or hist[-1]["dA"] < 0.05):
                d, nit, ok = _pcg(K, -R, lu.solve, eta, pcg_maxit)
            if not ok:
                lu = None  # release the old factors before computing new ones (memory)
                t1 = time.time()
                lu = splu(K.tocsc(), permc_spec="MMD_AT_PLUS_A", options={"SymmetricMode": True})
                t_fact += time.time() - t1
                nfact += 1
                d = -lu.solve(R)
                nit = 0
            elif nit > pcg_maxit // 2:
                lu = None  # converged, but slowly: refresh the preconditioner next time
            g0 = float(d @ R)
            alpha = 1.0
            R1 = self.residual(Af + d)
            g1 = float(d @ R1)
            nls = 0
            if g1 > 0 and g0 < 0 and len(self.nonlin):
                # energy minimum inside (0, 1): Illinois regula falsi on g(a) = d . R(A + a d)
                a0, ga, a1, gb = 0.0, g0, 1.0, g1
                side = 0
                for nls in range(1, 13):
                    a = (a0 * gb - a1 * ga) / (gb - ga)
                    Ra = self.residual(Af + a * d)
                    g = float(d @ Ra)
                    if abs(g) < 0.05 * abs(g0):
                        break
                    if g > 0:
                        a1, gb = a, g
                        if side == -1:
                            ga *= 0.5
                        side = -1
                    else:
                        a0, ga = a, g
                        if side == 1:
                            gb *= 0.5
                        side = 1
                alpha, R1 = a, Ra
            Af = Af + alpha * d
            R = R1
            rn = float(np.linalg.norm(R)) / Fn
            dn = float(np.max(np.abs(alpha * d))) / max(float(np.max(np.abs(Af))), 1e-300)
            hist.append({"it": it, "alpha": alpha, "line_search_evals": nls, "cg_iterations": nit,
                         "factorized": linear_solver == "lu" and nit == 0, "res": rn, "dA": dn})
            how = ("AMG-CG %3d" % nit) if linear_solver == "amg" else ("LU" if nit == 0 else "PCG %2d" % nit)
            self.log("  Newton %2d: %s, alpha %.3f (%d), |R|/|F| %.3e, max|dA|/max|A| %.3e"
                     % (it, how, alpha, nls, rn, dn))
            if not len(self.nonlin):
                converged = True
                break
            if rn <= rtol and dn <= atol:
                converged = True
                break
        self.timing["solve"] = time.time() - t0
        self.timing["factorize"] = t_fact
        self.timing["n_factorizations"] = nfact
        if not converged:
            self.log("WARNING: Newton did not converge in %d iterations" % maxit)
        self.A = self.full(Af)
        self.history = hist
        self.converged = converged
        return self.A

    def energy(self, A=None):
        """Stored magnetic energy (J) = 2 pi int w(|B|) r dr dz."""
        A = self.A if A is None else A
        W = 0.0
        for zi, m in enumerate(self.mats):
            el = np.where(self.mesh.zone == zi)[0]
            if len(el) == 0:
                continue
            Br, Bz = self._b_at_quad(A, el)
            W += float(np.sum(self.wr[el] * m.energy_density(np.sqrt(Br * Br + Bz * Bz))))
        return 2.0 * np.pi * W

    def zone_b_stats(self, A=None):
        """Per zone: volume-weighted mean and max |B| at the quadrature points (tesla)."""
        A = self.A if A is None else A
        out = []
        for zi, m in enumerate(self.mats):
            el = np.where(self.mesh.zone == zi)[0]
            if len(el) == 0:
                out.append(None)
                continue
            Br, Bz = self._b_at_quad(A, el)
            b = np.sqrt(Br * Br + Bz * Bz)
            out.append({"mean_T": float(np.sum(b * self.wr[el]) / np.sum(self.wr[el])), "max_T": float(b.max())})
        return out


def _pcg(K, b, msolve, rtol, maxiter):
    """Preconditioned conjugate gradients from x = 0; returns (x, iterations, converged)."""
    bn = float(np.linalg.norm(b))
    x = np.zeros_like(b)
    if bn == 0.0:
        return x, 0, True
    r = b.copy()
    z = msolve(r)
    p = z.copy()
    rz = float(r @ z)
    for k in range(1, maxiter + 1):
        Kp = K @ p
        pKp = float(p @ Kp)
        if pKp <= 0:
            return x, k, False
        al = rz / pKp
        x += al * p
        r -= al * Kp
        if float(np.linalg.norm(r)) <= rtol * bn:
            return x, k, True
        z = msolve(r)
        rz_new = float(r @ z)
        p = z + (rz_new / rz) * p
        rz = rz_new
    return x, maxiter, False


def _amg_cg(K, b, rtol, log=None):
    """CG with a smoothed-aggregation AMG preconditioner on the symmetrically diagonal-scaled system.

    An inexact step is acceptable inside Newton: if CG stalls above rtol but reduced the residual by at least
    1e3, the step is used (with a warning); otherwise an error is raised."""
    try:
        import pyamg
    except ImportError as e:  # pragma: no cover
        raise ImportError("linear_solver='amg' needs pyamg (pip install pyamg)") from e
    dg = 1.0 / np.sqrt(K.diagonal())
    D = sp.diags(dg)
    Ks = (D @ K @ D).tocsr()
    ml = pyamg.smoothed_aggregation_solver(Ks, symmetry="symmetric", strength=("evolution", {}), max_coarse=500)
    res = []
    y = ml.solve(b * dg, tol=rtol, accel="cg", maxiter=1000, residuals=res)
    red = res[-1] / res[0] if res and res[0] > 0 else 0.0
    if red > 10 * rtol:
        if red > 1e-3:
            raise RuntimeError("AMG-CG did not converge (relative residual %.1e)" % red)
        if log:
            log("    AMG-CG stopped at relative residual %.1e (asked %.1e); using the inexact step" % (red, rtol))
    return y * dg, len(res) - 1
