# SPDX-License-Identifier: LGPL-3.0-or-later
"""Self-checks of rzfield (python3 tests/test_rzfield.py, or pytest).

1. Air-core thick solenoid (r 100-104 cm, z +-100 cm, 1 MA-turn), open boundary at 40x the coil radius:
   Bz on the axis against the exact formula
       Bz(0, z) = mu0 J / 2 [F(z + L/2) - F(z - L/2)],  F(u) = u ln((a2 + sqrt(a2^2 + u^2)) / (a1 + sqrt(a1^2 + u^2))),
   plus Br = 0 on the axis and the stored energy against 1/2 int J A dV.
2. B-H curve: the POISSON 1010 table round trip (gamma, H, energy density, monotonicity).
3. GEMC map writer / reader round trip and header layout.
4. Map nodes on iron faces and the air-side option: a coil inside an iron cylinder with two faces on the 1 cm grid
   lines (r = 30, z = 25: the nodes on them must take the air-side field) and two between them (r = 40.6,
   z = -25.3: with air_side=True the bilinear interpolation of the map must follow the air field up to them).
5. Provenance of extract_rz_model.py: paths relative to the repository, ${DD4hepINSTALL}/... through a view of
   symbolic links, never a host directory; the geometry digest ignores comments, <info> and (for a top-level read
   with --no-field) the field includes, but not a changed dimension.
"""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from rzfield.bh import BHCurve, GAUSS, MU0, OERSTED  # noqa: E402
from rzfield.gemcmap import grid_axis, header, read_map, write_map  # noqa: E402

DECK = """Air-core test solenoid
&reg kprob=0, icylin=1, xjfact=1.0, mode=0, dx=1.0 &
&po x=0.0, y=-300.0 &
&po x=300.0, y=-300.0 &
&po x=300.0, y=300.0 &
&po x=0.0, y=300.0 &
&po x=0.0, y=-300.0 &
; the coil
&reg mat=1, cur=1.0e6 &
&po x=100.0, y=-100.0 &
&po x=104.0, y=-100.0 &
&po x=104.0, y=100.0 &
&po x=100.0, y=100.0 &
&po x=100.0, y=-100.0 &
"""


def exact_axis_bz(z, a1, a2, L, J):
    def F(u):
        return u * np.log((a2 + np.sqrt(a2 ** 2 + u ** 2)) / (a1 + np.sqrt(a1 ** 2 + u ** 2)))
    return 0.5 * MU0 * J * (F(z + L / 2) - F(z - L / 2))


def test_aircore_solenoid():
    try:
        import gmsh  # noqa: F401
    except ImportError:
        print("gmsh not available: air-core test skipped")
        return
    from rzfield.driver import solve_deck
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "aircore.am")
        with open(path, "w") as f:
            f.write(DECK)
        sol, deck = solve_deck(path, mode="open", open_factor=40.0, log=lambda *a, **k: None)
    a1, a2, L = 1.00, 1.04, 2.00                              # m
    J = 1.0e6 / ((a2 - a1) * L)                               # A/m^2
    z = np.linspace(-3.0, 3.0, 121)
    Br, Bz = sol.B(np.zeros_like(z), z)
    ex = exact_axis_bz(z, a1, a2, L, J)
    rel = np.max(np.abs(Bz - ex)) / np.max(np.abs(ex))
    print("air-core solenoid: Bz(0,0) %.6f T (exact %.6f), max |dBz| / Bz(0,0) on the axis = %.2e, max |Br| = %.1e T"
          % (Bz[60], ex[60], rel, np.max(np.abs(Br))))
    assert rel < 5e-5      # 1e-5 with the default mesh (1 cm elements along the axis)
    assert np.max(np.abs(Br)) == 0.0
    # energy = 1/2 int J A dV over the coil (linear problem)
    W = sol.info["energy_J"]
    rr, zz = np.meshgrid(np.linspace(a1, a2, 41), np.linspace(-L / 2, L / 2, 401), indexing="ij")
    e, l1, l2 = sol.locate(rr.ravel(), zz.ravel())
    from rzfield.fem import p2_basis
    N, _, _ = p2_basis(l1, l2)
    A = np.sum(sol.A[sol.tri[e]] * N, axis=1).reshape(rr.shape)
    integrand = 0.5 * J * A * 2 * np.pi * rr
    W2 = np.trapezoid(np.trapezoid(integrand, zz[0], axis=1), rr[:, 0]) if hasattr(np, "trapezoid") else \
        np.trapz(np.trapz(integrand, zz[0], axis=1), rr[:, 0])
    print("air-core solenoid: energy %.6f kJ, 1/2 int J A dV %.6f kJ (rel %.1e)" % (W / 1e3, W2 / 1e3, W / W2 - 1))
    assert abs(W / W2 - 1) < 1e-3


def test_bh_poisson1010():
    m = BHCurve.builtin("poisson1010")
    # table points are reproduced exactly: B = 1.89 kG at H = 1 Oe, B = 27.4772 kG at H = 6918.319 Oe
    for b_g, h_oe in ((1890.0, 1.0), (17336.2, 100.0), (27477.2, 6918.319)):
        assert abs(m.H(np.array([b_g * GAUSS]))[0] / OERSTED / h_oe - 1) < 1e-6
    # saturation beyond the table: dH/dB = 1/mu0
    b = np.array([3.0, 3.1])
    h = m.H(b)
    assert abs((h[1] - h[0]) / 0.1 * MU0 - 1) < 1e-9
    # monotonic H(B), positive differential reluctivity, energy = int H dB
    bb = np.linspace(0, 3.2, 20001)
    hh = m.H(bb)
    assert np.all(np.diff(hh) > 0)
    w = m.energy_density(bb)
    wnum = np.concatenate([[0], np.cumsum(0.5 * (hh[1:] + hh[:-1]) * np.diff(bb))])
    assert np.max(np.abs(w - wnum)) / w[-1] < 1e-6
    nu, dnub = m.nu(bb ** 2)
    assert np.all(nu + dnub * bb ** 2 > 0)
    print("B-H 1010 table: ok (mu_r max %.0f at %.2f T)" % (np.max(1 / (MU0 * nu)), bb[np.argmax(1 / nu)]))


def test_gemc_roundtrip():
    r = grid_axis(0, 5, 1)
    z = grid_axis(-3, 4, 1)
    Br = np.random.default_rng(1).normal(size=(len(r), len(z)))
    Bz = 1e4 * np.random.default_rng(2).normal(size=(len(r), len(z)))
    with tempfile.TemporaryDirectory() as td:
        fn = os.path.join(td, "m.dat")
        write_map(fn, "solenoid_test", r, z, Br, Bz)
        txt = open(fn).read().splitlines()
        d = read_map(fn)
    assert txt[5] == '                       <first  name="transverse" npoints="6"   min="0"  max="5"  units="cm"/>'
    assert txt[6] == '                       <second name="longitudinal"   npoints="8"  min="-3" max="4"  units="cm"/>'
    assert txt[13].split("\t")[:2] == ["0.00000", "-3.00000"]
    assert np.allclose(d["Br"], Br, rtol=1e-6, atol=1e-12) and np.allclose(d["Bz"], Bz, rtol=1e-6)
    # a free-text comment (solve_rz.py --map-comment) must not break the header attribute
    assert 'comment="rzfield \'open\' (test) + x"' in header("solenoid_test", r, z, comment='rzfield "open" <test> & x')
    print("GEMC map round trip: ok")


IRON_DECK = """Coil inside an iron cylinder: faces r = 30 and z = 25 on the 1 cm grid, r = 40.6 and z = -25.3 off it
&reg kprob=0, icylin=1, xjfact=1.0, mode=0, dx=1.0 &
&po x=0.0, y=-150.0 &
&po x=150.0, y=-150.0 &
&po x=150.0, y=150.0 &
&po x=0.0, y=150.0 &
&po x=0.0, y=-150.0 &
&reg mat=1, cur=2.0e5 &
&po x=20.0, y=-15.0 &
&po x=22.0, y=-15.0 &
&po x=22.0, y=15.0 &
&po x=20.0, y=15.0 &
&po x=20.0, y=-15.0 &
&reg mat=2 &
&po x=30.0, y=-25.3 &
&po x=40.6, y=-25.3 &
&po x=40.6, y=25.0 &
&po x=30.0, y=25.0 &
&po x=30.0, y=-25.3 &
"""


def test_map_faces_air_side():
    try:
        import gmsh  # noqa: F401
    except ImportError:
        print("gmsh not available: iron-face map test skipped")
        return
    from rzfield.driver import solve_deck
    from rzfield.gemcmap import map_from_solution
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "ring.am")
        with open(path, "w") as f:
            f.write(IRON_DECK)
        sol, deck = solve_deck(path, mode="open", refine_box=[0.0, 60.0, -50.0, 50.0], log=lambda *a, **k: None)
    r = grid_axis(0, 60, 1)
    z = grid_axis(-50, 50, 1)
    R, Z = np.meshgrid(r, z, indexing="ij")
    _, _, in_air = sol.B_cm_gauss(R, Z, return_air=True)
    inside = (R > 30.5) & (R < 40.5) & (Z > -25.0) & (Z < 24.5)
    on_face = ((R == 30.0) & (Z >= -25.3) & (Z <= 25.0)) | ((Z == 25.0) & (R >= 30.0) & (R <= 40.6))
    assert np.all(in_air[on_face]) and not np.any(in_air[inside])          # on a face: air side
    maps = {}
    for air_side in (False, True):
        info = {}
        maps[air_side] = map_from_solution(sol, r, z, air_side=air_side, info=info)
    assert info["n_air_side_nodes"] > 0
    # random air points within 1 cm of the faces between grid lines (r = 40.6, z = -25.3), against the solution
    rng = np.random.default_rng(3)
    pr = rng.uniform(29.0, 41.6, 40000)
    pz = rng.uniform(-26.3, 26.0, 40000)
    br, bz, air = sol.B_cm_gauss(pr, pz, return_air=True)
    near = air & ((pr > 40.6) | (pz < -25.3))
    err = {}
    for air_side, (mBr, mBz) in maps.items():
        u, v = pr - r[0], pz - z[0]
        i, j = np.minimum(u.astype(int), len(r) - 2), np.minimum(v.astype(int), len(z) - 2)
        a, b = u - i, v - j
        ib = [M[i, j] * (1 - a) * (1 - b) + M[i, j + 1] * (1 - a) * b + M[i + 1, j] * a * (1 - b)
              + M[i + 1, j + 1] * a * b for M in (mBr, mBz)]
        err[air_side] = np.hypot(ib[0] - br, ib[1] - bz)[near]
    B = np.hypot(br, bz)[near]
    p90 = {k: np.percentile(v, 90) for k, v in err.items()}
    print("iron cylinder: on-face nodes in air; air within 1 cm of the off-grid faces, map vs solution: default median "
          "%.1f G, 90%% %.0f G, max %.0f G; air-side median %.2f G, 90%% %.1f G (field there %.0f-%.0f G, %d air-side "
          "nodes)" % (np.median(err[False]), p90[False], err[False].max(), np.median(err[True]), p90[True], B.min(),
                      B.max(), info["n_air_side_nodes"]))
    assert np.median(err[True]) < 0.05 * np.median(err[False])
    assert p90[True] < 0.1 * p90[False]


def test_portable_paths():
    import extract_rz_model as x
    with tempfile.TemporaryDirectory() as td:
        repo = os.path.join(td, "checkout")
        real = os.path.join(td, "software", "dd4hep-1.38")          # the DD4hep package behind the view
        view = os.path.join(td, "local")                             # $DD4hepINSTALL = a view of symbolic links
        for d in (os.path.join(repo, "compact", "magnet"), os.path.join(real, "DDDetectors"), view):
            os.makedirs(d)
        os.symlink(os.path.join(real, "DDDetectors"), os.path.join(view, "DDDetectors"))
        mine = os.path.join(repo, "compact", "magnet", "solenoid_v4.xml")
        for f in (mine, os.path.join(real, "DDDetectors", "e.xml"), os.path.join(td, "other.xml")):
            open(f, "w").close()
        old = os.environ.get("DD4hepINSTALL")
        os.environ["DD4hepINSTALL"] = view
        try:
            assert x.portable_path(mine, repo) == "compact/magnet/solenoid_v4.xml"
            assert x.portable_path(os.path.join(view, "DDDetectors", "e.xml"), repo) == \
                "${DD4hepINSTALL}/DDDetectors/e.xml"
            assert x.portable_path(os.path.join(real, "DDDetectors", "e.xml"), repo) == ".../e.xml"
            assert x.portable_path(os.path.join(td, "other.xml"), repo) == ".../other.xml"
        finally:
            if old is None:
                os.environ.pop("DD4hepINSTALL")
            else:
                os.environ["DD4hepINSTALL"] = old
        top = os.path.join(repo, "top.xml")
        txt = ('<lccdd>\n  <info name="t" title="magnet + field_CLEOv9.xml"><comment>x</comment></info>\n'
               '  <include ref="${DETECTOR_PATH}/compact/fields/field_CLEOv9.xml"/>\n'
               '  <!-- the magnet -->\n  <constant name="r1" value="176.53*cm"/>\n</lccdd>\n')
        digests = []
        for t in (txt, txt.replace("CLEOv9.xml", "SoLID_v4c_2D.xml").replace("the magnet", "yoke v4c"),
                  txt.replace("176.53", "176.63")):
            with open(top, "w") as f:
                f.write(t)
            digests.append((x.geometry_digest(top, no_field=True), x.geometry_digest(top)))
        assert digests[0][0] == digests[1][0] != digests[2][0]      # field include + comment + title: same geometry
        assert digests[0][1] != digests[1][1]                       # read with its field: the include counts
    print("provenance paths and geometry digests: ok")


if __name__ == "__main__":
    test_bh_poisson1010()
    test_gemc_roundtrip()
    test_portable_paths()
    test_aircore_solenoid()
    test_map_faces_air_side()
