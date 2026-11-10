#!/usr/bin/env python3
"""evdisp -- prototype SoLID event display on top of solid_lineview (feasibility study, not a product).

Inputs
  geometry : Geant4 DAWNFILE .prim (geometry only) -> solid_lineview.parse_prim (cached as .pkl)
  event    : <run>_event.npz from extract_event.py (EDM4hep MCParticles / SimTrackerHits / SimCalorimeterHits
             + contributions, cellIDs decoded with the DD4hep readout, module placement map)
  steps    : steps.bin from the SolidStepRecorder DDG4 stepping action (trajectory points with global time)
             or --traj edm4hep: pseudo-trajectories from MCParticle vertex + SimTrackerHits (Hermite spline)
Outputs
  static   : SVG / PDF / PNG (solid_lineview writers)
  frames   : PNG frames (pycairo, cached geometry layers) -> GIF / WebP (Pillow), WebM (GStreamer) via encode.py
  svg-anim : one self-contained animated SVG (SMIL), tracks grow, hits and modules light up
Views
  side     : section through the beam axis at the azimuth of the calorimeter energy centroid; tracks and hits in
             rho-z projection (--proj ortho: orthographic), hit modules off the plane translucent (rotated in)
  sidezoom : the same around the hit calorimeter
  beam     : x-y section 2 cm behind the front of the hit EC; tracks/hits upstream of it, hits as rings
  cutaway  : 3-D hidden-line view, cut at the event plane, tracks hidden behind solids (depth buffer)
solid_lineview is imported from $SOLID_LINEVIEW or scripts/solid_lineview_v03.py (snapshot of bin/solid_lineview).
"""
import argparse
import math
import os
import pickle
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np

def _load_lineview():
    """solid_lineview as a module: $SOLID_LINEVIEW (path) or the snapshot of the repository's bin/solid_lineview
    (v0.3.0) next to this script.  Registered as 'solid_lineview' so that pickled scenes stay loadable."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.environ.get("SOLID_LINEVIEW") or os.path.join(here, "solid_lineview_v03.py")
    spec = importlib.util.spec_from_file_location("solid_lineview", path,
                                                  loader=importlib.machinery.SourceFileLoader("solid_lineview", path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["solid_lineview"] = mod
    spec.loader.exec_module(mod)
    return mod


import importlib.machinery  # noqa: E402
lv = _load_lineview()

STEP_DTYPE = np.dtype([("event", "<i4"), ("track", "<i4"), ("parent", "<i4"), ("pdg", "<i4"), ("flags", "<i4"),
                       ("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("t", "<f4"), ("ekin", "<f4"), ("edep", "<f4"),
                       ("vol", "<i4")])


def hex2rgb(h):
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (1, 3, 5))


# categorical slots 1-3 of the validated reference palette + neutral inks (dataviz skill)
STYLE = {
    "e-": dict(rgb=hex2rgb("#2a78d6"), dash=None, label="e−"),
    "e+": dict(rgb=hex2rgb("#eb6834"), dash=None, label="e+"),
    "gamma": dict(rgb=hex2rgb("#1baf7a"), dash=(2.0, 1.6), label="γ"),
    "hadron": dict(rgb=hex2rgb("#52514e"), dash=None, label="hadron / μ"),
    "neutral": dict(rgb=hex2rgb("#898781"), dash=(0.6, 1.4), label="n / other neutral"),
}
INK = hex2rgb("#0b0b0b")
INK2 = hex2rgb("#52514e")
# one-hue sequential ramp for energy (ColorBrewer Reds 9; the blue default ramp collides with e- and the
# steel-blue EC module fill)
REDS = ["#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"]
REDS_RGB = np.array([hex2rgb(h) for h in REDS])


def ptype(pdg):
    if pdg == 11:
        return "e-"
    if pdg == -11:
        return "e+"
    if pdg == 22:
        return "gamma"
    if abs(pdg) in (13, 211, 321, 2212) or (abs(pdg) > 1000000000):
        return "hadron"
    return "neutral"


def nice_ceil(x):
    p = 10.0 ** math.floor(math.log10(x))
    return next(m * p for m in (1, 2, 5, 10) if m * p >= x * 0.999)


def energy_color(E, emin, emax):
    """log scale, E in MeV -> rgb (from the 2nd ramp step so that tiny deposits stay visible)."""
    u = (math.log10(max(E, emin)) - math.log10(emin)) / (math.log10(emax) - math.log10(emin))
    u = min(max(u, 0.0), 1.0)
    x = 1 + u * (len(REDS_RGB) - 2)
    i = min(int(x), len(REDS_RGB) - 2)
    f = x - i
    return tuple(REDS_RGB[i] * (1 - f) + REDS_RGB[i + 1] * f)


# =====================================================================================================
#  inputs
# =====================================================================================================

def load_geometry(prim, cache=None):
    t0 = time.perf_counter()
    if cache and os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(prim):
        with open(cache, "rb") as fh:
            sc = pickle.load(fh)
        how = "cache"
    else:
        sc = lv.parse_prim(prim)
        how = "parse"
        if cache:
            with open(cache, "wb") as fh:
                pickle.dump(sc, fh, protocol=5)
    lv.log("geometry (%s) %.2f s: %d solids, %d facets" % (how, time.perf_counter() - t0, len(sc.solids), len(sc.F)))
    return sc


class Track:
    __slots__ = ("tid", "parent", "pdg", "P", "T", "ekin0", "kind")

    def __init__(self, tid, parent, pdg, P, T, ekin0):
        self.tid, self.parent, self.pdg, self.P, self.T, self.ekin0 = tid, parent, pdg, P, T, ekin0
        self.kind = ptype(pdg)


def load_steps(path, event):
    s = np.fromfile(path, STEP_DTYPE)
    s = s[s["event"] == event]
    order = np.lexsort((np.arange(len(s)), s["track"]))   # stable per track
    s = s[order]
    tracks = []
    bounds = np.nonzero(np.diff(s["track"]))[0] + 1
    for g in np.split(np.arange(len(s)), bounds):
        if len(g) == 0:
            continue
        r = s[g]
        P = np.stack([r["x"], r["y"], r["z"]], 1).astype(float)
        T = r["t"].astype(float)
        o = np.argsort(T, kind="stable")
        tracks.append(Track(int(r["track"][0]), int(r["parent"][0]), int(r["pdg"][0]), P[o], T[o],
                            float(r["ekin"][0])))
    return tracks


def load_event(npz, event):
    d = np.load(npz)
    ev = {"collections": [str(c) for c in d["collections"]]}
    for pre in ("mcp_", "th_", "ch_", "cc_"):
        m = d[pre + "event"] == event
        idx = np.nonzero(m)[0]
        for k in d.files:
            if k.startswith(pre):
                ev[k] = d[k][idx]
        if pre == "ch_":
            remap = -np.ones(len(m), np.int64)
            remap[idx] = np.arange(len(idx))
    ev["cc_hit"] = remap[ev["cc_hit"]]                      # contributions -> row in this event's ch_*
    for k in d.files:
        if k.startswith("mod_"):
            ev[k] = d[k]
    return ev


def module_solids(sc, ev):
    """(system, module) -> prim solid index, via the DD4hep placement name (exact).  Cross-check: the centre of
    the drawn polyhedron (bounding-box centre of the hexagonal prism) must coincide with the DD4hep placement."""
    by_name = {}
    for s in sc.solids:
        by_name[s.name.rsplit(":", 1)[0]] = s.index
    out = {}
    dmax, nmiss = 0.0, 0
    for sysid, mod, node, cx, cy, cz in zip(ev["mod_system"], ev["mod_module"], ev["mod_node"],
                                           ev["mod_cx"], ev["mod_cy"], ev["mod_cz"]):
        sid = by_name.get(str(node))
        if sid is None:
            nmiss += 1
            continue
        out[(int(sysid), int(mod))] = sid
        s = sc.solids[sid]
        dmax = max(dmax, float(np.linalg.norm(0.5 * (s.bmin + s.bmax) - np.array([cx, cy, cz]))))
    if not module_solids.reported:
        lv.log("module map: %d of %d DD4hep module placements matched to drawn solids by name, %d unmatched, "
               "max |polyhedron centre - DD4hep placement| = %.3g mm" % (len(out), len(ev["mod_node"]), nmiss, dmax))
        module_solids.reported = True
    return out


module_solids.reported = False


def rhoz_points(P, phi, sign=None):
    """rho-z projection into the plane through the beam axis at azimuth phi: every point keeps its radius and z and
    is rotated about the axis onto the half plane phi (sign +1) or phi + 180 deg (sign -1).  sign None: per point,
    from its own azimuth; a track uses one sign for all its points (no jumps between the half planes)."""
    P = np.asarray(P, float)
    if len(P) == 0:
        return P
    rho = np.hypot(P[:, 0], P[:, 1])
    if sign is None:
        sign = np.where(P[:, 0] * math.cos(phi) + P[:, 1] * math.sin(phi) >= 0, 1.0, -1.0)
    return np.stack([sign * rho * math.cos(phi), sign * rho * math.sin(phi), P[:, 2]], 1)


def rhoz_track(t, phi):
    """the side of a track: that of its radius-weighted mean direction (the outer part decides)."""
    proj = t.P[:, 0] * math.cos(phi) + t.P[:, 1] * math.sin(phi)
    sign = 1.0 if proj.sum() >= 0 else -1.0
    u = Track(t.tid, t.parent, t.pdg, rhoz_points(t.P, phi, sign), t.T, t.ekin0)
    return u


def module_deposits(ev, modmap):
    """prim solid -> (sorted contribution times, cumulative energies [MeV]); total = last value."""
    hit_sid = np.array([modmap.get((int(s), int(m)), -1) for s, m in zip(ev["ch_system"], ev["ch_module"])],
                       np.int64)
    sid = hit_sid[ev["cc_hit"]]
    t = ev["cc_time"]
    e = ev["cc_energy"] * 1000.0
    out = {}
    for s in np.unique(sid):
        if s < 0:
            continue
        m = sid == s
        o = np.argsort(t[m])
        out[int(s)] = (t[m][o], np.cumsum(e[m][o]))
    nmiss = int((hit_sid < 0).sum())
    if nmiss:
        lv.log("warning: %d calorimeter hits without a drawn module" % nmiss)
    return out


def tracks_from_edm4hep(ev, nseg=12):
    """Route B: pseudo-trajectories without a side file: MCParticle vertex -> SimTrackerHits (time ordered)
    -> endpoint, joined by cubic Hermite segments whose tangents are the momentum directions.
    Times: hit times; the endpoint time is extrapolated with v = c*beta along the chord."""
    c_mm_ns = 299.792458
    hits = defaultdict(list)
    for i, p in enumerate(ev["th_particle"]):
        hits[int(p)].append(i)
    tracks = []
    for j in range(len(ev["mcp_index"])):
        idx = int(ev["mcp_index"][j])
        mom = np.array([ev["mcp_px"][j], ev["mcp_py"][j], ev["mcp_pz"][j]])
        pm = float(np.linalg.norm(mom))
        if idx not in hits and ev["mcp_gen"][j] != 1:
            continue
        E = float(ev["mcp_energy"][j])
        beta = pm / E if E > 0 else 1.0
        pts = [np.array([ev["mcp_vx"][j], ev["mcp_vy"][j], ev["mcp_vz"][j]])]
        tan = [mom / max(pm, 1e-12)]
        ts = [float(ev["mcp_time"][j])]
        plast = pm
        for i in sorted(hits.get(idx, []), key=lambda i: ev["th_time"][i]):
            hm = np.array([ev["th_px"][i], ev["th_py"][i], ev["th_pz"][i]])
            hp = float(np.linalg.norm(hm))
            # hits of untracked secondaries are attributed to the kept ancestor (e.g. shower back-splash
            # into the SPD): drop hits whose momentum falls far below the particle's
            if hp < 0.3 * plast:
                continue
            plast = hp
            pts.append(np.array([ev["th_x"][i], ev["th_y"][i], ev["th_z"][i]]))
            tan.append(hm / max(hp, 1e-12))
            ts.append(float(ev["th_time"][i]))
        endp = np.array([ev["mcp_ex"][j], ev["mcp_ey"][j], ev["mcp_ez"][j]])
        if np.linalg.norm(endp - pts[-1]) > 1.0:
            me = np.array([ev["mcp_epx"][j], ev["mcp_epy"][j], ev["mcp_epz"][j]])
            pts.append(endp)
            tan.append(me / np.linalg.norm(me) if np.linalg.norm(me) > 0 else tan[-1])
            ts.append(ts[-1] + np.linalg.norm(endp - pts[-2]) / (c_mm_ns * max(beta, 1e-3)))
        P, T = [pts[0]], [ts[0]]
        for a in range(len(pts) - 1):
            p0, p1 = pts[a], pts[a + 1]
            L = float(np.linalg.norm(p1 - p0))
            m0, m1 = tan[a] * L, tan[a + 1] * L
            s = np.linspace(0, 1, nseg + 1)[1:, None]
            h00, h10, h01, h11 = 2 * s**3 - 3 * s**2 + 1, s**3 - 2 * s**2 + s, -2 * s**3 + 3 * s**2, s**3 - s**2
            P.extend(h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1)
            T.extend(ts[a] + s[:, 0] * (ts[a + 1] - ts[a]))
        parent = 0 if ev["mcp_gen"][j] == 1 else int(ev["mcp_parent"][j]) + 1    # 0 = primary, as in Geant4
        tracks.append(Track(idx, parent, int(ev["mcp_pdg"][j]), np.array(P), np.array(T),
                            1000.0 * (E - float(ev["mcp_mass"][j]))))                # GeV -> MeV
    return tracks


# =====================================================================================================
#  views and geometry layer
# =====================================================================================================

def lv_opts(argv):
    o = lv.build_parser().parse_args(["dummy.prim"] + argv)
    o.width, o.height = [int(v) for v in o.size.lower().split("x")]
    o.lw_outline = o.lw if o.lw_outline is None else o.lw_outline
    o.lw_cut = 1.4 * o.lw if o.lw_cut is None else o.lw_cut
    return o


def make_view(kind, ev, args):
    """side: section through the z axis and the azimuth of the calorimeter energy centroid (event plane);
    beam: x-y section at z (through the hit calorimeter)."""
    E = ev["ch_energy"]
    if len(E):
        cx, cy, cz = (np.average(ev["ch_" + a], weights=E) for a in "xyz")
    else:
        cx, cy, cz = 0.0, 1.0, 4600.0
    if kind in ("side", "sidezoom"):
        phi = math.atan2(cy, cx) if args.phi is None else math.radians(args.phi)
        n = np.array([-math.sin(phi), math.cos(phi), 0.0])
        up = np.array([math.cos(phi), math.sin(phi), 0.0])
        cam = lv.Camera(n, up, (0, 0, 0))
        zoom = None
        if kind == "sidezoom":
            rc = math.hypot(cx, cy)
            hw = args.zoom_half
            zoom = (cz - 1.3 * hw, cz + 1.3 * hw, rc - hw, rc + hw)
        what = ("tracks, hits and off-plane hit modules (translucent) in ρ-z projection" if args.proj == "rhoz"
                else "tracks and hits projected orthographically; off-plane hit modules (translucent) rotated into it")
        return dict(kind=kind, cam=cam, plane=(n, 0.0), phi=phi, zoom=zoom,
                    label="section through the beam axis at φ = %.0f° / %.0f° (event plane); %s" % (
                        math.degrees(phi) % 360, (math.degrees(phi) + 180) % 360, what))
    if kind == "cutaway":
        # 3-D hidden-line view, the half in front of the event plane removed (DAWN/dawncut style)
        phi = math.atan2(cy, cx) if args.phi is None else math.radians(args.phi)
        n = np.array([-math.sin(phi), math.cos(phi), 0.0])
        ur = np.array([math.cos(phi), math.sin(phi), 0.0])
        eye = lv._unit(0.75 * n + 0.48 * ur - 0.45 * np.array([0.0, 0.0, 1.0]))
        cam = lv.Camera(eye, ur, (0.0, 0.0, 700.0), args.perspective)
        keep = (-n, 0.0)                     # keep -n.P >= 0: the half behind the event plane
        return dict(kind=kind, cam=cam, plane=keep, cut=keep, phi=phi,
                    label="cut-away at the event plane (φ = %.0f°), hidden lines removed, perspective"
                    % (math.degrees(phi) % 360))
    # default: 2 cm behind the front face of the FAEC / LAEC shower modules (solid_JPsi_lite z, prototype)
    zsec = args.zsec if args.zsec is not None else (4463.0 if cz > 2000 else -557.0)
    cam = lv.Camera((0, 0, 1), (0, 1, 0), (0, 0, zsec))
    return dict(kind=kind, cam=cam, plane=(np.array([0.0, 0.0, 1.0]), -zsec), zsec=zsec,
                label="beam view, section at z = %.0f cm" % (zsec / 10))


def geometry_layer_hlr(sc, view, args, title):
    argv = ["--mode", "hlr", "--size", args.size, "--dpi", str(args.dpi), "--lw", str(args.lw),
            "--font-size", str(args.font_size), "--cap-fill", "color", "--fill-lighten", "0",
            "--perspective", "%g" % args.perspective, "--hlr-scale", "2"]
    if title:
        argv += ["--title", title]
    o = lv_opts(argv)
    saved = [s.color for s in sc.solids]
    for s_ in sc.solids:                         # lighten the cap colours like the section views do
        s_.color = tuple(0.35 + 0.65 * v for v in s_.color)
    try:
        rd = lv.Renderer(sc, o)
        rd._fitinfo = None
        orig_fit = rd._fit

        def fit_and_store(cam_, pts):
            r = orig_fit(cam_, pts)
            rd._fitinfo = r
            return r
        rd._fit = fit_and_store
        rd._cap_planes = [view["cut"]]
        dl = rd.render(view["cam"], [view["cut"]], [view["cut"]], "hlr")
    finally:
        for s_, c in zip(sc.solids, saved):
            s_.color = c
    if hasattr(lv, "draw_title"):            # solid_lineview >= 0.3 draws the title outside _decorations
        lv.draw_title(dl, o)
    return dl, rd, o


def geometry_layer(sc, view, args, title):
    if view["kind"] == "cutaway":
        return geometry_layer_hlr(sc, view, args, title)
    argv = ["--mode", "section", "--fill", "--size", args.size, "--dpi", str(args.dpi), "--scale-bar",
            "--lw", str(args.lw), "--font-size", str(args.font_size)]
    if title:
        argv += ["--title", title]
    if view["kind"] == "beam":
        r = args.beam_rmax
        argv += ["--fit-box", "%g" % -r, "%g" % r, "%g" % -r, "%g" % r]
    if view.get("zoom"):
        argv += ["--fit-box"] + ["%g" % v for v in view["zoom"]]
    o = lv_opts(argv)
    dl, rd = lv.render_section(sc, o, view["cam"], view["plane"], [])
    if hasattr(lv, "draw_title"):
        lv.draw_title(dl, o)
    return dl, rd, o


def section_rings(sc, rd, plane, sids):
    """closed section outlines (canvas px) of the given solids, like render_section does for fills."""
    mesh = rd.mesh
    fsel = np.zeros(len(sc.solids), bool)
    fsel[list(sids)] = True
    tsel = np.nonzero(fsel[sc.fsolid][mesh.tri_face])[0]
    rings = defaultdict(list)
    if len(tsel) == 0:
        return rings
    Q0, Q1, ti, k0, k1 = lv.section_segments(sc.V, mesh.tri[tsel], plane)
    sol = sc.fsolid[mesh.tri_face[tsel[ti]]]
    for s in np.unique(sol):
        m = sol == s
        for pts, closed in lv.chain_segments(k0[m], k1[m], Q0[m], Q1[m]):
            if closed and len(pts) >= 3:
                x, y, _ = rd.to_px(pts)
                rings[int(s)].append(np.stack([x, y], 1))
    return rings


def projected_hull(sc, rd, sid, phi_plane=None):
    """outline of the orthographic projection of a (convex) solid: convex hull of its projected vertices.
    With phi_plane the solid is first rotated about the beam axis into the section plane (rho-z projection:
    it keeps its radius, the upper or lower half chosen by its azimuth)."""
    s = sc.solids[sid]
    V = sc.V[s.v0:s.v1]
    if phi_plane is not None:
        c = V[:, :2].mean(0)
        phm = math.atan2(c[1], c[0])
        tgt = phi_plane if math.cos(phm - phi_plane) >= 0 else phi_plane + math.pi
        d = tgt - phm
        cd, sd = math.cos(d), math.sin(d)
        V = np.stack([cd * V[:, 0] - sd * V[:, 1], sd * V[:, 0] + cd * V[:, 1], V[:, 2]], 1)
    x, y, _ = rd.to_px(V)
    P = np.unique(np.round(np.stack([x, y], 1), 3), axis=0)
    if len(P) < 3:
        return None
    P = P[np.lexsort((P[:, 1], P[:, 0]))]

    def half(pts):
        h = []
        for p in pts:
            while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (p[1] - h[-2][1])
                                   - (h[-1][1] - h[-2][1]) * (p[0] - h[-2][0])) <= 0:
                h.pop()
            h.append(p)
        return h
    lo, hi = half(P), half(P[::-1])
    return np.array(lo[:-1] + hi[:-1])


# =====================================================================================================
#  overlays (into a solid_lineview DrawList)
# =====================================================================================================

def select_tracks(tracks, args):
    keep = []
    for t in tracks:
        if t.parent == 0:
            keep.append(t)
            continue
        emin = args.emin_gamma if t.kind == "gamma" else args.emin_neutral if t.kind == "neutral" else args.emin
        if t.ekin0 >= emin and len(t.P) >= 2:
            keep.append(t)
    return keep


def clip_track_time(P, T, tmax):
    """points with time <= tmax plus the interpolated head point."""
    if tmax is None or T[-1] <= tmax:
        return P
    k = int(np.searchsorted(T, tmax, side="right"))
    if k == 0:
        return P[:0]
    if k < len(T) and T[k] > T[k - 1]:
        f = (tmax - T[k - 1]) / (T[k] - T[k - 1])
        return np.vstack([P[:k], P[k - 1] + f * (P[k] - P[k - 1])])
    return P[:k]


def track_items(tracks, to_px, args, tmax=None):
    """[(kind, is_primary, polyline px)] in drawing order (photons first, primaries last)."""
    order = sorted(tracks, key=lambda t: (t.parent == 0, t.kind != "gamma", t.ekin0))
    out = []
    for t in order:
        P = clip_track_time(t.P, t.T, tmax)
        if len(P) < 2:
            continue
        x, y, _ = to_px(P)
        out.append((t.kind, t.parent == 0, np.stack([x, y], 1)))
    return out


def track_items_hlr(tracks, rd, args, tmax=None):
    """like track_items, but only the parts not hidden by the (cut) geometry, using the renderer's depth buffer."""
    hs = rd.o.hlr_scale
    order = sorted(tracks, key=lambda t: (t.parent == 0, t.kind != "gamma", t.ekin0))
    out = []
    # the cut faces (caps) must not hide the event: push them to infinite depth for this test
    cap = rd.zface == -2
    saved = rd.zA.copy(), rd.zB.copy(), rd.zC.copy()
    rd.zA[cap], rd.zB[cap], rd.zC[cap] = 0.0, 0.0, 1e30
    try:
        runs = []
        for t in order:
            P = clip_track_time(t.P, t.T, tmax)
            if len(P) < 2:
                continue
            xa, ya, za = rd.to_px(P[:-1])
            xb, yb, zb = rd.to_px(P[1:])
            m1 = np.full(len(xa), -9, np.int64)
            runs.append((t, xa, ya, xb, yb, rd._visible_runs(xa * hs, ya * hs, za, xb * hs, yb * hs, zb, m1, m1)))
    finally:
        rd.zA, rd.zB, rd.zC = saved
    for t, xa, ya, xb, yb, (e, t0, t1) in runs:
        if len(e) == 0:
            continue
        X0 = xa[e] + t0 * (xb[e] - xa[e]); Y0 = ya[e] + t0 * (yb[e] - ya[e])
        X1 = xa[e] + t1 * (xb[e] - xa[e]); Y1 = ya[e] + t1 * (yb[e] - ya[e])
        for pl in lv.chain_2d(X0, Y0, X1, Y1):
            out.append((t.kind, t.parent == 0, pl))
    return out


def add_tracks(dl, items, args):
    groups = defaultdict(list)
    for kind, prim, pl in items:
        groups[(kind, prim)].append(pl)
    for (kind, prim) in sorted(groups, key=lambda k: (k[1], k[0] != "gamma")):
        st = STYLE[kind]
        lw = args.traj_lw * (2.6 if prim else 1.0)
        dl.lines(groups[(kind, prim)], st["rgb"], lw, dash=st["dash"] if not prim else None)


def add_hits(dl, hitP, hitT, to_px, args, tmax=None, ring=False):
    """dots (white halo) or, in the beam view where they sit on coloured modules, hollow rings."""
    pt = dl.dpi / 72.0
    m = np.ones(len(hitT), bool) if tmax is None else hitT <= tmax
    if not m.any():
        return
    x, y, _ = to_px(hitP[m])
    r = args.hit_size * pt
    if ring:
        a = np.linspace(0, 2 * math.pi, 17)
        dl.lines([np.stack([xi + r * np.cos(a), yi + r * np.sin(a)], 1) for xi, yi in zip(x, y)], INK, 0.7)
        return
    for xi, yi in zip(x, y):
        dl.marker("circle", float(xi), float(yi), r + 0.6 * pt, (1, 1, 1))
        dl.marker("circle", float(xi), float(yi), r, INK)


def clip_track_z(t, zmax):
    """the part of a track upstream of z = zmax (cut at its first crossing, time interpolated)."""
    above = np.nonzero(t.P[:, 2] > zmax)[0]
    if len(above) == 0:
        return t
    k = int(above[0])
    if k == 0:
        return None
    f = (zmax - t.P[k - 1, 2]) / (t.P[k, 2] - t.P[k - 1, 2])
    P = np.vstack([t.P[:k], t.P[k - 1] + f * (t.P[k] - t.P[k - 1])])
    T = np.append(t.T[:k], t.T[k - 1] + f * (t.T[k] - t.T[k - 1]))
    return Track(t.tid, t.parent, t.pdg, P, T, t.ekin0)


def add_colorbar_legend(dl, args, emin, emax, x0, y0, with_legend=True, ring=False):
    """vertical log colour bar (MeV) + particle legend in the band right of the drawing."""
    pt = dl.dpi / 72.0
    fs = args.font_size
    bw, bh = 12 * pt, 150 * pt
    dl.text(x0, y0 - 6 * pt, "module energy deposit", fs, INK)
    n = 64
    for i in range(n):
        u0, u1 = i / n, (i + 1) / n
        E = 10 ** (math.log10(emin) + (u0 + u1) / 2 * (math.log10(emax) - math.log10(emin)))
        ya, yb = y0 + bh * (1 - u1), y0 + bh * (1 - u0)
        dl.fill([np.array([[x0, ya], [x0 + bw, ya], [x0 + bw, yb + 0.5], [x0, yb + 0.5]])],
                energy_color(E, emin, emax))
    dl.lines([np.array([[x0, y0], [x0 + bw, y0], [x0 + bw, y0 + bh], [x0, y0 + bh], [x0, y0]])], INK2, 0.5)
    for p in range(int(math.floor(math.log10(emin))), int(math.ceil(math.log10(emax))) + 1):
        for mlt in (1, 2, 5):
            E = mlt * 10.0 ** p
            if E < emin * 0.999 or E > emax * 1.001:
                continue
            u = (math.log10(E) - math.log10(emin)) / (math.log10(emax) - math.log10(emin))
            yy = y0 + bh * (1 - u)
            dl.lines([np.array([[x0 + bw, yy], [x0 + bw + 3 * pt, yy]])], INK2, 0.5)
            if mlt == 1 or (mlt == 5 and math.log10(emax / emin) < 2.5) or abs(E - emax) < 1e-6 * emax:
                dl.text(x0 + bw + 5 * pt, yy + 0.35 * fs * pt, ("%g MeV" % E), fs * 0.9, INK2)
    if not with_legend:
        return
    yl = y0 + bh + 28 * pt
    dl.text(x0, yl - 6 * pt, "particles", fs, INK)
    for i, k in enumerate(["e-", "e+", "gamma", "hadron", "neutral"]):
        st = STYLE[k]
        yy = yl + (i + 0.5) * 13 * pt
        dl.lines([np.array([[x0, yy], [x0 + 22 * pt, yy]])], st["rgb"], 1.2, dash=st["dash"])
        dl.text(x0 + 28 * pt, yy + 0.35 * fs * pt, st["label"], fs * 0.9, INK)
    yy = yl + 5.5 * 13 * pt
    if ring:
        a = np.linspace(0, 2 * math.pi, 17)
        r = args.hit_size * pt
        dl.lines([np.stack([x0 + 11 * pt + r * np.cos(a), yy + r * np.sin(a)], 1)], INK, 0.7)
    else:
        dl.marker("circle", x0 + 11 * pt, yy, args.hit_size * pt, INK)
    dl.text(x0 + 28 * pt, yy + 0.35 * fs * pt, "tracker hit", fs * 0.9, INK)
    yy += 13 * pt
    dl.lines([np.array([[x0, yy], [x0 + 22 * pt, yy]])], STYLE["e-"]["rgb"], args.traj_lw * 2.6)
    dl.text(x0 + 28 * pt, yy + 0.35 * fs * pt, "primary", fs * 0.9, INK)


class Scene:
    """Static layers of one view + everything needed to draw the event at any time t."""

    def __init__(self, sc, ev, tracks, view, args, title=None):
        t0 = time.perf_counter()
        self.sc, self.ev, self.view, self.args = sc, ev, view, args
        dl, rd, o = geometry_layer(sc, view, args, title)
        self.t_geometry = time.perf_counter() - t0
        self.rd, self.to_px = rd, rd.to_px
        pt = dl.dpi / 72.0
        self.band = int(150 * pt)                      # colour bar / legend band on the right
        self.W, self.H = dl.W + self.band, max(dl.H, int(300 * pt))
        self.dpi = dl.dpi
        first_line = next((i for i, it in enumerate(dl.items) if it[0] != "fill"), len(dl.items))
        self.fill_items = dl.items[:first_line]
        self.line_items = dl.items[first_line:]
        modmap = module_solids(sc, ev)
        self.deposits = module_deposits(ev, modmap)
        totals = [v[1][-1] for v in self.deposits.values()]
        self.emax = nice_ceil(max(totals)) if totals else 100.0
        self.emin = args.emin_cal
        sids = list(self.deposits)
        self.rings = section_rings(sc, rd, view["plane"], sids)
        self.ghosts = {}
        if view["kind"] in ("side", "sidezoom") and args.project_modules:
            for s in sids:
                if s not in self.rings:
                    h = projected_hull(sc, rd, s, view["phi"])
                    if h is not None:
                        self.ghosts[s] = h
        self.tracks = select_tracks(tracks, args)
        self.hitP = np.stack([ev["th_x"], ev["th_y"], ev["th_z"]], 1)
        self.hitT = ev["th_time"]
        if view["kind"] in ("side", "sidezoom") and args.proj == "rhoz":
            self.tracks = [rhoz_track(t, view["phi"]) for t in self.tracks]
            self.hitP = rhoz_points(self.hitP, view["phi"])
        self.ring_hits = view["kind"] == "beam"
        if view["kind"] == "beam":
            # section semantics: only what lies upstream of the section plane (the shower behind it would
            # cover the module colours)
            self.tracks = [u for u in (clip_track_z(t, view["zsec"]) for t in self.tracks) if u is not None]
            m = self.hitP[:, 2] <= view["zsec"]
            self.hitP, self.hitT = self.hitP[m], self.hitT[m]
        self.title = title
        self.title_band = int(math.ceil(3.0 * o.font_size * o.dpi / 72.0)) if title else 0
        lv.log("view %s: canvas %dx%d, %d hit modules (%d cut by the section, %d projected), %d/%d tracks drawn"
               % (view["kind"], self.W, self.H, len(sids), len(self.rings), len(self.ghosts), len(self.tracks),
                  len(tracks)))

    def energy_at(self, sid, tmax):
        ts, ce = self.deposits[sid]
        if tmax is None:
            return ce[-1]
        k = int(np.searchsorted(ts, tmax, side="right"))
        return ce[k - 1] if k > 0 else 0.0

    def module_fill_items(self, tmax=None):
        dl = lv.DrawList(self.W, self.H, self.dpi, None)
        pt = self.dpi / 72.0
        for s in sorted(self.ghosts, key=lambda s: self.energy_at(s, None)):
            E = self.energy_at(s, tmax)
            if E >= self.emin:
                dl.fill([self.ghosts[s]], energy_color(E, self.emin, self.emax), stroke=INK2, lw_pt=0.25,
                        alpha=0.55)
        for s, rings in self.rings.items():
            E = self.energy_at(s, tmax)
            if E >= self.emin:
                dl.fill(rings, energy_color(E, self.emin, self.emax), stroke=INK, lw_pt=0.3)
        return dl.items

    def overlay_items(self, tmax=None, clock=True):
        dl = lv.DrawList(self.W, self.H, self.dpi, None)
        pt = self.dpi / 72.0
        if self.view.get("zoom"):
            # zoomed views: geometry outside the window reaches into the title band and the legend band
            x1, band = self.W - self.band, self.title_band + int(8 * pt)
            dl.fill([np.array([[0, 0], [x1, 0], [x1, band], [0, band]], float),
                     np.array([[x1, 0], [self.W, 0], [self.W, self.H], [x1, self.H]], float)], (1, 1, 1))
        if self.view["kind"] == "cutaway":
            add_tracks(dl, track_items_hlr(self.tracks, self.rd, self.args, tmax), self.args)
        else:
            add_tracks(dl, track_items(self.tracks, self.to_px, self.args, tmax), self.args)
        add_hits(dl, self.hitP, self.hitT, self.to_px, self.args, tmax, ring=self.ring_hits)
        x1 = self.W - self.band
        dl.fill([np.array([[x1, 0], [self.W, 0], [self.W, self.H], [x1, self.H]], float)], (1, 1, 1))
        add_colorbar_legend(dl, self.args, self.emin, self.emax, self.W - self.band + 14 * pt, 40 * pt,
                            ring=self.ring_hits)
        fs = self.args.font_size
        if self.view.get("zoom") and self.title:
            dl.text((self.W - self.band) / 2, 1.6 * fs * pt, self.title, fs * 1.3, INK, "middle")
        dl.text(10 * pt, 1.6 * fs * pt + 14 * pt, self.view["label"], fs * 0.9, INK2)
        if tmax is not None and clock:
            dl.text(10 * pt, 1.6 * fs * pt + 30 * pt, "t = %5.1f ns" % tmax, fs * 1.4, INK)
        return dl.items

    def drawlist(self, tmax=None):
        dl = lv.DrawList(self.W, self.H, self.dpi, (1, 1, 1))
        dl.items = self.fill_items + self.module_fill_items(tmax) + self.line_items + self.overlay_items(tmax)
        return dl


# =====================================================================================================
#  raster frames with cached geometry layers (pycairo)
# =====================================================================================================

def cairo_draw(cr, items, dpi):
    import cairo
    pt = dpi / 72.0
    cr.set_line_join(cairo.LINE_JOIN_ROUND)
    for it in items:
        if it[0] == "lines":
            _, pls, rgb, lw, dash, cap = it
            cr.set_source_rgb(*rgb)
            cr.set_line_width(lw * pt)
            cr.set_line_cap(cairo.LINE_CAP_ROUND if cap == "round" else cairo.LINE_CAP_BUTT)
            cr.set_dash([v * pt for v in dash] if dash else [])
            for p in pls:
                cr.move_to(*p[0])
                for q in p[1:]:
                    cr.line_to(*q)
            cr.stroke()
        elif it[0] == "fill":
            _, rings, rgb, stroke, lw, alpha = it
            cr.set_fill_rule(cairo.FILL_RULE_EVEN_ODD)
            for r in rings:
                cr.move_to(*r[0])
                for q in r[1:]:
                    cr.line_to(*q)
                cr.close_path()
            cr.set_source_rgba(rgb[0], rgb[1], rgb[2], alpha)
            if stroke is not None and lw > 0:
                cr.fill_preserve()
                cr.set_source_rgb(*stroke)
                cr.set_line_width(lw * pt)
                cr.set_dash([])
                cr.stroke()
            else:
                cr.fill()
        elif it[0] == "marker":
            _, kind, x, y, r, rgb = it
            cr.set_source_rgb(*rgb)
            if kind == "circle":
                cr.new_sub_path()
                cr.arc(x, y, r, 0, 2 * math.pi)
            else:
                cr.rectangle(x - r, y - r, 2 * r, 2 * r)
            cr.fill()
        elif it[0] == "text":
            _, x, y, s, size, rgb, anchor, font = it
            cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            cr.set_font_size(size * pt)
            ext = cr.text_extents(s)
            dx = {"start": 0, "middle": -ext.x_advance / 2, "end": -ext.x_advance}[anchor]
            cr.set_source_rgb(*rgb)
            cr.move_to(x + dx, y)
            cr.show_text(s)


class FrameRenderer:
    def __init__(self, scene, scale=1.0):
        import cairo
        self.cairo = cairo
        self.scene, self.scale = scene, scale
        self.W, self.H = int(scene.W * scale) // 2 * 2, int(scene.H * scale) // 2 * 2   # even for video
        t0 = time.perf_counter()
        self.bg = self._layer(scene.fill_items, opaque=True)
        self.fg = self._layer(scene.line_items, opaque=False)
        self.t_cache = time.perf_counter() - t0

    def _layer(self, items, opaque):
        cairo = self.cairo
        surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.W, self.H)
        cr = cairo.Context(surf)
        cr.scale(self.scale, self.scale)
        if opaque:
            cr.set_source_rgb(1, 1, 1)
            cr.paint()
        cairo_draw(cr, items, self.scene.dpi)
        return surf

    def frame(self, tmax, path):
        cairo = self.cairo
        surf = cairo.ImageSurface(cairo.FORMAT_RGB24, self.W, self.H)
        cr = cairo.Context(surf)
        cr.set_source_surface(self.bg, 0, 0)
        cr.paint()
        cr.save()
        cr.scale(self.scale, self.scale)
        cairo_draw(cr, self.scene.module_fill_items(tmax), self.scene.dpi)
        cr.restore()
        cr.set_source_surface(self.fg, 0, 0)
        cr.paint()
        cr.scale(self.scale, self.scale)
        cairo_draw(cr, self.scene.overlay_items(tmax), self.scene.dpi)
        surf.write_to_png(path)


# =====================================================================================================
#  animated SVG (SMIL): static geometry + growing tracks + hits / modules switching on
# =====================================================================================================

def svg_items(items, dpi):
    out = []
    pt = dpi / 72.0
    for it in items:
        if it[0] == "lines":
            _, pls, rgb, lw, dash, cap = it
            d = "".join("M" + lv._fmt(p[0]) + "L" + " ".join(lv._fmt(q) for q in p[1:]) for p in pls)
            extra = ' stroke-dasharray="%s"' % ",".join("%g" % (v * pt) for v in dash) if dash else ""
            out.append('<path fill="none" stroke="%s" stroke-width="%.3f" stroke-linecap="%s" '
                       'stroke-linejoin="round"%s d="%s"/>' % (lv._rgbhex(rgb), lw * pt, cap, extra, d))
        elif it[0] == "fill":
            _, rings, rgb, stroke, lw, alpha = it
            d = "".join("M" + lv._fmt(r[0]) + "L" + " ".join(lv._fmt(q) for q in r[1:]) + "Z" for r in rings)
            st = ' stroke="%s" stroke-width="%.3f"' % (lv._rgbhex(stroke), lw * pt) \
                if stroke is not None and lw > 0 else ' stroke="none"'
            op = ' fill-opacity="%.3f"' % alpha if alpha < 1 else ""
            out.append('<path fill="%s" fill-rule="evenodd"%s%s d="%s"/>' % (lv._rgbhex(rgb), op, st, d))
        elif it[0] == "marker":
            _, kind, x, y, r, rgb = it
            out.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s"/>' % (x, y, r, lv._rgbhex(rgb)))
        elif it[0] == "text":
            _, x, y, s, size, rgb, anchor, font = it
            s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            out.append('<text x="%.2f" y="%.2f" font-family="Helvetica, Arial, sans-serif" font-size="%.2f" '
                       'fill="%s" text-anchor="%s">%s</text>' % (x, y, size * pt, lv._rgbhex(rgb), anchor, s))
    return out


def write_svg_animation(scene, path, t_end, duration, hold=2.0):
    """tracks: stroke-dashoffset animated with keyTimes following the step times (pathLength = 1);
    hits: opacity set at the hit time; modules: fill colour steps as the cumulative energy grows."""
    k = duration / t_end
    total = duration + hold
    dpi, pt = scene.dpi, scene.dpi / 72.0
    a = scene.args
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%.1fpt" height="%.1fpt" '
           'viewBox="0 0 %d %d" data-t-end-ns="%.4f" data-duration-s="%.4f" data-loop-s="%.4f">'
           % (scene.W * 72 / dpi, scene.H * 72 / dpi, scene.W, scene.H, t_end, duration, duration + hold),
           '<!-- SoLID event animation prototype: %.0f ns in %.0f s, loops; SMIL -->' % (t_end, duration),
           '<rect width="%d" height="%d" fill="#ffffff"/>' % (scene.W, scene.H)]
    out += svg_items(scene.fill_items, dpi)
    # modules: discrete colour steps
    loop = 'dur="%.3fs" repeatCount="indefinite"' % total
    for s, rings in [(s, [h]) for s, h in scene.ghosts.items()] + list(scene.rings.items()):   # ghosts below
        ts, ce = scene.deposits[s]
        if ce[-1] < scene.emin:
            continue
        ghost = s not in scene.rings
        d = "".join("M" + lv._fmt(r[0]) + "L" + " ".join(lv._fmt(q) for q in r[1:]) + "Z" for r in rings)
        # sample the cumulative energy at 12 times between the first and the last contribution
        tq = np.unique(np.concatenate([[0.0], np.linspace(ts[0], min(ts[-1], t_end), 12)]))
        vals, keys = [], []
        for tt in tq[1:]:
            E = scene.energy_at(s, tt)
            if E >= scene.emin:
                vals.append(lv._rgbhex(energy_color(E, scene.emin, scene.emax)))
                keys.append(min(tt * k / total, 1.0))
        if not vals:
            continue
        # hidden (opacity 0) until the first deposit above threshold, then the colour steps with the energy
        vals = [vals[0]] + vals
        keys = [0.0] + keys
        out.append('<path fill="%s" stroke="%s" stroke-width="%.2f"%s fill-rule="evenodd" opacity="0" d="%s">'
                   '<animate attributeName="opacity" calcMode="discrete" values="0;1" keyTimes="0;%.4f" %s/>'
                   '<animate attributeName="fill" calcMode="discrete" values="%s" keyTimes="%s" %s/></path>'
                   % (vals[0], "#52514e" if ghost else "#0b0b0b", (0.25 if ghost else 0.3) * pt,
                      ' fill-opacity="0.55"' if ghost else "", d, keys[1], loop, ";".join(vals),
                      ";".join("%.4f" % x for x in keys), loop))
    out += svg_items(scene.line_items, dpi)
    # tracks
    for t in sorted(scene.tracks, key=lambda t: (t.parent == 0, t.kind != "gamma", t.ekin0)):
        if len(t.P) < 2:
            continue
        x, y, _ = scene.to_px(t.P)
        seg = np.hypot(np.diff(x), np.diff(y))
        L = np.concatenate([[0], np.cumsum(seg)])
        if L[-1] < 0.5:
            continue
        frac = L / L[-1]
        st = STYLE[t.kind]
        prim = t.parent == 0
        lw = a.traj_lw * (2.6 if prim else 1.0) * pt
        tt = np.clip(t.T, 0, t_end)
        # thin the keyTimes to <= 24 entries
        sel = np.unique(np.concatenate([[0], np.linspace(0, len(tt) - 1, min(len(tt), 24)).astype(int)]))
        kt = [0.0] + list(np.clip(tt[sel] * k / total, 0, 1)) + [1.0]
        vals = [1.0] + list(1 - frac[sel]) + [1 - frac[sel][-1]]
        kt = np.maximum.accumulate(np.array(kt))
        dash = ""
        if st["dash"] and not prim:
            dash = ' stroke-dasharray="%s"' % ",".join("%g" % (v * pt) for v in st["dash"])
            # a dashed track cannot also use the dash trick: reveal it with an opacity switch instead
            out.append('<path fill="none" stroke="%s" stroke-width="%.2f"%s opacity="0" d="M%s">'
                       '<animate attributeName="opacity" calcMode="discrete" values="0;1;1" keyTimes="0;%.4f;1" %s/>'
                       '</path>' % (lv._rgbhex(st["rgb"]), lw, dash,
                                    "L".join("%.1f %.1f" % (xi, yi) for xi, yi in zip(x, y)), kt[1], loop))
            continue
        out.append('<path fill="none" stroke="%s" stroke-width="%.2f" stroke-linecap="round" pathLength="1" '
                   'stroke-dasharray="1 1" stroke-dashoffset="1" d="M%s"><animate attributeName="stroke-dashoffset" '
                   'values="%s" keyTimes="%s" %s/></path>'
                   % (lv._rgbhex(st["rgb"]), lw, "L".join("%.1f %.1f" % (xi, yi) for xi, yi in zip(x, y)),
                      ";".join("%.4f" % v for v in vals), ";".join("%.4f" % v for v in kt), loop))
    # hits
    x, y, _ = scene.to_px(scene.hitP)
    r = a.hit_size * pt
    for xi, yi, th in zip(x, y, scene.hitT):
        kt = min(th * k / total, 0.999)
        style = ('fill="none" stroke="#0b0b0b" stroke-width="%.2f"' % (0.7 * pt) if scene.ring_hits else
                 'fill="#0b0b0b" stroke="#ffffff" stroke-width="%.2f"' % (0.6 * pt))
        out.append('<circle cx="%.1f" cy="%.1f" r="%.2f" %s '
                   'opacity="0"><animate attributeName="opacity" calcMode="discrete" values="0;1" keyTimes="0;%.4f" %s/>'
                   '</circle>' % (xi, yi, r, style, kt, loop))
    static_dl = lv.DrawList(scene.W, scene.H, dpi, None)
    add_colorbar_legend(static_dl, a, scene.emin, scene.emax, scene.W - scene.band + 14 * pt, 40 * pt,
                        ring=scene.ring_hits)
    fs = a.font_size
    static_dl.text(10 * pt, 1.6 * fs * pt + 14 * pt, scene.view["label"], fs * 0.9, INK2)
    out += svg_items(static_dl.items, dpi)
    # clock: a set of text labels switched on in turn (every 1 ns)
    for i, t0 in enumerate(np.arange(0, t_end + 1e-9, 1.0)):
        k0 = t0 * k / total
        k1 = min((t0 + 1.0) * k / total, 1.0) if t0 + 1.0 <= t_end else 1.0
        out.append('<text x="%.1f" y="%.1f" font-family="Helvetica, Arial, sans-serif" font-size="%.1f" fill="#0b0b0b" '
                   'opacity="0">t = %4.0f ns<animate attributeName="opacity" calcMode="discrete" values="0;1;0" '
                   'keyTimes="0;%.4f;%.4f" %s/></text>' % (10 * pt, 1.6 * fs * pt + 30 * pt, fs * 1.4 * pt, t0,
                                                           max(k0, 1e-4), max(k1, k0 + 1e-4) if k1 < 1 else 0.9999,
                                                           loop))
    out.append("</svg>")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


# =====================================================================================================
#  CLI
# =====================================================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prim", required=True)
    ap.add_argument("--geo-cache", default=None)
    ap.add_argument("--event-npz", required=True)
    ap.add_argument("--steps", default=None, help="steps.bin of SolidStepRecorder (route A)")
    ap.add_argument("--traj", choices=["steps", "edm4hep", "both"], default="steps")
    ap.add_argument("-e", "--event", type=int, default=0)
    ap.add_argument("--views", default="side,beam")
    ap.add_argument("--phi", type=float, default=None, help="side view azimuth (deg; default: calo centroid)")
    ap.add_argument("--zsec", type=float, default=None, help="beam view section z [mm]")
    ap.add_argument("--beam-rmax", type=float, default=2600.0)
    ap.add_argument("--perspective", type=float, default=40000.0, help="cutaway camera distance [mm]")
    ap.add_argument("--zoom-half", type=float, default=600.0, help="sidezoom: half height of the window [mm]")
    ap.add_argument("--proj", choices=["rhoz", "ortho"], default="rhoz",
                    help="side views: rho-z projection of tracks and hits into the event plane (default) or "
                         "orthographic projection onto it")
    ap.add_argument("--no-project-modules", dest="project_modules", action="store_false",
                    help="side view: do not draw hit modules that are not cut by the section")
    ap.add_argument("--size", default="2400x1500")
    ap.add_argument("--dpi", type=float, default=200.0)
    ap.add_argument("--lw", type=float, default=0.35)
    ap.add_argument("--font-size", type=float, default=9.0)
    ap.add_argument("--traj-lw", type=float, default=0.45)
    ap.add_argument("--hit-size", type=float, default=1.6, help="tracker hit marker radius [pt]")
    ap.add_argument("--emin", type=float, default=5.0, help="charged secondaries: min vertex Ekin [MeV]")
    ap.add_argument("--emin-gamma", type=float, default=20.0)
    ap.add_argument("--emin-neutral", type=float, default=50.0)
    ap.add_argument("--emin-cal", type=float, default=0.1, help="colour scale minimum [MeV]")
    ap.add_argument("-o", "--outdir", default=".")
    ap.add_argument("--tag", default="event")
    ap.add_argument("--formats", default="svg,png,pdf")
    ap.add_argument("--animate", default=None, help="comma list of views to animate")
    ap.add_argument("--frames", type=int, default=120)
    ap.add_argument("--t-end", type=float, default=None, help="animation end time [ns]")
    ap.add_argument("--frame-scale", type=float, default=0.5)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--svg-anim", action="store_true")
    ap.add_argument("--svg-duration", type=float, default=10.0, help="animated SVG: seconds for t = 0..t_end")
    a = ap.parse_args()
    lv.VERBOSE = 1
    os.makedirs(a.outdir, exist_ok=True)
    timing = {}
    t0 = time.perf_counter()
    sc = load_geometry(a.prim, a.geo_cache)
    timing["geometry"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    ev = load_event(a.event_npz, a.event)
    tracks = []
    if a.traj in ("steps", "both"):
        tracks += load_steps(a.steps, a.event)
    if a.traj in ("edm4hep", "both"):
        tracks += tracks_from_edm4hep(ev)
    timing["event"] = time.perf_counter() - t0
    lv.log("event %d: %d tracks, %d tracker hits, %d calo hits, %d contributions (%.2f s)" % (
        a.event, len(tracks), len(ev["th_time"]), len(ev["ch_energy"]), len(ev["cc_energy"]), timing["event"]))
    scenes = {}
    for vk in [v for v in a.views.split(",") if v]:
        t0 = time.perf_counter()
        view = make_view(vk, ev, a)
        sc_ = Scene(sc, ev, tracks, view, a, title="SoLID-Lite J/ψ (solid_JPsi_lite), event %d" % a.event)
        scenes[vk] = sc_
        timing["scene_" + vk] = time.perf_counter() - t0
        dl = sc_.drawlist(None)
        for fmt in a.formats.split(","):
            p = os.path.join(a.outdir, "%s_%s.%s" % (a.tag, vk, fmt))
            t1 = time.perf_counter()
            lv.write_outputs(dl, [p], "cairo")
            timing["write_%s_%s" % (vk, fmt)] = time.perf_counter() - t1
    if a.animate:
        for vk in a.animate.split(","):
            scene = scenes[vk]
            t_end = a.t_end or float(max(t.T[-1] for t in scene.tracks if t.parent == 0)) + 2.0
            fr = FrameRenderer(scene, a.frame_scale)
            fdir = os.path.join(a.outdir, "frames_%s_%s" % (a.tag, vk))
            os.makedirs(fdir, exist_ok=True)
            t1 = time.perf_counter()
            for i in range(a.frames):
                fr.frame(t_end * i / (a.frames - 1), os.path.join(fdir, "f%04d.png" % i))
            dt = time.perf_counter() - t1
            timing["frames_%s" % vk] = dt
            lv.log("animation %s: %d frames %dx%d in %.1f s (%.3f s/frame, layer cache %.2f s), t_end %.1f ns" % (
                vk, a.frames, fr.W, fr.H, dt, dt / a.frames, fr.t_cache, t_end))
            if a.svg_anim:
                t1 = time.perf_counter()
                p = os.path.join(a.outdir, "%s_%s_anim.svg" % (a.tag, vk))
                write_svg_animation(scene, p, t_end, duration=a.svg_duration)
                timing["svg_anim_%s" % vk] = time.perf_counter() - t1
                lv.log("wrote %s (%.1f kB)" % (p, os.path.getsize(p) / 1024))
    for k, v in timing.items():
        print("timing %-22s %8.3f s" % (k, v))


if __name__ == "__main__":
    main()
