#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare the hall-frame extents of a SoLID DD4hep configuration with the corresponding solid_gemc gcard.

usage (inside eic-shell, after `source thissolid.sh`):
    compare_gemc.py solid_SIDIS_He3 --gemc /path/to/solid_gemc [--md out.md] [--json out.json] [--tol 0.1]

The configuration name selects the top-level compact ($DETECTOR_PATH/<config>.xml) and the GEMC reference gcard
(table CONFIG_GCARD below).  Each DD4hep top-level DetElement is matched to one or more GEMC detector files of the
gcard (table DETECTOR_MAP: regex on the GEMC file basename + optional regex on the GEMC volume name).

Extent conventions (identical on both sides, cm):
  * z range and r range of the material volumes (virtual/flux planes excluded on the GEMC side);
  * polyhedra: r_max is the outer corner, r_min the inner flat (GEMC analytic) / inner corner (DD4hep mesh), so a
    polyhedral inner radius differs by r_flat*(1/cos(pi/n)-1) between the two (shown as 'pgon' in the notes).
Only envelopes are compared here; internal positions (mirrors, modules, plates) are checked by the subsystem tests.
"""
import argparse
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gemc_txt  # noqa: E402

CONFIG_GCARD = {
    'solid': 'solid_SIDIS_He3_moved_full.gcard',
    'solid_SIDIS_He3': 'solid_SIDIS_He3_moved_full.gcard',
    'solid_SIDIS_NH3': 'solid_SIDIS_NH3_moved_full.gcard',
    'solid_JPsi': 'solid_JPsi_LH2_moved_full.gcard',
    'solid_DDVCS': 'solid_JPsi_DDVCS_LH2_moved_full.gcard',
    'solid_PVDIS': 'solid_PVDIS_LD2_moved_full.gcard',
    'solid_JPsi_lite': 'solid_JPsi_LH2_lite_full.gcard',
    'solid_DDVCS_lite': 'solid_JPsi_DDVCS_LH2_lite_full.gcard',
}

# (dd4hep DetElement name, regex on GEMC detector-file basename, regex on GEMC volume name or None)
DETECTOR_MAP = [
    ('MagnetYoke', r'^solid_solenoid$', None),
    ('MagnetCoil', r'^solid_solenoid_cryo$', None),
    ('Target', r'_target_', None),
    ('Beamline', r'_beamline', None),
    ('PVDISBaffle', r'_baffle_', None),
    ('GEMTrackerAssembly', r'^solid_(SIDIS|PVDIS)_gem$', None),        # holds GEMTracker (+ empty dummy, ACTS)
    ('GEMTrackerLiteAssembly', r'^solid_SIDIS_(LA|FA)_gem$', None),    # holds GEMTrackerLA + GEMTrackerFA
    ('GEMTracker', r'^solid_(SIDIS|PVDIS)_gem$', None),
    ('GEMTrackerLA', r'^solid_SIDIS_LA_gem$', None),
    ('GEMTrackerFA', r'^solid_SIDIS_FA_gem$', None),
    ('MuonGEMTracker', r'^solid_DDVCS_muon_gem_forwardangle$', None),
    ('LightGasCherenkov', r'^lg_cherenkov$', None),
    ('HeavyGasCherenkov', r'^solid_SIDIS_hgc$', None),
    ('LiteCherenkov', r'^solid_SIDIS_hgc$', None),
    ('FAECPreShower', r'_ec_forwardangle$', r'shield|presci|suppor'),
    ('FAECShower', r'_ec_forwardangle$', r'showe'),
    ('LAECPreShower', r'_ec_largeangle$', r'shield|presci|suppor'),
    ('LAECShower', r'_ec_largeangle$', r'showe'),
    ('FASPD', r'_spd_forwardangle$', None),
    ('LASPD', r'_spd_largeangle$', None),
    ('FAMRPC', r'_mrpc_forwardangle$', None),
    ('MuonIron', r'_muon_forwardangle$', r'(?i)iron|_[123][1-4]$'),
    ('MuonScint', r'_muon_forwardangle$', r'(?i)scint'),
]

VIRTUAL = re.compile(r'(?i)virtual|_virt|observer|_vp_|vplane')


# ------------------------------------------------------------------ GEMC side
def _ring(r, z, n=72, phi0=0.0):
    """n points on a circle (or the n corners of a regular polygon starting at phi0 [deg])."""
    p0 = math.radians(phi0)
    return [(r * math.cos(p0 + 2 * math.pi * i / n), r * math.sin(p0 + 2 * math.pi * i / n), z) for i in range(n)]


def local_points(v, vols, depth=0):
    """Point cloud on the surface envelope of a GEMC solid (local frame) + analytic coaxial r-range."""
    t = v.type
    if t.startswith('CopyOf'):
        orig = t.split()[1]
        return local_points(vols[orig], vols, depth + 1) if orig in vols else ([], None)
    if t.startswith('Operation'):
        expr = t.split(':', 1)[1].replace('@', ' ').replace('~', ' ').split()
        first = expr[0]
        if first not in vols:
            return [], None
        pts, rr = local_points(vols[first], vols, depth + 1)
        # a coaxial tube-shaped hole sets the inner radius (same rule as solid_gemc geometry_outline.py)
        if len(expr) >= 3 and expr[1] == '-' and expr[2] in vols and vols[expr[2]].type == 'Tube':
            hole = vols[expr[2]].dims[1]
            if rr:
                rr = (max(rr[0], hole), rr[1])
            pts = [p for p in pts if math.hypot(p[0], p[1]) >= hole - 1e-6] or pts
        return pts, rr
    d = v.dims
    pts, rr = [], None
    try:
        if t == 'Tube':
            rmin, rmax, hz = d[0], d[1], d[2]
            for z in (-hz, hz):
                pts += _ring(rmax, z) + (_ring(rmin, z) if rmin > 0 else [(0, 0, z)])
            rr = (rmin, rmax)
        elif t == 'Cons':
            rmin1, rmax1, rmin2, rmax2, hz = d[:5]
            pts += _ring(rmax1, -hz) + _ring(rmax2, hz) + _ring(max(rmin1, 1e-9), -hz) + _ring(max(rmin2, 1e-9), hz)
            rr = (min(rmin1, rmin2), max(rmax1, rmax2))
        elif t == 'Polycone':
            n = int(d[2]); rin, rout, z = d[3:3 + n], d[3 + n:3 + 2 * n], d[3 + 2 * n:3 + 3 * n]
            for a, b, c in zip(rin, rout, z):
                pts += _ring(b, c) + (_ring(a, c) if a > 0 else [(0, 0, c)])
            rr = (min(rin), max(rout))
        elif t == 'Pgon':
            ns, n = int(d[2]), int(d[3]); rin, rout, z = d[4:4 + n], d[4 + n:4 + 2 * n], d[4 + 2 * n:4 + 3 * n]
            k = 1 / math.cos(math.pi / ns)
            for a, b, c in zip(rin, rout, z):
                pts += _ring(b * k, c, ns, d[0]) + _ring(a * k, c, ns, d[0])
            rr = (min(rin), max(rout) * k)
        elif t in ('Box',):
            hx, hy, hz = d[:3]
            pts = [(sx * hx, sy * hy, sz * hz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
        elif t == 'Trd':
            hx1, hx2, hy1, hy2, hz = d[:5]
            pts = [(sx * hx1, sy * hy1, -hz) for sx in (-1, 1) for sy in (-1, 1)] + \
                  [(sx * hx2, sy * hy2, hz) for sx in (-1, 1) for sy in (-1, 1)]
        elif t == 'Sphere':
            rmax = d[1]
            pts = [(sx * rmax, sy * rmax, sz * rmax) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
        elif t == 'G4Trap':
            hz = d[0]; hy = max(d[3], d[7]); hx = max(d[4], d[5], d[8], d[9])
            pts = [(sx * hx, sy * hy, sz * hz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    except (IndexError, TypeError, ValueError):
        pass
    return pts, rr


def gemc_extents(gcard_path, gemc_script_dir):
    text = open(gcard_path).read().split('</gcard>')[0]
    dets = re.findall(r'<detector name="([^"]+)" factory="TEXT" variation="([^"]+)"', text)
    out = []
    for path, var in dets:
        f = os.path.normpath(os.path.join(gemc_script_dir, f'{path}__geometry_{var}.txt'))
        vols = gemc_txt.load(f)
        base = os.path.basename(path)
        for v in vols.values():
            if v.mother != 'root' or v.material == 'Component' or v.exist == '0':
                continue
            pts, rr = local_points(v, vols)
            if not pts:
                continue
            R, tt = gemc_txt.world_transform(vols, v.name)
            gp = [tuple(sum(R[i][k] * p[k] for k in range(3)) + tt[i] for i in range(3)) for p in pts]
            zs = [p[2] for p in gp]
            rs = [math.hypot(p[0], p[1]) for p in gp]
            coax = abs(tt[0]) < 1e-6 and abs(tt[1]) < 1e-6 and abs(abs(R[2][2]) - 1) < 1e-9
            rmin, rmax = (rr if (coax and rr) else (min(rs), max(rs)))
            out.append(dict(file=base, variation=var, volume=v.name, material=v.material,
                            virtual=bool(VIRTUAL.search(v.name)) or v.material == 'G4_Galactic',
                            zmin=min(zs), zmax=max(zs), rmin=rmin, rmax=rmax))
    return out


def aggregate(rows):
    if not rows:
        return None
    return dict(zmin=min(r['zmin'] for r in rows), zmax=max(r['zmax'] for r in rows),
                rmin=min(r['rmin'] for r in rows), rmax=max(r['rmax'] for r in rows), n=len(rows))


# ------------------------------------------------------------------ comparison
def compare(dd4hep_dets, gemc_rows, tol):
    table = []
    for det in dd4hep_dets:
        name = det['detector']
        e = det['extent']
        rule = [m for m in DETECTOR_MAP if m[0] == name]
        g = None
        if rule:
            _, frx, vrx = rule[0]
            sel = [r for r in gemc_rows if re.search(frx, r['file']) and not r['virtual']
                   and (vrx is None or re.search(vrx, r['volume']))]
            g = aggregate(sel)
        row = dict(detector=name, dd4hep=e, gemc=g)
        if g and e['zmin'] is not None:
            dz = max(abs(e['zmin'] - g['zmin']), abs(e['zmax'] - g['zmax']))
            dr = max(abs(e['rmin'] - g['rmin']), abs(e['rmax'] - g['rmax']))
            row.update(dz=dz, dr=dr, ok=(dz <= tol and dr <= tol))
        table.append(row)
    return table


def fmt_range(a, b):
    return '–' if a is None else f'{a:.2f} … {b:.2f}'


def to_markdown(config, gcard, table, tol):
    lines = [f'### {config}  (GEMC: `{gcard}`)', '',
             '| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |',
             '|---|---|---|---|---|---:|---:|---|']
    for r in table:
        e, g = r['dd4hep'], r['gemc']
        st = 'no GEMC match' if g is None else ('OK' if r.get('ok') else 'DIFF')
        lines.append(f"| {r['detector']} | {fmt_range(e['zmin'], e['zmax'])} | "
                     f"{fmt_range(g['zmin'], g['zmax']) if g else '–'} | {fmt_range(e['rmin'], e['rmax'])} | "
                     f"{fmt_range(g['rmin'], g['rmax']) if g else '–'} | "
                     f"{r.get('dz', float('nan')):.2f} | {r.get('dr', float('nan')):.2f} | {st} |")
    lines.append('')
    lines.append(f'tolerance {tol} cm; envelopes only (see the module docstring for conventions).')
    return '\n'.join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('config', help='configuration name, e.g. solid_SIDIS_He3')
    sibling = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'solid_gemc')
    ap.add_argument('--gemc', default=os.environ.get('SOLID_GEMC', os.path.normpath(sibling)),
                    help='path to a solid_gemc checkout (default: $SOLID_GEMC, else ../solid_gemc next to '
                         'this repository)')
    ap.add_argument('--compact', default=None, help='compact file (default $DETECTOR_PATH/<config>.xml)')
    ap.add_argument('--tol', type=float, default=0.1, help='tolerance in cm for OK (default 0.1)')
    ap.add_argument('--md', default=None)
    ap.add_argument('--json', default=None)
    a = ap.parse_args()

    gcard = CONFIG_GCARD.get(a.config)
    if not gcard:
        sys.exit(f'no GEMC gcard known for {a.config}')
    script_dir = os.path.join(a.gemc, 'script')
    gemc_rows = gemc_extents(os.path.join(script_dir, gcard), script_dir)

    compact = a.compact or os.path.join(os.environ.get('DETECTOR_PATH', '.'), a.config + '.xml')
    import ROOT  # noqa: F401
    import dd4hep
    import solid_geo_dump
    desc = dd4hep.Detector.getInstance()
    desc.fromXML(compact)
    dets = solid_geo_dump.dump_detectors(desc)

    table = compare(dets, gemc_rows, a.tol)
    md = to_markdown(a.config, gcard, table, a.tol)
    print(md)
    if a.md:
        open(a.md, 'w').write(md + '\n')
    if a.json:
        json.dump(dict(config=a.config, gcard=gcard, table=table, gemc=gemc_rows), open(a.json, 'w'), indent=1,
                  default=lambda o: o.item() if hasattr(o, 'item') else str(o))


if __name__ == '__main__':
    main()
