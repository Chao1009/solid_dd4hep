#!/usr/bin/env python3
"""Extract what an event display needs from a ddsim EDM4hep file, using DD4hep for the readout.

Runs inside eic_xl (podio/edm4hep + DD4hep python).  Writes one .npz with flat arrays:
  mcp_*    MCParticles               (event, pdg, status, charge, time, vertex, endpoint, momenta, parent)
  th_*     SimTrackerHits            (event, collection, cellID, edep, time, position, momentum, particle)
  ch_*     SimCalorimeterHits        (event, collection, cellID, energy, position, decoded ids, module key)
  cc_*     CaloHitContributions      (event, calo-hit row, pdg, energy, time, step position)
  mod_*    calorimeter module map    (detector, system, sector, module, TGeo node name, global centre [mm])
and checks the cellID -> position mapping (DD4hep CellIDPositionConverter vs the EDM4hep hit position).

usage: extract_event.py COMPACT.xml EDM4HEP.root OUT.npz
"""
import sys
import time
from array import array

import numpy as np

t0 = time.time()
compact, infile, outfile = sys.argv[1:4]

import ROOT  # noqa: E402
import dd4hep  # noqa: E402
import podio.root_io as rio  # noqa: E402

ROOT.gSystem.Load("libDDRec")
desc = dd4hep.Detector.getInstance()
desc.fromXML(compact)
conv = ROOT.dd4hep.rec.CellIDPositionConverter(desc)
t_geo = time.time() - t0

CALO = ["FAECShowerHits", "FAECPreShowerHits", "LAECShowerHits", "LAECPreShowerHits"]
CALO_DET = {"FAECShowerHits": "FAECShower", "FAECPreShowerHits": "FAECPreShower",
            "LAECShowerHits": "LAECShower", "LAECPreShowerHits": "LAECPreShower"}

# ---------------------------------------------------------------- readout descriptors
decoders = {}
for coll in CALO:
    ro = desc.readout(coll)
    idspec = ro.idSpec()
    dec = idspec.decoder()
    decoders[coll] = dec
    print("%-20s id spec: %s  segmentation: %s" % (coll, idspec.fieldDescription(), ro.segmentation().type()))

# ---------------------------------------------------------------- module map (placements with volIDs)
mod_rows = []
for coll in CALO:
    det = CALO_DET[coll]
    de = desc.detector(det)
    env = de.placement()
    vol = env.volume()
    wt = de.nominal().worldTransformation()
    sysid = int(de.id())
    n = vol.GetNdaughters()
    nmod = 0
    for i in range(n):
        node = vol.GetNode(i)
        pv = dd4hep.PlacedVolume(node)
        ids = {str(p.first): int(p.second) for p in pv.volIDs()}
        if "module" not in ids:
            continue
        tr = node.GetMatrix().GetTranslation()
        loc = array("d", [tr[0], tr[1], tr[2]])
        glo = array("d", [0.0, 0.0, 0.0])
        wt.LocalToMaster(loc, glo)
        # module full length along z: from the shape (TGeo cm)
        shape = node.GetVolume().GetShape()
        dz = shape.GetDZ() if hasattr(shape, "GetDZ") else 0.0
        mod_rows.append((det, sysid, ids.get("sector", 0), ids["module"], node.GetName(),
                         10 * glo[0], 10 * glo[1], 10 * glo[2], 10 * dz))
        nmod += 1
    print("%-14s system %3d: %d module placements (of %d daughters), e.g. %s ids=%s" % (
        det, sysid, nmod, n, mod_rows[-1][4], ids))

# ---------------------------------------------------------------- events
reader = rio.Reader(infile)
events = reader.get("events")
mcp, th, ch, cc = [], [], [], []
coll_names = []


def cid(name):
    if name not in coll_names:
        coll_names.append(name)
    return coll_names.index(name)


ncheck, maxdev = 0, 0.0
# VolumeManager cross-check: cellID -> DetElement 'module<id>' (its nominal world placement must be the module
# centre found above) and -> sensitive-slice context (world transform = cell centre = EDM4hep hit position)
vm = desc.volumeManager()
mod_centre = {(r[1], r[3]): np.array(r[5:8]) for r in mod_rows}
vm_n, vm_badname, vm_dmod, vm_dcell = 0, 0, 0.0, 0.0
vm_example = None
for ev, frame in enumerate(events):
    avail = set(frame.getAvailableCollections())
    for p in frame.get("MCParticles"):
        v, e, m, me = p.getVertex(), p.getEndpoint(), p.getMomentum(), p.getMomentumAtEndpoint()
        parents = p.getParents()
        par = parents[0].getObjectID().index if len(parents) else -1
        mcp.append((ev, p.getObjectID().index, p.getPDG(), p.getGeneratorStatus(), p.getSimulatorStatus(),
                    p.getCharge(), p.getTime(), p.getMass(), p.getEnergy(), v.x, v.y, v.z, e.x, e.y, e.z,
                    m.x, m.y, m.z, me.x, me.y, me.z, par))
    for name in sorted(avail):
        coll = frame.get(name)
        tname = type(coll).__name__
        if tname == "SimTrackerHitCollection":
            k = cid(name)
            for h in coll:
                x, mm = h.getPosition(), h.getMomentum()
                part = h.getParticle()
                pidx = part.getObjectID().index if part.isAvailable() else -1
                th.append((ev, k, h.getCellID(), h.getEDep(), h.getTime(), h.getPathLength(),
                           x.x, x.y, x.z, mm.x, mm.y, mm.z, pidx))
        elif tname == "SimCalorimeterHitCollection" and name in decoders:
            k = cid(name)
            dec = decoders[name]
            for h in coll:
                c = h.getCellID()
                x = h.getPosition()
                f = [dec.get(c, fld) for fld in ("system", "sector", "module", "layer", "slice")]
                row = len(ch)
                ch.append((ev, k, c, h.getEnergy(), x.x, x.y, x.z, *f))
                if ncheck < 2000:                      # DD4hep position of the cell vs EDM4hep hit position
                    pos = conv.position(c)
                    d = np.hypot(np.hypot(10 * pos.x() - x.x, 10 * pos.y() - x.y), 10 * pos.z() - x.z)
                    maxdev = max(maxdev, float(d))
                    ncheck += 1
                    de = vm.lookupDetElement(c)
                    if de.name() != "module%d" % f[2]:
                        vm_badname += 1
                    tm = de.nominal().worldTransformation().GetTranslation()
                    dm = np.linalg.norm(10 * np.array([tm[0], tm[1], tm[2]]) - mod_centre[(f[0], f[2])])
                    ctx = vm.lookupContext(c)
                    pw = ctx.localToWorld(ROOT.dd4hep.Position(0.0, 0.0, 0.0))
                    dc = np.linalg.norm(10 * np.array([pw.x(), pw.y(), pw.z()]) - np.array([x.x, x.y, x.z]))
                    vm_dmod, vm_dcell = max(vm_dmod, float(dm)), max(vm_dcell, float(dc))
                    vm_n += 1
                    if vm_example is None:
                        vm_example = ("%s cellID 0x%x = %s -> DetElement %s (centre %.1f, %.1f, %.1f mm), sensitive "
                                      "placement %s (centre %.1f, %.1f, %.1f mm)" % (
                                          name, c, dec.valueString(c), de.path(), 10 * tm[0], 10 * tm[1], 10 * tm[2],
                                          ctx.volumePlacement().name(), 10 * pw.x(), 10 * pw.y(), 10 * pw.z()))
                for con in h.getContributions():
                    s = con.getStepPosition()
                    part = con.getParticle()
                    cc.append((ev, row, con.getPDG(), con.getEnergy(), con.getTime(), s.x, s.y, s.z,
                               part.getObjectID().index if part.isAvailable() else -1))

print("cellID -> CellIDPositionConverter position vs EDM4hep hit position: max |d| = %.3g mm over %d hits"
      % (maxdev, ncheck))
print("VolumeManager: %d hits, %d DetElement name mismatches, max |module DetElement - placement map| = %.3g mm, "
      "max |sensitive-volume context - hit position| = %.3g mm" % (vm_n, vm_badname, vm_dmod, vm_dcell))
print("example:", vm_example)

f8, i4, i8, u8 = np.float64, np.int32, np.int64, np.uint64


def table(rows, cols):
    out = {}
    for j, (name, dt) in enumerate(cols):
        out[name] = np.array([r[j] for r in rows], dtype=dt)
    return out


arrs = {}
arrs.update({"mcp_" + k: v for k, v in table(mcp, [
    ("event", i4), ("index", i4), ("pdg", i4), ("gen", i4), ("sim", i4), ("charge", f8), ("time", f8),
    ("mass", f8), ("energy", f8), ("vx", f8), ("vy", f8), ("vz", f8), ("ex", f8), ("ey", f8), ("ez", f8),
    ("px", f8), ("py", f8), ("pz", f8), ("epx", f8), ("epy", f8), ("epz", f8), ("parent", i4)]).items()})
arrs.update({"th_" + k: v for k, v in table(th, [
    ("event", i4), ("coll", i4), ("cellID", u8), ("edep", f8), ("time", f8), ("path", f8),
    ("x", f8), ("y", f8), ("z", f8), ("px", f8), ("py", f8), ("pz", f8), ("particle", i4)]).items()})
arrs.update({"ch_" + k: v for k, v in table(ch, [
    ("event", i4), ("coll", i4), ("cellID", u8), ("energy", f8), ("x", f8), ("y", f8), ("z", f8),
    ("system", i4), ("sector", i4), ("module", i4), ("layer", i4), ("slice", i4)]).items()})
arrs.update({"cc_" + k: v for k, v in table(cc, [
    ("event", i4), ("hit", i8), ("pdg", i4), ("energy", f8), ("time", f8), ("x", f8), ("y", f8), ("z", f8),
    ("particle", i4)]).items()})
arrs.update({"mod_" + k: v for k, v in table(mod_rows, [
    ("det", "U16"), ("system", i4), ("sector", i4), ("module", i4), ("node", "U48"),
    ("cx", f8), ("cy", f8), ("cz", f8), ("dz", f8)]).items()})
arrs["collections"] = np.array(coll_names)
np.savez_compressed(outfile, **arrs)
print("events %d: %d MCParticles, %d tracker hits, %d calo hits, %d contributions, %d modules mapped"
      % (len(events), len(mcp), len(th), len(ch), len(cc), len(mod_rows)))
print("geometry load %.1f s, total %.1f s -> %s" % (t_geo, time.time() - t0, outfile))
