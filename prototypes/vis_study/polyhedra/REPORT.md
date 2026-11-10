> Verbatim copy of the polyhedra role report, which was returned inline to the workflow on 6 Nov 2026 and saved here by the final-revision step.
> Corrections in $SCRATCH/vis_study/FEASIBILITY.md section 7 (review log) take precedence over statements they name.

The prototype works. The Geant4 polyhedra route reproduces the DAWNFILE scenes without the driver, and I recommend it. The TGeo route only works with repairs, and one of its booleans stays wrong. The harness would not let me write `REPORT.md`, so the full report follows this summary.

Work directory (all scripts, outputs and images): `$SCRATCH/vis_study/polyhedra` (`$W` below).

**Summary**

1. **Geant4 route: works.** I walk the Geant4 volume tree and call `G4VSolid::GetPolyhedron()` with 72 rotation steps, the same call the vis drivers make.
   * In all 7 configurations it gives the same volumes as DAWNFILE, in the same order, with the same names, colours, depths and vertex/facet totals. All 628 unique meshes are watertight and none fails.
   * Rendered with `solid_lineview`, PVDIS is pixel-identical to the DAWNFILE scene in all 4 standard views. SoLID-Lite J/psi differs by 0.04–0.09 % of pixels. All of that comes from the section plane x = 0 passing exactly through mesh vertices; at x = 0.7 mm the images are identical.
   * It runs standalone (no physics; 11–15 s as a prototype, about 3–4 s once compiled) or inside ddsim at `/run/initialize`, where it adds 3.0 s.
   * It also exposes two DAWNFILE defects. The 330 PVDIS baffle blades are written without a volume name. No invisible-edge flags are written at all.
2. **Volume IDs work.** Each drawn volume carries its DetElement path, readout, and a partial volume ID with a bit mask.
   * These equal the IDs Geant4 assigns during tracking for 308/308 sensitive volumes.
   * 821/821 simulated hits (3 events) map to exactly one drawn volume, so calorimeter blocks can be coloured by energy.
3. **TGeo route: not robust.** Out of the box (ROOT's mesh export plus its CSG library, RootCsg) it has four problems:
   * cracked boolean meshes (1554–7623 open edges per configuration), giving spurious lines and unfilled sections;
   * spheres stuck at 20 divisions, which puts a Cherenkov mirror 19 mm off;
   * 72 steps on every small tube segment (PVDIS: 1.4 M vertices, 84 MB file);
   * extra seam faces on full tubes.

   After repairs SoLID-Lite J/psi matches Geant4 to ≤ 0.04 % of pixels. RootCsg still gets `MagnetYoke_EndcapInterface` wrong in all 5 full-SoLID configurations (+52 % volume), and the NH3 target support (+4.6 %).
4. **Recommendation.**
   * **Source:** Geant4 polyhedra.
   * **Interface:** an instanced `.npz` scene (unique meshes, placements, IDs): 0.13–0.33 MB against 5.5–14.9 MB for the DAWNFILE `.prim`.
   * **Reader:** a `load_scene()` in `solid_lineview`. The 70-line prototype gives the same images and parses 2–3× faster than the DAWNFILE file.
   * **Compatibility:** keep a `.prim` writer.
   * **Producer:** compiled C++, usable as a ddsim run action and as a standalone tool.
5. **Effort:** about **5.5 developer-days**. Making the TGeo route robust would add 5–8 days; not recommended.

Key images:
* `$W/img/g4_PVDIS/cutaway.png` (driver-free, identical to DAWNFILE)
* `$W/img/booleans/cmp3_collars_hlr.png`
* `$W/img/booleans/cmp3_cherenkov_hlr.png`
* `$W/img/pvdis/cmp_endcapif.png`

---

# REPORT: geometry tessellation without a vis driver (role: polyhedra)

**Setup.**
* Container: eic_xl 26.10.0-stable (DD4hep 1.38, Geant4 11.4.1, ROOT 6.40.04).
* Geometry: `scratchpad/verify_clone/install` (main).
* Load: all jobs ran under `nice -n 10`, one ddsim at a time; machine load was 3–6.
* Variables used in the commands:
  * `SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif`
  * `I=$SCRATCH/verify_clone/install`
  * `P=$SCRATCH/overview_job/prims`
  * `W` = the work directory above.

## 1. What DAWNFILE writes (reference)

* **Cost:** `ddsim --runType vis` with `macro/dawnfile.mac` on SoLID-Lite J/psi takes **15.8 s wall** and writes a 5.48 MB `.prim` (run kept in `$W/dawn_timing/`).
* **Content:** 3546 `#/PVName` blocks: 2946 `/Polyhedron`, plus 306 `/Tubs`, 270 `/Trd` and 24 `/Box`. `solid_lineview` tessellates the last three itself.
* **Conventions:**
  * vertices are in the local frame, with `/Origin` and `/BaseVector`;
  * the name tag is `name:copyNo`, indented by 5 + 5·depth spaces;
  * `/BoundingBox` is the bounding sphere of the visible volumes.
* **Defect 1: no invisible-edge flags.** All `/Facet` labels are positive: 0 of 62 802 (J/psi-lite) and 0 of 173 042 (PVDIS). The Geant4 polyhedra themselves flag 11 924 J/psi-lite facets.
* **Defect 2: anonymous extruded solids.** G4ExtrudedSolid volumes (G4TessellatedSolid-derived) get colour, transform and polyhedron but no `#/PVName` header.
  * All 330 PVDIS `PVDISBaffle_plate*_blade_*` become `<Polyhedron@line>`, with depth −1 and no hierarchy.
  * `--include`/`--exclude` cannot select them.
  * Their colour is whatever the previous volume set.

## 2. Geant4 route

### Implementation

**`solidvis_g4.h`** (C++ helper, compiled at run time by cling):
* **`Walker`** copies what `G4PhysicalVolumeModel` does:
  * depth-first walk with global transforms;
  * invisible volumes are not drawn, but their daughters are visited unless "daughters invisible" is set;
  * replicas and parameterised volumes are looped over (SoLID has none).
* **`tessellate()`**: `G4Polyhedron::SetNumberOfRotationSteps(72)`, then `GetPolyhedron()`, then `GetVertex`/`GetFacet` (index × edge flag).
* **`convert()`**: DD4hep `Geant4Converter` with no run manager.
* **`volume_id_at()`**: navigator touchable, then `Geant4VolumeManager::volumeID`; used only as a cross-check.

**`g4poly.py`** (Python driver):
* Builds one mesh per unique solid: 48 meshes for 3546 placements in J/psi-lite.
* Writes `.prim` and `.npz`.
* Maps each Geant4 volume to its TGeo node, through `g4Placements`, or `g4VolumeImprints` for assemblies.
* DetElement path: the deepest DetElement on that chain.
* Readout owner: the nearest DetElement that has a SensitiveDetector, so the GEM trackers inside `GEMTrackerLiteAssembly` resolve correctly.
* Volume ID: the placement `volIDs` are accumulated along the chain and encoded with the readout's `IDDescriptor`, giving a value and a mask.

**`polydump_steer.py`** (ddsim hook):
* It runs `steering/solid_sim.py`, then registers a function with `SIM.physics.setupUserPhysics(f)`.
* `f` adds a `Geant4PythonDetectorConstruction` whose sensitive-detector callback does the dump. That callback runs at `/run/initialize`, after sensitive detectors and field are set up.

Commands:

```bash
# standalone; --edge-flags both also writes .g4flags.prim
nice -n 10 singularity exec $SIF bash -c ". /etc/profile; source $I/bin/thissolid.sh solid_JPsi_lite; cd $W; \
  python3 g4poly.py -c \$DETECTOR_PATH/solid_JPsi_lite.xml -o out/g4_JPsi_lite \
  --bbox-from $P/solid_JPsi_lite/g4_0000.prim --edge-flags both"
# inside ddsim, 3 events, with the volume-ID cross-check
nice -n 10 singularity exec $SIF bash -c ". /etc/profile; source $I/bin/thissolid.sh solid_JPsi_lite; cd $W; \
  SOLID_OPTICAL=0 POLYDUMP_CHECK_IDS=1 POLYDUMP_OUT=out/hook_JPsi_lite ddsim --steeringFile polydump_steer.py \
  --compactFile \$DETECTOR_PATH/solid_JPsi_lite.xml -N 3 --outputFile out/hook_JPsi_lite.edm4hep.root"
./run_all.sh     # Geant4 and TGeo, all 7 configurations, plus comparisons -> out/all/
./render_views.sh <prim> <outdir>; python3 imgdiff.py A.png B.png diff.png
```

### Results versus DAWNFILE

Compared with `compare_prims.py`, which uses `solid_lineview`'s own parser:

| configuration | Geant4 volumes / drawn | unique meshes | vertices / facets | DAWNFILE `.prim` | dump `.prim` | `.npz` |
|---|---|---|---|---|---|---|
| SIDIS He3 | 4598 / 4519 | 68 | 115814 / 150538 | 9.0 MB | 7.7 MB | 0.20 MB |
| SIDIS NH3 | 4615 / 4537 | 89 | 126143 / 162947 | 9.7 | 8.4 | 0.30 |
| PVDIS | 8846 / 8604 | 231 | 171680 / 208814 | 14.9 | 11.6 | 0.33 |
| J/psi | 4587 / 4512 | 64 | 113582 / 148306 | 8.9 | 7.6 | 0.19 |
| DDVCS | 4879 / 4794 | 72 | 115838 / 149998 | 9.0 | 7.8 | 0.20 |
| SoLID-Lite J/psi | 3580 / 3546 | 48 | 53874 / 67950 | 5.5 | 3.6 | 0.13 |
| SoLID-Lite DDVCS | 3872 / 3828 | 56 | 56130 / 69642 | 5.6 | 3.8 | 0.14 |

* **Totals and names:** vertex and facet totals equal DAWNFILE in every configuration. Names, colours and depths are identical, except that the 330 PVDIS blades are now named.
* **Per mesh (J/psi-lite):**
  * 2906 polyhedra are identical in vertices (to 1e-6 mm), facets and edge flags;
  * 40 differ only in edge flags;
  * the 600 Box/Trd/Tubs volumes have the same counts but a different facet order.
* **Placement:** the largest bounding-box difference is 5e-6 mm.

Images (share of pixels that differ by more than 8/255):

| | section_yz | section_xz | cutaway | lines_yz |
|---|---|---|---|---|
| PVDIS, dump vs DAWNFILE | **0** | **0** | **0** | **0** |
| J/psi-lite, dump vs DAWNFILE | 0.064 % | 0.064 % | 0.039 % | 0.089 % |
| J/psi-lite, with Geant4 edge flags | 0.064 % | 0.064 % | 0.181 % | 0.112 % |

* **J/psi-lite residuals:** thin discs that DAWNFILE wrote as `/Tubs`, and the x = 0 section plane passing exactly through vertices. At x = 0.7 mm the section is identical (`img/diff_g4_section_yz.png`).
* **Edge flags:** keeping the Geant4 flags removes 3820 px of real outline edges (`img/diff_g4flags_cutaway.png`). Writing all edges visible, as DAWNFILE does, is correct.

**Timing (standalone prototype):** 10.7–14.7 s wall.

| step | time |
|---|---|
| cling compile of the Geant4/DDG4 headers | 5.5–6.6 s (0.65–11 s depending on cache and load) |
| DD4hep geometry load | 1.2–1.7 s |
| `Geant4Converter` | 1.1–1.5 s |
| tree walk | 0.1 s |
| all `GetPolyhedron` calls | 4–10 ms |
| Python loops (per volume, IDs) | 2.2–3.4 s |

* Compiled, the whole thing should take about 3–4 s.
* Inside ddsim the dump adds 3.0 s. The 3-event run took 19.0 s, against 15.8 s for the DAWNFILE run.

**IDs:**
* Coverage: 3546/3546 drawn volumes get a DetElement path; 3513 get a readout (beamline, magnet and target have none).
* Tracking cross-check: the encoded IDs equal the touchable IDs for 308/308 sensitive drawn volumes (8 Cherenkov tiles, 300 SPD pads).
* Hit mapping: `hitmap_check.py` (run with `python3 -I`) maps 821/821 hits in 3 events uniquely, each inside its volume's bounding box.
  * Collections: FAEC preshower and shower, LAEC shower, GEM FA/LA, FASPD, FAMRPC.
  * Calorimeter hits land on the drawn module, e.g. `/world/FAECShower/module457`.
* Caveat: `Geant4VolumeManager::volumeDescriptor(path)` returns "Insensitive" for the 300 assembly-imprinted SPD pads, although touchables work. Don't use the path API.
* Standalone mode has no sensitive detectors attached, so only the DD4hep-side encoding is available there.

**Vis semantics:** DDG4's `Geant4Converter` maps VisAttr *lineStyle* to Geant4 forced wireframe. The compact's `drawingStyle="wireframe"` (Cherenkov tanks, alpha 0.15) therefore never reaches Geant4 or DAWNFILE.

## 3. ROOT/TGeo route

**Implementation** (`solidvis_tgeo.h`, `tgeopoly.py`):
* Meshes: `TGeoShape::MakeBuffer3D()` → `RootCsg::ConvertToMesh`. Composite shapes go through `RootCsg::BuildUnion/Difference/Intersection` on transformed operands, the same way ROOT's web geometry viewer does it.
* Walk: global matrices and DD4hep VisAttr. Assemblies are never drawn and never stop the walk; depths are converted to Geant4 depths.
* IDs: same as the Geant4 route.
* Options: `--raw-tgeo --no-fix` gives plain ROOT behaviour.

```bash
nice -n 10 singularity exec $SIF bash -c ". /etc/profile; source $I/bin/thissolid.sh solid_PVDIS; cd $W; \
  python3 tgeopoly.py -c \$DETECTOR_PATH/solid_PVDIS.xml -o out/tgeoraw_PVDIS --raw-tgeo --no-fix \
  --bbox-from $P/solid_PVDIS/g4_0000.prim"      # without the two flags -> fixed version
python3 compare_npz.py out/g4_PVDIS.npz out/tgeofixed_PVDIS.npz; python3 topo_stats.py out/all/*.npz
```

| # | problem with raw ROOT | evidence | fix in the prototype |
|---|---|---|---|
| 1 | RootCsg output is not watertight (T-junction cracks) | 1554 / 1806 / 7623 open edges (J/psi-lite / PVDIS / NH3); crack lines in hlr; boolean sections unfilled (section_yz 1.26 % px off) | weld, insert T-junction vertices, re-fan faces; leaves 7–410 open edges |
| 2 | `TGeoSphere` ignores `SetNsegments` (stays at 20 divisions) | LiteCherenkov mirror: Geant4 1325.8 cm³, z_min 159 mm; raw RootCsg 1316.4 cm³, z_min 140 mm; Monte-Carlo truth via `Contains` (2e7 points) 1329.7 cm³ ±0.3 %, lowest inside point at z = 163 mm | set sphere divisions from nseg; after the fix +0.075 % volume, within 1.2 mm of Geant4 |
| 3 | tube and polycone segments get Nsegments over their own angle (a 6° pad gets 72 steps) | 2.6× (J/psi-lite) and 8× (PVDIS: 1.40 M vertices, 83.6 MB) more vertices than Geant4 | scale by Δφ/360; afterwards +10–26 % vertices relative to Geant4 |
| 4 | full 360° `TGeoTubeSeg` gets seam faces at φ = 0 | 22 / 106 non-manifold edges | tessellate as a full tube/cone |
| 5 | concave `TGeoXtru` caps (baffle blades) | 22 concave n-gons; `solid_lineview`'s fan triangulation would be wrong | ear clipping; blade volumes then equal Geant4 exactly |
| 6 | wireframe semantics differ from Geant4 | the tank turns transparent (`img/booleans/cmp_tank_hlr.png`) | copy DDG4's lineStyle mapping |
| 7 | **wrong RootCsg results that no fix repairs** | `MagnetYoke_EndcapInterface`: 2.07e6 cm³ vs 1.36e6 (Geant4) and 1.39e6 (TGeo's own Monte-Carlo capacity), section unfilled; NH3 `Target_magnet_support` +4.6 % | none; needs a different CSG library |

With the fixes:

| | section_yz | section_xz | cutaway | lines_yz |
|---|---|---|---|---|
| J/psi-lite, TGeo vs Geant4 | 0 | 0.009 % | 0.043 % | 0.021 % |
| PVDIS, TGeo vs Geant4 | 0.056 % | 0.23 % (EndcapInterface) | 0.070 % | 0.076 % |

* Colours: identical for all 8604 PVDIS volumes.
* Bounding boxes: agree within 2.8 mm (the mirrors).
* Run time: 8.3–16.6 s wall, mostly cling (5.4–8 s) and the Python T-junction repair (1.9–6.6 s).

## 4. Recommendation and interface

**Use Geant4 polyhedra as the source.**
* They are exactly what DAWNFILE and all Geant4 vis drivers draw.
* Geant4's boolean processor gave watertight, correct meshes for all 628 unique meshes.
* The same code covers geometry-only drawings (standalone) and event displays (inside the ddsim session that writes trajectories and hits).

**Scene file `solid_scene.npz` (v1)**, as the prototype writes it:
* **meta:** JSON with config, source, nseg, units (mm) and readout field specs.
* **unique meshes:** vertices, facets (N×4 int32, 1-based, sign = Geant4 edge flag), offsets, name and type. Vertices should be float64; the prototype's float32 changed 0.25 % of pixels at x = 0 and none at x = 0.7.
* **placements, in depth-first order:**
  * name, copy, depth, parent;
  * mesh index (−1 = not drawn);
  * 3×4 matrix;
  * RGBA, visible/wireframe flags (add DD4hep drawing style and alpha).
* **IDs:** DetElement path, subdetector, readout, TGeo path, volume ID and mask (uint64), sensitive flag.
* **Events (from the event-display roles):** per-event trajectory polylines, plus hits mapped to volumes by `(cellID & mask) == volid`, giving per-volume energy for colouring.

**`solid_lineview`:**
* `load_scene(path)` dispatches on `.prim` or `.npz` and returns the existing `PrimScene`.
* Each unique mesh is welded once, then instanced.
* The prototype `lineview_npz.py` (70 lines) gives identical images (0 px at x = 0.7).

Parse times:

| scene | DAWNFILE `.prim` | dump `.prim` | `.npz` |
|---|---|---|---|
| J/psi-lite | 1.25 s | 2.7 s | 0.39 s |
| PVDIS | 4.6 s | 19 s | 2.2 s |

The dump `.prim` is slow because every volume is an explicit polyhedron, while lineview caches DAWNFILE's primitives.

**Producer:** compiled C++ instead of cling, which removes 5–11 s of compilation and 2–12 s of Python loops.
* As a DDG4 run action, e.g. `ddsim … --action.run '["SoLIDSceneDump"]'` (a proposal; I did not build this) or the steering hook above.
* As a standalone `solid_scene_dump <compact> -o scene.npz [--prim]`.
* `make_solid_views` uses it by default and keeps `--dawnfile` as a fallback.

## 5. Effort to productise the Geant4 route

| item | days |
|---|---|
| C++ dump (walker, polyhedra, ID encoding, npz/prim writers), as DDG4 action and standalone tool, CMake | 2.5 |
| `solid_lineview`: scene loader, npz reader, float64, nudge exactly-degenerate section planes, per-volume colour override | 1.5 |
| `make_solid_views` switch (DAWNFILE as fallback), `docs/visualization.md` | 0.5 |
| CI: compare with DAWNFILE on 2 configurations (counts, names, image diff), hit-mapping test | 1.0 |
| **total** | **≈ 5.5** |

A production-quality TGeo route would need another 5–8 days and is not recommended.

## 6. Files (all in `$W`)

* **Producers:** `g4poly.py`, `solidvis_g4.h`, `polydump_steer.py`, `tgeopoly.py`, `solidvis_tgeo.h`
* **Interface prototype:** `lineview_npz.py`
* **Comparison and checks:** `render_views.sh`, `imgdiff.py`, `compare_prims.py`, `compare_detail.py`, `compare_npz.py`, `topo_stats.py`, `hitmap_check.py`, `mirror_mc.py`, `tgeo_unit.py`
* **Batch:** `run_all.sh`
* **Outputs:**
  * `out/all/{g4,tgeo}_<config>.{prim,npz,summary.json,log}`
  * `out/all/cmp_*.txt`, `out/all/topo_*.txt`
  * `out/hook_JPsi_lite.{prim,npz,edm4hep.root,log}`
  * `out/tgeoraw_PVDIS.prim` is 84 MB and can be deleted
* **Sample images (all under 1 MB):**
  * `img/g4_PVDIS/cutaway.png`
  * `img/g4_JPsi_lite/{section_yz,cutaway}.png`
  * `img/booleans/cmp3_{collars,cherenkov}_hlr.png` (Geant4, raw TGeo, fixed TGeo)
  * `img/pvdis/cmp_pvdis_beam96.png` (baffle slits at z = 96 cm)
  * `img/pvdis/cmp_endcapif.png`
  * diff maps (red = only in the first image, blue = only in the second): `img/diff_g4_section_yz.png`, `img/diff_g4flags_cutaway.png`, `img/diff_tgeo_section_yz.png`, `img/pvdis/diff_tgeofixed_section_xz.png`

**Environment note:** neither ffmpeg nor imageio is available on the host or in the container, which matters for the animation role.