> Verbatim copy of the geant4_native role report, which was returned inline to the workflow on 6 Nov 2026 and saved here by the final-revision step.
> Corrections in $SCRATCH/vis_study/FEASIBILITY.md section 7 (review log) take precedence over statements they name.

**geant4_native summary.** I could not write REPORT.md: the harness blocks report files from subagents, so the full report follows below the summary.

**Drivers that work headless** in ddsim with Geant4 11.4.1: TSG_OFFSCREEN (PNG/JPEG raster, gl2ps PDF/SVG/EPS vector), VRML2FILE, DAWNFILE and ASCIITree.
- **RayTracer crashes** (segfault in `G4TheMTRayTracer::StoreUserActions`: the MT ray tracer finds no MT run manager, since ddsim runs a sequential one).
- **OpenGL and the X11/Qt drivers need a display.** Qt's offscreen mode gets no OpenGL context, and there is no Xvfb in the container or on the host.
- Vtk is not built, and HepRepFile was removed in 11.4.

**Geometry drawings:**
- TSG colour renders of the x > 0 cutaway take 2.2–4 s at 3200×2000. Sections take under 1 s. No Boolean cutaway failures.
- DD4hep forces solid style and alpha < 1 on most volumes; override with `/vis/geometry/set/forceSolid all -1 false` and `/vis/geometry/set/colour all -1 R G B 1`.
- `sectionPlane` also moves the camera to the plane normal; only one clipper (section or cutaway) works at a time.
- **TSG has no real hidden-line mode:** wireframe plus hiddenEdge fills faces in the object colour, so a black drawing comes out as a solid black silhouette.
- **A workaround works:** render an opaque one-colour scene as gl2ps SVG, then recolour it with `gl2ps_svg_hlr.py` (white faces, black edges). The depth-sorted output hides the back lines and looks close to `solid_lineview`.
  - It takes 4.3 s and gives a 4.3 MB SVG; `solid_lineview` takes 2.9 s for 153 KB.
  - Curved surfaces get no outline lines, and some edges are broken or faint.

**Events:** trajectories (rich, by charge, with filters) work in raster and vector. `/vis/scene/add/hits` draws nothing, because DDG4 hits have no Draw method.

**Energy colouring:**
- **Touchables work** via a macro generated from EDM4hep hits, with a second pass on the same seed. The DD4hep module id is not the copy number, so the mapping goes through DetElements.
  - Each touchable costs 0.38 s (full search of 5.08 M volumes): 50 FAEC modules in 19 s.
  - The FAEC beam view reads very well.
- **Scoring meshes are crude.**
  - They need a cppyy hook, because ddsim never creates the scoring manager.
  - `/score/draw` is now `/score/drawProjection`.
  - The cylinder r-phi projection and `drawColumn` work; one box projection did not draw (not investigated).

**Animation works.**
- Time windows and the head-time label draw in TSG; the light front does not.
- The offscreen viewer ignores time-window changes on refresh, so every frame needs a full rebuild: 3.4 s/frame with a cutaway, 0.3 s without.
- 60 frames gave a 0.60 MB GIF and a 0.52 MB WebP. There is no ffmpeg, so no MP4.

**Recommendation:** keep DAWNFILE + `solid_lineview` for publication line drawings. Use TSG_OFFSCREEN for colour, event, energy and animation output, with the gl2ps hidden-line trick as a fallback that needs no DAWN.

**Effort:** about 10–14 developer-days for everything; standard views plus event displays alone take 4–5 days.

Key samples:
- $SCRATCH/vis_study/geant4_native/samples/02_vector_hlr_vs_lineview.png
- $SCRATCH/vis_study/geant4_native/samples/01_compare_cutaway_2x2.png
- $SCRATCH/vis_study/geant4_native/samples/05a_faec_energy_touchables_beamview.png
- $SCRATCH/vis_study/geant4_native/samples/07a_anim_640x400_60f.gif

---

# REPORT.md (full content)

## Geant4-native visualization for SoLID (role: geant4_native)

`$G` = $SCRATCH/vis_study/geant4_native

- Every run has its macro, `ddsim.log` and outputs in `$G/runs/NN_*`.
- Images under 1 MB are in `$G/samples/`.
- Configuration: `solid_JPsi_lite`, eic_xl 26.10.0, Geant4 11.4.1, DD4hep 1.38.
- Headless: `DISPLAY`, `WAYLAND_DISPLAY` and `XAUTHORITY` unset; one ddsim at a time under `nice -n 10`.

### 0. Verdict

| goal | headless, stock 11.4.1? | how / limit |
|---|---|---|
| colour surface drawings, cutaways, sections | yes | TSG_OFFSCREEN zb (PNG) and gl2ps (PDF/SVG/EPS); 1–5 s per image |
| black-line hidden-line drawing | not directly | TSG hlr fills faces in the object colour. Workaround: opaque gl2ps SVG + recolour script gives a depth-sorted vector drawing close to `solid_lineview`. Files 25–30× larger, no outlines of curved surfaces. |
| trajectories | yes | `/vis/scene/add/hits` draws nothing (DDG4 hits have no `Draw`) |
| calorimeter modules coloured by energy | yes | generated touchable macro, two passes with the same seed, 0.38 s per touchable |
| scoring-mesh overlay | partly | cppyy hook needed; cylinder r-phi and box `drawColumn` work; crude |
| animation | yes | generic 11.4 time windows; full rebuild per frame (offscreen viewer bug); GIF/WebP |
| RayTracer | no | segfault in `G4TheMTRayTracer::StoreUserActions` |
| OpenGL (true hidden lines) | no | needs X; Qt offscreen has no GL context; no Xvfb |

DAWNFILE is still built and shipped in 11.4 (last code change 2023). Only HepRepFile and Qt3D were removed.

### 1. Drivers (`/vis/list`, `$G/runs/01_list`)

Registered: ASCIITree, DAWNFILE, RayTracer, VRML2FILE, gMocrenFile, TSG_OFFSCREEN, OGLIQt/OGLSQt, OGLIX/OGLSX, RayTracerX, RayTracerQt, TSG_X11_GLES/ZB and TSG_QT_GLES/ZB (TSG_QT_GLES is the default "TSG"/"OGL"). `geant4-config --has-feature` reports vtk=no, freetype=no.

| driver | headless? | outputs | tested result |
|---|---|---|---|
| TSG_OFFSCREEN | yes | zb_png/jpeg/ps; gl2ps eps/ps/pdf/svg/tex/pgf | main candidate (sections 2–5) |
| VRML2FILE | yes | `.wrl` | see below |
| DAWNFILE | yes | `.prim` | the current reference pipeline |
| RayTracer | crashes | JPEG | `gdb`: `StoreUserActions` ← `CreateBitMap` ← `Trace`, via both `/vis/rayTracer/trace` and `/vis/viewer/flush` (`$G/runs/10_rt`) |
| gMocrenFile | yes | `.gdd` (medical format) | not relevant |
| ASCIITree | yes | text | whole-world dumps at verbosity ≥ 2 grow past 5 GB and were killed twice (rc 137, likely out of memory). The FAEC alone has 4 M volumes; dump subtrees instead. |
| OGL*/TSG_X11_*/TSG_QT_*/RTX/RTQt | need a display | | OGLSQt with `QT_QPA_PLATFORM=offscreen`: "QOpenGLWidget: Failed to create context", then gl2ps errors |
| Vtk | not built | | |

VRML2FILE took 2.3 s and wrote 2.06 MB for the cut geometry:
- 1833 face sets carrying volume names (in comments), colours and transparency.
- The cutaway is applied; polyhedron edge flags are lost.
- It could replace `.prim` as `solid_lineview` input or feed a web viewer. Only geometry was tested.

### 2. Geometry with TSG_OFFSCREEN

**Wrapper used for every vis run** (`$G/bin/run_vis.sh RUNDIR MACRO`; `STEER=` picks a wrapper steering file):

```
singularity exec $SIF bash -c '. /etc/profile >/dev/null 2>&1; unset DISPLAY WAYLAND_DISPLAY XAUTHORITY;
  source $INST/bin/thissolid.sh solid_JPsi_lite; export SOLID_OPTICAL=0;
  echo exit | nice -n 10 ddsim --steeringFile ${STEER:-$DETECTOR_PATH/steering/solid_sim.py} \
    --compactFile $DETECTOR_PATH/$DETECTOR_CONFIG.xml --runType vis --macroFile MACRO --outputFile sim.edm4hep.root'
```

**Core macro** (cutaway like `make_solid_views cutaway`):

```
/vis/open TSG_OFFSCREEN 1600x1000-0+0
/vis/tsg/offscreen/set/size 3200 2000        # render 2x and downsample: zb has no anti-aliasing
/vis/tsg/offscreen/set/format zb_png         # or gl2ps_svg | gl2ps_pdf | gl2ps_eps
/vis/tsg/offscreen/set/file cutaway.png      # one file per rebuild/refresh/flush
/vis/viewer/set/lineSegmentsPerCircle 72
/vis/drawVolume
/vis/viewer/set/background white
/vis/viewer/set/viewpointVector -0.75 0.48 -0.45
/vis/viewer/set/upVector 0 1 0
/vis/viewer/set/targetPoint 0 0 0.6 m
/vis/viewer/set/projection p 25.7 deg        # optional
/vis/viewer/zoomTo 1.75                      # scene extent is the 20 m world box
/vis/viewer/addCutawayPlane 0 0 0 m 1 0 0    # keeps x > 0
/vis/viewer/set/style surface
/vis/viewer/set/hiddenEdge true
/vis/viewer/set/edge true
/vis/viewer/rebuild
```

**Gotchas:**
- **DD4hep forces the drawing style.** `drawingStyle` maps to ForceSolid, so `style wireframe` is ignored for most volumes. Fix: `/vis/geometry/set/forceSolid all -1 false`.
- **DD4hep colours mostly have alpha < 1.** With `/vis/tsg/offscreen/set/transparency false` those faces are simply skipped, which leaves a coloured wireframe, or a blank image for sections (`$G/runs/05_section/s2`).
- **`sectionPlane` also sets the viewpoint to its normal.** Normal `1 0 0` views from +x and mirrors the side view (axes evidence in `$G/runs/05b_orient`); use `-1 0 0`.
- **Only one clipper at a time.** A section plus a cutaway means no clipping at all. A slab can be made from two cutaway planes in `cutawayMode intersection`.
- **The generic section is a thin box** (1e-5 × scene radius, about 0.17 mm). No Boolean-processor warnings appeared in any run.

**(a) Hidden-line drawing:**
- TSG's hlr mode (wireframe + hiddenEdge) draws filled triangles plus edges, both in the object colour (`G4ToolsSGSceneHandler::AddPrimitive(G4Polyhedron)`). Black lines therefore give a solid black silhouette (sample 08).
- Plain wireframe works and respects edge flags, but is cluttered by the EC module prisms and does not remove hidden lines.
  - zb PNG: 3.8 s.
  - gl2ps PDF: 3.3 MB unsorted, 15 MB depth-sorted.
  - gl2ps SVG: 42.8 MB, 294 k polylines.
- **Workaround (painter's algorithm):**
  - Make the scene opaque and uniform: `/vis/geometry/set/colour all -1 0.93 0.93 0.93 1` and `forceSolid all -1 false`, with surface + edges.
  - Write `gl2ps_svg` with depth sorting and options `SILENT|OCCLUSION_CULL|BEST_ROOT|DRAW_BACKGROUND|SIMPLE_LINE_OFFSET`.
  - Run `python3 $G/bin/gl2ps_svg_hlr.py in.svg out.svg [--shaded]`: faces become white (or keep their shading), geometry edges become black, trajectories keep their colours.
  - Cost: 4.3 s, 4.3 MB SVG (16.3 k triangles + 16.8 k polylines).
  - Quality: close to `solid_lineview` (sample 02). Curved surfaces have no outline lines and some edges are faint or broken.
  - OCCLUSION_CULL deletes hidden track segments, so drawing tracks on top needs it switched off.

**(b) Hidden surface with edges:**
- Colour zb: 2.4 s at 3200×2000, 1.3–1.7 MB raw, 0.4–0.6 MB after 2× downsampling (samples 01, 03). Faint triangle seams on large faces; edges come out in the face colour.
- Opaque grey "clay" with a headlight (`/vis/viewer/set/lightsVector` = viewpoint): 2.2 s, DAWN-like look. The default light renders most faces dark (`$G/runs/13_final/p3`).
- The same clay as gl2ps PDF: 1.55 MB, clean.
- gl2ps with the transparent DD4hep colours has depth-sorting artefacts (sample 09).
  - The SVG is 4.4 MB with 28 k primitives and takes 34 s to rasterize.

**(c) Section plane:**
- `sectionPlane on 0 0 0 m -1 0 0`: under 1 s per image.
- zb PNG: sample 10.
- gl2ps PDF with fills: 85 KB, sample 12. Fill shades vary with sort order; this one is mirrored (normal +x).
- Recoloured vector outlines with tracks: 3.5 MB, sample 04b. Comparable to `solid_lineview --mode section` without fill.
- Beam-eye section at z = 440 cm works but shows triangulation seams on the caps.

**(d) Trajectories:**
- Rich trajectories, drawByCharge, an IMag ≥ 20 MeV filter and `endOfEventAction accumulate`.
- A `rebuild` after `/run/beamOn` redraws the kept events in any view (samples 03, 04a, 04b).
- The offscreen viewer writes no file per event; you must rebuild after the run.
- The gun is limited to x > 0 by `$G/bin/steer_xpos.py` (phi, theta and seed from environment variables).

### 3. Animation (`$G/runs/07_anim`, `$G/runs/07b_animspeed`)

```
/vis/modeling/trajectories/drawByCharge-0/default/setTimeSliceInterval 0.05 ns
/vis/viewer/set/timeWindow/displayHeadTime true -0.85 0.85 24 0 0 0
/vis/viewer/set/timeWindow/fadeFactor 0.3
/run/beamOn 1
# for each frame, t = 0.5*i ns:
/vis/viewer/set/timeWindow/startTime 0 ns <t> ns
/vis/tsg/offscreen/set/file frames/frame_<i>.png
/vis/viewer/rebuild
```

- 60 frames at 1280×800 with the cutaway took 204 s: 3.4 s and about 450 KB per frame. Time windows and the head-time label work; the light front is OpenGL-only and is not drawn.
- **Bug:** `G4ToolsSGOffscreenViewer::DrawView` never calls `ProcessTransients()`. After a time-window change, `refresh` or `flush` writes identical stale frames (same file for 5, 15 and 25 ns), so every frame needs a full rebuild.
  - Per-frame cost: 2.5 s with a cutaway, 0.3 s without, under 1 s with a section.
  - A one-line upstream fix would bring cutaway frames to about 0.1 s.
- Encoding uses Pillow (the host `convert` is GraphicsMagick and has no `-layers`):
  - GIF 640×400: 0.60 MB.
  - WebP 800×500: 0.52 MB; 1280×800: 1.36 MB.
  - APNG: 6.3 MB.
  - No ffmpeg anywhere, so no MP4.

### 4. Energy colouring

**4a. Touchables** (`$G/runs/08_touch`):
1. **Pass 1, batch.** Run `$G/bin/run_batch.sh` with 3 events, seed 1001, theta 9–13°, phi −40..40°. It takes 16 s; FAEC shower hits are 306, 390 and 316 per event.
2. **Mapping, inside the container:**
   ```
   python3 $G/bin/hits_to_touchables.py $DETECTOR_PATH/$DETECTOR_CONFIG.xml pass1.edm4hep.root touch_all.mac all 0.0005
   ```
   - It runs in 4.7 s and sums energy per module from the cellID `module` field.
   - Each module id is mapped through its `DetElement` placement. The module id is not the copy number: module 953 is `FAECShower_module_631`.
   - Per module it writes:
     ```
     /vis/set/touchable world_volume_1 0 FAECShower_envelope_7 7 FAECShower_module_631 631
     /vis/touchable/set/colour r g b 1
     /vis/touchable/set/forceSolid true
     ```
   - Colours use a log heat scale; 50 modules are above 0.5 MeV.
3. **Pass 2, vis.** Same seed, so the events are identical.
   - Run `/control/execute touch_all.mac`, dim the unhit modules with `/vis/geometry/set/colour FAECShower_module 0 0.8 0.8 0.85 0.15`, then rebuild.
   - Timing: 50 touchables take 18.8 s, because each search walks the full 5.08 M-node tree; the rebuild takes 2.6 s.
   - Results:
     - Sample 05a: FAEC beam view from downstream, three coloured clusters and one track through the hole. Very readable.
     - Sample 05b: the iso view is tinted by the translucent disks in front.

**4b. Scoring** (`$G/runs/09_score`):
- ddsim never creates `G4ScoringManager`, so `/score/` does not exist. `$G/bin/steer_scoring.py` adds it through cppyy: `gInterpreter.Declare('#include "G4ScoringManager.hh"')`, then `G4ScoringManager.GetScoringManager()`.
- Box and cylinder meshes at the FAEC accumulate correctly; the CSV dump has 206 non-zero cells.
- In 11.4, `/score/draw` is renamed `/score/drawProjection`; the old name aborts the macro.
- Draw after the run, then `/vis/viewer/refresh` (not `rebuild`, which clears the score drawing).
- **What works:** cylinder r-phi projection with `logColorMap`, and box `drawColumn` (samples 06a, 06b).
- **Limits:**
  - Zero cells are drawn as an opaque plane.
  - The grid is regular, not the hexagonal modules.
  - The legend text is crude.
  - The box projection drew once (camera at +z) but not twice later (camera at −z); not investigated.
  - DDG4's EDM4hep output logs "ERROR no Geant4HitCollection: eDep" every event (harmless).

### 5. Speed and size

All times are after the ~7–9 s ddsim start-up.

| product | time | size |
|---|---|---|
| DAWNFILE + `solid_lineview` (reference) | ~13 s CPU for the prim + 2.9 s render | 5.5 MB prim, 153 KB SVG |
| TSG zb colour cutaway, 3200×2000 | 2.4–4 s | 1.3–2.6 MB (0.4–0.85 MB downsampled) |
| TSG zb clay | 2.2 s | 0.49 MB |
| TSG gl2ps SVG + recolour | 4.3 s + 0.2 s | 3.0–4.3 MB |
| gl2ps PDF clay / section | 5.5 s / 0.7 s | 1.55 MB / 85 KB |
| gl2ps wireframe | ~5 s | 3.3–43 MB |
| VRML2FILE | 2.3 s | 2.06 MB |
| animation frame (cutaway / none / section) | 3.4 / 0.3 / <1 s | ~450 KB PNG |
| touchable colouring | 0.38 s per module + rebuild | — |

### 6. Effort to productise (developer-days)

1. **`make_solid_views --backend tsg`: 3–4 d.**
   - Standard views (section_yz/xz, cutaway, lines, beam_z).
   - Supersampled PNG.
   - Vector hidden-line SVG/PDF via the recolour script.
   - Opacity and forced-style handling.
   - Camera presets matching `solid_lineview`, with zoom from the DD4hep envelope.
   - Checks on all 7 configurations.
2. **Event displays: 1–2 d.** Trajectories, filters, HepMC input, fixed seeds.
3. **Energy colouring: 2–3 d.**
   - Generic hit→touchable mapping for FAEC/LAEC/preshowers, later SPD/MRPC/Cherenkov PMTs.
   - Legend overlay with Pillow.
   - Two-pass driver.
   - An in-process DDG4 C++ action (no second pass, no per-touchable search) would add about 3 d.
4. **Animation: 2 d.** Frame and camera-path generator, GIF/WebP, ffmpeg in the container for MP4. Add 0.5 d plus Geant4 review time for the upstream `ProcessTransients` fix.
5. **Docs and a headless CI smoke test: 1–2 d.**

Total about 10–14 d; items 1 + 2 alone about 4–5 d.

Optional:
- A real hidden-line mode in `G4ToolsSGSceneHandler` (background-coloured faces + polygon offset, ~20 lines): 1–2 d plus Geant4 review.
- A VRML reader for `solid_lineview` as insurance against DAWNFILE going away: 1–2 d.

### 7. Files

**Scripts in `$G/bin/`:**
- `run_vis.sh`, `run_batch.sh`
- `steer_xpos.py`, `steer_scoring.py`
- `hits_to_touchables.py`, `gl2ps_svg_hlr.py`, `inspect_hits.py`

**Samples in `$G/samples/`:**
- 01_compare_cutaway_2x2.png
- 02_vector_hlr_vs_lineview.png
- 03_tsg_zb_cutaway_tracks.png
- 04a_tsg_gl2ps_hlr_cutaway_tracks.png
- 04b_tsg_gl2ps_hlr_section_tracks.png
- 05a_faec_energy_touchables_beamview.png
- 05b_faec_energy_touchables_iso.png
- 06a_scoring_box_proj_FAILED_cyl_rphi_ok.png
- 06b_scoring_box_column_log.png
- 07a_anim_640x400_60f.gif
- 07b_anim_800x500_60f.webp
- 07c_anim_frames_23ns_30ns.png
- 08_tsg_hlr_style_is_filled.png
- 09_tsg_gl2ps_svg_transparent_colour.png
- 10_tsg_zb_section_tracks.png
- 11_tsg_gl2ps_wireframe_nosort.png
- 12_tsg_gl2ps_pdf_section_colour.png

**Vector outputs:**
- `$G/runs/15_gl2ps_hlr/h1_hlr.svg`
- `$G/runs/15_gl2ps_hlr/h1_shaded_black.svg`
- `$G/runs/16_gl2ps_events/e1_hlr.svg`, `e1_shaded.svg`, `e2_hlr.svg`
- `$G/runs/13_final/p5_persp_clay.pdf`
- `$G/runs/11_vrml_gl2ps/g4_00.wrl`