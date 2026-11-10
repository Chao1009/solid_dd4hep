# SoLID visualization without DAWNFILE: feasibility study and plan

6 November 2026 · revision 2, after an independent review · eic_xl 26.10.0-stable (DD4hep 1.38, Geant4 11.4.1, ROOT 6.40.04), no display · nothing was committed to the repository

**Tags:** **[D]** demonstrated (an output file exists) · **[I]** inferred, estimated or read from source · **[N]** not tested.

**Sources.** Four prototype studies (polyhedra, geant4_native, events_animation, web_interactive), a synthesis, an independent review that re-ran the key claims, and this revision. The role reports are saved as `REPORT.md` in each role directory (§6). Section 7 lists every review finding and what changed. All work is in session scratch space under `/tmp`, so archiving it is the first task (Phase 0).

![overview]($SCRATCH/vis_study/final_revision/feasibility_overview_v2.png)

Panels (a)–(i) of this sheet are cited below. In this revision, panels (d) and (e) are drawn from the driver-free scene and panel (h) shows Geant4-native energy colouring.

## 0. Bottom line

| question | answer |
|---|---|
| 1. Clean line/surface drawings without DAWNFILE | **Yes [D].** DAWNFILE only writes out each visible volume's `G4Polyhedron`. A small dump calls `G4VSolid::GetPolyhedron()` itself and feeds the existing `solid_lineview`, with no vis driver. It gives the same volumes, names and colours as DAWNFILE in all 7 configurations. Images are identical for PVDIS and differ by at most 0.09 % of pixels for SoLID-Lite J/ψ. Geant4's own headless driver gives good colour renders, but its line drawings are only approximate. |
| 2. Event displays, calorimeter blocks coloured by energy | **Yes [D].** Trajectories, hits, and modules on a log energy scale, in four vector views, take about 5 s per event. This also works on the driver-free scene. |
| 3. Animation | **Yes [D].** A frame takes 0.04–0.2 s. GIF and animated SVG work in eic_xl and on the host. WebP and WebM were made on the host only. MP4 needs ffmpeg, which neither has. |

- **Recommendation:** option A, a Geant4 polyhedra scene + `solid_lineview` + an event layer. Keep DAWNFILE as a fallback and regression reference while Geant4 ships it. It is still in 11.4, so nothing forces a switch now.
- **Size of the job (revised):**
  - A prototype-quality first version takes about **7 developer-days**.
  - Production quality takes about **27 developer-days (22–35)**, about 5.5 weeks for one developer. That covers all 7 configurations, all readouts, multi-threaded runs, CI and docs, plus about 15 % for PR review.
  - Upkeep is 2–4 d per year.
  - The earlier figure of 15.5 d was the plain sum of two role estimates and left out several work items (§4, §7).
- **Decide now (§5):** D1 geometry source, D4 packaging, D7 MP4 route, and whether to spend 1 day on a spike of our own vis driver (D9).
- **Limits of the evidence:**
  - Events came from a 4 GeV e⁻ gun and 2 hand-written J/ψ-like HepMC3 events. No generator samples were run.
  - Each run had at most 5 events, on 4 shared cores.
  - Multi-threading was tested only for reproducibility of the simulated events, not with the event-display tools.

## 1. Findings

### 1.1 Drawings without DAWNFILE

**Recommended route: dump the Geant4 polyhedra (polyhedra study).**

- **How it works [D].**
  - A C++ helper repeats `G4PhysicalVolumeModel`'s depth-first walk and visibility rules, then calls `SetNumberOfRotationSteps(72)` and `GetPolyhedron()`. It keeps one mesh per unique solid; in the prototype, cling compiles it at run time.
  - It runs standalone after DD4hep's `Geant4Converter` (no physics), or inside ddsim at `/run/initialize`.
- **Fidelity [D].**
  - In all 7 configurations it writes the same drawn volumes as DAWNFILE, with the same order, names, colours, depths and vertex/facet totals. All 628 unique meshes are watertight.
  - Pixel differences in the 4 standard views:
    - **PVDIS:** 0 px.
    - **SoLID-Lite J/ψ:** 0.039–0.089 %, from two causes, both within the 0.1 % acceptance limit.
      - The x = 0 section plane passes exactly through mesh vertices. At x = 0.7 mm the section and cutaway views become identical.
      - DAWNFILE writes 600 volumes as DAWN primitives (`/Tubs`, `/Trd`, `/Box`), which `solid_lineview` tessellates itself. Because of this, lines_yz still differs by 0.067 % at x = 0.7 mm, on the thin forward SPD/MRPC plates.
  - **Independent re-run:** the `.prim` is byte-identical (md5 c2ff0423…). The largest bounding-box difference to DAWNFILE is 5e-6 mm.
- **More than DAWNFILE gives [D].**
  - It names the 330 PVDIS baffle blades, which DAWNFILE writes without a name.
  - Every drawn volume carries its DetElement path, readout and volume ID with a bit mask. For 308/308 sensitive volumes these IDs equal the ones Geant4 assigns during tracking.
- **Hit mapping [D], corrected.**
  - 821/821 hits in 3 events map to a drawn volume, and 809 of them match exactly one.
  - The 12 FAEC preshower hits also match the shared plate `FAECPreShower_Radiator_0`, whose mask covers only the system field (0xff). The most specific mask (their module, 0xffffff) resolves them.
  - The rule is therefore "most specific match", checked by a uniqueness test (§3).
- **Bounding box [D], to fix.**
  - The dump writes the world's bounding sphere (r = 17.3 m). DAWNFILE writes the sphere around the visible volumes (r = 8.02 m, centre z = +0.875 m).
  - Default views do not use it (0 px change), but `solid_lineview --mag` framing does. The dump should compute the visible-volume sphere, as `G4Scene` does.
- **Edge flags [D].** Write every edge as visible, as DAWNFILE does. Keeping Geant4's flags erased 3820 px of real outlines.
- **Cost [D].**

  | | polyhedra dump | DAWNFILE |
  |---|---|---|
  | scene file | 0.13–0.33 MB `.npz` (about 40× smaller) | 5.5–14.9 MB `.prim` |
  | parse time, J/ψ-lite / PVDIS | 0.39 s / 2.2 s | 1.25 s / 4.6 s |
  | wall time | 10.7–14.7 s standalone, 5.5–6.6 s of it cling (about 3–4 s once compiled [I]); +3.0 s inside ddsim | 15.8 s ddsim run |

- **Not covered [N].** SoLID has no reflected or replicated placements (every configuration logs "0 reflected placements"), so those code paths are untested. Phase 1b adds a unit test.

**Other routes:**

| route | result |
|---|---|
| Geant4 TSG_OFFSCREEN | [D] **Good colour and grey "clay" renders**, cutaways and sections at 3200×2000, in 2.2–4 s per image after 7–9 s of ddsim start-up (reproduced: 2.7 s). **No true hidden-line mode:** its hlr style fills faces in the line colour. A gl2ps SVG recoloured by a script gives **approximate hidden-line removal, with hidden-line leaks and seams**: back arcs leak above the endcap/Cherenkov box and at the end of the beam pipe, triangle seams cross the baffle faces, the yoke silhouette is dashed, and curved surfaces have no outlines. It is 4.3 MB against 153 KB and not publication grade (panel c). DD4hep forces solid style and alpha < 1; override with `/vis/geometry/set/forceSolid all -1 false` and `/vis/geometry/set/colour`. |
| OpenGL/Qt, RayTracer | [D] Unusable. OpenGL and Qt need X: there is no Xvfb, and Qt's offscreen mode gets no OpenGL context. RayTracer fails in both run modes. Sequential ddsim segfaults in `G4TheMTRayTracer::StoreUserActions`. With `-j 2`, the MT ray tracer issues `/run/beamOn 64000` and DDG4 aborts ("invalid object of type dd4hep::sim::Geant4Run"). This is a DDG4/RayTracer incompatibility, not worth pursuing. |
| VRML2FILE | [D] 2.3 s, 2.06 MB, names and colours, cutaway applied. It loses edge flags, so it is useful only as a fallback input format. |
| TGeo + RootCsg | [D] Not robust: cracked boolean meshes (1554–7623 open edges) and spheres stuck at 20 divisions. Even after repairs, `MagnetYoke_EndcapInterface` is 52 % too large in all 5 full-SoLID configurations. |
| JSROOT, Firebird, Phoenix | [D] Shaded WebGL only, no line art. JSROOT's SVG export sorts surfaces in the wrong order. |
| dusk/duskcut | [N] The earlier evaluation in `docs/visualization.md` found that `npdet_to_step` drops the polyhedra yoke and some booleans. |
| **Our own Geant4 vis driver** (new) | [N] `G4VGraphicsSystem`, `G4VSceneHandler` and `G4VViewer` subclasses, registered from a DDG4 plugin, would write `scene.npz` plus trajectory points with times. That makes it a DAWNFILE successor we own, built only on public Geant4 API [I]. It keeps vis-manager behaviour: culling, `/vis/geometry`, touchables, cutaways, trajectory filters, and times from rich trajectories. It could replace both the DDG4-internal walker and the step recorder; in multi-threaded runs, events are drawn in the vis sub-thread [I]. It needs a vis manager in the job, and volume IDs would still come from DDG4 or from names. Worth a 1-day spike (D9). |
| pyvista/VTK offscreen (new) | [N] Anti-aliased colour or clay renders and camera-path animations from `scene.npz`, with no ddsim session. Not in eic_xl [D]; it needs VTK wheels with OSMesa or EGL. Only worth doing if TSG colour renders are not enough. |
| glTF → Blender (new) | [N] Freestyle line art, Cycles renders and fly-throughs from the same meshes (Phase 4 glTF writer). Outreach quality; not needed for papers. |

**DAWNFILE status [D].** DAWNFILE is still built and shipped in Geant4 11.4 (last code change 2023); only HepRepFile and Qt3D were removed. The reasons to move are names for extruded solids, volume IDs, scenes about 40× smaller, and a code path we own.

### 1.2 Event displays with energy colouring

**Recommended route: an event layer on `solid_lineview` (events_animation study).**

- **Trajectories [D].** A 116-line DDG4 stepping action, loaded with `--action.step`, writes thinned points: one every 20 mm or 0.04 rad (sagitta ≤ 0.1 mm), each with global time and kinetic energy.
  - **Overhead (corrected).** DDSim's own timer gives 4.13 s without the recorder and 4.22 s with it, for 5 events (+2 %). These are single runs on a loaded machine and cannot resolve a 5 % effect. Phase 2a benchmarks at least 100 events.
  - **Size (corrected).** As run, the recorder wrote 2.2 MB per event for the e⁻ gun and 3.8 MB per event for the J/ψ-like events. A 5 MeV cut, applied by filtering the e⁻-gun file afterwards, gives 0.56–0.9 MB per event. Hadronic DIS events will be larger. The 5 MeV cut becomes the default.
- **Mapping [D].**
  - Calorimeter cellIDs are mapped to their module DetElement with DD4hep's VolumeManager: 0 mismatches over 3274 hits.
  - Modules are matched to drawn polygons by name: 2904/2904. Energy is summed per module.
- **Rendering [D].** 1034 lines on top of `solid_lineview` draw four views:
  - a section in the event plane;
  - a zoom on the hit calorimeter;
  - a beam-view section;
  - a hidden-line cutaway.

  Each view has tracks coloured by particle, hit markers, modules on a log energy scale, a colour bar and a legend. Four views in SVG, PNG and PDF take 4.8–6.0 s per event; each PDF is 40–155 kB. Panels (d)–(f).
- **On the driver-free scene [D] (reviewer's re-run).** Drawn on the polyhedra dump's `.prim`, the module map is 2904/2904. The side and cutaway images are pixel-identical to the DAWNFILE-based ones (0 of 2.9 M and 0 of 3.8 M px); panels (d) and (e). The `.npz` loader path has not been tried with the event layer yet [N]; it is now a Phase 2a check.
- **Multi-threading [D] (new).**
  - With `--random.enableEventSeed`, a 1-thread run (`G4RunManager`) and a 2-thread run (`G4MTRunManager`) gave identical events. All 4 events had the same primary momentum and MCParticle count, and the same hit count and energy sum in every collection.
  - Without that flag, the two run managers use different random streams.
  - The 2-thread run spent about 54 s in shutdown; not investigated.
- **Limits seen [D].**
  - The cutaway colours only the modules lying on the cut.
  - Hit markers are not depth-tested.
  - Beam-section z positions are hard-coded.
  - The recorder is single-threaded.

**Other routes:**
- **Geant4 native [D].**
  - Trajectories draw in both raster and vector output. `/vis/scene/add/hits` draws nothing, because DDG4 hits have no `Draw` method.
  - Energy colouring works through a generated touchable macro, run as a second pass with the same seed. It costs 0.38 s per module, because each lookup searches 5.08 M volumes, and it has no legend (samples 05a and 05b, panel h).
  - A single-pass version would remove both the second pass and the search [N]: an end-of-event action that calls `G4VVisManager::Draw` on module solids with energy colours, or a `Draw()` method for DDG4 calorimeter hits.
  - Scoring meshes need a cppyy hook and look crude.
- **Web [D].** JSROOT and the hosted Phoenix app show modules coloured by energy (panel i). Firebird's built-in hit painter ignores energy, so Firebird would need a SoLID painter [I].

### 1.3 Animation

- **`solid_lineview` frames [D].** In each frame, tracks stop at time t with an interpolated head, hits appear at their own time, and modules accumulate energy by contribution time. Geometry layers are cached, so a frame takes 0.04–0.07 s at 0.6 scale and 0.10–0.21 s at full size. Panel (g).

  | output, 120 frames | size | where it can be made |
  |---|---|---|
  | GIF (Pillow) | 222–357 kB | host and eic_xl [D] |
  | animated SVG (SMIL) | 0.40–0.85 MB; headless Firefox shows the same picture as the frames | host [D]; eic_xl [I] (plain text, no encoder) |
  | lossless WebP | 0.23–1.1 MB | host only: Pillow in eic_xl cannot write WebP [D] |
  | WebM (VP9) | 135–495 kB | host only (GStreamer): eic_xl has neither GStreamer nor ffmpeg [D] |
  | MP4 (H.264) | – | nowhere: no ffmpeg, and GStreamer has no H.264 encoder [D] |

- **Geant4 time windows [D].**
  - Time slices and the head-time label work; the light front does not.
  - `refresh` writes stale frames, so every frame needs a full rebuild: 2.5–3.4 s per frame with a cutaway, 0.3 s without.
  - From the source, `G4ToolsSGOffscreenViewer::DrawView` never calls `ProcessTransients()` [I]. The upstream fix would be one line.
- **Firebird [D].** It has a native time slider in ns, and tracks grow over time. Motion is coarse when a track has few points. 14 headless frames made a 333 kB GIF.
- **MP4 (D7).**
  - Adding ffmpeg to `containers/Dockerfile` only helps the optional solid_xl image, which has no build pipeline. `solid-shell` users and CI run eic_xl.
  - **Options:**
    - detect an ffmpeg binary at run time: on `PATH`, through `$SOLID_FFMPEG`, or a user-level static build such as the one in the imageio-ffmpeg wheel [N];
    - encode on the host, which has the FFmpeg 4.4 libraries and libx264 but no `ffmpeg` command [D];
    - ask for ffmpeg in eic_xl upstream.
  - x264 is GPL, which is acceptable for a separately installed tool [I].

## 2. Options compared

| option | publication figures | events + calorimeter energy | animation | headless / CI | interactive | effort (d) | maintenance risk |
|---|---|---|---|---|---|---|---|
| **A. Polyhedra scene + `solid_lineview` + event layer (recommended)** | vector line art, same as today [D] | yes, 4 vector views, also on the driver-free scene [D] | yes, 0.04–0.2 s/frame [D] | yes. numpy 2.2, scipy 1.17, pycairo 1.24, Pillow 12.2 and uproot 5.7 are in eic_xl [D]. | no (static or animated SVG) | first version 7; production ≈ 27 (22–35) including review | low–medium. About 2–2.5 k lines we own. The ID map uses DDG4 internals (`g4Placements`, `g4VolumeImprints`). |
| A2. As A, but scene and trajectories written by our own Geant4 vis driver | as A [I] | as A; trajectory times from rich trajectories [I] | as A [I] | yes [I] | no | spike 1; if adopted, roughly the same cost as A's walker + recorder [I] | low. Public Geant4 vis API; needs a vis manager in the job. |
| B. Keep DAWNFILE `.prim` + same event layer | as today [D] | yes [D] | yes [D] | yes | no | ≈ 18 (A without the geometry phases) [I] | medium. Unnamed extruded solids, no IDs, depends on a driver Geant4 may drop. |
| C. Geant4 TSG_OFFSCREEN | good colour/clay; line art only approximate [D] | tracks yes; energy via touchables, 0.38 s/module, no legend [D] | 2.5–3.4 s/frame [D] | yes; one ddsim per batch | no (needs X) | 10–14 | medium. Vis commands change (`/score/draw` renamed in 11.4); refresh bug; DD4hep forces styles. |
| D. JSROOT page | slides only | yes [D] | no | screenshots need a browser | rotate, clip, orthographic | 2–3 | low (ships with ROOT, MIT) |
| E. Firebird SoLID app + JSROOT fallback | slides only | tracks and hits [D]; energy painter to write [I] | yes, interactive [D] | browser | full | 10–12 (core 8–9) | medium–high. Beta, one main developer, licence unclear. |
| F. Phoenix app | slides only | yes [D] | intro only | browser | full, including VR/AR | 6–9 | medium–low (HSF) |
| G. TGeo/RootCsg meshes → `solid_lineview` | one boolean still wrong after repairs [D] | as A | as A | yes | no | A + 5–8 | high |
| H. REve, DDEve, OpenGL, RayTracer | – | – | – | blocked [D]. REve crashes inside civetweb when a browser connects; DDEve and OpenGL need X; RayTracer is incompatible with DDG4. | – | – | – |

A and B are re-based in this revision. The efforts for C–G are the investigators' own estimates and were not re-based. They share some of the gaps the review found in A's estimate: review cycles, all configurations, and validation.

## 3. Recommended architecture

```
ddsim + steering/solid_sim.py
 |- SOLID_SCENE=f.npz -> scene producer: G4 walk + G4Polyhedron + volume IDs --------> scene.npz
 |                       (same code standalone: bin/solid_scene_dump <compact>, no physics;
 |                        or our own vis driver if the D9 spike succeeds)
 |- SOLID_TRAJ=1      -> SoLIDTrajectoryRecorder (DDG4 stepping action, per-thread files) -> traj (x,y,z,t,Ekin)
 '- EDM4hep output: MCParticles, SimTracker/SimCalorimeterHits + contributions -----------> sim.edm4hep.root
bin/solid_evexport  (uproot or podio; hit -> most specific drawn volume) ----------------> event.npz
solid_lineview      (importable module + CLI; reads .prim or scene.npz)
bin/solid_evview    (views, legend, colour bar, event selection/overlay, --animate)
                     -> SVG/PDF/PNG, GIF, animated SVG; WebP/WebM/MP4 when an encoder is found
bin/make_solid_views: scene by default, --dawnfile fallback, --events N -> solid_evview
optional web/: glTF from scene.npz (or dd_web_display export); event.npz -> JSROOT/Firebird/Phoenix
```

- **Hit → drawn volume.**
  - Among the drawn volumes of the hit's readout with `(cellID & mask) == volID`, take the one with the largest mask, from sensitive or DetElement-level volumes only.
  - A second match with the same mask is an error, caught by a CI test.
  - For calorimeter modules, the VolumeManager route (`lookupDetElement(cellID)` → module → drawn volume by name) is an independent cross-check.
- **Scene schema v1**, as in the prototype plus two metadata fields:
  - unique meshes: float64 vertices; facets as 1-based indices, with the sign carrying the edge flag;
  - placements, in depth-first order: name, copy number, depth, parent, mesh index, 3×4 matrix, RGBA, visible and wireframe flags;
  - IDs: DetElement path, readout, volume ID, mask, sensitive flag;
  - metadata: configuration, software versions, number of segments, plus the **bounding sphere of the visible volumes** and a **geometry hash** for golden references.
- **New code:**
  - a small DDG4 plugin library (scene producer and recorder);
  - the `SOLID_SCENE` and `SOLID_TRAJ` switches in `solid_sim.py`;
  - `solid_scene_dump`, `solid_evexport` and `solid_evview`;
  - a fatal CI step;
  - docs.
- **Reused from `solid_lineview` 0.3.0:**
  - **Used as is:** sections, hidden-line removal, cutaways, camera presets, and SVG/PDF/EPS/PNG output.
  - **Small additions:** `load_scene()`, float64 vertices, a small shift of section planes that pass exactly through vertices, per-volume colour and alpha overrides, and hollow markers.

## 4. Plan and effort (revised)

**[first]** marks steps that are cheap and high value. Efforts are in developer-days.

**Phase 0 [first]: set-up, 0.5 d**
- Archive the prototypes (§6). Decide D1, D4 and D7, and whether to run spike S (D9).
- **Done when:** the archive is outside `/tmp` and D7 is settled.

**Phase 1a [first]: driver-free geometry, prototype quality, 2.5 d**
- **Deliverables:**
  - `bin/solid_scene_dump` from the prototype, still compiled by cling (11–15 s), writing the visible-volume bounding sphere;
  - an `.npz` loader in `solid_lineview`;
  - `make_solid_views --source scene|dawnfile`.
- **Done when,** for all 7 configurations:
  - volume counts, names, colours and depths equal DAWNFILE's;
  - the 4 views differ by at most 0.1 % of pixels, and by none for PVDIS;
  - `--mag` framing equals DAWNFILE's;
  - `docs/images` can be regenerated;
  - every hit in 3 events of 2 configurations has exactly one most-specific volume.

**Spike S: our own vis driver (D9), 1 d, right after 1a**
- A minimal scene handler, registered from a steering hook, that writes `scene.npz` and trajectory points with times.
- **Done when:** its scene equals the walker's for J/ψ-lite and its points carry times, or there is a written reason why not. The result picks the producer for Phases 1b and 2a.

**Phase 1b: geometry, production, 4 d**
- **Deliverables:**
  - the producer compiled ahead of time as a DDG4 plugin, with the `SOLID_SCENE` hook, a standalone tool and an optional `.prim` writer (2 d);
  - a unit test with a reflected and a replicated placement (0.5 d);
  - a separate, fatal CI step (1 d) with:
    - a comparison against DAWNFILE on 2 configurations (counts, names, image difference);
    - the hit-uniqueness test;
    - golden scene hashes per geometry version, so regression tests survive DAWNFILE's removal;
    - diff images uploaded as artefacts.
  - docs (0.5 d).
- **Done when:**
  - the dump takes at most 5 s standalone and at most 3 s extra inside ddsim;
  - the scene from the hook equals the standalone scene;
  - the fatal step passes for all 7 configurations.
- **Risk:** the ID map uses DDG4 internals. The weekly CI against new eic_xl tags would catch changes.

**Packaging P (D4): 1.5 d**
- Make `solid_lineview` importable, as a package or as a module with a thin CLI, with no importlib or pickle tricks, and add unit tests.
- **Done when:** `solid_evview` imports it and the CLI output is unchanged.

**Phase 2a [first]: static event figures, prototype quality, 4 d**
- **Deliverables:**
  - the recorder (single-threaded), with `SOLID_TRAJ` and a 5 MeV default cut;
  - an exporter for the four EC readouts, using the most-specific rule;
  - side and beam views with legend and colour bar, drawn from `scene.npz`;
  - `make_solid_views -n` switched to the recorder.
- **Done when:**
  - module energy sums equal the hit energy sums;
  - 100 % of EC hits map uniquely over at least 100 events;
  - figures from `scene.npz` equal the figures from the DAWNFILE `.prim`;
  - the recorder costs at most 5 % of DDSim's "Processing and Init" time over at least 100 events;
  - the size per event at the default cut is measured;
  - you sign off 3 sample events: FAEC, LAEC and J/ψ.

**Phase 2b: event displays, production, 6 d**
- **Deliverables:**
  - a thread-safe recorder (per-thread files) and multi-thread checks (1 d);
  - the other readouts: SPD, MRPC, GEM, muon (DDVCS), and Cherenkov PMTs with optical physics (1 d);
  - generalised zoom and cutaway views, with section positions taken from the scene (1 d);
  - all 7 configurations (1 d);
  - event selection, for example "J/ψ → e⁺e⁻ in acceptance", and a multi-event overlay (1 d);
  - a CI smoke test and docs (1 d).
- **Done when:**
  - with `--random.enableEventSeed`, `-j1` and `-j2` runs give identical figures (the simulated events are already identical, §1.2);
  - the per-thread files are complete (record counts and energy sums agree with EDM4hep);
  - 4 views take at most 6 s per event;
  - the smoke test (1 event, 2 views) passes.

**Validation V: physics samples, 1.5 d**
- Real J/ψ generator events, DIS with hadrons, and background, at least 100 events each. Measure recorder overhead and size, and review the figures.

**Phase 3: animation, 2.5 d**
- **Deliverables:**
  - `solid_evview --animate`;
  - GIF and animated SVG always, tested in CI;
  - WebP where Pillow supports it;
  - WebM and MP4 through a detected ffmpeg, with GStreamer as the WebM fallback;
  - docs.
- **Done when:**
  - 120 frames take at most 15 s per view at 0.6 scale;
  - the GIF is at most 1 MB;
  - the last frame matches the static view;
  - where ffmpeg exists, the MP4 plays in PowerPoint and Keynote.

**Optional:**
- **Phase 4, web display:**
  - a JSROOT page (2–3 d) or a Firebird app (8–12 d);
  - a glTF writer for `scene.npz` (about 1 d [I]), so that Phoenix, three.js and Blender use our exact Geant4 meshes instead of their own CSG code. This would settle Q4. Firebird's glTF import was not checked [N].
- **Moving camera:** +1 d. It needs hidden-line removal on every frame, about 2 s per frame [I].
- **Phase 5:** a native Geant4 colour/clay backend and the upstream reports (D8), 3–4 d plus review.

**Totals:**

| scope | developer-days |
|---|---|
| [first] path: 0 + 1a + 2a, plus the prototype's GIF output | **7**, prototype quality: SoLID-Lite and EC only, single-threaded, cling start-up, CI not gated |
| production: 0, 1a, S, 1b, P, 2a, 2b, V, 3 | 23.5 (range 19.5–29) |
| PR review and iterations with JLab upstream, about 15 % (15–20 % for the range) | 3.5 |
| **total, production** | **≈ 27 (22–35)**, about 5.5 weeks |
| upkeep | 2–4 per year: DD4hep/Geant4 drift (notably the DDG4 internals), new eic_xl tags |
| optional | JSROOT +2–3 · Firebird +8–12 · glTF +1 · moving camera +1 · Phase 5 +3–4 plus review |

**Why this is higher than 15.5 d.** The earlier figure was 5.5 d (polyhedra) + 10 d (events_animation), with no overlap removed. Thread-safety, the other readouts, the CI smoke test and an importable renderer were already in those role estimates. Items no role had budgeted, now added:
- multi-thread checks;
- Cherenkov PMTs;
- checks in all 7 configurations;
- packaging;
- physics-sample validation;
- golden references and a fatal CI step;
- the reflection/replica test;
- event selection;
- the spike;
- review cycles.

## 5. Decisions and open questions

**Decisions:**
- **D1. Default geometry source.**
  - Make the polyhedra scene the default, and keep DAWNFILE as fallback and regression reference while Geant4 ships it.
  - Keep a `.prim` writer for dawn/dawncut users.
  - *Recommended: yes to both.* Option B saves about 9 d, but keeps the unnamed solids and the dependence on the driver.
- **D2. Scene format.** `.npz` v1 now. Add a glTF export later if a web viewer or Blender should use our meshes.
- **D3. Trajectory storage.** Options:
  - a sidecar file, as in the prototype;
  - an EDM4hep collection (+1.5 d), which travels with the simulation output;
  - Firebird's DDG4 writer `FirebirdTrajectoryWriterEventAction` (Geant4 trajectories with times), downloaded but not built because the permission policy blocked it; evaluate it before writing our own;
  - our own vis driver (D9).
- **D4. Packaging.** Keep a single-file CLI that `solid_evview` imports by path, or make it an installed Python package (P, 1.5 d)? *Recommended: package.*
- **D5. Scope.**
  - Detectors: EC first, then SPD, MRPC, GEM and muon. PMTs need optical physics.
  - Configurations: SoLID-Lite first?
  - Physics samples: J/ψ generator events, DIS, background mixing?
  - Event selection and overlay.
- **D6. Interactive web display.** None, a JSROOT page, Firebird or Phoenix? Hosted where, for example JeffersonLab GitHub Pages?
- **D7. MP4 route.** Settle it in Phase 0.
  - Options: (a) detect ffmpeg at run time (`PATH`, `$SOLID_FFMPEG`, a user-level static binary); (b) ask for ffmpeg in eic_xl upstream; (c) encode on the host.
  - *Recommended: (a) and (b).* GIF and animated SVG are the formats CI tests; MP4 is optional.
- **D8. Upstream reports.**
  - **Geant4:** the offscreen `ProcessTransients` bug; a TSG hidden-line mode; DAWNFILE writing extruded solids without a name.
  - **DD4hep:** VisAttr drawingStyle is not passed to Geant4; `Geant4VolumeManager::volumeDescriptor(path)` returns "Insensitive" for the 300 SPD pads placed through assemblies.
  - **eic_xl:** the REve/civetweb crash; pyrobird 2026.03rc1 writing empty EDM4hep events; ffmpeg (D7).
- **D9. Scene producer.** The proven walker, which uses DDG4 internals for IDs, or our own vis driver, which uses public Geant4 API and may also replace the recorder? *Recommended: the 1-day spike S, then decide.*

**Open questions:**
- **Q1. Wireframe drawing style.** DDG4 maps a VisAttr's *lineStyle* to Geant4's forced solid or wireframe and ignores *drawingStyle*. The Cherenkov tanks' `drawingStyle="wireframe"` therefore never reaches DAWNFILE [D]. Should the dump honour drawingStyle, which changes today's pictures, or stay identical to DAWNFILE?
- **Q2. MCParticle records [D].** FAEC showers keep 5–11 MCParticles, LAEC showers about 1000 or more; the multi-thread test shows the same pattern.
  - **Likely cause [I, not verified]:** ddsim's `Geant4TCUserParticleHandler` with `tracker_region_zmax = 3.5 m` and `rmax = 2.6 m`. The FAEC (z ≈ 4.4 m) lies outside that region, the LAEC (z ≈ −0.6 m) inside.
  - **To decide:** is this intended for production? The recorder route does not depend on it.
- **Q3. Multi-threading.** Will event displays be made from multi-threaded runs? With `--random.enableEventSeed` they can be checked against 1-thread runs [D].
- **Q4. Booleans in web viewers [N].** JSROOT and Firebird use their own CSG code. A glTF export of our meshes would make this question moot.
- **Q5. Firebird licence.** The repository says GPL-3.0, the npm package says MIT. Which applies?
- **Q6. dusk.** Revisit it once `npdet_to_step` handles polyhedra and booleans [N]?

## 6. Files

**Role reports**, saved verbatim in this revision from the inline reports in the workflow journal:
- `$SCRATCH/vis_study/polyhedra/REPORT.md`
- `$SCRATCH/vis_study/geant4_native/REPORT.md`
- `$SCRATCH/vis_study/events_animation/REPORT.md`
- `$SCRATCH/vis_study/web_interactive/REPORT.md`

**Prototype code to start from:**
- **Scene dump:**
  - `$SCRATCH/vis_study/polyhedra/solidvis_g4.h`
  - `$SCRATCH/vis_study/polyhedra/g4poly.py`
  - `$SCRATCH/vis_study/polyhedra/polydump_steer.py`
  - `$SCRATCH/vis_study/polyhedra/lineview_npz.py`
  - Scenes for all 7 configurations: `$SCRATCH/vis_study/polyhedra/out/all/`
- **Hit-mapping check:** `$SCRATCH/vis_study/reviewer/scripts/nmatch.py`. It tests uniqueness; the prototype's `hitmap_check.py` only picks the deepest match.
- **Events and animation:**
  - `$SCRATCH/vis_study/events_animation/plugin/SolidStepRecorder.cpp`
  - `$SCRATCH/vis_study/events_animation/scripts/extract_event.py`
  - `$SCRATCH/vis_study/events_animation/scripts/evdisp.py`
  - `$SCRATCH/vis_study/events_animation/scripts/encode.py`
  - `$SCRATCH/vis_study/events_animation/scripts/run_all.sh`
- **Geant4 native:**
  - `$SCRATCH/vis_study/geant4_native/bin/run_vis.sh`
  - `$SCRATCH/vis_study/geant4_native/bin/gl2ps_svg_hlr.py`
  - `$SCRATCH/vis_study/geant4_native/bin/hits_to_touchables.py`
- **Web:**
  - `$SCRATCH/vis_study/web_interactive/scripts/edm4hep2web.py`
  - `$SCRATCH/vis_study/web_interactive/scripts/make_display_geo.py`
  - `$SCRATCH/vis_study/web_interactive/pages/jsroot_event.html`
  - Firebird plugin build files, prepared but not run: `$SCRATCH/vis_study/web_interactive/scripts/fbplugin/`

**Best images** (all under 1 MB):
- Overview sheet: `$SCRATCH/vis_study/final_revision/feasibility_overview_v2.png`
- PVDIS cutaway without DAWNFILE, identical to the DAWNFILE render: `$SCRATCH/vis_study/polyhedra/img/g4_PVDIS/cutaway.png`
- Event on the driver-free scene: `$SCRATCH/vis_study/reviewer/ev/jpsi1_g4poly_side.png`
- The same event as a cutaway: `$SCRATCH/vis_study/reviewer/ev/jpsi1_g4poly_cutaway.png`
- Beam view with module energies: `$SCRATCH/vis_study/events_animation/out/jpsi0_beam.png`
- Animation, 120 frames: `$SCRATCH/vis_study/events_animation/out/jpsi1_side_anim.gif`. The `.webm` and `.svg` next to it are the same animation.
- Geant4-native touchable energy colouring: `$SCRATCH/vis_study/geant4_native/samples/05a_faec_energy_touchables_beamview.png`
- TSG colour render: `$SCRATCH/vis_study/geant4_native/samples/03_tsg_zb_cutaway_tracks.png`
- Web:
  - `$SCRATCH/vis_study/web_interactive/shots2/jsroot_ev2_side_ortho.png`
  - `$SCRATCH/vis_study/web_interactive/shots2/firebird_ev0_anim.gif`

**Archive for Phase 0.** Files under 5 MB, without downloads and frames: 849 files, 250 MB before compression. Code and text alone are 148 files, 2.8 MB, small enough for a `vis-prototypes` branch on the fork.

```bash
cd $SCRATCH
find vis_study -type f -size -5120k ! -path '*/dl_*' ! -path '*/jsroot_root6.40/*' ! -path '*/pyrobird_2026_10/*' \
  ! -path '*/frames*' ! -path '*/__pycache__/*' ! -name '*.pkl' | tar czf <persistent-dir>/vis_study_2026-11-06.tar.gz -T -
```

## 7. Review log

The review re-ran the core claims independently:
- The standalone dump re-run gave a byte-identical `.prim`; all 7 configurations gave 628 meshes with 0 open edges.
- PVDIS showed 0 px difference in all 4 views.
- TSG_OFFSCREEN worked with `DISPLAY` unset: a cutaway in 2.7 s, and a gl2ps SVG with 2 events in 15.1 s (6.2 MB).

Its evidence is in `$SCRATCH/vis_study/reviewer/`. This revision added a multi-thread test and an eic_xl encoder check, both in `$SCRATCH/vis_study/final_revision/`.

| # | finding (severity) | verdict | what changed, with evidence |
|---|---|---|---|
| 1 | No deliverable on disk (high) | right | This file and the four role `REPORT.md` files now exist (§6). Archiving is Phase 0, with the command in §6. |
| 2 | 15.5 d was a plain sum and left scope unbudgeted (medium) | right | Re-estimated bottom-up: ≈ 27 d (22–35) including review; the 7 d [first] path is kept, labelled prototype quality; upkeep is 2–4 d per year; a reflection/replica unit test is added (all configurations log "0 reflected placements"). Several listed items were already in the role estimates, so the increase is 11.5 d rather than the whole list (§4). |
| 3 | "2-thread = 1-thread" cannot pass (medium) | right for default seeding; criterion kept with a flag | DDSim uses `G4RunManager` for `-j1` and `G4MTRunManager` for `-j>1` (`DD4hepSimulation.py`, lines 576–585). With `--random.enableEventSeed`, `-j1` and `-j2` gave 4 of 4 identical events (`final_revision/mt/mt_compare_j1_j2.txt`). The criterion now requires that flag and adds per-thread file checks. |
| 4 | ffmpeg in the solid_xl Dockerfile misses users and CI (medium) | right, with one correction | D7 is now settled in Phase 0 and has three options (§1.3). Correction to the suggested fix: CI on eic_xl cannot make WebM (no ffmpeg, no GStreamer) or WebP (Pillow fails with "unknown file extension: .webp") [D]. CI therefore tests GIF and animated SVG; an animated GIF wrote fine in eic_xl. |
| 5 | Hit mapping is not unique for the FAEC preshower (medium) | right | Re-run: all 12 FAECPreShowerHits match 2 drawn volumes; every other collection matches 1. The rule is now the most specific match plus a uniqueness test, and the result is restated as 821 mapped, 809 unique and 12 resolved (§1.1, §3). |
| 6 | Residual at x = 0.7 mm in lines_yz (low) | right | Both causes are named (§1.1); the ≤ 0.1 % limit holds. Evidence: `reviewer/poly/diff_lines_yz_x07.png`, `lines_yz_x07_residual_zoom.png`. |
| 7 | Event layer never run on the driver-free scene (low) | right, now resolved | The reviewer's run gave 0 px difference. It is cited in §1.2, panels (d) and (e) now use it, and a `scene.npz` check is added to the Phase 2a done-when. |
| 8 | gl2ps hidden-line output described too kindly; panel (h) mis-cited (low) | right | It is now "approximate HLR with hidden-line leaks and seams", not publication grade (`reviewer/tsg_hlr_artifacts_zoom.png`). Panel (h) is now touchable sample 05a; touchables are cited as 05a and 05b. |
| 9 | Recorder overhead and size figures (low) | right | DDSim timer: 4.13 → 4.22 s (+2 %), from single runs. Sizes as run: 2.2 and 3.8 MB per event. 0.56–0.9 MB per event applies only with the 5 MeV cut, which becomes the default. A benchmark of at least 100 events is in Phase 2a. |
| 10 | Bounding box taken from the world volume (low) | right | The visible-volume sphere goes into the scene metadata, and a `--mag` framing check is in Phase 1a. Only the `/BoundingBox` line differs between the hook and standalone `.prim` files. |
| 11 | CI views step is non-fatal (low) | right | A separate, fatal step with diff artefacts is in Phase 1b (`ci.yml`: `make_solid_views ... \|\| echo "views failed (non-fatal)"`). |
| 12 | RayTracer cause incomplete (low) | right | Recorded as a DDG4/RayTracer incompatibility in both run modes (`reviewer/rt/ddsim.log`: "invalid object of type dd4hep::sim::Geant4Run"). |
| 13 | Independent reproduction (no defect) | — | Cited at the top of this section and in §1.1. |
| M | Missing options | added, briefly | Our own vis driver: A2, D9 and spike S. Single-pass native energy colouring (§1.2). Firebird's DDG4 trajectory writer (D3). MP4 without a container rebuild (D7). glTF export (Phase 4, D2, Q4). pyvista/VTK (§1.1). Golden references (1b). Event selection and overlay (2b). |

**Multi-thread test (this revision).** Runs took 14.8 s wall for `-j1` and 68.6 s for `-j2`. The comparison script is `final_revision/scripts/mt_compare.py`.

```bash
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
INST=$SCRATCH/verify_clone/install
F=$SCRATCH/vis_study/final_revision
for j in 1 2; do singularity exec $SIF bash -c ". /etc/profile; source $INST/bin/thissolid.sh solid_JPsi_lite; export SOLID_OPTICAL=0; cd $F/mt;
  nice -n 10 ddsim --steeringFile \$DETECTOR_PATH/steering/solid_sim.py --compactFile \$DETECTOR_PATH/\$DETECTOR_CONFIG.xml \
    -N 4 -j $j --random.seed 4242 --random.enableEventSeed --outputFile seed_j$j.edm4hep.root"; done
singularity exec $SIF bash -c ". /etc/profile; python3 -I $F/scripts/mt_compare.py $F/mt/seed_j1.edm4hep.root $F/mt/seed_j2.edm4hep.root"
# -> identical events: 4 of 4
```
