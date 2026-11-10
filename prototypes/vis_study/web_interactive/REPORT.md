> Verbatim copy of the web_interactive role report, which was returned inline to the workflow on 6 Nov 2026 and saved here by the final-revision step.
> Corrections in $SCRATCH/vis_study/FEASIBILITY.md section 7 (review log) take precedence over statements they name.

I couldn't write REPORT.md because the harness doesn't let subagents write report files, so the full report is included below the summary.

## Summary (web_interactive)

Tested on solid_JPsi_lite with 3 ddsim events (4 GeV e-), using headless Firefox 159 on the host (software WebGL2).

- **Clean line drawings without DAWNFILE:** not possible with any web display. They all draw shaded triangles.
  - JSROOT's SVG output sorts surfaces in the wrong depth order.
  - Its `wire` option shows triangulation edges, not clean outlines.
  - The cleanest image is a JSROOT orthographic cut view with axes. It works for slides, not as line art.
- **Geometry export:** works with both DD4hep tools.
  - `dd_web_display --export` takes 2.4 s and gives 424 kB with 5.08 M nodes; `geoConverter` takes 0.7 s.
  - Pruning the ECal layers leaves 10,371 nodes in 385 kB.
  - Both files load in JSROOT (0.7 s pruned, 2.6 s full) and Firebird (ready after 2.5 s / 5.7 s).
- **Firebird 2026.10** is the best interactive option. Its frontend comes from the PyPI pyrobird wheel and runs from a plain static server, with no backend.
  - It loads SoLID geometry, our event file, and the ddsim EDM4hep file directly in the browser (tracker hits and MC trajectories, but no calorimeter hits).
  - It has a built-in time animation in ns: 14 frames took 6.5 s and made a 333 kB GIF. It also has WebM recording and glTF export, which I did not test.
  - Gap: its hit painter ignores energy and size, so colouring ECal blocks by energy needs a small SoLID add-on. Firebird has an extension API for this.
  - The container's pyrobird 2026.03rc1 writes empty events for EDM4hep files.
- **JSROOT** (ships with ROOT): the prototype already shows tracks, hits and ECal modules coloured by summed energy, with cut planes and orthographic views. It has no animation.
- **Phoenix:** the unmodified hosted app shows our SoLID event with energy-coloured ECal cells. Showing SoLID geometry would need our own app. It has no real-time animation.
- **REve:** the scene builds in 0.2 s, but ROOT crashes when the browser connects over websocket (inside its civetweb web server), and the long-poll fallback returns 404. This happens with both a C++ macro and PyROOT, so REve is unusable in this container.
- **DDEve:** not run. It needs an X display, and it can only read DDG4 ROOT files, not EDM4hep.
- **Blocked by the permission policy, never run:** the pyrobird 2026.10 command line and a build of Firebird's DD4hep trajectory plugin (Geant4 trajectories with times). I only assessed both by reading their source; running them needs your approval.

**Effort estimate:**

| Option | Developer-days |
|---|---|
| SoLID Firebird app plus a JSROOT fallback page | about 10–12 (core 8–9) |
| JSROOT page only | 2–3 |
| Our own Phoenix app | 6–9 |
| REve | blocked until a ROOT/civetweb build with working websockets |

**Licences and risk:**
- Firebird: GPL-3.0 repo, though its npm library says MIT. It is beta with about one main developer, so we should pin versions.
- Phoenix: Apache-2.0, maintained by HSF.
- JSROOT: MIT. ROOT/REve: LGPL.

**Sample images** (in `$SCRATCH/vis_study/web_interactive/shots2/`):
- `jsroot_ev2_side_ortho.png`
- `jsroot_ev0_cutaway.png`
- `firebird_ev0_clip.png`
- `firebird_ev0_time_montage.png`
- `firebird_ev0_anim.gif`
- `firebird_edm4hep_inbrowser_ev2_panel.png`
- `phoenix_playground_ev0_crop.png`
- `reve_ev0.png`

---

# Full report: interactive and web event displays for SoLID

Role: web_interactive. Test configuration: solid_JPsi_lite, 3 ddsim events (4 GeV e- gun, `SOLID_OPTICAL=0`).
Container: eic_xl-26.10.0-stable (DD4hep 1.38, Geant4 11.4.1, ROOT 6.40.04 with JSROOT 7.11.1, pyrobird 2026.03rc1).
Host: Firefox 159 (snap) with geckodriver, run headless; WebGL2 is software-rendered (llvmpipe).

```bash
SIF=$HOME/Projects/solid/containers/eic_xl-26.10.0-stable.sif
W=$SCRATCH/vis_study/web_interactive
# $W/scripts/env.sh: . /etc/profile; source <install>/bin/thissolid.sh solid_JPsi_lite; export SOLID_OPTICAL=0
```

## 1. Verdict

| Question | Answer from the web and interactive tools |
|---|---|
| DAWN-like line drawings without DAWNFILE? | **No.** These tools draw shaded triangles with WebGL. JSROOT's SVG mode sorts surfaces in the wrong depth order (`$W/shots/jsroot_ev0_side_svgrender_raster.png`). Its `wire` option draws triangulation edges without hiding back lines (`$W/shots2/jsroot_ev1_wire.png`). The best static figure is a JSROOT orthographic cut view with axes (`$W/shots2/jsroot_ev2_side_ortho.png`): good for slides, not line art. Firebird's glTF export into Blender Freestyle is a possible route; I did not test it. |
| Tracks plus calorimeter blocks coloured by energy | **Yes.** JSROOT does it now (prototype). Phoenix does it with per-cell RGB `IrregularCaloCells`. Firebird needs a custom painter for the colouring. |
| Animation | **Yes, in Firebird.** It has a native time animation in ns (play/step; hits appear at their time; tracks grow), a WebM recorder and glTF export. I captured frames headless and made a GIF. |

**Recommendation:**
- Build a SoLID Firebird app for interactive viewing, animation and public web pages.
- Keep a JSROOT static page as the low-risk fallback for energy-coloured snapshots.
- Do not use REve (broken in this container) or DDEve (needs X11; no EDM4hep reader).

## 2. What I ran (exact commands)

### 2.1 Geometry export (container)

```bash
singularity exec $SIF bash -c "source $W/scripts/env.sh; cd $W/geo2
  dd_web_display --export -o solid_JPsi_lite.root \$DETECTOR_PATH/solid_JPsi_lite.xml                         # 2.4 s, 320 MB
  geoConverter -compact2tgeo -input \$DETECTOR_PATH/solid_JPsi_lite.xml -output solid_JPsi_lite_geoconv.root    # 0.7 s
  python3 $W/scripts/make_display_geo.py \$DETECTOR_PATH/solid_JPsi_lite.xml solid_JPsi_lite_display.root"    # 4.0 s
```

- The full file is 424 kB and has 5,083,659 physical nodes. 5.07 M of these are the absorber/scintillator layers inside the ECal modules.
- `make_display_geo.py` does two things:
  - Clears the daughters of the 4 `*Shower_module` volumes, leaving 10,371 nodes in 385 kB.
  - Maps DD4hep's run-time colour indices (1000 and above) to standard ROOT colours.
- Both files load in JSROOT and Firebird, which cap how many nodes they draw.

| Viewer | Pruned geometry | Full geometry |
|---|---|---|
| JSROOT, time to draw | 0.73 s | 2.6 s |
| Firebird, time to ready | 2.5 s | 5.7 s |

### 2.2 ddsim (one job, run with `nice`)

```bash
singularity exec $SIF bash -c "source $W/scripts/env.sh; cd $W/sim2
  nice -n 10 ddsim --steeringFile \$DETECTOR_PATH/steering/solid_sim.py --compactFile \$DETECTOR_PATH/solid_JPsi_lite.xml \
       -N 3 --random.seed 4242 --outputFile e4GeV_3evt.edm4hep.root"
```

- 13.5 s wall: 4.9 s startup plus 1.36 s per event; 925 MB peak memory; 573 kB output.
- Events 0 and 1 hit the large-angle ECal (LAEC); event 2 hits the forward ECal (FAEC).

### 2.3 Event conversion

```bash
singularity exec $SIF bash -c "source $W/scripts/env.sh; cd $W/events2
  python3 $W/scripts/edm4hep2web.py $W/sim2/e4GeV_3evt.edm4hep.root $W/geo2/solid_JPsi_lite_display.root e4GeV_3evt   # 6.5 s
  pyrobird convert -e 0-2 -o pyrobird0326.firebird.json $W/sim2/e4GeV_3evt.edm4hep.root"
python3 -I $W/scripts/dex2phoenix.py $W/events2/e4GeV_3evt.firebird.json $W/site2/e4GeV_3evt.phoenix.json   # 53 kB
python3 -I $W/scripts/dex2txt.py $W/events2/e4GeV_3evt.firebird.json 0 $W/reve2/ev0.txt
```

`edm4hep2web.py` writes two files.

**`e4GeV_3evt.firebird.json`** (221 kB, Firebird event format "DEX" 1.0). It contains these pieces:
- tracker hits as `BoxHit`
- ECal layer hits as `BoxHit` (394–455 per event)
- `*ShowerHitsModules` as `BoxHit`: energy summed per module, with each module's centre and size taken from the geometry (19–28 modules per event)
- `McHitTrajectories` as `PointTrajectory`: vertex, then tracker hits, then the first ECal contribution, each with a time
- `NeutralParticles`: straight lines sampled every 0.2 ns

**`e4GeV_3evt.jsroot.root`** (13 kB). It contains ROOT objects that JSROOT draws directly:
- `TPolyLine3D` tracks and `TPolyMarker3D` hits
- `TEveGeoShapeExtract` copies of the hit ECal modules, coloured by summed energy on a log scale from 1 MeV to 1 GeV

**Container pyrobird 2026.03rc1** only handles EDM4eic. On our EDM4hep file it exits without error but writes 3 empty events (old format 0.04, `"groups": []`).

### 2.4 Serving and headless screenshots (host)

Two static servers:

```bash
python3 -I $W/scripts/serve.py 8811 /jsroot/=$W/jsroot_root6.40 /pages/=$W/pages /solid/=$W/site2
python3 -I $W/scripts/serve.py 8812 /=$W/pyrobird_2026_10/pyrobird/server/static /solid/=$W/site2 --spa
```

- `serve.py` is a small static server with HTTP Range, CORS and an SPA fallback. No pyrobird backend is needed.
- `jsroot_root6.40` is identical to the container's `$(root-config --prefix)/share/root/js`.
- `pyrobird_2026_10` is the unpacked PyPI wheel. Only its static frontend was used (Firebird 2026.10.0-beta.1, 11 MB without sample data).
- GitHub Pages returns `206 Partial Content` with `Access-Control-Allow-Origin: *` for Range requests (curl against eic.github.io), so this kind of static hosting works there.

Screenshots use `wd.py`, my small WebDriver client for Firefox. Each shot takes 6–9 s including browser start.

```bash
python3 -I $W/scripts/wd.py shot "http://127.0.0.1:8811/pages/jsroot_event.html?geo=/solid/solid_JPsi_lite_display.root&ev=/solid/e4GeV_3evt.jsroot.root&item=ev0;1&opt=clipx" \
   $W/shots2/jsroot_ev0_cutaway.png --w 1200 --h 750 --wait-js "return window.__done === true"
#   item=ev2;1&opt=clipx,orthoZOY,axis -> jsroot_ev2_side_ortho.png ; item=ev1;1&opt=wire,clipx,zoom150 -> jsroot_ev1_wire.png
python3 -I $W/scripts/wd.py shot "http://127.0.0.1:8812/?geometry=/solid/solid_JPsi_lite_display.root&dex=/solid/e4GeV_3evt.firebird.json&event=0&config.geometry.clippingEnabled=true" \
   $W/shots2/firebird_ev0_clip.png --w 1400 --h 860 --sleep 3 --wait-js "return !!(window.firebird && window.firebird.ready)"
python3 -I $W/scripts/fb_frames.py "<same Firebird URL>" $W/shots2/fb_ev0 --frames 14 --w 1100 --h 700     # clicks "Step" (1 ns)
convert -delay 35 -loop 0 $W/shots2/fb_ev0/f0*.png -crop 640x330+120+190 +repage -colors 64 $W/shots2/firebird_ev0_anim.gif
python3 -I $W/scripts/wd.py shot "http://127.0.0.1:8812/?geometry=http://127.0.0.1:8812/solid/solid_JPsi_lite_display.root&dex=http://127.0.0.1:8812/solid/e4GeV_3evt.edm4hep.root&config.events.rootEventRange=0-2&event=0" \
   $W/shots2/firebird_edm4hep_inbrowser_ev0.png --w 1400 --h 860 --sleep 3 --wait-js "return !!(window.firebird && window.firebird.ready)"
python3 -I $W/scripts/wd.py shot "https://hepsoftwarefoundation.org/phoenix/playground?file=http://127.0.0.1:8811/solid/e4GeV_3evt.phoenix.json&type=json" \
   $W/shots2/phoenix_playground_ev0.png --w 1400 --h 860 --sleep 6 --prefs '{"network.lna.blocking": false, "network.lna.enabled": false}' \
   --wait-js "return !/Loading/.test(document.body.innerText||'') && document.querySelectorAll('canvas').length > 0" \
   --js "<12 wheel events on the canvas to zoom out>"
```

- Firebird's in-browser EDM4hep conversion needs an absolute URL. A relative `/solid/...` path makes it ask for a pyrobird backend.
- Firefox 159's local-network protection blocks the HTTPS Phoenix page from reading `http://127.0.0.1`. Without the prefs above, the request never returns (checked with `$W/scripts/probe_fetch.py`). Real hosting on a public server avoids this.

## 3. Results per tool

### 3.1 EIC Firebird and pyrobird

**What works:**
- **Geometry:** the SoLID geometry file loads from a URL (`?geometry=`) and is drawn plain grey because there is no SoLID colour theme yet.
  - Images: `firebird_ev0_clip.png`, and `firebird_geo_solid_JPsi_lite.root.png` for the full 5 M-node file.
- **Events:** our event file loads with `?dex=…&event=N`.
- **In-browser EDM4hep conversion** reads the ddsim file directly, with no server: ready in 3.5 s.
  - For event 2 it produced FAMRPCHits 14, FASPDHits 3, GEMTrackerFAHits 5, McHitTrajectories 2 and MCParticles 5. The tracker counts match my converter (`firebird_edm4hep_inbrowser_ev2_panel.png`).
  - It produces **no calorimeter hits**. The same is true of `pyrobird convert -t edm4hep` 2026.10: its source and help text list only tracker_hits, mc_trajectories and mc_particles.
- **Time animation:**
  - Controls: a 0–200 ns slider; each "Step" click advances 1 ns; play/rewind and a "Collision" intro.
  - Hits appear at their time and tracks grow (`firebird_ev0_track_zoom.png`, `firebird_ev0_time_montage.png`).
  - Speed: 14 frames took 6.5 s after about 10 s of startup; the GIF is 333 kB.
  - The animation is coarse when a track has only a few points. Two things would smooth it: Firebird's DD4hep plugin, which records Geant4 trajectory points with times, and `pyrobird smooth`.
- **Automation:** the page exposes `window.firebird.{ready, errors, …}`, and deep links accept `?geometry=`, `?dex=`, `?event=` and `?config.<key>=`.
- **Other features:** clipping (on by default), z-clipping, orthographic camera, camera presets, a Recording button (WebM), Export to GLTF and Make Screenshot. I did not test recording or export.

**What does not work out of the box:**
- **ECal energy colouring.** The built-in hit painter draws fixed 10 mm cubes, ignores each hit's size and energy, and gives each collection a random colour.
- Firebird 2026.10 is built as a library (`@dexvis/firebird-ng`, an Angular package) with an extension API: `provideFirebird`, `withLazyPainter`, `withCameraPreset`, `withGeometryTheme`, `withRootGeometryRules`.
- A SoLID add-on could replace the built-in hit painter with one that colours by energy, but that means building our own small app.

**How ePIC-specific it is:**
- The core works with any geometry file and event file.
- The ePIC-specific parts sit in one add-on package (`@dexvis/firebird-epic`): the default geometry URL, colour themes, presets, and the default node-removal rules used by `pyrobird geo process` (which `-r` can replace).
- The public app has moved from eic.github.io/firebird (now documentation only) to https://seeeic.org, which I did not test.

**Not run:**
- The DD4hep plugin (`FirebirdTrajectoryWriterEventAction`, 572 lines) writes Geant4 trajectory points with their times into the event file. I prepared `$W/scripts/fbplugin/CMakeLists.txt` and `$W/scripts/solid_firebird_steering.py`.
- The permission policy blocked compiling this downloaded C++ and running the downloaded pyrobird 2026.10 command line.

### 3.2 HSF Phoenix 4.2.x (Apache-2.0)

- **What I tested:** the unmodified hosted app loaded the SoLID event file. It showed tracks, hits, and `IrregularCaloCells` coloured by energy (`phoenix_playground_ev0_crop.png`).
- **Event format:** it also supports `CaloCells` and `PlanarCaloCells` (box length scales with energy), jets and vertices. Tracks can be point lists or track parameters that Phoenix propagates itself.
- **Geometry:** glTF, OBJ or ROOT files through app code or the import menu. The hosted app has no URL option for geometry, so SoLID would need its own app. I did not try loading SoLID geometry through the UI.
- **Animation:** only a "collision" intro and camera presets; no real-time animation.
- **Other features:** VR/AR, screenshot mode, shareable links. No EDM4hep reader.
- **Relation to Firebird:** the Firebird 2026.10 bundle no longer includes the Phoenix library; earlier Firebird versions were built on it.

### 3.3 JSROOT 7.11.1 (bundled with ROOT, MIT)

- **Best quick result.** It draws the geometry file directly. Useful draw options: `clipx/y/z`, `orthoZOY` and other orthographic views, `axis`, `projx/y/z`, `transpNN`, `camx…`, `wire`, `outline`.
- **Event overlay:** ECal modules coloured by summed energy, plus an HTML legend.
- **Timings:** geometry in 0.84 s, with the event 1.26 s, 6.3 s including browser start.
- **Images:** `jsroot_ev0_cutaway.png` (3D cutaway) and `jsroot_ev2_side_ortho.png` (side section with axes in cm).
- **Limits:**
  - No time animation (would have to be scripted) and no physics-object UI.
  - Cut solids are left open rather than capped.
  - SVG export is broken for 3D, and `wire` shows triangulation edges.

### 3.4 ROOT EVE-7 (REve)

- **Scene:** `$W/scripts/solid_reve.C` builds it in 0.2 s:
  - 73 geometry shapes
  - one track propagated through a uniform 1.5 T solenoid field
  - a hit set
  - 24 ECal boxes coloured by energy
- **Server:** I ran it with `root -l -q` (with `-b` no server starts). It printed `New web window: http://localhost:8795/win1/?key=…`.
- **Crash:** the browser loads the REve interface, but when its websocket connects, ROOT segfaults inside `libcivetweb.so.1.16.0` (a null function pointer in the websocket handler). This happens with both C++ and PyROOT (`$W/reve2/reve_cpp2.log`, `$W/reve/reve.log`).
- **Long-poll fallback:** with `&ws=longpoll` the server survives, but `win1/root.longpoll` returns 404 and the page shows "Client Disconnected" (`reve_ev0.png`, `reve_ev0_longpoll.png`).
- **Likely cause (not verified):** ROOT 6.40.04 here is linked to an external spack civetweb `1.16-213-g5864b55a`, so a version mismatch is likely. If so, every ROOT 7 web widget is affected in this container: RBrowser, RGeomViewer and REve.
- **Hosting:** REve needs a running ROOT process, so it cannot be hosted as static files.

### 3.5 DD4hep DDEve

- `ddeve -config <xml> -data <file>` needs an X display with OpenGL, and neither the host nor the container has Xvfb.
- `libDDEvePlugins.components` lists only `DDG4EventHandler`, which reads DDG4 ROOT files; there is no EDM4hep reader.
- It is legacy desktop software. Not recommended.

## 4. Capability matrix

| | Firebird 2026.10 | Phoenix 4.2 | JSROOT 7.11 | REve | DDEve |
|---|---|---|---|---|---|
| SoLID geometry | works from URL | needs own app | works | works (server side) | reads compact XML |
| EDM4hep input | tracker and MC in browser; calo via our converter | via our JSON | via our ROOT overlay | via our code | no (DDG4 ROOT or LCIO only) |
| ECal coloured by energy | needs custom painter | **yes** | **yes** | yes | yes |
| Time animation | **yes** (ns slider, WebM) | intro only | would need scripting | no | no |
| Cut planes / orthographic | yes / yes | yes / yes | yes / yes, with axes | not tested | yes |
| Headless screenshot here | **yes** | **yes** | **yes** | no (crash) | no (no X) |
| Static hosting | yes; GitHub Pages works | yes | yes | no | no |
| Licence | GPL-3.0 (repo, pyrobird); npm library says MIT | Apache-2.0 | MIT | LGPL-2.1+ | LGPL-3.0 |
| Maintenance risk | medium-high: beta, about one developer, renames and site moves | medium-low: HSF, multi-experiment, active | low | medium; broken here | high (legacy) |

## 5. Effort to productise (developer-days)

| Task | Days |
|---|---|
| Geometry export for all 7 configurations in CI (export, pruning, colour map) | 1 |
| EDM4hep-to-Firebird converter including calorimeter modules: harden the prototype, look up modules by cell ID, tests; optionally contribute upstream | 2 |
| SoLID Firebird app: energy-coloured hit painter with legend, camera presets, subsystem colour theme, per-configuration presets | 3–4 |
| GitHub Pages deployment | 1 |
| Batch screenshot/animation tool (MP4 needs ffmpeg, which is not on the host or in the container) | 1 |
| Optional: Geant4 trajectories with times via Firebird's DD4hep plugin | 1–2 |
| Documentation and buffer | 1 |
| **Total** | **about 10–12 (core 8–9)** |

Alternatives: JSROOT page only, 2–3 days. Our own Phoenix app, 6–9 days. REve, at least 3 days once a ROOT/civetweb build with working websockets exists. DDEve is not recommended.

## 6. Files

Everything is under `$W`, which is `$SCRATCH/vis_study/web_interactive`.

**Scripts** (`$W/scripts/`):
- `make_display_geo.py`, `edm4hep2web.py`, `dex2phoenix.py`, `dex2txt.py`
- `serve.py`, `wd.py`, `fb_frames.py`, `probe_fetch.py`
- `solid_reve.C`, `run_reve.sh`, `env.sh`
- prepared but not run: `fbplugin/`, `solid_firebird_steering.py`

The JSROOT page is `$W/pages/jsroot_event.html`.

**Data** (`$W/site2/`): `solid_JPsi_lite_display.root`, `solid_JPsi_lite.root`, `e4GeV_3evt.edm4hep.root`, `e4GeV_3evt.firebird.json`, `e4GeV_3evt.jsroot.root`, `e4GeV_3evt.phoenix.json`.

**Logs:** `$W/geo2/*.log`, `$W/sim2/ddsim.log`, `$W/events2/*.log`, `$W/reve2/*.log`.

**Sample images** (`$W/shots2/`, all under 0.5 MB):
- `jsroot_ev2_side_ortho.png`
- `jsroot_ev0_cutaway.png`
- `jsroot_ev1_wire.png`
- `firebird_ev0_clip.png`
- `firebird_ev0_time_montage.png`
- `firebird_ev0_track_zoom.png`
- `firebird_ev0_anim.gif`
- `firebird_edm4hep_inbrowser_ev2_panel.png`
- `phoenix_playground_ev0_crop.png`
- `reve_ev0.png`

**Earlier interrupted attempt** (14:16–15:01, same pipeline and same event content): `$W/shots`, `$W/geo`, `$W/sim`, `$W/events`, `$W/site_data`. Useful extra images there: `jsroot_ev0_side_svgrender_raster.png`, `jsroot_ortho_pair.png`, `firebird_front_clip.png`.

**Downloads** (`$W/dl_*`, `$W/pyrobird_2026_10`): reference material only. Of the downloaded code, only the Firebird frontend was run, inside the browser sandbox.