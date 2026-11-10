<!-- SPDX-License-Identifier: LGPL-3.0-or-later -->
# r-z magnet models (POISSON decks) from the DD4hep geometry

`../extract_rz_model.py` turns the DD4hep CLEO-II magnet into 2D axisymmetric POISSON Superfish decks. It reads
the DetElements `MagnetYoke` and `MagnetCoil` (plugin `SoLID_CLEOMagnet`). The decks follow
`solid_gemc/field/CLEOv10.am`, so the same file runs in the solver of this directory and in POISSON.

All iron and coil outlines come from the TGeo shapes, and no dimension is copied by hand. Only the distribution
of the current comes from a parameter file (`cleo2_coils.json`), because the DD4hep coil is a passive volume.

## Files

| file | magnet | compact read | regions (iron) | iron mass: r-z model / DD4hep |
|---|---|---|---|---|
| `solid_v4c.am` | full SoLID, yoke v4c (`solenoid_v4.xml`) | `test/test_magnet_v4.xml` | 9 (2) | 620.3 t / 642.1 t |
| `solid_lite_stage1.am` | SoLID-Lite stage 1, unmodified CLEO-II iron (`solenoid_nomod.xml`) | `test/test_magnet_nomod.xml` | 8 (2) | 458.8 t / 488.0 t |
| `solid_lite_stage2.am` | SoLID-Lite stage 2, **provisional** (`solenoid_lite_sidis.xml`) | `test/test_magnet_lite_sidis.xml` | 11 (5) | 583.3 t / 612.3 t |
| `<name>.json` | sidecar for each deck: provenance (sha256 of every compact file read), options, every DD4hep piece with its r-z outline, volumes, sampling checks, snaps and merges | | | |
| `cleo2_coils.json` | current-distribution parameter set (main coil, correction coils) | | | |
| `equal_area/<name>_equal_area.am` | the same magnets with polygonal iron at equal-area radii (`--convention equal-area`), for sensitivity studies | | | v4c 642.1 t, stage 1 488.0 t, stage 2 612.5 t (r-z) |

All three decks:
- were generated on 2026-11-06 from solid_dd4hep commit `aa23999` (branch `main`; no magnet file has changed since
  `48360d8`, from which the first decks with identical regions were made), with `--no-field`;
- read compact files with no uncommitted changes, using DD4hep 1.38 and ROOT 6.40.04 in `eic_xl-26.10.0-stable`;
- use a problem box of r 0-500 cm, with z -500..800 cm (v4c, stage 2) or -500..700 cm (stage 1, the CLEOv10.am box);
- carry the main coil at 1282 turns x 3266 A = 4 187 012 A-turns.

The full configurations give identical regions: `solid_SIDIS_He3.xml` for v4c and `solid_JPsi_lite.xml` for stage 1.
Each deck header records the compact path, the commit and the sha256 of every magnet fragment read; the sidecar
records, for every compact file read, its sha256 and a geometry digest (the sha256 without comments, `<info>` and,
for the top-level read with `--no-field`, the field includes), plus the sha256 of the plugin source
`src/CLEOMagnet_geo.cpp`. Paths are relative to the repository root, and files of the container's DD4hep appear as
`${DD4hepINSTALL}/...`. `../extract_rz_model.py check --repo <checkout> <name>.json` compares them with a checkout
(no DD4hep needed); `../make_maps.sh -n` runs it and stops if the geometry digest of a magnet file
(`compact/magnet/`, `test/test_magnet_*`) has changed since the deck was extracted.

## Regenerating

`../make_maps.sh` runs `extract` for the three decks (with `--no-field --no-sampling`) and then the solver; the
commands below are the individual steps. `extract` needs DD4hep python (the eic_xl container). `compare` and
`plot-fraction` need only numpy and matplotlib. Without `--no-field`, every field map the compact refers to must be
resolvable: either in `./fieldmaps/`, or in `${DETECTOR_PATH}/fieldmaps/`, or downloadable.

```bash
REPO=/path/to/solid_dd4hep            # a git checkout (test/ is not installed)
GEMC=/path/to/solid_gemc/field        # reference decks and the Stage 1/2 html
MF=$REPO/magnet_field
# inside eic_xl, after:  source <install>/bin/thissolid.sh ; export DETECTOR_PATH=$REPO
python3 $MF/extract_rz_model.py extract $REPO/test/test_magnet_v4.xml -o $MF/models/solid_v4c.am \
    --name solid_v4c --no-field --sample-out work/solid_v4c_fraction.npz
python3 $MF/extract_rz_model.py extract $REPO/test/test_magnet_nomod.xml -o $MF/models/solid_lite_stage1.am \
    --name solid_lite_stage1 --no-field --sample-out work/solid_lite_stage1_fraction.npz
python3 $MF/extract_rz_model.py extract $REPO/test/test_magnet_lite_sidis.xml -o $MF/models/solid_lite_stage2.am \
    --name solid_lite_stage2 --no-field --commented-addon cleov11_correction_coil \
    --sample-out work/solid_lite_stage2_fraction.npz
# equal-area variants (likewise nomod -> solid_lite_stage1, lite_sidis -> solid_lite_stage2)
python3 $MF/extract_rz_model.py extract $REPO/test/test_magnet_v4.xml -o work/solid_v4c_equal_area.am \
    --name solid_v4c_equal_area --convention equal-area --no-field --no-sampling
# any python with numpy + matplotlib
python3 $MF/extract_rz_model.py compare $MF/models/solid_lite_stage1.json $GEMC/CLEOv10.am \
    --report work/compare_stage1_vs_CLEOv10.md --plot $MF/validation/geometry_stage1_vs_CLEOv10.png
python3 $MF/extract_rz_model.py compare $MF/models/solid_v4c.json $GEMC/CLEOv9.am \
    --report work/compare_v4c_vs_CLEOv9.md --plot $MF/validation/geometry_v4c_vs_CLEOv9.png \
    --zoom -290 -190 20 100 --zoom 140 270 150 290 --zoom 440 590 0 120
python3 $MF/extract_rz_model.py compare $MF/models/solid_lite_stage2.json "$GEMC/CLEO-II Stage 1 and 2.html" \
    --ref-key geometry --report work/compare_stage2_vs_CLEOv11.md --plot $MF/validation/geometry_stage2_vs_CLEOv11.png \
    --zoom -295 -195 25 105 --zoom 180 330 0 120 --zoom 240 660 0 300
python3 $MF/extract_rz_model.py plot-fraction $MF/models/solid_v4c.json work/solid_v4c_fraction.npz \
    -o $MF/validation/geometry_v4c_fraction.png          # likewise stage1, stage2
```

Each `extract` takes about 35 s, mostly the sampling (about 20 s with `--no-sampling`). Useful options:

| option | effect |
|---|---|
| `--convention equal-area` | equal-area radii for polygonal surfaces |
| `--no-merge` | one region per DD4hep piece, e.g. for force bookkeeping |
| `--snap 0` | do not close small gaps |
| `--bars full` | full-length spacer ring |
| `--aux-ni CorrectionCoil=<A-turns>` | current of the v4c correction coil |
| `--addon-coil NAME[:NI=..|:J=..]` | add a coil from the parameter set |
| `--box-z ZMIN ZMAX` | problem box in z |
| `--no-field` | load the compact without its field fragments, so no field map is needed |
| `--no-sampling` | skip the sampling check of the analytic pieces (the bars are still sampled) |

The equal-area variants of the three decks (`equal_area/`, decks only) were made with the `--convention equal-area`
command above, for sensitivity studies: solved like the main decks, they show the size of the 2D yoke convention
(../README.md, Limitations). The maps are made from the main decks.

## Conventions

**Frame and deck format.**
- Hall frame, origin at the coil centre; x = r and y = z, in cm.
- Region 1 is the problem box. The main coil follows as six `mat=1` regions with `cur=` (A-turns of the region), in CLEOv10.am order.
- Next come auxiliary coils, if any, then the `mat=2` iron regions, sorted by z.
- `mode=0` with no MT table, i.e. POISSON's internal iron curve. `rhogam=0.02`, because the Stage 1/2 study found 0.08 can oscillate in saturated steel.
- Dirichlet boundary on all four edges, like CLEOv10.am; open boundaries are handled by the solver.
- Regions never overlap and contain no holes. Every polygon is closed, and no line is longer than 100 characters.
- Radial mesh: the CLEOv10.am zones (`xreg/kreg`, 0.1 cm at the coil), refined automatically where iron is thinner than two cells. Only v4c needs this: r 48-54 cm at 0.2 cm, for the 0.635 cm steel corrector ring. z is uniform at 1 cm.

**What is iron.** A material counts as iron if its Fe mass fraction is at least 0.9 (`Magnet_YokeSteel` 0.9935, G4 `Iron`).
- The cryostat walls and flanges (`Magnet_StainlessSteel`, Fe 0.68) are air.
- So are the Al coil shell and the conductor.
- Each deck header lists these non-magnetic volumes.

**Polygonal (TGeoPgon) surfaces** are placed at the apothem ("flat" distance) by default, as in CLEOv8/9/10.am:
- octagon flats 176.53 / 212.62 / 221.51 / 257.61 cm;
- collar outside 180.34 cm;
- the inner side of the EndcapInterface 257.61 cm.

TGeoPgon radii are apothems; the sampled volumes confirm it (they match n tan(pi/n) a^2 to within 0.04 %).

With `--convention equal-area` every polygonal radius is multiplied by sqrt(n tan(pi/n)/pi) = 1.02704 (n = 8). That keeps the cross-section area of each piece, for example 257.61 -> 264.58 cm.

Consequences of the flat convention:
- Octagon-octagon annuli carry pi/(8 tan 22.5 deg) = 0.948 of the real iron.
- An octagonal *inner* surface with round iron outside (EndcapInterface) gains iron: 14.3 t against 10.8 t.
- Net effect: -29 t for stage 1 and stage 2, -22 t for v4c.
- The equal-area variants reproduce the DD4hep masses: 488.0 / 642.1 / 612.5 t.

**Coaxial booleans** (round bore in an octagonal collar, octagonal bore in a round interface) are computed exactly in r-z. Each TGeoCompositeShape leaf is converted in its own convention first.

**Geometric check of every piece.** Each piece is sampled with `TGeoShape::Contains_v`:
- annular cells of 0.25 x 0.5 cm (r x z), 360 phi points each;
- low-discrepancy r and z positions inside each cell;
- a golden-ratio phi offset per radius, against aliasing with the 8-fold structures.

The sampled volumes agree with the analytic ones (or, for the bars, with the finer ring sampling) to within 0.09 % for all 35 pieces of the three decks. The sampled azimuthal fraction map is plotted in `../validation/geometry_*_fraction.png`. It shows the octagon corners as partial fill outside the flats, and shows that the deck outlines sit at the flats.

**Discrete bars** (any off-axis volume, here the 8 spacer bars per end, 149.86 x 8.89 x 34.04 cm) become an equivalent full ring:
1. Ray casting (`DistFromOutside` / `DistFromInside`) through the centroid of each copy gives the normal radial and axial extents: r 212.62-221.51 cm, z -266.70..-232.66 cm. These agree for all 8 copies.
2. The iron volume comes from (r, phi, z) sampling: 362 894 cm3 per end, against a TGeo capacity of 362 800 cm3 (+0.026 %). The sampling also finds no material outside the ray-cast extents.
3. Azimuthal coverage = sampled volume / volume of the full model ring = **0.8793**.
4. As in the decks, the ring keeps its full radial extent and its length is scaled: 34.04 -> **29.93 cm**, fixed at the end with the larger |z|. The rings are z -266.70..-236.77 and +236.77..+266.70 (stage 1 and 2), and +179.30..+209.23 downstream in v4c.

Comparison with the decks: they used "90 % of the circumference", 34 -> 30.6 cm. That is close to the bar width over the circumference at the inner face: 8 x 149.86 / (2 pi 212.62) = 0.8975.

In the equal-area variant the ring sits between the scaled flats, and the coverage becomes 0.8336, the fraction of the real octagonal gap (28.38 cm).

**Small gaps and offsets** below `--snap` (0.5 cm) between neighbouring iron pieces are closed:
- coordinates snap to the earlier piece;
- a piece never loses one of its own dimensions;
- every snap is listed in the deck header.

Only stage 2 needs snaps:

| piece | snap | reason |
|---|---|---|
| UpstreamShield1-3, FrontPiece | r 144.00 -> 144.14 | the DD4hep pieces leave a 1.4 mm radial air gap to the collar bore (v0 numbers 144.0 against the real bore 144.14) |
| UpstreamShield3 | z -266.50 -> -266.70 | flush with the yoke face |
| EndcapBottomInner, EndcapBottomOuter | z 530.01 -> 530.00, 545.01 -> 545.00 | the 0.1 mm gaps Geant4 needs between the donut and the plates |

After snapping no iron feature is thinner than 0.5 cm. The thinnest is the v4c corrector ring at 0.635 cm; the 0.4 cm coil sheets are deliberate.

**Merging.** Touching iron pieces are merged greedily: largest contact length first, and among equal contacts the smaller merged area first. A merge is made only while the union stays one simple polygon without holes. The air gap between the two yoke layers, closed by the spacer rings at both ends, therefore keeps the layers in separate regions:

| deck | iron regions |
|---|---|
| stage 1 | [InnerOctagon + both collars], [OuterOctagon + both spacer rings] |
| v4c | [InnerOctagon + collars + plug + corrector ring], [OuterOctagon + interface + spacer rings + endcap cylinder + nose + plates] |
| stage 2 | stage-1 groups (FrontPiece and Shields 1-3 join the inner group), Shield4, Shield5, [donut + bottom plates] |

Each region's comment lists its pieces. `--no-merge` keeps one region per piece.

## Coil current

The DD4hep `CoilWinding` is r 151.70-154.90 cm, z +-173.75 cm: 2 layers x 1.6 cm of Al-stabilised conductor. The POISSON decks put the current in sheets 152.3-152.7 and 153.9-154.3 cm. These are 0.4 cm thick and centred on the two layer mid-radii, 152.5 and 154.1 cm, where the Cu/NbTi insert in the middle of the conductor carries the DC current.

The mapping is therefore:
1. split the winding into `len(layer_names)` = 2 equal layers;
2. put each layer's current in a sheet of `sheet_thickness_cm` = 0.4 cm centred on the layer mid-radius (`null` spreads it over the whole layer);
3. grade it axially as CLEOv10.am: 166 / 309 / 166 turns per layer;
4. map the section edges (+-85.5 of +-173.8 cm in the deck) linearly onto the DD4hep winding, giving +-85.475 cm;
5. use the current per turn `current_A` = 3266 A.

The result reproduces the CLEOv10.am coil to within 0.05 cm (ends +-173.75 against +-173.80), with identical currents: 542 156 and 1 009 194 A-turns per region, 4 187 012 in total. NIM quotes 1300 turns x 3300 A; to use that, edit `cleo2_coils.json`.

## Correction coils

**v4c.** `solenoid_v4.xml` has a `CorrectionCoil` (Copper, r 40-49 cm, z -277..-268 cm) and a steel `CorrectorRing` (r 50.80-51.435 cm, z -278..-266.70 cm).
- The ring is iron.
- The coil is written as a `mat=1, cur=0.0` region: its current in the 2022 TOSCA v4 model is not available (`solenoid_v4.dat` returns 404).
- Set the current with `--aux-ni CorrectionCoil=<A-turns>`. Bz at the target is linear in it (Stage 1/2 html), so two solves give the cancelling value.

**Stage 2.** `solenoid_lite_sidis.xml` contains **no** correction coil, and no plug, nose or bridge either. The CLEOv11 coil of "CLEO-II Stage 1 and 2.html" (2026-10-31) is in `cleo2_coils.json` as `cleov11_correction_coil`:
- z -286.5..-266.5 cm, r 40.0-47.5 cm;
- NI = -16 240 A-turns, J = -108.27 A/cm2, air cooled;
- it sits on an aluminium ring, r 38.94-39.57 cm, which is non-magnetic and therefore not modelled;
- `cleov12_correction_coil` (-16 166 A-turns) is also in the file.

In the DD4hep stage-2 iron this coil would overlap UpstreamShield4 and 5, and 0.2 cm of the snapped Shield3, by 54 cm2. The deck therefore carries it only as a commented-out region, and the map from `solid_lite_stage2.am` has no correction coil. `--addon-coil cleov11_correction_coil[:J=<A/cm2>|:NI=<A-turns>]` adds it; while it overlaps iron it is refused unless `--allow-overlap` is given.

## Validation (`../validation/geometry_*.png`; full tables in the compare reports)

### (a) Stage 1: `solenoid_nomod.xml` vs `CLEOv10.am`, agreement within 0.25 cm except the spacer length

Iron: 458.8 t against 457.8 t. Of that, 1.5 t is only in the DD4hep model and 0.5 t only in CLEOv10. The Hausdorff distance of the iron outlines is 0.87 cm.

| CLEOv10 region | DD4hep piece | outer-edge offsets, DD4hep minus deck (cm) |
|---|---|---|
| Barrel yoke Inner | InnerOctagon | ends z +-266.70 vs +-266.50: +0.20; bore 176.53 vs 176.6: +0.07; outside 212.62 vs 212.6: +0.02 |
| Barrel yoke Outer | OuterOctagon | ends +0.20; inside 221.51 vs 221.5: -0.01; outside 257.61 vs 257.5: +0.11 |
| Slab spacer Up/Down | SpacerBarsUp/Down (ring) | outer end +0.20; inner end -236.77 vs -235.90: **-0.87** (coverage 0.879 vs the deck's 0.90) |
| Coil collar Up/Down | Up/DownstreamCoilCollar | bore 144.14 vs 144.0: -0.14; inner face z +-189.23 vs +-189.0: -0.23; outer end +0.20 |
| 6 coil regions | main coil | <= 0.05 cm, A-turns identical |

The collar/yoke boundary itself differs: the DD4hep collars end at the 180.34 cm counterbore, CLEOv10's at 177.9 cm. This boundary lies inside the merged iron and changes no material.

### (b) Full SoLID v4c: `solenoid_v4.xml` vs `CLEOv9.am`, expected real differences

Iron: 620.3 t (DD4hep 642.1 t) against 603.6 t; 577.6 t is common. The barrel and the upstream collar agree as in (a). The coil is identical.

The real differences: CLEOv9 is the 2019 v0-like long endcap, whereas v4c is the 2022 Benesch v4 with the octagon flats corrected.

| where | v4c (DD4hep) | CLEOv9 |
|---|---|---|
| upstream bore | UpstreamPlug z -266.70..-214.25, 25 deg bore 38.94 -> 63.32 cm, out to 144.14 (+8.3 t) | front piece z -237..-207 (r 55.6/70.2-144) and shields 1-3 inside the bore |
| upstream of the yoke | steel CorrectorRing r 50.80-51.44 and the Cu CorrectionCoil (current unknown) | shields 4 and 5 at z -274.5..-270.5 and -280.5..-277.5 to r 257.5 (-6.4 and -4.8 t) |
| barrel downstream end | inner octagon to z 224.79 (bevel 212.62 -> 200.17), outer to 209.23, spacer ring z 179.30..209.23 | barrel to 189.0, then 0.2 cm of air before the collar at 189.2; spacer z 159..189 (-1.9 t) |
| downstream collar | 144 -> 164.45 cm cone bore (z 199.23-224.79) + EndcapInterface r 257.61-285.24, z 170.67-209.23 (+7.0 t) | "Coil collar Downstream 1", one piece r 144-285, z 189.2-224.8 with a sloped face (-2.4 t) |
| endcap cylinder | r 270-285.24, z 209.23-529.59 | donut r 270-285, z 209-530 (0.24 cm thinner) |
| back plates | two 15.24 cm plates r 85-285.24 (z 529.59-560.07) + cone plate z 560.07-575.31 (+18.6 and +2.8 t) | bottom inner r 30-285 (z 530-545) + bottom outer r 30/45-185/170 (z 545-560) |
| nose | 6.5 deg cone from z -330: bore 19 -> 28.31 -> 35, outside 59.16 -> 85, z 189.23-575.31 | bore 20 -> 30, outside 61 -> 69 -> 90, z 189.2-530 (+2.7 / -7.3 t) |

### (c) Stage 2: `solenoid_lite_sidis.xml` vs CLEOv11 of the Stage 1/2 html, different designs

`solenoid_lite_sidis.xml` digitises the **2026-09-30 workshop** CLEOV11.AM (C068 p11, the GEMC v0 pieces): front piece, 5 shields, donut and 2 bottom plates, no nose. Its own header says the later CLEOv11 is not built. The html's CLEOv11 (2026-10-31, stage 2A) is a different, later design.

Only the stage-1 iron agrees, with the offsets of (a). Iron: 583.3 t (DD4hep model) against 724.7 t, of which 556.0 t is common.

| where | DD4hep stage 2 | html CLEOv11 |
|---|---|---|
| upstream, inside the bore | FrontPiece z -237..-207 (r 55.6/70-144.14, 12.4 t) + Shields 1-3 | plug z -266.5..-240.375, 25 deg bore 38.94 -> 51.12, out to 144 (6.5 t not covered by the shields) |
| upstream of the yoke | Shield4/5 to r 257.5 (6.4 + 4.8 t) | Al ring + correction coil -16 240 A-turns (the coil is not in DD4hep) |
| downstream nose | none | nose z 189.2-225, r 20-60 (2.8 t) + bridge r 50-60, z 225-530 (8.3 t) |
| donut | z 266.70-530 | z 260-540 |
| back plates | r 30-285 (z 530-545) + r 30/45-185/170 (z 545-560) | three 36 cm plates r 40-267 at z 530-566, 570-606, 610-646 (25 + 62 + 62 t only in the html) |

The html data block itself lacks the "Upstream coil Inner" region, so its coil total reads 3 628 616 A-turns. The geometry is symmetric, and the other five regions match within 0.05 cm.

The 2026-09-30 design was not force balanced: -59.7 t on the coil according to the workshop. A stage-2 map from `solid_lite_stage2.am` is therefore the DD4hep provisional design, not CLEOv11 2A. To map the html design, `solenoid_lite_sidis.xml` must first be updated to it; then rerun `extract` with `--addon-coil cleov11_correction_coil`.

## Known limitations

- 2D model of 3D iron. The octagons and bars follow the conventions above. The octagon corners (up to +8.2 % in radius) appear only in the equal-area mass bookkeeping, and the azimuthal variation of the field is not modelled.
- The v4c correction-coil current is unknown (0 A-turns). The stage-2 correction coil is not in the DD4hep geometry.
- The solver resolves the 0.4 cm coil sheets and the 0.635 cm corrector ring (thin layers get 0.7 times their thickness).
- Iron is single-material `mat=2` with POISSON's built-in curve, the same for all steels.
