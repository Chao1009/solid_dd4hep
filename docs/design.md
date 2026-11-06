# SoLID DD4hep simulation: design

Status: v0.3 (November 2026). This document describes how the SoLID detector description in this repository is
organized, the conventions every subsystem follows, and how it is validated against `solid_gemc` and the reference
documents (2019 preCDR, 2022 white paper, Oct 2026 High-Luminosity workshop).

## 1. Goals

* One DD4hep description that serves **all SoLID configurations**: full SoLID (SIDIS 3He, SIDIS NH3, PVDIS, J/psi,
  DDVCS) and **SoLID-Lite** (J/psi-lite and DDVCS-lite now, SIDIS-3He-lite once its magnet stage is defined).
* Geometry that matches the current GEMC baseline (`*_moved` = longer endcap, `*_separated`/`*_ext` = lite) unless a
  newer or more authoritative source exists (documents, real hardware such as the CLEO-II magnet). Every deviation
  from GEMC is listed in `docs/geometry_comparison.md`.
* Simulation output in EDM4hep (ddsim/npsim), readouts ready for EICrecon (or Gaudi) reconstruction and ACTS tracking.
* Swappable components for SoLID-Lite, whose layout is still a set of options (EC type, forward-angle segmentation,
  magnet stage, Cherenkov sector list, MRPC on/off, muon system on/off).

## 2. Coordinate frame and units

The **GEMC hall frame** is kept so that positions, field maps and plots compare 1:1 with `solid_gemc`:
right-handed, z along the beam, origin at the center of the CLEO-II coil.

| configuration | target center z | target |
|---|---|---|
| SIDIS 3He | -350 cm | 40 cm polarized 3He glass cell |
| SIDIS NH3 | -350 cm | 3 cm polarized NH3 (+ target magnet) |
| J/psi, DDVCS (full and lite) | -315 cm | 15 cm LH2 |
| PVDIS | +10 cm | 40 cm LD2 (LH2 option) |

Compact files always carry explicit units. GEMC rotations are passive; all rotations in this repository are
**active** (DD4hep convention): a GEMC rotation "a b c" becomes `RotationX(-a)*RotationY(-b)*RotationZ(-c)`.

## 3. Repository layout

```
solid_<CONFIG>.xml            top-level configurations (solid.xml -> solid_SIDIS_He3.xml)
compact/
  definitions.xml             world, system-ID table, global constants
  materials.xml optics.xml    shared materials and optical tables
  colors.xml display.xml      named colors and common vis attributes
  configs/                    detector sets shared by several configurations (full_common, lite_common)
  fields/                     field maps (CLEOv9: full SoLID, CLEOv10: lite stage 1)
  magnet/ target/ beamline/ baffle/ tracking/ cherenkov/ ecal/ timing/ muon/   subsystem fragments
src/                          DD4hep plugins (SoLID_* detector constructors, field map, file loader)
steering/solid_sim.py         generic ddsim/npsim steering (optical photons, gun per configuration)
scripts/geometry/             check_geometry.py (load + overlaps), solid_geo_dump.py (extents), gemc_txt.py,
                              compare_gemc.py (dd4hep vs GEMC extents)
bin/                          thissolid.sh (generated), solid_lineview, make_solid_views
docs/                         design.md (this file), geometry_comparison.md, visualization.md, reconstruction.md
containers/                   solid-shell launcher and the solid_xl image recipe
```

## 4. Configurations

| top-level | detectors | field | GEMC reference gcard |
|---|---|---|---|
| `solid_SIDIS_He3.xml` (default) | full_common + He3 target + He3 beamline | CLEOv9 | solid_SIDIS_He3_moved_full |
| `solid_SIDIS_NH3.xml` | full_common + NH3 target + NH3 beamline | CLEOv9 | solid_SIDIS_NH3_moved_full |
| `solid_JPsi.xml` | full_common + LH2 15 cm + J/psi beamline | CLEOv9 | solid_JPsi_LH2_moved_full |
| `solid_DDVCS.xml` | J/psi + muon system | CLEOv9 | solid_JPsi_DDVCS_LH2_moved_full |
| `solid_PVDIS.xml` | v4 magnet, PVDIS target, beamline, baffles, PVDIS GEM, PVDIS LGC, PVDIS FAEC | CLEOv9 | solid_PVDIS_LD2_moved_full |
| `solid_JPsi_lite.xml` | lite_common + LH2 15 cm + J/psi beamline | CLEOv10 | solid_JPsi_LH2_lite_full |
| `solid_DDVCS_lite.xml` | J/psi-lite + muon system | CLEOv10 | solid_JPsi_DDVCS_LH2_lite_full |

* `full_common` = magnet v4 (longer endcap) + cryostat/coil, 6-plane GEM, LGC (SIDIS), HGC, LAEC, FAEC (95-240 cm),
  LA/FA SPD, FA MRPC.
* `lite_common` = magnet "nomod" (CLEO-II iron as built) + cryostat/coil, LA GEM (4 planes) + FA GEM (5 planes), one
  CO2 Cherenkov (HGC-based tank, sectors 4,5,11,12,19,20,26,27 = "4seg24degFA"), LAEC, FAEC (95-220 cm), SPD, FA MRPC.
* A provisional `solid_SIDIS_He3_lite.xml` (stage-2 magnet with front piece + endcap without nose) is planned; the
  iron outline exists only as a 2D POISSON test (Oct 2026) and there is no 3D field map yet.

Components are swapped by including a different fragment (e.g. `ecal/ec_lite_FA.xml` vs `ecal/ec_SIDIS_FA.xml`).
A variant that differs only in parameters is either a copy of the fragment with different constants, or is
selected by an **optional switch constant** that the fragment itself never defines and the top-level file may set,
e.g. `Target_NH3_Variant` (SM2023 / Oxford2020), `Target_He3_Collimators`, `Target_PVDIS_Fill` (LD2 / LH2),
`Baffle_Enable`, `Baffle_Mode` (exact / gemc staircase), `Baffle_Chirality`, `Baffle_Observers_Enable`,
`FAEC_PhiMask_enable` (+ `_nseg`, `_dphi`, `_centers`), `LiteCherenkov_Sectors`. The plugins fall back to the
documented default when the switch is absent. Validation-only copies reproduce GEMC exactly
(`magnet/solenoid_v4_gemc.xml`, `solenoid_nomod_gemc.xml`, `cleo_coil_cryostat_gemc.xml`). Do not redefine a fragment's constant in the top-level file to "override" it:
DD4hep 1.38 keeps the first definition in the Detector but the last one in the expression evaluator, so the result
depends on how the constant is used.

## 5. Conventions for fragments and plugins

### 5.1 System IDs and names
System IDs are constants in `compact/definitions.xml` (blocks: 1-9 passive, 10-19 tracking, 20-29 Cherenkov,
30-39 calorimetry, 40-49 timing, 50-59 muon, 90-99 virtual planes, 100+ beam tests).

| detector (DetElement name) | id constant | readout | plugin |
|---|---|---|---|
| Beamline | Beamline_ID | – | SoLID_Beamline |
| Target | Target_ID | – | SoLID_PolHe3Target / SoLID_PolNH3Target / SoLID_CryoTarget |
| MagnetYoke, MagnetCoil | MagnetYoke_ID, MagnetCoil_ID | – | SoLID_CLEOMagnet (alias SoLID_Solenoid) |
| PVDISBaffle | Baffle_ID | (BaffleObserverHits, optional) | SoLID_PVDISBaffle |
| GEMTracker / GEMTrackerLA / GEMTrackerFA / MuonGEMTracker | GEM_ID / GEMLA_ID / GEMFA_ID / MuonGEM_ID | `<name>Hits` | SoLID_GEMTracker |
| LightGasCherenkov | LGC_ID | LightGasCherenkovHits | SoLID_LightGasCherenkov |
| HeavyGasCherenkov, LiteCherenkov | HGC_ID, LiteCherenkov_ID | `<name>Hits` | SoLID_SectorGasCherenkov |
| FAECPreShower, FAECShower, LAECPreShower, LAECShower | *_ID | `<name>Hits` | SoLID_HexShashlykEMCal |
| FASPD, LASPD | SPDFA_ID, SPDLA_ID | FASPDHits, LASPDHits | SoLID_ScintPadDisk |
| FAMRPC | MRPCFA_ID | FAMRPCHits | SoLID_MRPCDisk |
| MuonIron, MuonScint | MuonIron_ID, MuonScint_ID | MuonScintHits | SoLID_MuonDetector |

Rules:
1. Every constant, vis attribute, optical surface and material defined in a fragment carries the subsystem prefix
   (DD4hep throws on duplicate vis names and **silently ignores** duplicate constants and materials).
2. Detectors reference their id by constant (`id="LGC_ID"`), never by number.
3. Readouts: `<DetElementName>Hits`, first field `system:8`, then the subsystem fields (sector, layer, module, ...).
   Field names follow GEMC identifiers where they exist so hits can be compared.
4. Sensitive types: `tracker` (GEM, SPD, MRPC, scintillators, Cherenkov photocathodes with the optical-photon filter),
   `calorimeter` (EC). The Cherenkov readouts use `Geant4OpticalTrackerAction` (see `steering/solid_sim.py`).
5. Plugins are named `SoLID_<Name>`, live in `src/<Name>_geo.cpp`, read all numbers from XML (no hard-coded
   geometry), build DetElements per sector/module so the VolumeManager can resolve every cell, and use
   `printout()` for diagnostics.
6. A fragment documents the constants it requires from the top level (normally only `TargetCenter_z`).
7. Inside a fragment, a nested `<include>` is processed after the fragment's own `<detectors>`. Materials, optical
   tables and surfaces a detector needs must therefore sit directly in the fragment, or come through
   `<includes><file ref="..."/></includes>`, which is processed first.
8. Readouts without segmentation need an explicit `<segmentation type="NoSegmentation"/>`. Otherwise DD4hep 1.38's
   `CellIDPositionConverter` throws.

### 5.2 Materials and optics
`compact/materials.xml` holds generic materials. Subsystem-specific materials (C4F8 at 1.7 atm, GE180, NH3+LHe mix,
borated polyethylene, CFRP mirror substrate, ...) are defined in the subsystem fragment with a prefixed name. Optical
property tables follow GEMC tables, quoted with their source.

### 5.3 Field
`SoLID_FieldMapB` reads the GEMC ASCII maps including their `<mfield>` header (grid, units, symmetry), converts the
header unit (gauss) to tesla itself (`scale` is a plain extra multiplier, normally 1), is thread-safe, and can
translate/rotate a map (cylindrical-z solenoid maps, cylindrical-x maps such as the NH3 target magnet). Maps:
`solenoid_CLEOv9.dat` (full SoLID; computed for an older long-endcap iron, the 3D map for the v4 yoke is not
available), `solenoid_CLEOv10.dat` (lite stage 1, unmodified CLEO-II iron), `oxford_ptarget.dat` (NH3 target, 5 T).
The `FileLoader` downloads them once from pinned `solid_gemc` GitHub raw URLs into `$DETECTOR_PATH/fieldmaps`
(offline: copy `solid_gemc/field/*.dat` there).

## 6. Validation

Every configuration must pass, in the eic-shell container:
1. `scripts/geometry/check_geometry.py solid_<CONFIG>.xml` (loads, TGeo overlap check at 1 um);
2. `checkOverlaps` (Geant4) before productions;
3. `scripts/geometry/compare_gemc.py` — extents of every subsystem vs the GEMC gcard outline
   (`solid_gemc/script/geometry_outline.py`); differences are either fixed or explained in
   `docs/geometry_comparison.md`;
4. a short ddsim run with the default gun (hits in every sensitive detector);
5. line-view drawings (`bin/make_solid_views`) reviewed by eye.

## 7. Visualization

`bin/solid_lineview` renders Geant4 DAWNFILE `.prim` files (from `ddsim`/`npsim` with `macro/dawnfile.mac`) into
clean vector line drawings (SVG/PDF/PNG): exact 2D sections, cutaway views with hidden-line removal, optional
trajectories. It replaces the dawn/dawncut/ps2pdf chain. See `docs/visualization.md`.

## 8. Reconstruction and environment

See `docs/reconstruction.md` for the framework recommendation (EICrecon out-of-tree plugins over a
framework-neutral algorithm library, with optional Gaudi/k4FWCore wrappers) and `containers/` for the eic-shell based
`solid-shell` environment.
