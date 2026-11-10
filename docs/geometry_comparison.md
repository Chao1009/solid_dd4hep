# SoLID geometry: DD4hep vs GEMC vs reference documents

This page compares the geometry in this repository with `solid_gemc` and with the reference documents:
- DD4hep: this repository, v0.3, Nov 2026;
- GEMC: `JeffersonLab/solid_gemc` at 6826f84 (2026-11-02);
- the updated preliminary CDR (preCDR, Nov 2019);
- the SoLID white paper (WP, arXiv:2209.13357);
- the JLab High-Luminosity workshop (Oct 2026, Indico 1089), whose SoLID-Lite talks are cited as C0xx;
- for the existing CLEO-II magnet, the hardware publication (NIM A320 (1992) 66) and the CAD models.

The source-priority rule:
- GEMC "moved" (longer endcap) is the full-SoLID baseline, and GEMC "lite" (Oct 2026) is the SoLID-Lite baseline.
  Both are newer than the preCDR. The preCDR tables mostly describe the older short-endcap layout, so preCDR/GEMC
  differences in z are expected.
- Where GEMC is wrong or outdated, the documents or the real hardware win. Those cases are flagged as deviations in
  §3.

All lengths in cm unless noted, hall frame (z along the beam, origin at the coil center).

## 1. Envelope comparison per configuration (generated)

`scripts/geometry/compare_gemc.py` builds every DD4hep configuration and compares the z and r extents of each
top-level detector with the volumes the matching GEMC gcard loads (virtual/flux planes excluded). It runs in CI for
every configuration. "DIFF" rows are explained in §3. A polyhedral r_min can also differ for convention reasons: the
DD4hep side uses mesh corners, the GEMC side the inner flat.

#### solid_SIDIS_He3  (GEMC: `solid_SIDIS_He3_moved_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -575.00 … 550.00 | -575.00 … 550.00 | 0.00 … 19.00 | 0.00 … 19.00 | 0.00 | 0.00 | OK |
| Target | -420.00 … -280.00 | -420.00 … -280.00 | 0.00 … 79.10 | 0.00 … 79.10 | 0.00 | 0.00 | OK |
| MagnetYoke | -278.00 … 575.31 | -278.00 … 575.31 | 19.00 … 285.24 | 19.00 … 308.74 | 0.00 | 23.50 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| GEMTrackerAssembly | -175.80 … 92.80 | -175.80 … 92.80 | 20.89 … 135.74 | 21.00 … 135.00 | 0.00 | 0.74 | DIFF |
| LightGasCherenkov | 97.00 … 312.00 | 97.00 … 312.00 | 50.00 … 263.00 | 50.00 … 263.00 | 0.00 | 0.00 | OK |
| HeavyGasCherenkov | 316.00 … 426.00 | 316.00 … 426.00 | 77.00 … 240.00 | 77.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |

*tolerance 0.1 cm*

#### solid_SIDIS_NH3  (GEMC: `solid_SIDIS_NH3_moved_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -300.00 … 550.00 | -300.00 … 550.00 | 0.00 … 19.00 | 0.00 … 19.00 | 0.00 | 0.00 | OK |
| Target | -397.94 … -302.06 | -415.00 … -285.00 | 43.18 … 64.52 | 91.92 … 91.92 | 17.06 | 48.74 | DIFF |
| MagnetYoke | -278.00 … 575.31 | -278.00 … 575.31 | 19.00 … 285.24 | 19.00 … 308.74 | 0.00 | 23.50 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| GEMTrackerAssembly | -175.80 … 92.80 | -175.80 … 92.80 | 20.89 … 135.74 | 21.00 … 135.00 | 0.00 | 0.74 | DIFF |
| LightGasCherenkov | 97.00 … 312.00 | 97.00 … 312.00 | 50.00 … 263.00 | 50.00 … 263.00 | 0.00 | 0.00 | OK |
| HeavyGasCherenkov | 316.00 … 426.00 | 316.00 … 426.00 | 77.00 … 240.00 | 77.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |

*tolerance 0.1 cm*

#### solid_PVDIS  (GEMC: `solid_PVDIS_LD2_moved_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -687.00 … 550.00 | -687.00 … 550.00 | 0.00 … 30.00 | 0.00 … 30.00 | 0.00 | 0.00 | OK |
| Target | -27.00 … 47.00 | -27.00 … 47.00 | 0.00 … 45.00 | 0.00 … 45.00 | 0.00 | 0.00 | OK |
| MagnetYoke | -278.00 … 575.31 | -278.00 … 575.31 | 19.00 … 285.24 | 19.00 … 308.74 | 0.00 | 23.50 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| PVDISBaffle | 35.50 … 330.90 | 49.50 … 332.50 | 3.89 … 200.00 | 5.56 … 200.00 | 14.00 | 1.67 | DIFF |
| GEMTrackerAssembly | 156.70 … 325.80 | 156.70 … 325.80 | 50.72 … 245.34 | 51.00 … 244.00 | 0.00 | 1.34 | DIFF |
| LightGasCherenkov | 194.00 … 312.00 | 194.00 … 312.00 | 63.20 … 263.00 | 63.20 … 263.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 332.18 … 337.30 | 332.18 … 337.30 | 94.21 … 269.90 | 94.22 … 269.89 | 0.00 | 0.01 | OK |
| FAECShower | 337.30 … 382.70 | 337.30 … 382.70 | 94.21 … 269.90 | 94.22 … 269.89 | 0.00 | 0.01 | OK |

*tolerance 0.1 cm*

#### solid_JPsi  (GEMC: `solid_JPsi_LH2_moved_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -627.50 … 802.50 | -627.50 … 802.50 | 0.00 … 29.00 | 0.00 … 29.00 | 0.00 | 0.00 | OK |
| Target | -322.51 … -307.49 | -322.51 … -307.49 | 0.00 … 2.00 | 0.00 … 2.00 | 0.00 | 0.00 | OK |
| MagnetYoke | -278.00 … 575.31 | -278.00 … 575.31 | 19.00 … 285.24 | 19.00 … 308.74 | 0.00 | 23.50 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| GEMTrackerAssembly | -175.80 … 92.80 | -175.80 … 92.80 | 20.89 … 135.74 | 21.00 … 135.00 | 0.00 | 0.74 | DIFF |
| LightGasCherenkov | 97.00 … 312.00 | 97.00 … 312.00 | 50.00 … 263.00 | 50.00 … 263.00 | 0.00 | 0.00 | OK |
| HeavyGasCherenkov | 316.00 … 426.00 | 316.00 … 426.00 | 77.00 … 240.00 | 77.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |

*tolerance 0.1 cm*

#### solid_DDVCS  (GEMC: `solid_JPsi_DDVCS_LH2_moved_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -627.50 … 802.50 | -627.50 … 802.50 | 0.00 … 29.00 | 0.00 … 29.00 | 0.00 | 0.00 | OK |
| Target | -322.51 … -307.49 | -322.51 … -307.49 | 0.00 … 2.00 | 0.00 … 2.00 | 0.00 | 0.00 | OK |
| MagnetYoke | -278.00 … 575.31 | -278.00 … 575.31 | 19.00 … 285.24 | 19.00 … 308.74 | 0.00 | 23.50 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| MuonIron | 620.00 … 748.00 | 620.00 … 748.00 | 97.24 … 393.75 | 97.24 … 393.75 | 0.00 | 0.00 | OK |
| MuonGEMTracker | 657.70 … 751.30 | 657.70 … 751.30 | 99.45 … 286.57 | 100.00 … 285.00 | 0.00 | 1.57 | DIFF |
| GEMTrackerAssembly | -175.80 … 92.80 | -175.80 … 92.80 | 20.89 … 135.74 | 21.00 … 135.00 | 0.00 | 0.74 | DIFF |
| LightGasCherenkov | 97.00 … 312.00 | 97.00 … 312.00 | 50.00 … 263.00 | 50.00 … 263.00 | 0.00 | 0.00 | OK |
| HeavyGasCherenkov | 316.00 … 426.00 | 316.00 … 426.00 | 77.00 … 240.00 | 77.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 246.13 | 89.26 … 246.12 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |
| MuonScint | 660.50 … 757.50 | 660.50 … 757.50 | 100.00 … 285.00 | 100.00 … 285.00 | 0.00 | 0.00 | OK |
| MuonObserver | -266.70 … 758.01 | – | 100.00 … 286.10 | – | nan | nan | no GEMC match |

*tolerance 0.1 cm*

#### solid_JPsi_lite  (GEMC: `solid_JPsi_LH2_lite_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -627.50 … 802.50 | -627.50 … 802.50 | 0.00 … 29.00 | 0.00 … 29.00 | 0.00 | 0.00 | OK |
| Target | -322.51 … -307.49 | -322.51 … -307.49 | 0.00 … 2.00 | 0.00 … 2.00 | 0.00 | 0.00 | OK |
| MagnetYoke | -266.70 … 266.70 | -266.70 … 266.70 | 144.14 … 278.84 | 144.00 … 301.81 | 0.00 | 22.97 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| GEMTrackerLiteAssembly | -175.80 … 180.80 | -175.80 … 180.80 | 35.80 … 133.73 | 36.00 … 133.00 | 0.00 | 0.73 | DIFF |
| LiteCherenkov | 205.00 … 426.00 | 205.00 … 426.00 | 60.00 … 240.00 | 60.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 225.77 | 89.26 … 225.76 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 225.77 | 89.26 … 225.76 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |

*tolerance 0.1 cm*

#### solid_DDVCS_lite  (GEMC: `solid_JPsi_DDVCS_LH2_lite_full.gcard`)

| DD4hep detector | DD4hep z [cm] | GEMC z [cm] | DD4hep r [cm] | GEMC r [cm] | max Δz | max Δr | status |
|---|---|---|---|---|---:|---:|---|
| Beamline | -627.50 … 802.50 | -627.50 … 802.50 | 0.00 … 29.00 | 0.00 … 29.00 | 0.00 | 0.00 | OK |
| Target | -322.51 … -307.49 | -322.51 … -307.49 | 0.00 … 2.00 | 0.00 … 2.00 | 0.00 | 0.00 | OK |
| MagnetYoke | -266.70 … 266.70 | -266.70 … 266.70 | 144.14 … 278.84 | 144.00 … 301.81 | 0.00 | 22.97 | DIFF |
| MagnetCoil | -189.00 … 189.00 | -189.00 … 189.00 | 144.00 … 174.00 | 144.00 … 176.59 | 0.00 | 2.59 | DIFF |
| MuonIron | 620.00 … 748.00 | 620.00 … 748.00 | 97.24 … 393.75 | 97.24 … 393.75 | 0.00 | 0.00 | OK |
| MuonGEMTracker | 657.70 … 751.30 | 657.70 … 751.30 | 99.45 … 286.57 | 100.00 … 285.00 | 0.00 | 1.57 | DIFF |
| GEMTrackerLiteAssembly | -175.80 … 180.80 | -175.80 … 180.80 | 35.80 … 133.73 | 36.00 … 133.00 | 0.00 | 0.73 | DIFF |
| LiteCherenkov | 205.00 … 426.00 | 205.00 … 426.00 | 60.00 … 240.00 | 60.00 … 240.00 | 0.00 | 0.00 | OK |
| FAECPreShower | 439.18 … 444.30 | 439.18 … 444.30 | 89.25 … 225.77 | 89.26 … 225.76 | 0.00 | 0.01 | OK |
| FAECShower | 444.30 … 489.70 | 444.30 … 489.70 | 89.25 … 225.77 | 89.26 … 225.76 | 0.00 | 0.01 | OK |
| LAECPreShower | -62.82 … -57.70 | -62.82 … -57.70 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| LAECShower | -57.70 … -12.30 | -57.70 … -12.30 | 81.18 … 142.40 | 81.19 … 142.39 | 0.00 | 0.01 | OK |
| FASPD | 427.70 … 428.30 | 427.75 … 428.25 | 93.00 … 216.00 | 93.00 … 216.00 | 0.05 | 0.00 | OK |
| LASPD | -66.00 … -64.00 | -66.00 … -64.00 | 80.00 … 135.00 | 80.00 … 135.00 | 0.00 | 0.00 | OK |
| FAMRPC | 433.58 … 436.42 | 433.58 … 436.42 | 93.00 … 216.00 | 93.00 … 216.00 | 0.00 | 0.00 | OK |
| MuonScint | 660.50 … 757.50 | 660.50 … 757.50 | 100.00 … 285.00 | 100.00 … 285.00 | 0.00 | 0.00 | OK |
| MuonObserver | -266.70 … 758.01 | – | 100.00 … 286.10 | – | nan | nan | no GEMC match |

*tolerance 0.1 cm*

**What the DIFF rows mean.** All of them are deliberate; details are in §2 and §3.

| detector | difference | reason |
|---|---|---|
| MagnetYoke | r_max 285.2 vs 308.7 (full), 278.8 vs 301.8 (lite) | real CLEO-II octagon: GEMC entered the corner distances as flat distances, making the iron ~8 % too large (NIM A320; CAD). GEMC-parity copies: `magnet/solenoid_*_gemc.xml` |
| MagnetCoil | cryostat r_max 174.0 vs 176.6 | cryostat walls 10/12 mm and outer radius 1.74 m from the CLEO-II paper (GEMC: 1 mm walls, radius copied from the 2D field model); Al coil support shell added |
| GEMTracker(Lite)Assembly, MuonGEMTracker | r differs by <= 1.6 cm | planar trapezoid modules covering GEMC's annular sectors, needed for ACTS/EICrecon. GEMC-shaped option: `tracking/gem_SIDIS_annular.xml` |
| Target (NH3) | not a like-for-like envelope | default is the real SM2023 target magnet; the GEMC gcard still loads Oxford2020 (selectable with `Target_NH3_Variant`). Its chamber axis is vertical, which the envelope method only approximates; placements were checked against GEMC to 1e-9 cm in the subsystem test |
| PVDISBaffle | z starts at 35.5 vs 49.5, ends at 330.9 vs 332.5 | plate 1 belongs to the baffle here; GEMC builds it inside the target enclosure. The photon blocker moved 1.6 cm upstream because the GEMC position overlaps the FAEC |
| MuonObserver | no GEMC match | virtual planes (GEMC flux planes are excluded from the comparison) |

Not visible as DIFF (below the 0.1 cm tolerance) but deliberate: the FASPD is 0.6 cm thick (2026 prototype) instead of
GEMC's 0.5 cm (GEMC-thickness test: `test/test_timing_muon_gemc.xml`).

## 2. Key parameters per subsystem

Columns: value implemented here | GEMC latest variant | reference document (page) | note. Each implementer verified
the numbers against the GEMC geometry files with Python, placements to between 1e-9 cm and 0.1 mm depending on the
subsystem. The validation details are in the subsystem test files under `test/`.

### 2.1 CLEO-II solenoid and field maps

Plugin `SoLID_CLEOMagnet`; fragments `compact/magnet/*`, `compact/fields/*`. Reference "NIM" = CLEO-II detector paper, NIM A320 (1992) 66.

| Parameter (cm) | dd4hep | GEMC latest | reference | note |
|---|---|---|---|---|
| inner octagon flats | 176.53–212.62 | 191.075–230.142 (v4/nomod) | NIM A320 p70: 36 cm layers, 9 cm gaps | GEMC entered corners as apothems |
| outer octagon flats | 221.51–257.61 | 239.763–278.832 | NIM p70 | |
| octagon orientation | flats ⊥ x, y | same | NIM p70 (22.5° joints) | |
| barrel z | −266.70…209.23 (full), ±266.70 (lite) | same | 210" slabs | |
| spacers | 8 bars/end, 149.86×8.89×34.04 | 25.5 cm ring, down from 182.73 | preCDR p66 | ring = option |
| upstream collar | round 144.14 / oct 180.34 | Pgon 144.145–191.075 | NIM p70–71 | |
| downstream collar bore | 144→164.45 | 144→164.447 | preCDR p66 (PVDIS 35°) | |
| plug bore | 38.94→63.32 | 38.94→66.19 | 25° from −350 (CAD v4) | |
| endcap interface | z 170.67, round 285.24 | z 183.83, corners 308.7 | CAD v4 | |
| endcap cylinder, plates, nose | as GEMC v4 | v4 | preCDR Fig 46/50; WP Fig 29 | |
| coil winding | 151.7–154.9, ±173.75 | same, G4_Al | preCDR p65; NIM Table 1 | conductor 2.92 g/cm³ |
| coil shell (Al) | 155.0–159.5 | none | NIM: 7000 kg cold mass | estimate |
| cryostat walls | 144.0–145.0 / 172.8–174.0 | 144.0–144.1 / 176.48–176.59 | NIM §2.8: 10/12 mm, r 1.74 m | |
| flanges | 2 cm, ±187–189 | 1.99 cm | NIM 20 mm; preCDR 3.8 m | |
| yoke mass v4 / lite (t) | 642.1 / 488.0 | 715.8 / 585.1 | NIM 800 t (3 layers + poles) | |
| coil + cryostat (t) | 14.9 | 4.6 | NIM 7 t + 20 t | |
| bending integral ∫B⊥ dl along a straight 8° line from z −350 to 200 (T·m), v9 / v10 | −0.9629 / −0.5651 | same | workshop C068 p12 'signed integral': −0.9629 / −0.5634 | B⊥ = component normal to the line in the r-z plane; the talk used an older v10 map |

### 2.2 Targets and beamlines

Plugins `SoLID_PolHe3Target`, `SoLID_PolNH3Target`, `SoLID_CryoTarget`, `SoLID_Beamline`.

| Parameter | dd4hep | GEMC latest | reference | note |
|---|---|---|---|---|
| He3 center / length | −350 / 40 cm | same | preCDR p53, p218 | |
| He3 window / wall | GE180 120 µm / 1 mm | same | – | |
| He3 density | 10 amg (3He isotope) | same | preCDR p71; WP p7 | |
| He3 collimators | W powder, 7–23° / 7–15° | collimator 2022 | preCDR p71, p202 | 2018 option |
| Helmholtz pairs | 2 (R 75.8, 66.7 cm) | 2 | preCDR p71: 3 | 3rd pair disabled |
| NH3 cup | r 1.36144 × 2.82702 cm | SM2023 | preCDR p53: 3 cm | |
| NH3 chamber | r 47.94, exit ±28°/±22 cm | SM perl | preCDR p74: ±25° | |
| NH3 field | Oxford 5 T overlay | map in box | preCDR p72 | |
| J/psi cell | LH2 15 cm, Ø3.774 cm, Al 178/102 µm | same | preCDR p61; WP p20 | z −315 |
| PVDIS cell | LD2 40 cm × r 2.5 cm, Al 0.4 mm / 120 µm | same | preCDR p74–75: 40 cm LH2/LD2, 120 µm windows (radius and side wall from GEMC only) | z +10 |
| PVDIS enclosure | Al r 45, z −27…47 | same | preCDR p75: qualitative (two-section IVC), dimensions from GEMC only | |
| Plate-1 cavity | r 3.89–44.3, z 35.498–44.502 | plate 1 inside TACV | – | slot for baffle plate 1 (§3) |
| He3 exit pipe | Al 1.05→19 cm, −325…550 | beamline_moved | – | |
| J/psi exit pipe | 1.8→29 cm, to 802.5 | same | – | |
| PVDIS pipes | r 30 at −687…−27; 5→19 at 47…550 | same | – | |
| Kill disk | vacuum + LimitSet | Kryptonite | – | |

### 2.3 PVDIS baffle

Plugin `SoLID_PVDISBaffle`.

| Parameter | dd4hep | GEMC baffle_moved (Enclosure) | reference | note |
|---|---|---|---|---|
| plates × slits | 11 × 30 | 11 × 30 | preCDR p50, T7 p52 | |
| plate thickness | 9 cm Pb | 9 cm | preCDR p50 | |
| plate z centres | 40…180, pitch 14 | same | none (preCDR Fig.54 shows the plate shapes only) | |
| plate 1 r | 3.89/37.0/44.3 | same | preCDR Fig.54 | |
| plate 11 r | 58.8/60.4/132/140 | same | WP Fig.30 | |
| plate 11 kink r1 | 102.328 | same | – | |
| sector-1 rotation | 99° (96+3) | rot −99° (passive) | – | active RotZ |
| lead edges | exact polygon | 20-step staircase | – | ≤0.07 mm vs ≤3.5 mm |
| min slit opening | 4.24° (pl. 3, r0) | same | preCDR ">4°" | |
| open fraction | 0.502–0.518 | same | – | |
| Pb volume | 1.0253e6 cm³ | 1.0253e6 | – | 17.3 t incl. rings |
| photon blocker z | 325.9–330.9 | 327.5–332.5 | none (preCDR T7 gives no z) | moved 1.6 cm: GEMC position overlaps the FAEC |
| photon blocker | 30 × 2.5°, r 110–200, 5 cm | same | preCDR T7 p52: 5 cm, r 110–200 (30 × 2.5° blocks from GEMC only) | |
| BorPoly | 70% PE + 30% B, 1.19 g/cm³ | same | none | |
| inner absorber rmin | ≥ pipe + 1 mm | 5.56 (44 µm clearance) | – | |
| outer absorber z end | 175.498 | 188.99 (zero thick) | – | |
| observers | 22 × 20 µm, off | 22 (flux) | – | |
| plate 1 mother | slot in TACH/TACV cavity | TACV | – | §3 |

### 2.4 GEM trackers

Plugin `SoLID_GEMTracker`.

| Parameter | dd4hep | GEMC latest | reference doc | note |
|---|---|---|---|---|
| SIDIS planes, z (cm) | 6: −175/−150/−119/−68/5/92 | same (V2 2025) | preCDR T11 p86 | stack centre |
| SIDIS R (cm) | 39-87, 21-98, 25-112, 32-135, 42-99, 53-122 | same | T11 "needed" 36-87…55-123 | trapezoid ±0.55 % |
| PVDIS z (cm) | 157.5/185.5/190/315/325 | same (V3) | T10 306/315 (short endcap) | longer endcap |
| PVDIS R (cm) | 51-115, 61-135, 63-138, 104-237, 107-244 | same | T10 50-118…110-228 | — |
| PVDIS offrot (deg) | +0.5/0/0/−0.5/−0.5 | same | — | — |
| lite LA z / R | −175…−68 / 43-82…70-132 | same (gem_separated, 2026-10-25) | none (C068 discusses the lite GEM only qualitatively) | 14-25° |
| lite FA z / R | 0…180 / 36-88…55-133 | same (gem_separated) | none (C068 p17: qualitative only) | 6-14° |
| muon GEM | z 658.5/704.5/750.5, R 100-285 | same (V8) | WP p22 | µRWELL planned |
| sectors | 30 × 12°, s at 96+12(s−1)° | same | preCDR p83 | — |
| module shape | covering trapezoid | annular sector | "trapezoidal modules" p83 | annular variant available |
| envelope r (SIDIS 1) | 38.786-87.479 | 39-87 | — | trapezoids covering the GEMC annulus |
| thickness | 16.005 mm, 25 layers | same | — | — |
| X/X0 per plane | 0.915 % | 0.915 % | — | legacy dd4hep 0.219 % |
| gaps | 3/2/2/2 mm | same | Fig. 57 p83 | — |
| gas | Ar/CO2 mass 0.68/0.32, 1.8407 mg/cm³ | same | — | beam test 75:25 |
| sensitive volume | drift gap (slice 7) | slices 6,7,8,10,11,12,20 | p152 | option: `sensitive="true"` |
| resolution emulation | σ_φ 49.8 µm, σ_r 474 µm | true hits | 70 µm per view, WP p26 | — |
| ids | (sys, layer, module, slice) | 10⁶+n·10⁵+s·10³+l | — | — |

### 2.5 Light gas Cherenkov (LGC)

Plugin `SoLID_LightGasCherenkov`.

| Parameter | dd4hep | GEMC lgc_moved | reference | note |
|---|---|---|---|---|
| sectors, sector 1 φ | 30, 102 | 30, 102 | preCDR p96 | |
| tank z SIDIS/PVDIS | 97-312 / 194-312 | same | p49 2 m / 1 m | |
| tank r (SIDIS) | 50-263 | same | p96 r_out 265 | |
| gas | CO2 20 °C, n-1=4.18e-4 | same | p96 | old dd4hep: 0 °C index |
| windows | PVF 0.05 mm front, 0.1 mm back, full faces | SIDIS: front only (58.1–126.9, extruding the tank); PVDIS: front + back | p96 | SIDIS exit window added |
| M1 R; centre | 362.535; (0,326.619,-24.579) | same | Tab.12: 277.51 | |
| M1 coverage (z=-350) | 7.54-14.89 | same | 8-15 (p97, WP p25) | |
| M2 R; centre | 155.987; (0,208.119,147.66) | same | Tab.12: 157.99 | |
| M2 coverage (z=+10) | 32.1-37.5 | 32.05-37.52 | 22-35 total | |
| PVDIS M1 centre after tilt | (0,288.030,-57.275) | (0,288.033,-57.279) | ≈8° (p97) | Rx(-8°) about (0,73.537,235) |
| mirror | CFRP 3 mm | SL_RCFP (33% H) | CF (p97) | |
| PMT frame | (0,245,250), Rx63 + Rz1.8 | same | p97 | |
| PMT array | 3×3 at 5.2; 8×8 pixels, 6.48 at 6.5 mm | same | H12700 (WP p26) | |
| Winston cone | r 7.8→21, L 30 | same | p97 | |
| shield | z -11.52→7.34 | →8.48, overlaps cone | 30 cm (p97) | |
| p.e. per electron | 30.8 / 16.1 | — | 25-30 / 20-40 (Figs 75-76) | |

### 2.6 Heavy gas Cherenkov (HGC) and SoLID-Lite CO2 Cherenkov

Plugin `SoLID_SectorGasCherenkov`.

| Parameter | dd4hep | GEMC latest (variant) | reference doc | note |
|---|---|---|---|---|
| tank z | 316/326–426 | same (hgc_moved) | preCDR T13 p103: 306–406 | preCDR = short endcap |
| tank r | 77→87 / 240, back to 225 | same | 90 / 260 | |
| front window | 1.016 mm Al | same | WP p27: thin Al | |
| gas | C4F8 1.7 atm, 14.99 mg/cm³ | C4F8O 1.5 atm, 13.47 | preCDR p103, WP p27 | GEMC option kept |
| n−1 at 2.04 / 6.2 eV | 2.05e-3 / 2.34e-3 | same (SL_HGC_C4F8O table) | none (preCDR p109 gives no index values) | |
| sectors / φ₁ | 30 / 96° | same | 30 (p103) | |
| mirror R / thickness | 210 / 0.3 CFRP | same | ~210 (p103) | |
| mirror centre (r, z) | 199.226, 210.120 | same | — | |
| mirror φ span | 12°; s1/s16 11.49° | 12° | — | fix |
| image point / tilt | 215.479, 343.741 / 39° | same | — | |
| funnel (full sizes) | Fe, opening 22.3 × 22.3 → 44.82 × 32.0 cm, wall 0.5 cm, 16.18 cm long | same | Lexan on steel (p103) | Trd half sizes 11.15 → 22.41 × 16.0 |
| shield | Fe 23.3² box, 19.9 cm long | 20 cm long | 2 low-C steel layers (p108) | fix |
| PMT plane | 21.3 cm × 40 µm | same | 4×4 2" MaPMT (p103) | 4×4 option |
| channels | 30 × 32×32 cells | 30 planes | 480 (T13) | |
| lite gas | CO2 1 atm, 1.842 mg/cm³ | same (ext) | workshop C068 | |
| lite front section | z 205–270, r 60/70–140 | same | — | |
| lite sectors | 4,5,11,12,19,20,26,27 = segments centred at 42/138/222/318° | 4seg2sec (same) | talks: segments at 0/90/180/270° | **conflict** (45° apart): with 12° tank sectors a 2-sector segment can only be centred at 6° (mod 12°); kept GEMC, list is a parameter (`LiteCherenkov_Sectors`) — decision needed |
| extshift | main chamber 277–386 | −40 cm gcard shift | C068 | |

### 2.7 Electromagnetic calorimeters (FAEC, LAEC)

Plugin `SoLID_HexShashlykEMCal`.

| Parameter | dd4hep | GEMC latest | reference doc | note |
|---|---|---|---|---|
| hexagon apothem / side | 5.413 / 6.25 cm | 5.413 (moved) | preCDR Fig.86 p113/PDF119 "6.25-cm-side" | old plugin 7.19 |
| orientation | sides ∥ y | same (rot 30°) | ANL 200318 | |
| pitch / row pitch | 11.176 / 9.6787 | same | Fig.94 p119 | |
| lattice positions | 1748 (54 rows) | 1748 | p119: 1764/1800 | |
| N SIDIS FA / lite / PVDIS / LA | 1413/1145/1697/307 | same | T15 p111: ≈1360/–/1830/400 | |
| layer stack | 0.005+0.05+0.005+0.012+0.15+0.012 cm | same | T17 p112: 0.224 cm (no paint) | |
| shower length | 45.396 cm (194 layers) | 45.396 | 43.4 cm (T17) | paint = parameter |
| PS Pb / scint / Al | 1.12 / 2.0 / 2.0 cm | same | T16, p118 | old: scint 6.25 cm, no Al |
| rods | 6 × Ø6 mm SS, r 5.6255 | same | p118 "six rods" | |
| FA SIDIS z (Pb front / shower) | 439.182 / 444.302–489.698 | same | T15 (415,475) | |
| PVDIS z | 332.182 / 337.302–382.698 | same | T15 (320,380) | |
| LAEC z | −62.818 / −57.698…−12.302 | same | T15 (−65,−5) | |
| FA r (hex edges) | 89.26–246.12 (lite −225.76) | same | T15 (98,230) | |
| PVDIS r | 94.22–269.89 | same | (110,265) | 1.1 mm to endcap r 270 |
| LAEC r | 81.19–142.39 | same | (83,140) | 1.6 cm to cryostat 144 |
| scintillator | polystyrene 1.06, kB 0.126 mm/MeV | G4_POLYSTYRENE | T17 | |
| shower sampling fraction | 0.203–0.208 (Birks) | 0.253 (code comment, EM) | – | |

### 2.8 Scintillator pads, MRPC and the DDVCS muon system

Plugins `SoLID_ScintPadDisk`, `SoLID_MRPCDisk`, `SoLID_MuonDetector`.

| Parameter | dd4hep | GEMC latest | reference | note |
|---|---|---|---|---|
| FASPD z / R | 428 / 93-216 | 428 / 93-216 (spd_moved) | preCDR p131: qualitative (between HGC and EC), no z/R | |
| FASPD thickness | 0.6 cm | 0.5 cm | preCDR 5 mm; 2026 slide 6 mm | spec choice |
| FASPD segmentation | 60 phi x 4 rings (93/115/135/165/216) | none | preCDR p131: 240 pads | boundaries from GEMC analysis |
| LASPD z / R / t | -65 / 80-135 / 2.0 | same | preCDR p132 | |
| LASPD segmentation | 60 x 6 deg | none | preCDR p132 | sector 1 = 90-96 deg |
| MRPC z / R | 435 / 93-216 | same | preCDR p135: no z; 'inner diameter ~1 m, outer ~2 m' (inconsistent with 93-216 either way) | values from GEMC |
| MRPC stack | preCDR Fig.110, 28.3 mm, 2 x 5 x 0.25 mm gaps | 35 layers, 28.3 mm | preCDR p136 | GEMC list = variant 1 |
| MRPC gas | 90/5/5 vol, 4.26 mg/cm3 | SF6/iC4H10 fractions swapped | preCDR p135 | fixed |
| MRPC strips | PolarGridRPhi 2.8 cm x 7.2 deg | none | preCDR 1650 strips | |
| Iron layer z | 620-656/666-702/712-748 | same | cleo_iron_layout.md | |
| Iron blocks | 2 x 250x534x36 at x +-164.4; 2 x 78.8x177.8x36 at y +-177.8 | same | real 251.5x533.4x35.56 | |
| Scintillators | z 663/709/755, r 100-285, 5 cm, 60 paddles | unsegmented | C040: 60 segments | Nphi=1 -> GEMC |
| LA barrel VP | octagon flat 258.61, z -266.7..170 | flat 280, z -266..182 | magnet spec 4.8 | |
| LA endcap VP | r 286, z 170..619 | r 286, z 209..620 | magnet spec 4.8 | |

## 3. Deliberate deviations from GEMC

Every DIFF in §1 and every "fix" in §2 is one of these.

**CLEO-II solenoid and field maps**

- The real CLEO flats and the NIM coil/cryostat replace the GEMC values (§2.1).
- The v4 fixes are applied: plug 63.32, downstream ring from 183.73, interface from 170.67 with a round outside, and discrete spacer bars.
- The provisional lite-SIDIS Shield5 is moved +0.6 cm to −279.9…−276.9 (CLEOV11 accuracy is ±0.5 cm). Without the move it overlaps the He3 target coils, which end at z −280.0, by 5 mm.

**Targets and beamlines**

1. **Kill disk.** It is a vacuum disk with the `Beamline_KillLimits` limitset (track_length_max=0). ddsim/npsim always register G4StepLimiterPhysics, which provides the needed G4UserSpecialCuts. One difference: G4UserSpecialCuts uses fStopButAlive, so a stopped e+ can still annihilate.
2. **Bore profile.** Windows and kill disks follow the bore profile, which makes them +15 µm larger in r. The J/psi kill disk upstream face is r 28.895 instead of 28.8.
3. **J/psi cell.** It is built as an Al cylinder with an LH2 daughter; the material map is the same. BMD1 (vacuum inside vacuum) is omitted.
4. **NH3 field box.** It is not built. The field is a global overlay, not a replacement inside the box.

**PVDIS baffle**

1. **Slit edges.** Exact polygons replace the GEMC staircase. The blades are inset 30 µm from the rings, and the default sampling is Δr 4 cm and 1° arcs. This costs 10.3 s per 200 events versus 9.7 s for the staircase, and 22.5 s at the spec's 1 cm / 0.25° sampling.
2. **Photon blocker.** Moved from z 330 to 328.4 because it overlapped the FAEC.
3. **Inner absorber.** Its rmin follows the exit pipe plus 1 mm.
4. **Outer absorber.** The zero-thickness tail is dropped.
5. **Plate 1.** It goes in a cavity slot instead of inside TACV.

**GEM trackers**

- Planar trapezoid modules, needed by ACTS/EICrecon. The covered radius is −0.55 % at the inner flat and +0.55 % at the corners.
- Only the drift gap is sensitive.
- The ±1 um ACTS placement offset is compensated, so there is no net shift.
- A 1 mm `step_length_max` limitset is attached to the drift gap. It acts only if a step-limiter physics constructor is enabled.
- Flux virtual planes are not built.

**Light gas Cherenkov (LGC)**

- Windows cover the full faces, and SIDIS gets an exit window. The PVDIS "LEM" placeholder is dropped.
- The shield is shortened by 1.14 cm, which removes the GEMC overlap with the cone.
- Material fixes: PVF, CFRP and UV glass compositions; the photocathode is K2CsSb with the CO2 RINDEX.
- The last ABSLENGTH point is 0.1 mm instead of 0.
- The PVDIS pivot is exactly on the sphere. The tilted M1 wedge shrinks by up to 0.28° per side (no overlap); `clip="true"` is optional.

**Heavy gas Cherenkov (HGC) and SoLID-Lite CO2 Cherenkov**

1. Mirrors 1 and 16 are trimmed by 0.51° at the 0.5° Al wedges. In GEMC they overlap the wedges by 4.3 % of the mirror volume.
2. Shield shortened by 1 mm at its back end. The GEMC shield extrudes 80 µm.
3. The subtracted solids of the booleans are ε = 0.01 cm longer, with the wall planes unchanged.
4. An Al absorber replaces "Kryptonite".
5. The lite sensitive layer is CO2. GEMC kept C4F8O there.
6. Duplicate lite polycone planes are dropped.
7. The flux virtual planes are not built.

**Electromagnetic calorimeters (FAEC, LAEC)**

- ANL id 1748 uses x = 55.88 cm, not GEMC's typo 55.8 cm.
- Only the scintillator slices are sensitive. GEMC makes the whole module its sensitive mother.
- The Mylar gaps use the shared `Mylar` (1.39 g/cm³, GEMC 1.40).
- The PS scintillator and the Al tile share one hexagonal stack per module.
- No virtual planes.

## 4. Open questions and TODO

These are for the collaboration and the subsystem groups.

**CLEO-II solenoid and field maps**

- **Interfaces:** the real cryostat bore is r = 144.0, so there is **no clash**.
  - LAEC envelope / module corner (r 142.39 / 142.6 at z −62.8…−12.3): 1.4–1.6 cm clearance.
  - GEM: r ≤ 135.7 (SIDIS) and 138.8 (PVDIS, z 190).
  - FAEC to nose: 4.2 cm.
  - LGC and HGC fit inside the endcap (r 270).
- **Open issues:**
  - The 3D TOSCA map `solenoid_v4.dat` is missing. v4c runs with a 2D map of its own iron (`magnet_field/`, `field_SoLID_v4c_2D.xml`), with the octagon at its flats; an equal-area octagon changes the field by 0.2–1.3 %.
  - The stage-2 map (`field_SoLIDlite_stage2_2D.xml`) is as provisional as the stage-2 iron.
  - Still unconfirmed: shell thickness, the cryostat diameter (11.8 ft), and bars vs ring.
  - The v4c map reaches z = 800 cm, over the muon system; the muon iron itself is not magnetized.

**Targets and beamlines**

- No field map exists for the SM magnet.
- Still undefined: the J/psi scattering chamber, the PVDIS insulating vacuum chamber (IVC) and the NH3 beam chicane.
- The dimensions of the third Helmholtz pair are unknown.
- Values to confirm: NH3 and LD2 densities, GE180 composition.
- `BafflePlate1_cavity_halfz`, which the baffle now exports, could replace the cavity half-length I currently use (`BafflePlate_halfz` + 20 µm).

**PVDIS baffle**

- The photon-blocker z is a constant (328.4 cm). It cannot reference the EC constants because the baffle is built first, so it needs an update if the FAEC moves.
- The preCDR hot-edge trim (0.6° on even plates) is available but off.
- The BorPoly shielding has no document reference.
- The real PVDIS IVC is not yet defined.

**GEM trackers**

- Write a U/V strip digitizer for EICrecon: pitch 0.4 mm, efficiency 97 %, VMM3 or APV25.
- Optional frames, spacers and staggering (V7 / preCDR) need UVa dimensions.
- Check the PVDIS sector offset against the baffle ribs once the baffle exists.
- Decide the NH3 sector mask.
- The GEMC lite GEM is full 2π. A 4-segment FA version is only a `sectors` copy.
- The `SurfaceManager` built from Python showed empty maps. Surfaces were verified per DetElement instead (`volSurfaceList`); confirm with the `InstallSurfaceManager` plugin.

**Light gas Cherenkov (LGC)**

- Collaboration questions: the 2022 extZmid_7.5/modNose optics (Q1), the cut-cone apex offset (Q3, parameter `apex_z`), and the sector phase (Q4, parameter `LGC_phi_sec1`).
- Not modelled: tank walls and supports, the preCDR steel shield, the real 6.08 mm H12700 pixels, and the virtual planes.
- Physics checks still to do: a run with the full-SoLID field, a comparison with GEMC using the same gun, and digitization (QE, packing, trigger).

**Heavy gas Cherenkov (HGC) and SoLID-Lite CO2 Cherenkov**

- No field in the tests; the physical-config integration (FAEC/SPD/MRPC near z 426) is still to be checked.
- Mirror misalignment (pivot (220.48, 418.92)) and the virtual plane are not implemented.
- The QE and pixel summing belong in digitization.
- There is no measured C4F8 dispersion or absorption, and the CFRP composition is questionable.
- The lite optics were designed for z = −350.

**Electromagnetic calorimeters (FAEC, LAEC)**

- Paint vs paper thickness, rod diameter and Mylar vs paper are still open. They are parameters: `*_paint_thick`, `*_rods_enable`, `*_rod_radius`.
- Not modelled: WLS fibres, end plates, the 4 cm rear plate, and the inter-module gap material.
- The PVDIS photon blocker overlap with the Pb radiator depends on the baffle fragment (z_c 328.40).
- No Geant4 `checkOverlaps` was run, and no check was done against the magnet geometry. The clearances above are geometric.

**Scintillator pads, MRPC and the DDVCS muon system**

- Enable the `muon_gem_DDVCS.xml` include once it exists (GEM gaps 657.7-659.3 etc. are kept free).
- I could not check overlaps against the magnet, Cherenkov or EC fragments because they are not in the stage. Clearances: FASPD front 427.7 vs HGC end 426.0; LA VPs vs corrected yoke.
- Still unknown: SPD radial boundaries and wrapping, LASPD PMTs, WP thin-gap MRPC, real MRPC glass, real CLEO slab steel.
- Some zero-deposit secondary hits on the vacuum planes carry t = 0.
