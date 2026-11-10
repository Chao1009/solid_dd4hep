# Container

SoLID runs in the EIC software image `eic_xl` (DD4hep, Geant4, ddsim/npsim, EICrecon, ACTS). Nothing else needs to
be installed on the host.

## Install (Linux, JLab ifarm)

Needs Apptainer or Singularity. From the repository root:

```bash
containers/install.sh --build                     # get eic_xl:26.10.0-stable and build this repository in it
containers/solid-shell -c solid_JPsi_lite          # shell with the SoLID environment for one configuration
containers/solid-shell -- ddsim --help             # or run a single command
```

`install.sh` uses the image from CVMFS when it is mounted (ifarm), so nothing is downloaded. Otherwise it downloads
the image (about 5 GB) into `containers/images/`. It saves its choice in `containers/solid-shell.local` (not tracked
by git), where `solid-shell` finds it. `--build` installs the geometry into `install/`; run it again after you change
the code. Useful options (`containers/install.sh --help` lists all of them):

| option | effect |
|---|---|
| `-t TAG` | another eic_xl release, e.g. `26.11.0-stable` |
| `-d DIR` | store the image elsewhere (the download needs about 15 GB free there) |
| `-i IMAGE` | use an image you already have: a `.sif` file, a CVMFS path or `docker://...` |
| `--no-cvmfs` | download even when CVMFS is mounted |

## Docker (macOS)

```bash
docker run -it --rm -v "$PWD:$PWD" -w "$PWD" eicweb/eic_xl:26.10.0-stable
```

Then build and set up as in "Quick start" in the top-level README. The image is amd64 only; on Apple silicon add
`--platform linux/amd64` (it runs under emulation, slowly).

## solid_xl image

`solid_xl` is our own image: eic_xl plus what eic_xl does not ship, with this repository preinstalled in
`/opt/solid`. It needs nothing from the eic_xl or ifarm maintainers.

| added to eic_xl | used by |
|---|---|
| ffmpeg | MP4 output of the event-display animations |
| gmsh (and libXcursor, libXinerama) | the field solver, `magnet_field/solve_rz.py` |
| solid_dd4hep in `/opt/solid` | `ddsim` works without building anything |
| Gaudi and k4FWCore, only with `--build-arg WITH_GAUDI=1` | reconstruction studies ([docs/reconstruction.md](../docs/reconstruction.md)) |

Get it, on a laptop or on ifarm:

```bash
containers/install.sh --solid-xl                  # ghcr.io/chao1009/solid_xl:latest -> containers/images/
containers/install.sh --solid-xl -t v0.4.0 -d /volatile/halla/solid/$USER/images    # a release, stored elsewhere
containers/solid-shell -c solid_JPsi_lite          # the preinstalled geometry
```

- The image is not on CVMFS, so it is always downloaded (about 5 GB, and about 15 GB free space while it is
  converted). On ifarm, put it on `/work` or `/volatile` with `-d`; home directories are too small.
- With an install of your own checkout in `install/` (`containers/install.sh --solid-xl --build`), `solid-shell`
  uses that instead of `/opt/solid`. `-p PREFIX` names another install.
- `--solid-xl OWNER` takes the image from `ghcr.io/OWNER/solid_xl`, for example a fork's or JeffersonLab's.

### Publishing

`.github/workflows/image.yml` builds the image from `containers/Dockerfile`, tests it (ffmpeg, gmsh, loading
`solid_JPsi_lite.xml` and its field) and pushes it to the container registry of the repository's owner:

| trigger | published as |
|---|---|
| a pushed tag `v*` | `ghcr.io/<owner>/solid_xl:<tag>` and `:latest` |
| "Run workflow" on the Actions page, or `gh workflow run image.yml -f image_tag=dev` | `ghcr.io/<owner>/solid_xl:<image_tag>`; `-f latest=true` also moves `:latest`; `-f eic_tag=26.11.0-stable` builds on another eic_xl release |

- The workflow runs in the repository it is started in, so started in a fork it publishes under the fork owner's
  account. The images are published from the fork Chao1009/solid_dd4hep for now.
- A new package is private. Make it public once: GitHub profile, Packages, `solid_xl`, Package settings, Change
  visibility. Apptainer cannot pull a private image without a login.
- Each image is built on one eic_xl release (label `org.opencontainers.image.base.name`). To move to a new DD4hep
  or Geant4, change `EIC_TAG` in the Dockerfile and `EIC_PLATFORM` in `ci.yml` together.

To build it by hand, see the commands in the Dockerfile header.
