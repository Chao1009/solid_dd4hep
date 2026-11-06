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

`Dockerfile` builds `solid_xl`: eic_xl with this repository preinstalled in `/opt/solid`, and optionally Gaudi and
k4FWCore (`--build-arg WITH_GAUDI=1`). Build commands are in the Dockerfile header. `solid-shell -i <image>` uses the
preinstalled geometry automatically.
