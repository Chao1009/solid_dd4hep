#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# install.sh: set up the container for containers/solid-shell (the EIC software image eic_xl, run with
# Apptainer/Singularity), and optionally build this repository inside it.
#
#   containers/install.sh [options]
#
#   -t, --tag TAG      image release (default: 26.10.0-stable for eic_xl, the one this repository is tested with;
#                      latest for solid_xl)
#       --solid-xl [OWNER]  get solid_xl instead: eic_xl plus ffmpeg, gmsh and a preinstalled SoLID geometry,
#                      from ghcr.io/OWNER/solid_xl (default OWNER: chao1009). It is not on CVMFS, so it is
#                      always downloaded, also on ifarm (choose a /work or /volatile directory with -d)
#   -d, --dir DIR      where to store the .sif image (default: <repo>/containers/images). The image is about 5 GB;
#                      the download needs about 15 GB free in DIR while it runs
#   -i, --image IMAGE  use an image you already have instead: a .sif file, a /cvmfs/... directory or docker://...
#       --no-cvmfs     download even when the image is available on CVMFS (default: use CVMFS when mounted)
#       --build        also build solid_dd4hep and install it into <repo>/install, inside the container
#   -f, --force        download again even if the .sif file exists; skip the free-space check
#   -h, --help         show this help
#
# The chosen image is saved in containers/solid-shell.local (untracked), which containers/solid-shell reads.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/.." && pwd)"
tag=""
dir="${repo}/containers/images"
image=""
source="docker://eicweb/eic_xl"   # --source: another image repository (for testing)
use_cvmfs=1
build=0
force=0
solid_xl=0
owner="chao1009"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--tag)    tag="$2"; shift 2 ;;
    -d|--dir)    dir="$2"; shift 2 ;;
    -i|--image)  image="$2"; shift 2 ;;
    --source)    source="$2"; shift 2 ;;
    --solid-xl)  solid_xl=1
                 if [[ $# -gt 1 && "$2" != -* ]]; then owner="$2"; shift 2; else shift; fi ;;
    --no-cvmfs)  use_cvmfs=0; shift ;;
    --build)     build=1; shift ;;
    -f|--force)  force=1; shift ;;
    -h|--help)   sed -n '3,21p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)           echo "install.sh: unknown option '$1' (see --help)" >&2; exit 2 ;;
  esac
done

if [[ ${solid_xl} -eq 1 ]]; then
  source="docker://ghcr.io/${owner,,}/solid_xl"
  tag="${tag:-latest}"
fi
tag="${tag:-26.10.0-stable}"

say() { echo "install.sh: $*"; }

runner="$(command -v apptainer || command -v singularity || true)"
if [[ -z "${runner}" ]]; then
  echo "install.sh: Apptainer/Singularity not found. Install Apptainer" \
       "(https://apptainer.org/docs/admin/main/installation.html) or, with Docker, see containers/README.md." >&2
  exit 1
fi

if [[ -z "${image}" ]]; then
  cvmfs_image="/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:${tag}"
  if [[ ${use_cvmfs} -eq 1 && "${source}" == "docker://eicweb/eic_xl" && -d "${cvmfs_image}" ]]; then
    image="${cvmfs_image}"
    say "using ${image} from CVMFS (nothing to download)"
  else
    image="${dir}/$(basename "${source}")-${tag}.sif"
    if [[ -f "${image}" && ${force} -eq 0 ]]; then
      say "${image} already exists (--force downloads it again)"
    else
      mkdir -p "${dir}"
      avail_kb="$(df -Pk "${dir}" 2>/dev/null | awk 'NR == 2 {print $4}')" || avail_kb=""
      if [[ -n "${avail_kb}" && ${force} -eq 0 && $(( avail_kb / 1024 / 1024 )) -lt 15 ]]; then
        echo "install.sh: only $(( avail_kb / 1024 / 1024 )) GB free in ${dir}; the download needs about 15 GB." \
             "Choose another place with -d DIR (or skip this check with --force)." >&2
        exit 1
      fi
      # keep the layer cache and temporary files next to the image, not in a small \$HOME or /tmp
      tmp="$(mktemp -d "${dir}/.pull.XXXXXX")"
      trap 'rm -rf "${tmp}"' EXIT
      say "downloading ${source}:${tag} into ${image} (about 5 GB, this takes a while)"
      APPTAINER_TMPDIR="${tmp}" SINGULARITY_TMPDIR="${tmp}" \
      APPTAINER_CACHEDIR="${tmp}/cache" SINGULARITY_CACHEDIR="${tmp}/cache" \
        "${runner}" pull "${tmp}/image.sif" "${source}:${tag}"
      mv -f "${tmp}/image.sif" "${image}"
    fi
  fi
fi

# the EIC images set their environment in /etc/profile.d; solid-shell sources /etc/profile the same way
if "${runner}" exec "${image}" bash -c '. /etc/profile >/dev/null 2>&1; command -v ddsim' >/dev/null 2>&1; then
  say "image OK (ddsim found)"
else
  say "warning: no ddsim in ${image}; is it an EIC image?"
fi

conf="${here}/solid-shell.local"
{
  echo "# written by containers/install.sh on $(date +%Y-%m-%d); read by containers/solid-shell"
  printf 'SOLID_SHELL_IMAGE=%q\n' "${image}"
} > "${conf}"
say "saved the image choice in ${conf}"

if [[ ${build} -eq 1 ]]; then
  say "building solid_dd4hep into ${repo}/install"
  cd "${repo}"
  cmd='cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install -DCMAKE_BUILD_TYPE=Release'
  cmd+=' && cmake --build build -j"$(nproc)" && cmake --install build'
  "${here}/solid-shell" -- "${cmd}"
fi

cat <<EOF

Done. Open a shell with the SoLID environment (the default configuration is solid = SIDIS 3He):
  containers/solid-shell
  containers/solid-shell -c solid_JPsi_lite
EOF
if [[ ${build} -eq 0 ]]; then
  echo "To build the geometry run containers/install.sh --build, or see \"Quick start\" in README.md."
fi
