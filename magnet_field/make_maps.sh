#!/usr/bin/env bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# make_maps.sh: SoLID solenoid field maps from the DD4hep magnet geometry, in one command.
#
#   magnet_field/make_maps.sh [options] [model ...]
#
#   model    v4c     full SoLID, yoke v4c        compact/magnet/solenoid_v4.xml        -> solenoid_SoLID_v4c_2D.dat
#            stage1  SoLID-Lite stage 1          compact/magnet/solenoid_nomod.xml     -> solenoid_SoLIDlite_stage1_2D.dat
#            stage2  SoLID-Lite stage 2 (SIDIS)  compact/magnet/solenoid_lite_sidis.xml -> solenoid_SoLIDlite_stage2_2D.dat
#            (default: all three)
#
#   For each model:
#     1. extract_rz_model.py turns the DD4hep magnet (test/test_magnet_*.xml, without its field) into an r-z deck.
#        It needs DD4hep python and the solid_dd4hep plugins, so it runs in the eic_xl container through
#        containers/solid-shell (or directly when dd4hep can be imported).
#     2. solve_rz.py solves the deck with an open boundary and writes the GEMC map (R 0..500 cm, Z -600..800 cm,
#        1 cm, gauss) that compact/fields/field_<name>.xml reads. It needs python3 with numpy, scipy and gmsh
#        (requirements.txt); the eic_xl container has no gmsh, so this step runs on the host.
#
#   -o, --out DIR        output directory (default: magnet_field/work/maps): decks/, solutions/, logs/, the maps
#   -f, --fieldmaps DIR  also copy the maps into DIR, e.g. $DETECTOR_PATH/fieldmaps (where SoLID_FieldMapB looks)
#   -p, --prefix DIR     solid_dd4hep install prefix for step 1 (default: solid-shell's, $SOLID_PREFIX or
#                        <repo>/install; build it with containers/install.sh --build)
#   -i, --image IMAGE    container image for step 1 (default: solid-shell's choice)
#   -n, --no-extract     skip step 1 and solve the decks committed in magnet_field/models/ (no container needed).
#                        First the repository files each deck was extracted from (sha256 in models/<deck>.json;
#                        comments, <info> and field includes ignored) are compared with the checkout: if a magnet
#                        file (compact/magnet/, test/test_magnet_*) changed, the deck is stale and the script
#                        stops (run without -n); other changed files give a warning.
#   --allow-stale        with -n: only warn when the magnet files changed
#   --air-side           maps with air-side values at the iron nodes next to air (write_map.py --air-side):
#                        accurate in air up to the iron faces, wrong in the iron within about 2 cm of them
#   --python PY          python for step 2 (default: $MAGNET_FIELD_PYTHON, else python3)
#   -j N                 run N solves at the same time (default 1; each takes about 2 GB and 3 threads)
#   -h, --help           this text
#
# Runtime on a 4-core host: about 6 min for the three models (step 1 about 30 s, step 2 70-130 s per model).
set -euo pipefail

MF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${MF}/.." && pwd)"
OUT="${MF}/work/maps"
FIELDMAPS=""
PREFIX=""
IMAGE=""
EXTRACT=1
ALLOW_STALE=0
AIR_SIDE=0
PY="${MAGNET_FIELD_PYTHON:-python3}"
JOBS=1
MODELS=()

usage() { sed -n '3,37p' "$0" | sed 's/^# \{0,1\}//'; }
die() { echo "make_maps.sh: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out)       OUT="$2"; shift 2 ;;
    -f|--fieldmaps) FIELDMAPS="$2"; shift 2 ;;
    -p|--prefix)    PREFIX="$2"; shift 2 ;;
    -i|--image)     IMAGE="$2"; shift 2 ;;
    -n|--no-extract) EXTRACT=0; shift ;;
    --allow-stale)  ALLOW_STALE=1; shift ;;
    --air-side)     AIR_SIDE=1; shift ;;
    --python)       PY="$2"; shift 2 ;;
    -j)             JOBS="$2"; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
    v4c|stage1|stage2) MODELS+=("$1"); shift ;;
    *)              usage >&2; die "unknown argument '$1'" ;;
  esac
done
[[ ${#MODELS[@]} -gt 0 ]] || MODELS=(v4c stage1 stage2)
[[ "${JOBS}" =~ ^[1-9][0-9]*$ ]] || die "-j needs a positive integer"
# relative paths refer to the current directory (the extraction runs from the output directory)
abspath() { if [[ -e "$1" ]]; then (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); else echo "$1"; fi; }
[[ -n "${PREFIX}" ]] && PREFIX="$(abspath "${PREFIX}")"
[[ -n "${IMAGE}" ]] && IMAGE="$(abspath "${IMAGE}")"
if [[ -n "${FIELDMAPS}" ]]; then mkdir -p "${FIELDMAPS}"; FIELDMAPS="$(cd "${FIELDMAPS}" && pwd)"; fi

# model -> COMPACT (relative to the repository), DECK name, MAP name, extra extraction options XOPT
model_info() {
  XOPT=""
  case "$1" in
    v4c)    COMPACT=test/test_magnet_v4.xml;         DECK=solid_v4c;         MAP=solenoid_SoLID_v4c_2D ;;
    stage1) COMPACT=test/test_magnet_nomod.xml;      DECK=solid_lite_stage1; MAP=solenoid_SoLIDlite_stage1_2D ;;
    stage2) COMPACT=test/test_magnet_lite_sidis.xml; DECK=solid_lite_stage2; MAP=solenoid_SoLIDlite_stage2_2D
            XOPT="--commented-addon cleov11_correction_coil" ;;
  esac
}
MAP_R="0:500:1"
MAP_Z="-600:800:1"
POINTS=(--point 0,-315 --point 0,-350 --point 215,300 --point 215,344)

mkdir -p "${OUT}/decks" "${OUT}/solutions" "${OUT}/logs" "${OUT}/run"
OUT="$(cd "${OUT}" && pwd)"
T0=${SECONDS}

"${PY}" -c "import numpy, scipy, gmsh" 2>/dev/null ||
  die "${PY} cannot import numpy, scipy and gmsh: pip install -r ${MF}/requirements.txt (in a venv), or --python PY"

# ---------------------------------------------------------------------------------------------- 1. extraction
if [[ ${EXTRACT} -eq 1 ]]; then
  if python3 -c "import dd4hep" 2>/dev/null; then
    run_dd4hep() { (cd "${OUT}/run" && bash -c "$1"); }
    echo "== extraction with the DD4hep python of this environment"
  else
    shell="${REPO}/containers/solid-shell"
    [[ -x "${shell}" ]] || die "no DD4hep python here and no ${shell}: run inside eic_xl, or use --no-extract"
    pfx="${PREFIX:-${SOLID_PREFIX:-${REPO}/install}}"
    [[ -f "${pfx}/bin/thissolid.sh" ]] ||
      die "no ${pfx}/bin/thissolid.sh: build solid_dd4hep first (containers/install.sh --build), pass -p PREFIX," \
          "or use --no-extract"
    sargs=(-p "${pfx}")
    [[ -n "${IMAGE}" ]] && sargs+=(-i "${IMAGE}")
    run_dd4hep() { (cd "${OUT}/run" && "${shell}" "${sargs[@]}" -- "$1"); }
    echo "== extraction in the container (containers/solid-shell, prefix ${pfx})"
  fi
  for m in "${MODELS[@]}"; do
    model_info "${m}"
    t=${SECONDS}
    deck="${OUT}/decks/${DECK}.am"
    log="${OUT}/logs/${DECK}_extract.log"
    # DETECTOR_PATH = this checkout, so the compact files of the checkout are read (test/ is not installed)
    if ! run_dd4hep "export DETECTOR_PATH='${REPO}'; python3 '${MF}/extract_rz_model.py' extract \
           '${REPO}/${COMPACT}' -o '${deck}' --name ${DECK} --no-field --no-sampling ${XOPT}" > "${log}" 2>&1; then
      tail -20 "${log}" >&2
      die "extraction of ${m} failed (log: ${log})"
    fi
    ref="${MF}/models/${DECK}.am"
    if [[ -f "${ref}" ]] && diff -q <(grep -v '^;' "${deck}" | tail -n +2) <(grep -v '^;' "${ref}" | tail -n +2) \
         > /dev/null; then
      same="same regions as models/${DECK}.am"
    else
      same="DIFFERENT from models/${DECK}.am: the DD4hep magnet changed; review ${deck} and update models/"
    fi
    printf '  %-7s %s -> %s (%d s; %s)\n' "${m}" "${COMPACT}" "decks/${DECK}.am" $((SECONDS - t)) "${same}"
  done
  DECKDIR="${OUT}/decks"
else
  DECKDIR="${MF}/models"
  echo "== no extraction: decks from ${DECKDIR}; checking them against the magnet files of ${REPO}"
  decks=()
  for m in "${MODELS[@]}"; do model_info "${m}"; decks+=("${MF}/models/${DECK}.json"); done
  rc=0
  "${PY}" "${MF}/extract_rz_model.py" check --repo "${REPO}" "${decks[@]}" || rc=$?
  if [[ ${rc} -eq 2 ]]; then
    if [[ ${ALLOW_STALE} -eq 1 ]]; then
      echo "  WARNING: solving stale decks (--allow-stale): the maps do not describe this checkout's magnet" >&2
    else
      die "the committed decks no longer match the DD4hep magnet of this checkout: run without -n (re-extracts" \
          "the decks in the eic_xl container; then review them and update models/), or --allow-stale"
    fi
  elif [[ ${rc} -ne 0 && ${rc} -ne 1 ]]; then
    die "the deck check failed (exit ${rc})"
  fi
fi

# --------------------------------------------------------------------------------------------- 2. solve + map
rev="$(git --no-optional-locks -C "${REPO}" rev-parse --short=12 HEAD 2>/dev/null || echo "no git")"
[[ -n "$(git --no-optional-locks -C "${REPO}" status --porcelain -- compact 2>/dev/null)" ]] &&
  rev="${rev}, modified compact files"
[[ ${EXTRACT} -eq 0 && ${rc:-0} -eq 2 ]] && rev="${rev}, STALE deck"
MAPOPT=()
[[ ${AIR_SIDE} -eq 1 ]] && MAPOPT=(--map-air-side)
solve_one() {
  local m=$1 t=${SECONDS} rc=0
  model_info "${m}"
  local log="${OUT}/logs/${DECK}_solve.log" src="committed deck, solid_dd4hep ${rev}"
  [[ ${EXTRACT} -eq 1 ]] && src="extracted from solid_dd4hep ${COMPACT}, ${rev}"
  "${PY}" -u "${MF}/solve_rz.py" "${DECKDIR}/${DECK}.am" --mode open \
    -o "${OUT}/solutions/${DECK}_open.npz" --map "${OUT}/${MAP}.dat" --map-name "${MAP}" \
    --map-comment "rzfield open-boundary solution of ${DECK}.am (${src}), $(date +%F)" \
    --map-r "${MAP_R}" --map-z "${MAP_Z}" ${MAPOPT[@]+"${MAPOPT[@]}"} "${POINTS[@]}" > "${log}" 2>&1 || rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "  ${m}: solve FAILED (exit ${rc}, log ${log})" >&2
    return 1
  fi
  printf '  %-7s %s.dat (%d s)\n' "${m}" "${MAP}" $((SECONDS - t))
}
echo "== solving (open boundary) and writing the maps, ${JOBS} at a time"
fail=0
pids=()
for m in "${MODELS[@]}"; do
  solve_one "${m}" &
  pids+=($!)
  if [[ ${#pids[@]} -ge ${JOBS} ]]; then
    for p in "${pids[@]}"; do wait "${p}" || fail=1; done
    pids=()
  fi
done
for p in ${pids[@]+"${pids[@]}"}; do wait "${p}" || fail=1; done
[[ ${fail} -eq 0 ]] || die "at least one solve failed"

# ------------------------------------------------------------------------------------------------- summary
echo "== maps in ${OUT} (R ${MAP_R}, Z ${MAP_Z} cm, gauss)"
sols=()
for m in "${MODELS[@]}"; do model_info "${m}"; sols+=("${OUT}/solutions/${DECK}_open.npz"); done
MF_DIR="${MF}" "${PY}" - "${sols[@]}" <<'PYEOF'
import os, sys
sys.path.insert(0, os.environ["MF_DIR"])
import numpy as np
from rzfield.field import RZSolution
print("  %-20s %10s %10s %10s %12s %12s %9s %7s" % ("deck", "Bz(0,0)", "Bz(0,-315)", "Bz(0,-350)",
      "|B|(215,344)", "int Bz dz", "energy", "Newton"))
for path in sys.argv[1:]:
    sol = RZSolution.load(path)
    z = np.arange(-600.0, 800.5, 1.0)
    br, bz = sol.B_cm_gauss(np.array([0.0, 0.0, 0.0, 215.0]), np.array([0.0, -315.0, -350.0, 344.0]))
    _, bza = sol.B_cm_gauss(np.zeros_like(z), z)
    integ = float(np.sum(0.5 * (bza[1:] + bza[:-1]) * np.diff(z))) * 1e-6
    print("  %-20s %8.1f G %8.2f G %8.2f G %10.2f G %8.4f T m %6.2f MJ %7d" % (
        os.path.basename(path).replace("_open.npz", ""), bz[0], bz[1], bz[2], np.hypot(br[3], bz[3]), integ,
        sol.info["energy_J"] / 1e6, len(sol.info["newton"])))
PYEOF
for m in "${MODELS[@]}"; do
  model_info "${m}"
  f="${OUT}/${MAP}.dat"
  printf '  %-34s %10d bytes  sha256 %s\n' "$(basename "${f}")" "$(wc -c < "${f}" | tr -d ' ')" \
    "$( { sha256sum "${f}" 2>/dev/null || shasum -a 256 "${f}"; } | cut -c1-16)"
  if [[ -n "${FIELDMAPS}" ]]; then
    install -m 644 "${f}" "${FIELDMAPS}/"
  fi
done
[[ -n "${FIELDMAPS}" ]] && echo "== copied to ${FIELDMAPS}"
echo "== done in $((SECONDS - T0)) s"
