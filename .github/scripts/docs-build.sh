#!/bin/bash
# Documentation build, run on the self-hosted runner.
#
# The collection map uses the externally retained conversion references.
# Doxygen builds the site with the pinned key4hep toolchain.
#
# Overrides: DELPHI_CI_REFS (reference store).
set -e

REPO="${GITHUB_WORKSPACE:-$PWD}"
REFS="${DELPHI_CI_REFS:-/home/delphi-ci/refs}"
MAP="${REPO}/docs/collection-map"

PIN="${REPO}/.github/key4hep-production-release"
K4="$(tr -d '[:space:]' < "${PIN}" 2>/dev/null || true)"
[ -n "${K4}" ] || { echo "no key4hep release in ${PIN}" >&2; exit 1; }

# `source setup.sh` with no arguments passes ours through, so clear them first.
set --
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${K4}" >/dev/null 2>&1 || true
[ -n "${KEY4HEP_STACK:-}" ] || { echo "key4hep ${K4} did not set up" >&2; exit 1; }

cd "${REPO}"
python3 "${MAP}/extract.py" --refs "${REFS}" --out "${MAP}/collection_map.json"
python3 "${MAP}/render.py"  --map "${MAP}/collection_map.json" \
                            --out "${MAP}/collection_map.html"
doxygen docs/Doxyfile
