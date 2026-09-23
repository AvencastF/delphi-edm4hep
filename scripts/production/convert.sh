#!/bin/bash
# Container-only longDST pass 1. No fullDST pass is implied.
set -eo pipefail
data=0
if [ "${1:-}" = --data ]; then data=1; shift; fi
if { [ "$data" = 0 ] && [ "$#" -ne 4 ]; } ||
   { [ "$data" = 1 ] && [ "$#" -ne 3 ] && [ "$#" -ne 5 ]; }; then
  echo 'usage: convert.sh BINARY INPUT.ldst OUTPUT.root AUDIT_DIRECTORY
       convert.sh --data BINARY INPUT.al OUTPUT.root [-n MAX_EVENTS]' >&2
  exit 2
fi
here=$(cd "$(dirname "$0")" && pwd)
export DELPHI_INSTALL_DIR=/delphi
if [ "${DELPHI_CONVERTER_ENV_READY:-}" != 1 ]; then
  source /delphi/setup.sh
  source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
fi
unset CXXFLAGS CFLAGS LDFLAGS
if [ "$data" = 1 ]; then
  # ponytail: reuse pass 1 and its alignment audit; data have no generator ledger.
  binary=$1 input=$2 output=$3
  shift 3
  if [ "$#" -gt 0 ] && { [ "$1" != -n ] || [[ ! "$2" =~ ^[1-9][0-9]*$ ]]; }; then
    echo 'expected -n followed by a positive event count' >&2; exit 2
  fi
  if [ -e "$output" ] || [ -L "$output" ]; then echo "output already exists: $output" >&2; exit 2; fi
  "$binary" "$input" "$output" "$@"
else
  python3 "$here/prepare_weights.py"
  "$1" "$2" "$3" --mc-weights mc-weights.txt
  python3 "$here/audit_root.py" "$3" my_events.fadgen "$4"
  output=$3
fi
python3 "$here/align_audit.py" "$output"
