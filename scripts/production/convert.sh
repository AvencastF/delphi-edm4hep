#!/bin/bash
# Container-only longDST pass 1. No fullDST pass is implied.
set -eo pipefail
if [ "$#" -ne 4 ]; then echo 'usage: convert.sh BINARY INPUT.ldst OUTPUT.root AUDIT_DIRECTORY' >&2; exit 2; fi
export DELPHI_INSTALL_DIR=/delphi
if [ "${DELPHI_CONVERTER_ENV_READY:-}" != 1 ]; then
  source /delphi/setup.sh
  source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
fi
unset CXXFLAGS CFLAGS LDFLAGS
python3 "$(dirname "$0")/prepare_weights.py"
"$1" "$2" "$3" --mc-weights mc-weights.txt
python3 "$(dirname "$0")/audit_root.py" "$3" my_events.fadgen "$4"
python3 "$(dirname "$0")/../../delphi_edm4hep/tests/align_audit.py" "$3"
