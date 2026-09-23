# LongDST production conversion boundary

`convert.sh BINARY INPUT.ldst OUTPUT.root AUDIT_DIRECTORY` runs inside the pinned
DELPHI OCI image. Its working directory must contain `my_events.fadgen` from the
same generation task, plus `task.json`, `generation.json`, and `event-weights.csv`.
Weights and normalization provenance are embedded in the EDM4hep frames; see
[WEIGHTS.md](WEIGHTS.md). Source DELPHI first, then the pinned key4hep release; do not
reuse the simulation shell for conversion. The orchestrator provides isolated subprocess environments, mounts, logging,
timing, counts and publication. In the node-container runtime it captures this
wrapper's DELPHI-then-key4hep environment once per node and supplies
`DELPHI_CONVERTER_ENV_READY=1`; standalone calls still source both environments.

Pass 1 is the verified longDST conversion path. Output is an EDM4hep ROOT event
file; no second fullDST pass or flat ntuple is implied. The wrapper runs the
collection-alignment audit plus the event/content audit inherited from the
successful 1994 test. The orchestrator additionally requires exact requested
counts, unique run/event IDs and FADGEN truth momentum correspondence.

The audit does not certify all physics fields: known truth-vertex scale/origin,
HPC energyError and VD coordinate issues remain visible and are not corrected by
these operational scripts. Do not claim an empty collection is necessarily a
converter failure. Detector mappings are unchanged; the converter now supports
an optional production-weight input and metadata hook.

Wrapper and weight output passed the 2026-09-23 all-process NERSC test:
8,000 converted events and 80 independently checked ROOT metadata frames.
See [the report](validation/20260923/REPORT.md) for scope and remaining issues.
