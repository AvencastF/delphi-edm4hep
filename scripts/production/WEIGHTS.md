# Production weights in EDM4hep

Added 2026-09-23 for future production, not retrofitting test samples.

The production wrapper requires `task.json`, `generation.json` and
`event-weights.csv` from the same chunk. It validates complete sequential
generator identities, original trial identities when provided, finite signed
weights, cross sections and generator sumw/sumw2. The converter joins by the
expected negative DELSIM run and event number, never by ROOT row position.
Missing, duplicate or wrong-run weight identities fail conversion. Data and
standalone conversions without `--mc-weights` keep their existing behavior.

## Event frame

| Parameter | Meaning |
|---|---|
| `mc_gen_weight` | Nominal generator weight, double precision, including its sign. No luminosity, efficiency, spin or theory correction applied. |
| `mc_generator_trial` | Original BHWIDE trial ID; -1 when the generator input does not provide one. |

The existing `sDST_EVT_runNumber` and `sDST_EVT_eventNumber`, together with the
chunk provenance, identify the event. Nominal-only is explicit: the current
generators do not provide PDF/scale/spin variation-weight vectors here.

## Single file metadata frame

- `mc_weight_schema=1`, `mc_weight_names=[nominal]`, `mc_run_number`.
- `mc_cross_section_pb`, `mc_cross_section_error_pb`: this chunk's generator
  estimate and reported integration statistical uncertainty, not a theory error
  or an automatically combined campaign estimate.
- `mc_production_json`: complete generation summary and production task,
  including process, energy point, seed, run, physics card and card hash;
  BHWIDE full-trial-batch statistics and output sampling when supplied;
  source CSV checksum and normalization conventions.
- Parallel `mc_generated_event_ids`, `mc_generated_trial_ids`,
  `mc_generated_weights`: the complete generated chunk ledger, including records
  absent from the converted output. This makes bookkeeping independent of CSVs.
- `mc_converted_event_ids`, `mc_missing_ids_before_last_output`,
  `mc_ids_after_last_output`: distinguish converted records, internal gaps and
  the remaining tail. The last-output boundary does not prove which tail records
  were attempted; no simulated-efficiency interpretation is assigned to it.
- Each prefix `mc_generated_`, `mc_converted_`, `mc_missing_`,
  `mc_after_last_output_` has `count`, `sumw`, `sumw2`, `sumabsw`,
  `negative_count`. Signed sumw is the normalization quantity, not sumabsw.

## Analysis normalization

Pool compatible chunks/tags at the same process, energy and physics settings;
retain one metadata record per distinct chunk when merging. Do not assign full
luminosity independently to every chunk or sum cross sections across replicas.
Use an explicitly chosen compatible cross-section estimate and a normalization
population defined before physics cuts. Account for the deliberate simulation
subset and known failures; blindly normalizing surviving events can hide bias.
The complete ledger permits those choices without inventing lost-event weights.

The usual representative-sample formula is
`weight_nominal = luminosity * cross_section * gen_weight / sum_gen_weights`.
No separate filter efficiency or branching fraction is automatically multiplied:
first establish whether it is already included in the supplied cross section.
No k-factor is presumed. Luminosity, detector/PID efficiency scale factors,
trigger corrections, calibration variations, and deliberate spin reweights are
analysis products and must not overwrite `mc_gen_weight`. If a new generator
supplies named systematic weights, extend the schema to preserve those explicitly.

## Validation and deployment

Local tests exercise negative weights, non-unit weights, skipped IDs, the tail,
wrong runs, duplicate/missing IDs, nonfinite weights and inconsistent summaries.
The production ROOT audit reads the written event weights and generated ledger
back and checks them against production inputs, plus converted sums and IDs.

The full converter must be rebuilt in the approved NERSC compute environment
and the podio/ROOT round trip smoke-tested before deployment. Existing reusable
build identity checks reject the old binary after this source change. No old
samples, frozen campaigns, build installations or jobs are changed by this edit.
Scope is the production longDST pass-1 path; a future pass-2 or file-merging
workflow must explicitly carry the weight metadata as well as event parameters.
