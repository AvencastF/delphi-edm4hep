# All-process longDST → EDM4hep weight validation

Date: 2026-09-23 (project timezone). **PASS for weight persistence and bounded
chain/content checks.** This does not establish data/MC agreement or spin closure.

## Execution and Git provenance

- Full-chain interactive allocation **58779930**, node **nid004186**, account
  m5019, one CPU node, 60-minute limit. Completed exit 0 in **8m15s**, including
  rebuilding the converter. One container/node with 80 independent workers.
- Independent metadata read-back: shared-interactive allocation **58780041**,
  node **nid004164**, completed exit 0 in **23 seconds**.
- Local commits were pushed to `feature/nersc-long-production`, then pulled with
  `git pull --ff-only` into clean NERSC checkouts. Existing untracked production
  work was preserved.
- Executed simulation revision: `87eac4bf567b48729c945b4a664523668df6428b`.
- Executed converter revision: `13def39d95d439e95eddf5c842e1a6d25f3213a3`.
- Independent metadata checker: `ebdee7c` (after correcting singleton parameter
  handling in the checker; no converter or generated-file correction was needed).
- Same pinned DELPHI OCI image and key4hep 2026-04-08; PYTHIA 8.315 and the
  previously verified BHWIDE 1.05 installation. Detector v94c, ECM=91.202 GeV.
- Build tests: six passed, three sample-dependent tests skipped. Every generated
  output subsequently passed the runtime content and collection-alignment audits.
- Local preparation tests: two weight tests and 27 production-backend tests passed.

## Sample and results

Each process used ten independent chunks, each with 112 generated records and
100 requested converted records: **8,960 generated and 8,000 converted** in total.
All eight catalog processes were included, including diagnostic Z_ll. Different
physics processes are not combined into a single prediction by this test.

| Process | Generated | Converted | Converted sumw | Converted sumw² |
|---|---:|---:|---:|---:|
| Z_qq | 1120 | 1000 | 1000 | 1000 |
| Z_ll | 1120 | 1000 | 1000 | 1000 |
| Z_ee | 1120 | 1000 | 1000 | 1000 |
| Z_mumu | 1120 | 1000 | 1000 | 1000 |
| Z_tautau | 1120 | 1000 | 1000 | 1000 |
| yy_qq | 1120 | 1000 | 1000 | 1000 |
| yy_ll | 1120 | 1000 | 1000 | 1000 |
| Bhabha (BHWIDE) | 1120 | 1000 | 873.977468617526 | 1619.268192693067 |

All PYTHIA nominal weights were 1. BHWIDE converted weights ranged from
9.76667132108803e-14 to 3.3661326600064947. Its generated sumw was
985.6493342668147, distinct from the converted sum. None of the converted events
in this small test had negative weights; negative-weight support was exercised
by synthetic sampler/bookkeeping tests, not by a negative event in this detector
sample. Do not claim this run specifically validates negative-event detector transport.

## Verified directly from the ROOT output

- Every event's double-precision nominal weight and original trial ID matched
  its CSV record exactly under run/event identity matching.
- All 80 files passed independent read-back of cross sections/errors in pb,
  embedded task/generator JSON, generated IDs and original trial-ID ledger.
- Generated, converted, missing and tail count/sumw/sumw²/sumabsw/negative-count
  metadata matched independent sums from the source records.
- Converted IDs and the missing/tail partition matched the output event lists.
- All event identities/counts, FADGEN truth-order/PDG checks and collection/relation
  alignment checks passed. Truth momenta matched FADGEN exactly in these files.
- The first independent-checker attempt failed on a singleton integer-vector
  parameter: podio's convenience getter unwraps it to an integer. The checker
  was fixed locally, committed, pushed and pulled, then all 80 metadata frames
  passed. This was a checker API mismatch, not lost weights.

## Event losses and existing physics limitations

Z_ee chunks 000024 and 000026 skipped detector events 66 and 82, respectively
(stored runs -118716 and -118718). Both match exact ZTELUS out-of-memory log
messages under the owner's accepted policy. Their identities and generator
weights remain in the embedded ledger. Each chunk still produced 100 events
using the predeclared reserve; losses were not silently relabelled or corrected.

The known zero-momentum charged-object issue persists: 580 objects across all
files, **none with selection flag zero**. Counts by process are 40, 91, 158, 8,
79, 1, 26 and 177 in the table's order. These are missing/unmatched measurements,
not physical zero-momentum tracks. Truth relative-vertex differences also persist,
up to 17.42 mm overall and 13.98 mm for Z_tautau. Weight changes do not repair or
certify these detector/truth mappings. Small per-chunk cross-section estimates
are noisy and do not replace the higher-statistics campaign estimates.

## Artifacts and reproduction

Remote campaign:
`/pscratch/sd/a/avencast/Ztautau/validation/weights-all-20260923-v1`

Build:
`/pscratch/sd/a/avencast/Ztautau/validation/weights-all-20260923-v1-build`

Per-chunk products: `work/<six-digit-id>/events.edm4hep.root`, `simana.ldst`,
generator files, input weight CSV/JSON, and `audit/root-audit.json`.
The frozen source tree, plan, card hashes, build manifest, stage logs/resources
and Slurm identity remain on scratch. No old samples were modified and no CFS
publication or cleanup was performed. The default production build path was
not changed; future production must explicitly select the verified new build
or prepare a compatible one.

To repeat, use the committed simulation `production/weight_validation.py prepare`
with a new output/tag, its `production/weight-validation-config`, the converter
checkout and existing registry controller. Preparation reserves fresh seeds/runs.
Then execute its `run` action inside an explicitly approved interactive allocation;
it rebuilds, runs the full chain and writes `weight-validation.json`. Run
`scripts/production/check_weight_metadata.py CAMPAIGN` inside the pinned key4hep
container on compute for independent metadata verification.

Small local evidence: [process summary](weight-validation.json),
[metadata read-back](metadata-readback.json), [loss review](loss-review.json),
[allocation](allocation.json), [source provenance](test-provenance.json).
