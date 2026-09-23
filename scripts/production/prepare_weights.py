"""Validate a production chunk and prepare its converter weight input.

The line protocol is internal: magic, numeric header, one JSON provenance line,
then (detector event, original generator trial or -1, nominal weight) rows.
"""
import csv
import hashlib
import json
import math
from pathlib import Path


def prepare(directory):
    directory = Path(directory)
    task = json.loads((directory / 'task.json').read_text())
    generation = json.loads((directory / 'generation.json').read_text())
    path = directory / 'event-weights.csv'
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    records = []
    for row in rows:
        event = int(row['hepmc_event'])
        if 'fadgen_record' in row and int(row['fadgen_record']) != event:
            raise ValueError('FADGEN/HepMC identity mismatch')
        weight = float(row['weight'])
        trial = int(row.get('generator_trial', -1))
        if not math.isfinite(weight) or (trial != -1 and trial < 1):
            raise ValueError('Invalid event weight/trial')
        records.append((event, trial, weight))
    records.sort()
    count = generation['generated']
    if count != task['generated_requested'] or [r[0] for r in records] != list(range(1, count + 1)):
        raise ValueError('Incomplete or duplicated generator identities')
    trials = [r[1] for r in records if r[1] != -1]
    if len(set(trials)) != len(trials):
        raise ValueError('Duplicate generator trial')
    for key, value in [('sumw', math.fsum(r[2] for r in records)),
                       ('sumw2', math.fsum(r[2] ** 2 for r in records))]:
        if not math.isfinite(generation[key]) or not math.isclose(value, generation[key], rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError('Generator summary mismatch: ' + key)
    sigma, error = generation['sigma_mb'] * 1e9, generation['sigma_error_mb'] * 1e9
    if not math.isfinite(sigma) or sigma <= 0 or not math.isfinite(error) or error < 0:
        raise ValueError('Invalid cross section')
    if not 0 < task['simulated_requested'] <= count or not 0 < task['run'] <= 999999:
        raise ValueError('Invalid detector run/count')
    provenance = dict(schema=1, task=task, generation=generation,
                      weight_csv_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      cross_section_scope='Generator estimate for the configured process and cuts; see task.card and generation.',
                      normalization='Pool compatible chunks before analysis cuts. Generated and converted sums are distinct; account for the simulation subset and losses explicitly. No luminosity, filter-efficiency, loss correction or theory correction has been applied.',
                      variation_weights='Only nominal supplied by current generators; no PDF/scale/spin variation weights are implied.')
    lines = ['DELPHI_MC_WEIGHTS_V1', '{} {} {:.17g} {:.17g}'.format(-task['run'], count, sigma, error),
             json.dumps(provenance, separators=(',', ':'), allow_nan=False)]
    lines.extend('{} {} {:.17g}'.format(*r) for r in records)
    (directory / 'mc-weights.txt').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    prepare(Path.cwd())
