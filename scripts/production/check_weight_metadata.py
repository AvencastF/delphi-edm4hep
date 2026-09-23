"""Independent metadata read-back for a completed weight-validation campaign."""
import csv
import json
import math
from pathlib import Path
import sys
from podio import root_io


def check(work):
    generation=json.loads((work/'generation.json').read_text())
    task=json.loads((work/'task.json').read_text())
    with (work/'event-weights.csv').open() as stream:
        records={int(r['hepmc_event']):r for r in csv.DictReader(stream)}
    reader=root_io.Reader(str(work/'events.edm4hep.root'))
    frames=list(reader.get('metadata'))
    if len(frames)!=1: raise ValueError('Expected one metadata frame')
    meta=frames[0]
    def get(key):return meta.get_parameter(key)
    def equal(key,expected):
        actual=get(key)
        if isinstance(expected,float):
            if not math.isclose(float(actual),expected,rel_tol=1e-12,abs_tol=1e-12): raise ValueError(key)
        elif actual!=expected: raise ValueError(key)
    equal('mc_weight_schema',1)
    equal('mc_run_number',-task['run'])
    equal('mc_cross_section_pb',generation['sigma_mb']*1e9)
    equal('mc_cross_section_error_pb',generation['sigma_error_mb']*1e9)
    provenance=json.loads(str(get('mc_production_json')))
    if provenance['task']!=task or provenance['generation']!=generation: raise ValueError('Embedded provenance mismatch')
    generated=list(get('mc_generated_event_ids'))
    if generated!=sorted(records): raise ValueError('Generated IDs mismatch')
    if list(get('mc_generated_trial_ids'))!=[int(records[e].get('generator_trial',-1)) for e in generated]: raise ValueError('Trial ledger mismatch')
    converted=list(get('mc_converted_event_ids'));last=max(converted)
    missing=[e for e in generated if e<=last and e not in converted]
    tail=[e for e in generated if e>last]
    if list(get('mc_missing_ids_before_last_output'))!=missing or list(get('mc_ids_after_last_output'))!=tail: raise ValueError('Loss/tail partition mismatch')
    for prefix,ids in [('generated',generated),('converted',converted),('missing',missing),('after_last_output',tail)]:
        ws=[float(records[e]['weight']) for e in ids]
        for field,value in [('count',len(ids)),('sumw',math.fsum(ws)),('sumw2',math.fsum(w*w for w in ws)),
                            ('sumabsw',math.fsum(abs(w) for w in ws)),('negative_count',sum(w<0 for w in ws))]:
            equal('mc_'+prefix+'_'+field,value)
    return dict(task=task['id'],process=task['process'],converted=len(converted),passed=True)


if __name__=='__main__':
    root=Path(sys.argv[1]);plan=json.loads((root/'plan.json').read_text())
    rows=[check(root/'work'/task['id']) for task in plan['tasks']]
    report=dict(passed=all(row['passed'] for row in rows),files=len(rows),results=rows)
    (root/'metadata-readback.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],files=len(rows))))
