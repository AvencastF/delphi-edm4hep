"""Small routing check; mocks the unavailable DELPHI runtime, not physics."""
import json
import os
from pathlib import Path
import subprocess
import tempfile

wrapper = Path(__file__).resolve().parents[1] / 'convert.sh'
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    mock = '#!/usr/bin/env python3\nimport json,os,sys\nwith open(os.environ["CALLS"],"a") as f: f.write(json.dumps(sys.argv)+"\\n")\nsys.exit(int(os.environ.get("FAIL",0)))\n'
    # Absolute interpreter avoids recursion through the mock python3.
    import sys
    mock = mock.replace('/usr/bin/env python3', sys.executable)
    for name in ('python3', 'converter'):
        p = work / name
        p.write_text(mock)
        p.chmod(0o755)
    calls = work / 'calls'
    env = dict(os.environ, PATH=str(work) + ':' + os.environ['PATH'],
               CALLS=str(calls), DELPHI_CONVERTER_ENV_READY='1')

    def run(args, code=0, **extra):
        calls.write_text('')
        result = subprocess.run(['bash', str(wrapper), *args], cwd=work,
                                env=dict(env, **extra), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
        assert result.returncode == code, result.stderr
        return [json.loads(line) for line in calls.read_text().splitlines()]

    data = ['--data', str(work / 'converter'), 'input.al', 'output.root']
    for limit in ([], ['-n', '2']):
        rows = run(data + limit)
        assert rows[0][1:] == ['input.al', 'output.root', *limit]
        assert len(rows) == 2 and rows[1][1].endswith('/align_audit.py')
    rows = run([str(work / 'converter'), 'input.ldst', 'output.root', 'audit'])
    assert len(rows) == 4
    assert rows[0][1].endswith('/prepare_weights.py')
    assert rows[1][1:] == ['input.ldst', 'output.root', '--mc-weights', 'mc-weights.txt']
    assert rows[2][1].endswith('/audit_root.py')
    assert len(run(data, code=7, FAIL='7')) == 1  # no audit after failure
    assert not run(data + ['--mc-weights', 'bad'], code=2)
    assert not run(data + ['-n', '0'], code=2)
    (work / 'output.root').touch()
    assert not run(data, code=2)
print('PASS: data/MC routing, event limit, failure propagation, overwrite guard')
