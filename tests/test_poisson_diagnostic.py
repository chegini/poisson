import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class PoissonDiagnostic(unittest.TestCase):
    def test_matched_references_failures_precision_and_cost_accounting(self):
        runner = Path(__file__).with_name('run_poisson_quantization_diagnostic.py')
        with tempfile.TemporaryDirectory(prefix='poisson-diagnostic-control-') as temporary:
            repo = Path(temporary)
            model = repo/'bddc'
            model.write_text('#!/bin/sh\nexit 0\n')
            model.chmod(0o755)
            launcher = repo/'fake-mpi'
            launcher.write_text('''#!/usr/bin/env python3
import csv,json,os,sys
from pathlib import Path
a=sys.argv[1:]
def value(key): return a[a.index(key)+1]
with open('capture.jsonl','a') as stream: stream.write(json.dumps(a)+'\\n')
out=Path(value('--prefix'))
direct=value('--direct')=='1'
compressed=value('--compression')=='1'
full=value('--huffman')=='1'
bad=compressed and value('--compressionBits')=='8' and value('--restrictionPasses')=='1'
relative=0.2 if bad else 7e-8
calls=0 if direct else 3 if bad else 10
solver='direct' if direct else 'cg' if value('--useCg')=='1' else 'sd'
if not direct:
 with (out/'convergence.csv').open('w') as stream:
  w=csv.writer(stream)
  w.writerow(['solver_calls','accepted_updates','relative_residual','interface_jump'])
  w.writerow([calls,calls,relative,0])
 print('  raw sent bytes: '+str(1200 if compressed else 1000))
 print('  transmitted sent bytes: '+str(500 if full else 1250 if compressed else 1000))
 print('Poisson solver wall seconds: 1.5')
with (out/'solution.csv').open('w') as stream:
 w=csv.writer(stream); w.writerow(['global_dof','value']); w.writerows([[0,1],[1,2]])
print(f"Poisson result: status={'breakdown' if bad else 'converged'} solver={solver} calls={calls} rel={relative}")
sys.exit(2 if bad or (direct and os.environ.get('FAIL_DIRECT')=='1') else 0)
''')
            launcher.chmod(0o755)
            common = [sys.executable,str(runner),'--repo',str(repo),
                      '--launcher',str(launcher),'--mpi-args=--allow-run-as-root',
                      '--bits','8','16','--full-compression','--ranks','8']
            result = subprocess.run(common+['--output','sweep'],cwd=repo,capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            rows = json.loads((repo/'sweep/summary.json').read_text())
            self.assertEqual([r['status'] for r in rows],['PASS','PASS','FAIL','FAIL','PASS','PASS'])
            calls = [json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()]
            self.assertEqual(calls[0][calls[0].index('-np')+1],'1')
            self.assertEqual(calls[1][calls[1].index('-np')+1],'8')
            for a in calls:
                for flag, expected in [('--subres','4'),('--elmres','8'),('--sigma','1.0'),
                                       ('--relativeTolerance','1e-07'),('--useCg','1'),
                                       ('--graphLifting','0')]:
                    self.assertEqual(a[a.index(flag)+1],expected)
            for a in calls[2:]:
                self.assertEqual(a[a.index('--restrictCompressionBits')+1],a[a.index('--compressionBits')+1])
                self.assertEqual(a[a.index('--prolongateCompressionBits')+1],a[a.index('--compressionBits')+1])
            self.assertAlmostEqual(rows[-1]['net_saved_percent'],50)
            self.assertEqual(rows[-1]['comparison_to_direct']['relative_l2'],0)
            self.assertEqual(rows[2]['calls'],3)
            self.assertNotIn('comparison_to_direct',rows[2])
            refined = subprocess.run(common+['--output','passes2','--restriction-passes','2'],
                                     cwd=repo,capture_output=True,text=True)
            self.assertEqual(refined.returncode,0,refined.stdout+refined.stderr)
            refined_calls = [json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-6:]
            self.assertEqual([a[a.index('--restrictionPasses')+1] for a in refined_calls],
                             ['1','1','2','2','2','2'])
            invalid = subprocess.run(common+['--output','invalid','--bits','48'],cwd=repo,
                                     capture_output=True,text=True)
            self.assertEqual(invalid.returncode,2)
            self.assertFalse((repo/'invalid').exists())
            failed = subprocess.run(common+['--output','failed-reference'],cwd=repo,
                                    env=dict(os.environ,FAIL_DIRECT='1'),capture_output=True,text=True)
            self.assertEqual(failed.returncode,1)
            self.assertFalse((repo/'failed-reference/uncompressed-cg').exists())


if __name__=='__main__':
    unittest.main()
