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
if os.environ.get('KASKADE_FROZEN_TRANSFER_DIR'):
 capture=Path(os.environ['KASKADE_FROZEN_TRANSFER_DIR'])
 (capture/'mock.buffer').write_text('fake MPI input')
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
order='uncompressed' if not compressed else 'transform-first' if value('--graphLifting')=='1' else 'quantize-first' if value('--graphLiftingAfterQuantization')=='1' else 'quantization-only'
print('Poisson compression order: '+order)
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
            selected = [arg for arg in common if arg!='--full-compression']
            paired = subprocess.run(selected+['--bits','16','--transform-orders','before','after',
                                              '--output','orders'],cwd=repo,capture_output=True,text=True)
            self.assertEqual(paired.returncode,0,paired.stdout+paired.stderr)
            order_rows=json.loads((repo/'orders/summary.json').read_text())
            self.assertEqual([row['case'] for row in order_rows],
                             ['direct','uncompressed-cg','quantization-16','before-16','after-16'])
            self.assertEqual(order_rows[-1]['comparison_to_before']['relative_l2'],0)
            order_commands=[json.loads(line) for line in (repo/'capture.jsonl').read_text().splitlines()][-2:]
            for command,before,after in zip(order_commands,['1','0'],['0','1']):
                self.assertEqual(command[command.index('--graphLifting')+1],before)
                self.assertEqual(command[command.index('--graphLiftingAfterQuantization')+1],after)
                self.assertEqual(command[command.index('--huffman')+1],'1')
            # Exercise opt-in capture routing and the observer gate using a fake
            # replay executable. Production codec behavior is tested by C++ fixtures.
            replay=repo/'fake-replay'
            replay.write_text("""#!/usr/bin/env python3
from pathlib import Path
import sys
source, dest=map(Path,sys.argv[1:3])
assert (source/'mock.buffer').read_text()=='fake MPI input'
dest.mkdir()
(dest/'packets.csv').write_text('mode,phase,call,producer,peer,raw_bytes,wire_bytes,grid_bytes,framing_bytes,huffman_bytes,tail_bytes,padding_bytes,exact_bytes,slot_bytes,plan_bytes,squared_error,squared_norm,max_absolute_error\\nquantization,restrict,1,0,1,8,16,8,0,0,0,0,0,8,0,0,1,0\\n')
(dest/'stages.csv').write_text('mode,phase,call,producer,peer,stage,position,value\\nquantization,restrict,1,0,1,symbols,0,1\\n')
(dest/'codebooks.csv').write_text('mode,phase,call,training_symbols,estimated_table_bits\\nquantization,restrict,1,0,0\\n')
""")
            replay.chmod(0o755)
            frozen=subprocess.run(selected+['--bits','16','--transform-orders','before','after',
                '--frozen-replay-executable',str(replay),'--output','frozen'],
                cwd=repo,capture_output=True,text=True)
            self.assertEqual(frozen.returncode,0,frozen.stdout+frozen.stderr)
            frozen_root=repo/'frozen'
            self.assertTrue(json.loads((frozen_root/'capture-16/observer-parity.json').read_text())['passed'])
            for name in ('direct','uncompressed-cg','quantization-16','before-16','after-16'):
                self.assertNotIn('KASKADE_FROZEN_TRANSFER_DIR',
                    json.loads((frozen_root/name/'run-environment.json').read_text()))
            self.assertIn('KASKADE_FROZEN_TRANSFER_DIR',
                json.loads((frozen_root/'capture-16/run-environment.json').read_text()))
            self.assertTrue((frozen_root/'frozen-comparison/summary.json').is_file())
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
