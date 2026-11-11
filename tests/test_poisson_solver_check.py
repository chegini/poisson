import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import run_poisson_solver_check as check


class SolverCampaign(unittest.TestCase):
    def test_matrix_and_launch_layout(self):
        a=check.arguments(['--bits','8','12','--split','both','--transform-orders','before','after','both',
                           '--check-routing','--check-local-bypass'])
        modes=check.plan(a)
        self.assertEqual(len(modes),115)
        self.assertEqual(len({m['name'] for m in modes}),115)
        for m in modes:
            cmd=check.command(a,m,Path('/tmp/poisson-case'))
            def get(k): return cmd[cmd.index(k)+1]
            self.assertEqual(get('--originalResidual'),'1')
            self.assertEqual(get('--constantSplit'),str(int(m['split'])))
            self.assertEqual(get('--graphLiftingBothOrders'),str(int(m['order']=='both')))
            self.assertEqual(get('--localCompressionBypass'),str(int(m['route']=='local-bypass')))
        a=check.arguments(['--slurm-container','/tmp/test.sif','--ranks','8','--expected-nodes','4'])
        cmd=check.command(a,check.plan(a)[1],Path('/tmp/result'))
        self.assertIn('--nodes=4',cmd);self.assertIn('--ntasks=8',cmd)
        self.assertNotIn('--cleanenv',cmd)
        self.assertIn('OMPI_MCA_btl_vader_single_copy_mechanism=none',cmd[cmd.index('-c')+1])
        direct=check.command(a,check.plan(a)[0],Path('/tmp/result'))
        self.assertIn('--nodes=1',direct);self.assertIn('--ntasks=1',direct)

    def test_full_compression_only_3d_plan(self):
        a=check.arguments(['--dimension','3','--subres','2','--elmres','4','--bits','12',
            '--transform-orders','before','after','both','--check-routing','--check-local-bypass',
            '--no-check-routing','--no-check-local-bypass'])
        modes=check.plan(a)
        self.assertEqual(len(modes),13)
        self.assertTrue(all(m['route'] in ('direct','destination') for m in modes))
        self.assertEqual(sum(m['bits']==12 for m in modes),9)
        cmd=check.command(a,modes[-1],Path('/tmp/three-dimensional'))
        self.assertEqual(cmd[cmd.index('--dimension')+1],'3')

    def test_quantization_only_fcg_split(self):
        a=check.arguments(['--dimension','2','--subres','16','--elmres','32',
            '--bits',*[str(b) for b in range(4,25)],'--solvers','fcg',
            '--split','on','--transform-orders','quantization-only',
            '--ranks','16','--expected-nodes','8',
            '--no-check-routing','--no-check-local-bypass'])
        modes=check.plan(a)
        self.assertEqual(len(modes),23)
        self.assertEqual(modes[1]['name'],'fcg-uncompressed-destination')
        for mode in modes[2:]:
            self.assertIn('quantization-only',mode['name'])
            cmd=check.command(a,mode,Path('/tmp/quantization-only'))
            def value(flag): return cmd[cmd.index(flag)+1]
            for flag in ['--graphLifting','--graphLiftingAfterQuantization',
                         '--graphLiftingBothOrders','--huffman','--bitlength','--zigzag',
                         '--localCompressionBypass']:
                self.assertEqual(value(flag),'0')
            for flag in ['--compression','--constantSplit','--flexibleCg','--mpiRouting']:
                self.assertEqual(value(flag),'1')
            self.assertEqual(value('--compressionBits'),str(mode['bits']))

    def test_fcg_vtk_plan(self):
        a=check.arguments(['--dimension','3','--subres','2','--elmres','4','--bits','12',
            '--solvers','fcg','--transform-orders','before','--vtk','1',
            '--no-check-routing','--no-check-local-bypass'])
        modes=check.plan(a)
        self.assertEqual([m['name'] for m in modes],
            ['direct','fcg-uncompressed-destination','fcg-before-12-destination'])
        for mode in modes:
            cmd=check.command(a,mode,Path('/tmp/vtk-case'))
            self.assertEqual(cmd[cmd.index('--vtk')+1],'1')
        a=check.arguments([])
        cmd=check.command(a,check.plan(a)[0],Path('/tmp/no-vtk'))
        self.assertEqual(cmd[cmd.index('--vtk')+1],'0')

    def test_topology_and_traffic_rejections(self):
        text='BDDC MPI routing: destination; shared-memory domains=1; rank leaders= 0 0; compression policy=unchanged\n'
        text+='\n'.join(f'BDDC routed payload: phase={p} route={r} bytes={b} packets={int(b>0)}'
            for p in ('restrict','prolongate') for r,b in [('same-rank',8),('same-node',16),('inter-node',0)])
        self.assertEqual(check.routing_counts(text,1,2)['same-node'],32)
        with self.assertRaises(ValueError): check.routing_counts(text,1,4)
        for invalid,nodes in [(text,2),(text+'\n'+text,1),(text.replace('inter-node bytes=0','inter-node bytes=8'),1)]:
            with self.assertRaises(ValueError): check.routing_counts(invalid,nodes)

    def test_full_campaign_and_failed_references(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); exe=root/'bddc';exe.write_text('fixture'); launcher=root/'mpi'
            launcher.write_text('''#!/usr/bin/env python3
import csv,os,sys
from pathlib import Path
a=sys.argv[1:]
def v(k):return a[a.index(k)+1]
out=Path(v('--prefix')); direct=v('--direct')=='1'; compressed=v('--compression')=='1'
bypass=v('--localCompressionBypass')=='1'; routed=v('--mpiRouting')=='1'
solver='direct' if direct else 'fgmres' if v('--fgmres')=='1' else 'fcg' if v('--flexibleCg')=='1' else 'cg'
bad=os.environ.get('BAD_REFERENCE') and not direct and not compressed
rel=0.1 if bad else 1e-9
calls=0 if direct else 4
with (out/'solution.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['global_dof','value'])
 dim=int(v('--dimension'));axis=int(v('--subres'))*int(v('--elmres'))
 count=(axis+1)**dim+axis**dim+(3*axis*axis*(axis+1) if dim==3 else 0)
 w.writerows((i,i+1) for i in range(count))
if not direct:
 with (out/'convergence.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['solver_calls','relative_residual','interface_jump']);w.writerow([calls,rel,0])
 size=50 if compressed and not bypass else 100
 print('  raw sent bytes: 100')
 print('  transmitted sent bytes:',size)
 print('Poisson solver wall seconds: 0.5')
 print('Poisson outer solver logical payload bytes: 0 collectives=0')
 if bypass: print('BDDC local compression bypass: on; mixed-owner DOFs exact on all copies; nodes=1')
 if routed:
  print('BDDC MPI routing: destination; shared-memory domains=1; rank leaders= 0 0; compression policy=unchanged')
  for p in ('restrict','prolongate'):
   for r,b in [('same-rank',0),('same-node',size//2),('inter-node',0)]:
    print(f'BDDC routed payload: phase={p} route={r} bytes={b} packets={int(b>0)}')
order='uncompressed' if not compressed else 'both' if v('--graphLiftingBothOrders')=='1' else 'transform-first' if v('--graphLifting')=='1' else 'quantize-first' if v('--graphLiftingAfterQuantization')=='1' else 'quantization-only'
print('Poisson compression order:',order)
print(f"Poisson geometry: dimension={os.environ.get('BAD_DIMENSION',v('--dimension'))} subdomains={int(v('--subres'))**int(v('--dimension'))}")
print(f'Poisson result: status=converged solver={solver} calls={calls} rel={rel}')
''');launcher.chmod(0o755)
            command=[sys.executable,str(Path(check.__file__)), '--repo',str(root),'--launcher',str(launcher),
                     '--split','both','--transform-orders','before','after','both','--check-routing','--check-local-bypass']
            result=subprocess.run(command+['--output',str(root/'good')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            rows=json.loads((root/'good/summary.json').read_text())
            self.assertEqual(len(rows),61)
            self.assertTrue(all(r['status']=='PASS' for r in rows))
            self.assertTrue(any(r.get('combined_parity') for r in rows))
            self.assertTrue(all(r['net_saved_percent']==50 for r in rows if r['configuration']['bits'] and r['configuration']['route']!='local-bypass'))
            result=subprocess.run(command+['--transform-orders','quantization-only',
                '--solvers','fcg','--split','on','--no-check-routing','--no-check-local-bypass',
                '--output',str(root/'quant')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            quant_rows=json.loads((root/'quant/summary.json').read_text())
            self.assertEqual(len(quant_rows),3)
            self.assertTrue(all(r['status']=='PASS' for r in quant_rows))
            result=subprocess.run(command+['--output',str(root/'bad')],env=dict(os.environ,BAD_REFERENCE='1'),capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            rows=json.loads((root/'bad/summary.json').read_text());self.assertEqual(len(rows),2)
            self.assertEqual(rows[-1]['status'],'FAIL')
            result=subprocess.run(command+['--dimension','3','--subres','2','--elmres','2',
                '--output',str(root/'wrong-dimension')],env=dict(os.environ,BAD_DIMENSION='2'),capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stdout+result.stderr)
            rows=json.loads((root/'wrong-dimension/summary.json').read_text())
            self.assertEqual(len(rows),1)
            self.assertIn('dimension',rows[0]['reason'])


if __name__=='__main__': unittest.main()
