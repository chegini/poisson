#!/usr/bin/env python3
"""Stationary Poisson solver/compression/routing validation; no timing-speedup claim."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import time
import xml.etree.ElementTree as ET

from run_poisson_quantization_diagnostic import compare, read_solution


def require(condition, message):
    if not condition:
        raise ValueError(message)


def arguments(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--executable', type=Path, help='Optional separately built executable')
    p.add_argument('--kaskade', type=Path, default=Path(os.environ.get('KASKADE7', '/kaskade7')))
    p.add_argument('--vtk', type=int, choices=[0,1], default=0, help='Export the converged global solution.vtu for ParaView')
    p.add_argument('--dimension', type=int, choices=[2,3], default=2)
    p.add_argument('--subres', type=int, default=4)
    p.add_argument('--elmres', type=int, default=8)
    p.add_argument('--sigma', type=float, default=1.)
    p.add_argument('--ranks', type=int, default=2)
    p.add_argument('--expected-nodes', type=int, default=1)
    p.add_argument('--bits', type=int, nargs='+', default=[16])
    p.add_argument('--solvers', nargs='+', choices=['cg','fcg','fgmres'], default=['cg','fcg','fgmres'])
    p.add_argument('--split', choices=['off','on','both'], default='off')
    p.add_argument('--transform-orders', nargs='+', choices=['before','after','both','quantization-only'], default=['before'],
                   help='quantization-only disables lifting and entropy coding; splitting is controlled separately')
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--tolerance', type=float, default=1e-7)
    p.add_argument('--fgmres-restart', type=int, default=30)
    p.add_argument('--check-routing', action=argparse.BooleanOptionalAction, default=False)
    p.add_argument('--check-local-bypass', action=argparse.BooleanOptionalAction, default=False)
    p.add_argument('--launcher', default=os.environ.get('MPIEXEC','mpirun'))
    p.add_argument('--mpi-args', default='')
    p.add_argument('--slurm-container', type=Path)
    p.add_argument('--slurm-btl', choices=['self,vader,tcp','self,tcp','self,vader'], default='self,vader,tcp')
    p.add_argument('--slurm-vader-single-copy', choices=['default','none'], default='none')
    p.add_argument('--timeout', type=int, default=600, help='Bound each solve in seconds')
    p.add_argument('--output', type=Path)
    p.add_argument('--plan-only', action='store_true')
    a = p.parse_args(argv)
    if min(a.subres,a.elmres,a.ranks,a.expected_nodes,a.iterations,a.fgmres_restart,a.timeout)<1:
        p.error('Sizes, ranks, nodes and limits must be positive')
    if a.expected_nodes>a.ranks or a.ranks%a.expected_nodes:
        p.error('Require equal ranks per node and ranks >= nodes')
    if not all(math.isfinite(x) and x>0 for x in (a.sigma,a.tolerance)):
        p.error('Sigma and tolerance must be finite and positive')
    if any(b<1 or b>32 for b in a.bits) or any(len(x)!=len(set(x)) for x in (a.bits,a.solvers,a.transform_orders)):
        p.error('Use distinct choices and bits in 1..32')
    if a.expected_nodes>1 and a.slurm_btl=='self,vader':
        p.error('self,vader cannot communicate between nodes')
    a.repo=a.repo.resolve(); a.kaskade=a.kaskade.resolve()
    a.executable=(a.executable or a.repo/'bddc').resolve()
    if a.slurm_container: a.slurm_container=a.slurm_container.resolve()
    return a


def plan(a):
    modes=[dict(name='direct',solver='cg',bits=None,split=False,order=None,route='direct')]
    routes=['legacy','destination'] if a.check_routing else ['destination']
    splits=[False,True] if a.split=='both' else [a.split=='on']
    for solver in a.solvers:
        for route in routes:
            modes.append(dict(name=f'{solver}-uncompressed-{route}',solver=solver,bits=None,split=False,order=None,route=route))
        for bits in a.bits:
            for order in a.transform_orders:
                for split in splits:
                    for route in routes+(['local-bypass'] if a.check_local_bypass else []):
                        name=f'{solver}-{order}-{bits}'+('-split' if split else '')+'-'+route
                        modes.append(dict(name=name,solver=solver,bits=bits,split=split,order=order,route=route))
    return modes


CONTAINER_ENV='''set -eu
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export OMPI_MCA_pml=ob1
exec "$@"
'''


def command(a, mode, folder):
    direct=mode['route']=='direct'; solver=mode['solver']; bits=mode['bits']
    before=mode['order'] in ('before','both'); after=mode['order'] in ('after','both')
    entropy = mode['order'] != 'quantization-only'
    payload=[str(a.executable),'--mpi','1','--threads','1','--vtk',str(a.vtk),'--histogram','0','--timing','0',
        '--dimension',str(a.dimension),'--subres',str(a.subres),'--elmres',str(a.elmres),'--sigma',str(a.sigma),
        '--crisscross','1','--interfacetypes','7','--iter',str(a.iterations),
        '--relativeTolerance',str(a.tolerance),'--verbose','0','--direct',str(int(direct)),
        '--useCg','1','--originalResidual','1','--flexibleCg',str(int(solver=='fcg' and not direct)),
        '--fgmres',str(int(solver=='fgmres' and not direct)), '--fgmresRestart',str(a.fgmres_restart),
        '--compression',str(int(bits is not None)), '--compressionBits',str(bits or 32),
        '--constantSplit',str(int(mode['split'])), '--restrictionPasses','1',
        '--graphLifting',str(int(before)), '--graphLiftingAfterQuantization',str(int(after)),
        '--graphLiftingBothOrders',str(int(mode['order']=='both')),
        '--huffman',str(int(entropy)),'--bitlength',str(int(entropy)),'--zigzag',str(int(entropy)),'--bddcExcludeCornersFromCompression','1',
        '--mpiRouting',str(int(mode['route'] in ('destination','local-bypass'))),
        '--localCompressionBypass',str(int(mode['route']=='local-bypass')),'--prefix',str(folder)]
    ranks,nodes=(1,1) if direct else (a.ranks,a.expected_nodes)
    if not a.slurm_container:
        launcher=shlex.split(a.launcher)+shlex.split(a.mpi_args)
        require(bool(shlex.split(a.launcher)), 'Empty MPI launcher')
        return launcher+['-np',str(ranks)]+payload
    env=CONTAINER_ENV.replace('exec "$@"','export OMPI_MCA_btl='+a.slurm_btl+'\nexec "$@"')
    if a.slurm_vader_single_copy=='none':
        env=env.replace('exec "$@"','export OMPI_MCA_btl_vader_single_copy_mechanism=none\nexec "$@"')
    return ['srun','--mpi=pmix_v5','--nodes='+str(nodes),'--ntasks='+str(ranks),
        '--ntasks-per-node='+str(ranks//nodes),'--cpus-per-task=1','--distribution=block:block',
        '--kill-on-bad-exit=1','--export=ALL','--time='+str(max(1,math.ceil(a.timeout/60))),
        'apptainer','exec','--bind',f'{a.repo.parent}:{a.repo.parent}',
        '--bind',f'{a.kaskade}:/kaskade7','--pwd',str(a.repo),str(a.slurm_container),
        '/bin/bash','--noprofile','--norc','-c',env,'poisson-container']+payload


def routing_counts(text, nodes, ranks=None):
    require(re.findall(r'^BDDC MPI routing: destination; shared-memory domains=(\d+);',text,re.M)==[str(nodes)],
            'Missing/wrong node topology')
    leaders=re.findall(r'^BDDC MPI routing: destination; shared-memory domains=\d+; rank leaders=([^;]+);',text,re.M)
    require(len(leaders)==1 and (ranks is None or len(leaders[0].split())==ranks),
            'Wrong actual MPI rank count')
    rows=re.findall(r'^BDDC routed payload: phase=(restrict|prolongate) route=(same-rank|same-node|inter-node) bytes=(\d+) packets=(\d+)$',text,re.M)
    require(len(rows)==6 and len({(p,r) for p,r,_,_ in rows})==6,'Missing/duplicate routing counters')
    result={r:sum(int(b) for _,route,b,_ in rows if route==r) for r in ('same-rank','same-node','inter-node')}
    if nodes==1: require(result['inter-node']==0,'Single-node run reports remote bytes')
    else: require(all(int(b)>0 and int(n)>0 for _,r,b,n in rows if r=='inter-node'),
                  'Both inter-node directions must be exercised; use more subdomains')
    return result


def parse_run(text, folder, mode, a):
    matches=re.findall(r'^Poisson result: status=(\S+) solver=(\S+) calls=(\d+)(?: accepted_updates=\d+)? rel=(\S+)$',text,re.M)
    require(len(matches)==1,'Missing/duplicate final result')
    status,solver,calls,rel=matches[0]; rel=float(rel)
    require(status=='converged' and math.isfinite(rel) and 0<=rel<=a.tolerance,'True residual/completion check failed')
    require(solver==('direct' if mode['route']=='direct' else mode['solver']),'Wrong solver')
    require(re.findall(r'^Poisson geometry: dimension=(\d+) subdomains=(\d+)$',text,re.M)==
            [(str(a.dimension),str(a.subres**a.dimension))], 'Wrong/missing spatial dimension or subdomain count')
    row=dict(iterations=int(calls),relative_residual=rel)
    solution=read_solution(folder/'solution.csv')
    # Symmetric cuboid grids: corner nodes + cell centers (+ face centers in 3D).
    axis=a.subres*a.elmres
    expected=(axis+1)**a.dimension+axis**a.dimension
    if a.dimension==3: expected+=3*axis*axis*(axis+1)
    require(len(solution)==expected,'Unexpected global grid DOF count')
    if a.vtk:
        file=folder/'solution.vtu'
        require(re.findall(r'^Poisson VTK file: (.*)$',text,re.M)==[str(file)], 'Missing VTK export marker')
        document=ET.parse(file);piece=document.find('./UnstructuredGrid/Piece')
        require(piece is not None and int(piece.get('NumberOfPoints'))==len(solution), 'Invalid VTK point count')
        u=piece.find('./PointData/DataArray[@Name="u"]')
        ids=piece.find('./PointData/DataArray[@Name="global_dof"]')
        require(u is not None and ids is not None,'Missing VTK solution/DOF arrays')
        require(list(map(float,(u.text or '').split()))==solution,'VTK and CSV solutions differ')
        require(list(map(int,(ids.text or '').split()))==list(range(len(solution))), 'VTK global DOF numbering differs')
        require(int(piece.get('NumberOfCells'))==(24 if a.dimension==3 else 4)*axis**a.dimension, 'Unexpected VTK cell count')
        row['vtk_file']=str(file)
    if mode['route']=='direct': return row,solution
    with (folder/'convergence.csv').open() as stream: trace=list(csv.DictReader(stream))
    require(bool(trace) and int(trace[-1]['solver_calls'])==int(calls),'Trace iteration mismatch')
    require(float(trace[-1]['relative_residual'])==rel,'Trace residual mismatch')
    require(all(math.isfinite(float(t['interface_jump'])) and 0<=float(t['interface_jump'])<=a.tolerance for t in trace),
            'Shared-copy continuity failed')
    for label,key,convert in [('raw sent bytes','raw_bytes',int),('transmitted sent bytes','bytes',int),
                              ('Poisson solver wall seconds','solver_seconds',float)]:
        values=re.findall(r'^\s*'+label+r':\s*(\S+)',text,re.M)
        require(len(values)==1,'Missing/duplicate '+label)
        row[key]=convert(values[0]); require(math.isfinite(row[key]) and row[key]>=0,'Invalid '+label)
    outer=re.findall(r'^Poisson outer solver logical payload bytes: (\d+) collectives=(\d+)$',text,re.M)
    require(len(outer)==1,'Missing outer solver traffic')
    row['outer_solver_bytes']=int(outer[0][0])
    expected='uncompressed' if mode['bits'] is None else {'before':'transform-first','after':'quantize-first','both':'both','quantization-only':'quantization-only'}[mode['order']]
    require(re.findall(r'^Poisson compression order: (.*)$',text,re.M)==[expected],'Wrong transform order')
    if mode['route']=='local-bypass':
        require(re.findall(r'^BDDC local compression bypass: on; mixed-owner DOFs exact on all copies; nodes=(\d+)$',text,re.M)==[str(a.expected_nodes)],
                'Missing/wrong local bypass configuration')
    if mode['route'] in ('destination','local-bypass'):
        row['route_bytes']=routing_counts(text,a.expected_nodes,a.ranks)
        require(sum(row['route_bytes'].values())==row['bytes'],'Routing and logical sent bytes disagree')
    return row,solution


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def parity(first,second):
    # Timings are deliberately excluded from numerical parity.
    return all((first/name).read_bytes()==(second/name).read_bytes()
               for name in ('solution.csv','convergence.csv'))


def main(argv=None):
    a=arguments(argv); modes=plan(a)
    print('Planned runs: '+', '.join(m['name'] for m in modes),flush=True)
    if a.plan_only: return 0
    if a.slurm_container:
        require(bool(os.environ.get('SLURM_JOB_ID')),'Run inside an active allocation')
        require(int(os.environ.get('SLURM_JOB_NUM_NODES','0'))>=a.expected_nodes,'Allocation has too few nodes')
        require(a.slurm_container.is_file(),'Container image missing')
    require(a.executable.is_file(),'Build the Poisson executable first')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=(a.output or a.repo/'output'/f'poisson-solver-check-{stamp}-{os.getpid()}').resolve()
    root.mkdir(parents=True,exist_ok=False)
    sources=[a.repo/name for name in ('bddc.cpp','poisson_bddc_options.hh','poisson_validation.hh','poisson_fgmres.hh','poisson_fgmres_mpi.hh','poisson_geometry.hh','poisson_vtk.hh')]
    sources += [Path(__file__).resolve()]+[a.kaskade/'mg'/name for name in ('bddc.hh','bddc.hpp','bddcMpi.hh','bddcMpiRouting.hh','bddcSpaceTransfer.hh','bddcSpaceTransfer.hpp')]
    fingerprint=digest(a.executable)
    manifest=dict(vars(a),executable_sha256=fingerprint,source_sha256={str(p):digest(p) for p in sources if p.is_file()},
        timing_scope='Solver includes setup/validation/CSV; application is launcher wall time. Not a performance benchmark.',
        byte_scope='Logical R/P payloads; excludes framing/setup/coarse/codebooks/validation; outer FGMRES traffic reported separately.',
        slurm={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_JOB_NODELIST')})
    (root/'environment.json').write_text(json.dumps(manifest,default=str,indent=2)+'\n')
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    for k in ('KASKADE_FROZEN_TRANSFER_DIR','KASKADE_RESTRICTION_CAPTURE','KASKADE_TRANSFER_AUDIT_DIR'):
        env.pop(k,None)
    rows=[]; solutions={}; completed={}
    def save():
        (root/'summary.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
        lines=['| Case | Status | Iterations | Relative error/direct | Bytes saved % | Solver s |',
               '|---|---|---:|---:|---:|---:|']
        for row in rows:
            def number(key):
                return format(row[key],'.6g') if key in row else '--'
            lines.append('| '+ ' | '.join([row['case'],row['status'],number('iterations'),
                number('l2_vs_direct'),number('net_saved_percent'),number('solver_seconds')])+' |')
        lines += ['', 'Logical R/P payloads only; outer FGMRES communication is separate.',
                  'Errors have no scientific acceptance threshold. Validation times include setup and output.']
        (root/'summary.md').write_text('\n'.join(lines)+'\n')
        columns=('case','status','iterations','relative_residual','l2_vs_direct','l2_vs_uncompressed','bytes',
                 'net_saved_percent','inter_node_saved_percent','solver_seconds','application_seconds','outer_solver_bytes')
        with (root/'summary.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print('Output: '+str(root),flush=True)
    for mode in modes:
        folder=root/mode['name'];folder.mkdir()
        cmd=command(a,mode,folder);(folder/'command.sh').write_text(shlex.join(cmd)+'\n')
        row=dict(case=mode['name'],configuration=mode,status='FAIL');rows.append(row)
        print('Running '+mode['name']+' ...',flush=True)
        try:
            require(digest(a.executable)==fingerprint,'Executable changed during campaign')
            start=time.monotonic()
            with (folder/'run.log').open('w') as f:
                proc=subprocess.run(cmd,cwd=a.repo,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=a.timeout)
            row['application_seconds']=time.monotonic()-start
            log=(folder/'run.log').read_text(errors='replace')
            final=re.findall(r'^Poisson result: status=(\S+) solver=(\S+) calls=(\d+)(?: accepted_updates=\d+)? rel=(\S+)$',log,re.M)
            if len(final)==1:
                status,_,calls,relative=final[0]
                row['solver_status']=status; row['iterations']=int(calls)
                if math.isfinite(float(relative)): row['relative_residual']=float(relative)
            row['exit_code']=proc.returncode
            require(proc.returncode==0,'MPI/application exit code '+str(proc.returncode))
            metrics,x=parse_run(log,folder,mode,a);row.update(metrics)
            if mode['route']!='direct':
                row['error_vs_direct']=compare(x,solutions['direct']);row['l2_vs_direct']=row['error_vs_direct']['relative_l2']
                baseline=f"{mode['solver']}-uncompressed-destination"
                if mode['bits'] is not None:
                    require(baseline in completed,'Matched uncompressed reference failed')
                    row['error_vs_uncompressed']=compare(x,solutions[baseline]);row['l2_vs_uncompressed']=row['error_vs_uncompressed']['relative_l2']
                    base=completed[baseline]
                    if base['bytes']: row['net_saved_percent']=100*(1-row['bytes']/base['bytes'])
                    if 'route_bytes' in row and base['route_bytes']['inter-node']:
                        row['inter_node_saved_percent']=100*(1-row['route_bytes']['inter-node']/base['route_bytes']['inter-node'])
                else:
                    row['net_saved_percent']=0.
                if a.check_routing and mode['route']=='destination':
                    peer=mode['name'].removesuffix('destination')+'legacy'
                    require(peer in completed and parity(root/peer,folder),'Legacy/destination solution or trace mismatch')
                    require(row['bytes']==completed[peer]['bytes'],'Routing changed logical bytes')
                    row['routing_parity']='EXACT_SOLUTION_AND_TRACE'
                if mode['route']=='local-bypass' and a.expected_nodes==1:
                    require(parity(root/baseline,folder),'Single-node bypass differs from matched uncompressed solution/trace')
                    row['bypass_parity']='EXACT_SOLUTION_AND_TRACE'
                if mode['order']=='both':
                    peer=mode['name'].replace('-both-','-before-',1)
                    if peer in completed:
                        require(parity(root/peer,folder),'Combined/before solution or trace mismatch')
                        row['combined_parity']='EXACT_SOLUTION_AND_TRACE'
            row['status']='PASS';solutions[mode['name']]=x;completed[mode['name']]=row
        except (ValueError,OSError,KeyError,ET.ParseError,subprocess.TimeoutExpired) as e:
            row['reason']=str(e)
        save();print('  '+row['status']+(': '+row['reason'] if row['status']=='FAIL' else ''),flush=True)
        if row['status']=='FAIL' and (mode['route']=='direct' or mode['bits'] is None):
            print('Stopped: reference or routing baseline failed.');break
    print((root/'summary.md').read_text())
    print('Results: '+str(root))
    return int(any(r['status']!='PASS' for r in rows))


if __name__=='__main__':
    raise SystemExit(main())
