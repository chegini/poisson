#!/usr/bin/env python3
"""Fixed 2-D Poisson direct/CG comparison, with quantization and optional full coding."""
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
import sys
import time


def read_solution(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    ids = [int(row['global_dof']) for row in rows]
    values = [float(row['value']) for row in rows]
    if not rows or ids != list(range(len(rows))) or not all(map(math.isfinite, values)):
        raise ValueError('Invalid/nonfinite global solution: '+str(path))
    return values


def compare(first, reference):
    if len(first) != len(reference):
        raise ValueError('Global solution sizes differ')
    error = [a-b for a, b in zip(first, reference)]
    absolute = math.hypot(*error)
    norm = math.hypot(*reference)
    result = {'relative_l2': absolute/norm if norm else absolute,
              'max_absolute': max(map(abs, error))}
    if not all(map(math.isfinite,result.values())):
        raise ValueError('Nonfinite solution comparison')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--bits', type=int, nargs='+', default=[8,12,16,24,32])
    parser.add_argument('--restriction-passes', type=int, default=1)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--subres', type=int, default=4)
    parser.add_argument('--elmres', type=int, default=8)
    parser.add_argument('--sigma', type=float, default=1.0)
    parser.add_argument('--tolerance', type=float, default=1e-7)
    parser.add_argument('--iterations', type=int, default=300)
    parser.add_argument('--use-cg', type=int, choices=[0,1], default=1)
    parser.add_argument('--full-compression', action='store_true')
    parser.add_argument('--launcher', default=os.environ.get('MPIEXEC','mpirun'))
    parser.add_argument('--mpi-args', default='')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not args.bits or len(set(args.bits))!=len(args.bits) or any(b<1 or b>32 for b in args.bits):
        parser.error('--bits must be distinct integers in 1..32; 48/64 are unsupported')
    if not 1<=args.restriction_passes<=16:
        parser.error('--restriction-passes must be in 1..16')
    if min(args.ranks,args.subres,args.elmres,args.iterations)<1:
        parser.error('ranks, subres, elmres and iterations must be positive')
    if not all(math.isfinite(v) and v>0 for v in [args.sigma,args.tolerance]):
        parser.error('sigma and tolerance must be finite and positive')
    launcher = shlex.split(args.launcher)+shlex.split(args.mpi_args)
    if not shlex.split(args.launcher):
        parser.error('MPI launcher is empty')
    repo = args.repo.resolve()
    executable = repo/'bddc'
    if not executable.is_file() or not os.access(executable,os.X_OK):
        parser.error('Build the updated MPI executable first: '+str(executable))
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root = args.output or Path('output')/f'poisson-quantization-diagnostic-{stamp}-{os.getpid()}'
    if not root.is_absolute(): root = repo/root
    root.mkdir(parents=True,exist_ok=False)
    digest = hashlib.sha256()
    with executable.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024),b''): digest.update(chunk)
    metadata = vars(args).copy()
    metadata.update(executable_sha256=digest.hexdigest(),direct_ranks=1,
                    byte_scope='logical transfer payloads; excludes diagnostic collectives and MPI/codebook overhead',
                    solution_comparison='raw global coefficients; no mean offset removed',
                    timings='validation runs include true-residual checks and CSV output; not pure communication timings')
    (root/'environment.json').write_text(json.dumps(metadata,default=str,indent=2)+'\n')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    common = [str(executable),'--mpi','1','--threads','1','--vtk','0','--timing','0',
              '--subres',str(args.subres),'--elmres',str(args.elmres),'--sigma',str(args.sigma),
              '--crisscross','1','--interfacetypes','7','--iter',str(args.iterations),
              '--relativeTolerance',str(args.tolerance),'--useCg',str(args.use_cg),
              '--verbose','0','--histogram','0','--bddcExcludeCornersFromCompression','1']
    records = []
    references = {}

    def save():
        (root/'summary.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
        lines = ['case                 status passes calls      true rel  L2 vs direct       bytes   net saved   solver s']
        for row in records:
            def number(key, fmt):
                value = row.get(key)
                return format(value,fmt) if value is not None else '--'
            lines.append(f"{row['case']:20s} {row['status']:6s} {str(row.get('passes','--')):>6s} "
                         f"{str(row.get('calls','--')):>5s} {number('relative_residual','.3e'):>13s} "
                         f"{number('l2_vs_direct','.3e'):>13s} {str(row.get('bytes','--')):>11s} "
                         f"{number('net_saved_percent','.2f'):>10s} {number('solver_seconds','.3f'):>10s}")
        lines += ['', 'Same generated mesh, sigma and residual tolerance in every case; direct uses one rank.',
                  'Both transfer directions use the listed bit depth. Extra restriction passes keep the solution fixed.',
                  'PASS checks the true relative residual, shared-copy continuity and finite solution comparisons.',
                  'Raw solution differences have no hidden acceptance threshold and require scientific review.',
                  'Net savings include all iterations/passes, relative to uncompressed BDDC.',
                  'Byte counts include intra-rank transfers; exclude validation reductions, MPI framing and codebook training.',
                  'Validation timings include residual checks and CSV output; repeat separately for performance claims.']
        (root/'summary.txt').write_text('\n'.join(lines)+'\n')

    def run(name, bits=None, direct=False, full=False):
        folder = root/name
        folder.mkdir()
        passes = args.restriction_passes if bits is not None else 1
        coding = '1' if full else '0'
        command = launcher+['-np',str(1 if direct else args.ranks)]+common+[
            '--direct',str(int(direct)),'--compression',str(int(bits is not None)),
            '--compressionBits',str(bits if bits is not None else 32),
            '--restrictCompressionBits',str(bits if bits is not None else 32),
            '--prolongateCompressionBits',str(bits if bits is not None else 32),
            '--restrictionPasses',str(passes),'--graphLifting','0',
            '--graphLiftingAfterQuantization',coding,'--huffman',coding,
            '--bitlength',coding,'--zigzag',coding,'--prefix',str(folder)]
        (folder/'command.sh').write_text(shlex.join(command)+'\n')
        print('Running '+name+' ...',flush=True)
        started = time.monotonic()
        with (folder/'run.log').open('w') as stream:
            result = subprocess.run(command,cwd=repo,env=env,stdout=stream,stderr=subprocess.STDOUT)
        row = dict(case=name,status='FAIL',bits=bits,passes=passes,return_code=result.returncode,
                   elapsed_seconds=time.monotonic()-started,log=str(folder/'run.log'))
        text = (folder/'run.log').read_text(errors='replace')
        try:
            matches = re.findall(r'Poisson result: status=(\S+) solver=(\S+) calls=(\d+)(?: accepted_updates=\d+)? rel=(\S+)',text)
            if len(matches)!=1: raise ValueError('Missing/duplicate final result')
            status, solver, calls, relative = matches[0]
            relative = float(relative)
            row['calls'] = int(calls)
            if math.isfinite(relative): row['relative_residual'] = relative
            if not direct:
                with (folder/'convergence.csv').open() as stream: samples = list(csv.DictReader(stream))
                final = samples[-1]
                if int(final['solver_calls'])!=int(calls): raise ValueError('Solver-call count disagrees with history')
                jump = float(final['interface_jump'])
                if not math.isfinite(jump) or jump>args.tolerance: raise ValueError('Shared solution copies disagree')
                if float(final['relative_residual'])!=relative: raise ValueError('Residual disagrees with history')
                for label, key, convert in [('transmitted sent bytes','bytes',int),
                                             ('raw sent bytes','raw_bytes',int),
                                             ('Poisson solver wall seconds','solver_seconds',float)]:
                    values = re.findall(r'(?:^|\n)\s*'+label+r':\s*(\S+)',text)
                    if len(values)==1:
                        value = convert(values[0])
                        if not math.isfinite(value): raise ValueError('Nonfinite counter: '+label)
                        row[key] = value
            if result.returncode or status!='converged' or not math.isfinite(relative) or relative>args.tolerance:
                raise ValueError('Solver did not converge: '+status)
            expected = 'direct' if direct else 'cg' if args.use_cg else 'sd'
            if solver!=expected: raise ValueError('Unexpected solver type')
            solution = read_solution(folder/'solution.csv')
            if direct:
                references['direct'] = solution
            else:
                if not all(k in row for k in ['bytes','raw_bytes','solver_seconds']):
                    raise ValueError('Missing byte/time counters')
                if (min(row['bytes'],row['raw_bytes'])<0 or
                        not math.isfinite(row['solver_seconds']) or row['solver_seconds']<0):
                    raise ValueError('Invalid byte/time counters')
                row['comparison_to_direct'] = compare(solution,references['direct'])
                row['l2_vs_direct'] = row['comparison_to_direct']['relative_l2']
                if bits is None:
                    references['uncompressed'] = solution
                    references['bytes'] = row['bytes']
                else:
                    row['comparison_to_uncompressed'] = compare(solution,references['uncompressed'])
                if references['bytes']:
                    row['net_saved_percent'] = 100*(1-row['bytes']/references['bytes'])
            row['status'] = 'PASS'
        except (ValueError,OSError,KeyError,IndexError) as error:
            row['reason'] = str(error)
        records.append(row)
        (folder/'metrics.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
        save()
        print('  '+row['status']+(': '+row['reason'] if row['status']=='FAIL' else ''),flush=True)
        return row['status']=='PASS'

    print('Diagnostic output: '+str(root),flush=True)
    if not run('direct',direct=True) or not run('uncompressed-cg' if args.use_cg else 'uncompressed-sd'):
        print('Stopped: inspect the reference/baseline logs.'); return 1
    for bits in args.bits:
        run(f'quantization-{bits}',bits=bits)
        if args.full_compression: run(f'full-{bits}',bits=bits,full=True)
    print('\n'+(root/'summary.txt').read_text())
    print('Summary and logs: '+str(root))
    return int(any(row['status']!='PASS' for row in records))


if __name__=='__main__':
    sys.exit(main())
