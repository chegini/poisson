"""Summarize production-codec frozen replay; these are sampled logical payloads."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def verify_observer(control, observed):
    """Require matched counts/bytes and tight numerical parity; report exactness."""
    a, b = [json.loads((p/'metrics.json').read_text()) for p in (control, observed)]
    for key in ('status', 'calls', 'bytes', 'raw_bytes'):
        if a[key] != b[key] or (key == 'status' and a[key] != 'PASS'):
            raise ValueError('Capture/control mismatch: '+key)
    exact = True
    for file in ('solution.csv', 'convergence.csv'):
        left, right = rows(control/file), rows(observed/file)
        if len(left) != len(right) or not left:
            raise ValueError('Capture/control row mismatch: '+file)
        for l, r in zip(left, right):
            if l.keys() != r.keys():
                raise ValueError('Capture/control columns differ: '+file)
            for key in l:
                x, y = float(l[key]), float(r[key])
                if not math.isfinite(x) or not math.isfinite(y):
                    raise ValueError('Nonfinite observer comparison')
                exact &= x == y
                if abs(x-y) > 1e-12*max(abs(x), abs(y), 1e-300):
                    raise ValueError('Observer parity exceeds 1e-12 relative: '+file+'/'+key)
    return {'passed': True, 'numerically_identical': exact,
            'relative_tolerance': 1e-12, 'calls_and_bytes_identical': True}


def summarize(root):
    root = Path(root)
    packets = rows(root/'packets.csv')
    if not packets:
        raise ValueError('No frozen packets')
    totals = {}
    sums = ('raw_bytes', 'wire_bytes', 'grid_bytes', 'framing_bytes', 'huffman_bytes',
            'tail_bytes', 'padding_bytes', 'exact_bytes', 'slot_bytes')
    plans = {}
    for p in packets:
        key = p['phase'], p['mode']
        a = totals.setdefault(key, dict.fromkeys(sums, 0))
        for field in sums:
            a[field] += int(p[field])
        a['squared_error'] = a.get('squared_error', 0)+float(p['squared_error'])
        a['squared_norm'] = a.get('squared_norm', 0)+float(p['squared_norm'])
        a['max_absolute_error'] = max(a.get('max_absolute_error', 0),
                                      float(p['max_absolute_error']))
        plan_key = key+(int(p['producer']), int(p['peer']))
        size = int(p['plan_bytes'])
        if plans.setdefault(plan_key, size) != size:
            raise ValueError('Lifting plan size changed during capture')
    for key, a in totals.items():
        a['net_saved_percent'] = 100*(1-a['wire_bytes']/a['raw_bytes']) if a['raw_bytes'] else None
        a['absolute_l2_error'] = math.sqrt(a['squared_error'])
        a['relative_l2_error'] = (math.sqrt(a['squared_error']/a['squared_norm'])
                                  if a['squared_norm'] else None)
        a['plan_setup_bytes_once'] = sum(size for pk, size in plans.items() if pk[:2] == key)
        # Same plans serve R/P; never add these two per-direction copies together.
    for b in rows(root/'codebooks.csv'):
        key = b['phase'], b['mode']
        totals[key]['codebook_bits_sum_over_calls'] = (
            totals[key].get('codebook_bits_sum_over_calls', 0)+int(b['estimated_table_bits']))

    frequencies = defaultdict(Counter)
    reals = defaultdict(list)
    with (root/'stages.csv').open() as stream:
        for row in csv.DictReader(stream):
            key = row['phase'], row['mode'], row['stage']
            value = float(row['value'])
            if not math.isfinite(value):
                raise ValueError('Nonfinite stage value')
            if row['stage'] in ('symbols', 'bitlength', 'quantized', 'lifted'):
                frequencies[key][int(value)] += 1
            else:
                reals[key].append(value)
    with (root/'symbol_histograms.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['phase', 'mode', 'stage', 'symbol', 'count', 'probability'])
        for key, counts in sorted(frequencies.items()):
            n = sum(counts.values())
            for value, count in sorted(counts.items()):
                writer.writerow([*key, value, count, count/n])
            if key[2] in ('symbols', 'bitlength'):
                totals[key[:2]][key[2]+'_entropy_bits'] = -sum(
                    c/n*math.log2(c/n) for c in counts.values())
                totals[key[:2]][key[2]+'_zero_fraction'] = counts[0]/n
                totals[key[:2]][key[2]+'_mean'] = sum(v*c for v, c in counts.items())/n
    # Common bin edges across modes within each phase/stage; real units retained.
    ranges = {}
    for (phase, mode, stage), values in reals.items():
        low, high = ranges.get((phase, stage), (math.inf, -math.inf))
        ranges[phase, stage] = min(low, min(values)), max(high, max(values))
    with (root/'coefficient_histograms.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['phase', 'mode', 'stage', 'bin_left', 'bin_right', 'count'])
        for (phase, mode, stage), values in sorted(reals.items()):
            low, high = ranges[phase, stage]
            width = (high-low)/64 if high > low else 1
            counts = Counter(min(63, int((v-low)/width)) for v in values)
            for i in range(64):
                writer.writerow([phase, mode, stage, low+i*width, low+(i+1)*width, counts[i]])
    output = {
        'scope': 'First at most 32 captured calls per direction on one quantization-only trajectory; '
                 'logical directed-interface packets including intra-rank copies; not measured network traffic.',
        'baseline': '8 bytes per captured scalar on the same packets/calls, including exact corners; '
                    'not total bytes of a separately converged uncompressed solve.',
        'error': 'Decoded versus frozen raw values; shared entries counted once per directed packet. '
                 'This is transfer error, not global solution error.',
        'setup': 'Plan bytes are shown once per directed interface per mode; R/P use the same plans. '
                 'Codebook estimates are per-call tables, not measured traffic. Neither is in wire bytes. '
                 'Actual production codebooks are rebuilt per exchange; one-time training is not assumed.',
        'histograms': 'Packet-weighted empirical counts without synthetic alphabet seeds. '
                      'Floating quantization-input values are coefficients for before, nodal values otherwise. '
                      'Integer lifted values are stored modulo 2^32; see stages.csv. '
                      'Modes use the existing production training populations, which can differ '
                      'between union-based and per-neighbor collection.',
        'results': [dict(phase=k[0], mode=k[1], **a) for k, a in sorted(totals.items())],
        'sha256': {name: hashlib.sha256((root/name).read_bytes()).hexdigest()
                   for name in ('packets.csv', 'stages.csv', 'codebooks.csv')}
    }
    (root/'summary.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    lines = ['Frozen identical-input comparison (logical payloads, setup excluded)',
             'phase       mode              bytes   saved %    relative L2     max error']
    for row in output['results']:
        rel = row['relative_l2_error']
        saving = row['net_saved_percent']
        lines.append(f"{row['phase']:11s} {row['mode']:12s} {row['wire_bytes']:10d} "
                     f"{saving if saving is not None else 0:9.2f} "
                     f"{rel if rel is not None else 0:14.5e} {row['max_absolute_error']:13.5e}")
    lines += ['', output['baseline'], output['setup']]
    (root/'summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
    return output


if __name__ == '__main__':
    summarize(Path(sys.argv[1]))
