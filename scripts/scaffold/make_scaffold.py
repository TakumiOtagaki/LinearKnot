#!/usr/bin/env python3
"""Build LinearKnot scaffold G using the paper's LinearPartition rule.

Default: LinearPartition-V, beam 200, dangles 2; select all BPP pairs with
p >= 0.7, then keep stacked helices of at least two pairs. No MEA step.
--from-mea restores the old MEA-filter mode (gamma 3, threshold 0.5,
minimum helix 1, beam 100), with optional parameter overrides.

Input: FASTA or one sequence per line. Output: >name / sequence / G records,
preceded by a # comment with settings and LinearPartition provenance.
LinearKnot ignores # comments. --metadata also writes a JSON sidecar;
--no-metadata-header emits only the three-line records.

    scripts/scaffold/build_linearpartition.sh
    scripts/scaffold/make_scaffold.py seqs.fa > scaffolds.txt
    scripts/scaffold/make_scaffold.py --from-mea seqs.fa > old_scaffolds.txt
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LP = ROOT / 'build/third_party/LinearPartition/linearpartition'
CANONICAL = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}


def read_sequences(handle):
    name, chunks, count = None, [], 0
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('>'):
            if chunks:
                yield name, ''.join(chunks)
            name, chunks = line[1:].strip(), []
            continue
        if name is None:
            count += 1
            yield 'seq%d' % count, line
            continue
        chunks.append(line)
    if chunks:
        yield name, ''.join(chunks)


def pairs_of(structure):
    stack, pairs = [], {}
    for pos, char in enumerate(structure, start=1):
        if char == '(':
            stack.append(pos)
        elif char == ')':
            if not stack:
                raise ValueError('unbalanced MEA structure')
            pairs[stack.pop()] = pos
        elif char != '.':
            raise ValueError('unsupported MEA structure character')
    if stack:
        raise ValueError('unbalanced MEA structure')
    return pairs


def read_bpp(handle, n):
    bpp = {}
    for number, line in enumerate(handle, 1):
        fields = line.split()
        if not fields or fields[0].startswith(('#', '>')):
            continue
        if len(fields) != 3:
            raise ValueError('malformed BPP row %d' % number)
        i, j, p = int(fields[0]), int(fields[1]), float(fields[2])
        # Match the paper reader's tolerance for LinearPartition float output.
        if not 1 <= i < j <= n or not math.isfinite(p) or not 0 <= p <= 1.00001:
            raise ValueError('invalid BPP coordinates or probability at row %d' % number)
        if (i, j) in bpp:
            raise ValueError('duplicate BPP pair at row %d' % number)
        bpp[i, j] = p
    return bpp


def run_linearpartition(binary, sequence, beam, gamma=None, threshold=0.7, dangles=2):
    with tempfile.TemporaryDirectory() as tmp:
        bpp_path = os.path.join(tmp, 'bpp.txt')
        cmd = [str(binary), '-V', '-b', str(beam), '--dangles', str(dangles)]
        if gamma is None:
            cmd += ['-c', str(threshold)]
        else:
            # The old path needs MEA and its original complete BPP output.
            cmd += ['-M', '--gamma', str(gamma)]
        cmd += ['-r', bpp_path]
        result = subprocess.run(cmd, input=sequence + '\n', capture_output=True,
                                text=True, check=True)
        mea = None
        if gamma is not None:
            lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
            mea = next((line for line in reversed(lines)
                        if len(line) == len(sequence) and set(line) <= set('().')), None)
            if mea is None:
                raise ValueError('LinearPartition did not return an MEA structure')
        with open(bpp_path) as handle:
            bpp = read_bpp(handle, len(sequence))
    return mea, bpp


def scaffold(sequence, mea, bpp, threshold=0.7, min_helix=2):
    """Select from all BPP pairs, or filter an explicit MEA when supplied."""
    if mea is not None and len(mea) != len(sequence):
        raise ValueError('MEA and sequence lengths differ')
    candidates = bpp.keys() if mea is None else pairs_of(mea).items()
    selected = [(i, j) for i, j in candidates if bpp.get((i, j), 0.) >= threshold]
    positions = []
    for i, j in selected:
        if not 1 <= i < j <= len(sequence):
            raise ValueError('selected pair is outside the sequence')
        if sequence[i - 1] + sequence[j - 1] not in CANONICAL:
            raise ValueError('selected pair is noncanonical')
        positions.extend((i, j))
    # Do not silently choose between incompatible high-probability pairs.
    if len(positions) != len(set(positions)):
        raise ValueError('selected BPP pairs share an endpoint; raise the threshold or use --from-mea')
    kept = dict(selected)
    raw = ['.'] * len(sequence)
    for i, j in kept.items():
        raw[i - 1], raw[j - 1] = '(', ')'
    if pairs_of(''.join(raw)) != kept:
        raise ValueError('selected BPP pairs cross; raise the threshold or use --from-mea')
    filtered = {}
    for i, j in kept.items():
        if kept.get(i - 1) == j + 1:
            continue
        run = [(i, j)]
        while kept.get(run[-1][0] + 1) == run[-1][1] - 1:
            run.append((run[-1][0] + 1, run[-1][1] - 1))
        if len(run) >= min_helix:
            filtered.update(run)
    out = ['.'] * len(sequence)
    for i, j in filtered.items():
        out[i - 1], out[j - 1] = '(', ')'
    return ''.join(out)


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def linearpartition_provenance(binary):
    path = Path(binary).resolve()
    info = {'wrapper': str(path), 'wrapper_sha256': file_hash(path),
            'commit': None, 'tracked_source_modified': None, 'vienna_binary_sha256': None}
    native = path.parent / 'bin/linearpartition_v'
    if native.is_file():
        info['vienna_binary_sha256'] = file_hash(native)
    # A standalone wrapper inside this checkout must not inherit LinearKnot's
    # revision and report it as a LinearPartition commit.
    try:
        top = subprocess.run(['git', '-C', str(path.parent), 'rev-parse', '--show-toplevel'],
                             capture_output=True, text=True, check=True).stdout.strip()
        if Path(top).resolve() == path.parent:
            subprocess.run(['git', '-C', top, 'ls-files', '--error-unmatch', path.name],
                           capture_output=True, text=True, check=True)
            info['commit'] = subprocess.run(['git', '-C', top, 'rev-parse', 'HEAD'],
                                           capture_output=True, text=True, check=True).stdout.strip()
            status = subprocess.run(['git', '-C', top, 'status', '--porcelain', '--untracked-files=no'],
                                    capture_output=True, text=True, check=True,
                                    env=dict(os.environ, GIT_OPTIONAL_LOCKS='0')).stdout
            info['tracked_source_modified'] = bool(status)
    except (OSError, subprocess.CalledProcessError):
        pass  # Unknown revision is explicit; executable hashes are still recorded.
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', nargs='?', default='-', help='FASTA file (default: stdin)')
    parser.add_argument('--from-mea', action='store_true', help='use the previous MEA-filter rule')
    parser.add_argument('--threshold', type=float, help='minimum BPP (default: 0.7; --from-mea: 0.5)')
    parser.add_argument('--min-helix', type=int, help='minimum stacked pairs (default: 2; --from-mea: 1)')
    parser.add_argument('--beam', type=int, help='LinearPartition beam (default: 200; --from-mea: 100)')
    parser.add_argument('--gamma', type=float, help='MEA gamma; requires --from-mea (default: 3)')
    parser.add_argument('--dangles', type=int, choices=(0, 2), default=2, help='LinearPartition dangles (default: 2)')
    parser.add_argument('--linearpartition', type=Path, default=DEFAULT_LP, help='path to the LinearPartition wrapper')
    parser.add_argument('--metadata', type=Path, help='also write settings and record hashes to this JSON file')
    parser.add_argument('--no-metadata-header', action='store_true', help='omit the # metadata comment from stdout')
    args = parser.parse_args()
    if args.gamma is not None and not args.from_mea:
        parser.error('--gamma requires --from-mea')
    threshold = args.threshold if args.threshold is not None else (.5 if args.from_mea else .7)
    min_helix = args.min_helix if args.min_helix is not None else (1 if args.from_mea else 2)
    beam = args.beam if args.beam is not None else (100 if args.from_mea else 200)
    gamma = (args.gamma if args.gamma is not None else 3.) if args.from_mea else None
    if not math.isfinite(threshold) or not 0 <= threshold <= 1:
        parser.error('--threshold must be finite and between 0 and 1')
    if min_helix < 1 or beam < 1:
        parser.error('--min-helix and --beam must be positive')
    if gamma is not None and (not math.isfinite(gamma) or gamma < 0):
        parser.error('--gamma must be finite and nonnegative')
    if not args.linearpartition.is_file():
        parser.error('linearpartition not found; run scripts/scaffold/build_linearpartition.sh')
    try:
        metadata = {'schema_version': 1, 'mode': 'mea_filter' if args.from_mea else 'bpp_threshold',
                    'threshold': threshold, 'min_helix': min_helix, 'beam': beam, 'dangles': args.dangles,
                    'gamma': gamma, 'linearpartition': linearpartition_provenance(args.linearpartition)}
        if not args.no_metadata_header:
            print('# linearknot_scaffold_metadata=' + json.dumps(metadata, sort_keys=True))
        records = []
        handle = sys.stdin if args.input == '-' else open(args.input)
        try:
            for name, raw in read_sequences(handle):
                sequence = raw.upper().replace('T', 'U')
                if not sequence or set(sequence) - set('ACGU'):
                    raise ValueError('record %s: sequence must contain only A/C/G/U (T is normalized to U)' % name)
                mea, bpp = run_linearpartition(args.linearpartition, sequence, beam, gamma, threshold, args.dangles)
                g = scaffold(sequence, mea, bpp, threshold, min_helix)
                print('>%s\n%s\n%s' % (name, sequence, g))
                records.append({'name': name, 'length': len(sequence),
                                'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                                'scaffold_sha256': hashlib.sha256(g.encode()).hexdigest()})
        finally:
            if handle is not sys.stdin:
                handle.close()
        if args.metadata:
            args.metadata.write_text(json.dumps(dict(metadata, records=records), indent=2) + '\n')
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(2, 'make_scaffold: %s\n' % error)


if __name__ == '__main__':
    main()
