"""Short exact enumeration, conditional identities and finite-difference checks.

Uses only Python's standard library. The enumeration scorer walks a structure
rather than using LinearKnot's deductions, schedule or reverse pass.
"""
import argparse
import json
import math
import subprocess


def nested(sequence):
    def visit(a, b):
        if a > b:
            return ['']
        out = ['.' + x for x in visit(a + 1, b)]
        for j in range(a + 4, b + 1):
            if sequence[a] + sequence[j] in ('AU', 'UA', 'CG', 'GC', 'GU', 'UG'):
                out += ['(' + x + ')' + y for x in visit(a + 1, j - 1) for y in visit(j + 1, b)]
        return out
    return visit(0, len(sequence) - 1)


def pairs(structure):
    stack, result = [], []
    for i, c in enumerate(structure, 1):
        if c == '(':
            stack.append(i)
        elif c == ')':
            result.append((stack.pop(), i))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--probe', required=True)
    parser.add_argument('--score', required=True)
    args = parser.parse_args()
    process = subprocess.Popen([args.probe], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    def send(message):
        process.stdin.write(message + '\n')
        process.stdin.flush()
        line = process.stdout.readline()
        if not line:
            raise RuntimeError('probe closed unexpectedly')
        response = json.loads(line)
        if 'error' in response:
            raise RuntimeError(response)
        return response

    def evaluate(cid, displacement):
        return send('E ' + cid + ' ' + ' '.join(map(str, displacement)))

    def matrix(result):
        return {(i, j): p for i, j, p in result['bpp']}

    def dot(result, adjoint):
        return math.fsum(p * adjoint.get((i, j), 0.) for i, j, p in result['bpp'])

    max_z = max_bpp = max_vjp = max_likelihood = 0.
    zero = [0.] * 12
    try:
        metadata = send('M metadata')
        assert metadata['backend'] == 'contrafold-pk-v2'
        assert abs(metadata['ordinary_branch_log_score'] + .9253883752) < 1e-12
        for k, seq in enumerate(['GGGGGAACCAACCC', 'GCGAAACGAAACGC', 'GGAAACGAAACCC',
                                 'ACGUACGUACGUACGUAC', 'GCGAAACGCGAAACGC']):
            structures = nested(seq)
            scored = subprocess.run([args.score], input=''.join(seq + ' ' + s + '\n' for s in structures),
                                    text=True, capture_output=True, check=True)
            scores = list(map(float, scored.stdout.splitlines()))
            assert len(scores) == len(structures)
            hi = max(scores)
            z = hi + math.log(math.fsum(math.exp(x - hi) for x in scores))
            expected = {}
            for structure, score in zip(structures, scores):
                for pair in pairs(structure):
                    expected[pair] = expected.get(pair, 0.) + math.exp(score - z)
            cid = 'free' + str(k)
            send('R ' + cid + ' ' + seq + ' ' + '.' * len(seq))
            result = evaluate(cid, zero)
            actual = matrix(result)
            max_z = max(max_z, abs(z - result['logZ']))
            max_bpp = max(max_bpp, max((abs(actual.get(p, 0.) - expected.get(p, 0.))
                                      for p in actual.keys() | expected.keys()), default=0.))
            moved = evaluate(cid, [.123] * 12)
            assert moved['logZ'] == result['logZ'] and moved['bpp'] == result['bpp']
            assert max(map(abs, moved['expected_feature'])) < 1e-10
            send('D ' + cid)
        assert max_z < 3e-6 and max_bpp < 3e-6, (max_z, max_bpp)
        seq, scaffold = 'GGGGGAACCAACCC', '((.....)).....'
        send('R pk ' + seq + ' ' + scaffold)
        result = evaluate('pk', zero)
        actual = matrix(result)
        for p in pairs(scaffold):
            assert abs(actual[p] - 1.) < 1e-12
        # At least one generated pair crosses a fixed pair in this witness.
        assert any(p > 1e-5 and any(i < k < j < l or k < i < l < j for k, l in pairs(scaffold))
                   for i, j, p in result['bpp'])
        incident = [0.] * (len(seq) + 1)
        for i, j, p in result['bpp']:
            assert 0. <= p <= 1. + 1e-12
            incident[i] += p
            incident[j] += p
        assert max(incident) <= 1. + 1e-10
        adjoint = {(i, j): ((11 * i + 7 * j) % 19 - 9) / 10 for i, j, p in result['bpp']}
        request = ' '.join(f'{i} {j} {a}' for (i, j), a in adjoint.items())
        analytic = send(f'V pk {result["version"]} {len(adjoint)} ' + request)['gradient']
        step = 1e-4
        for k in range(12):
            d = zero.copy(); d[k] = step
            plus = evaluate('pk', d)
            d[k] = -step
            minus = evaluate('pk', d)
            numerical = (dot(plus, adjoint) - dot(minus, adjoint)) / (2 * step)
            max_vjp = max(max_vjp, abs(analytic[k] - numerical))
        assert max_vjp < 2e-6, max_vjp
        # Ordinary branch is fixed even when the legacy DP09 coordinate moves.
        d = zero.copy(); d[6] = 17.3
        moved = evaluate('pk', d)
        assert moved['logZ'] == result['logZ'] and moved['bpp'] == result['bpp']
        assert abs(moved['expected_feature'][6]) < 1e-12
        evaluate('pk', zero)
        likelihood = send('L pk (([[[..))..]]]')
        assert likelihood['reachable'] and likelihood['log_probability'] <= 1e-10
        for k in range(12):
            d = zero.copy(); d[k] = step
            evaluate('pk', d)
            plus = send('L pk (([[[..))..]]]')['log_probability']
            d[k] = -step
            evaluate('pk', d)
            minus = send('L pk (([[[..))..]]]')['log_probability']
            max_likelihood = max(max_likelihood, abs((plus - minus) / (2 * step)
                                                     - likelihood['gradient_log_probability'][k]))
        assert max_likelihood < 2e-6, max_likelihood
        impossible = send('L pk ..............')
        assert not impossible['reachable'] and impossible['log_probability'] is None
        print(json.dumps(dict(passed=True, enumeration_sequences=5, derivative_coordinates=12,
                              max_logZ_error=max_z, max_bpp_error=max_bpp,
                              max_vjp_error=max_vjp, max_likelihood_gradient_error=max_likelihood), indent=2))
    finally:
        process.stdin.close()
        process.wait(timeout=10)
        process.stdout.close()
    assert process.returncode == 0


if __name__ == '__main__':
    main()
