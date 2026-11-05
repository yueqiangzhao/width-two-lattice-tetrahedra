#!/usr/bin/env python3
"""Independent affine-map cross-check of bounded Q024 enumeration."""
import argparse
import ast
import hashlib
import importlib.util
import itertools
import json
import re
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('q024_enum', Path(__file__).with_name('enumerate.py'))
enum = importlib.util.module_from_spec(spec)
spec.loader.exec_module(enum)


def direct_equivalent(p, q):
    """Search all vertex bijections and solve U E = F directly.

    Any affine equivalence induces one of these 24 vertex permutations.
    A bijection is integral unimodular iff F adj(E)/det(E) is integral
    and |det(F)|=|det(E)|. This does not use Hermite normal forms.
    """
    e = enum.edge_matrix(p)
    de = enum.determinant(e)
    aq = enum.adjugate(e)
    if abs(de) != abs(enum.determinant(enum.edge_matrix(q))):
        return False
    for order in enum.PERMS:
        f = enum.edge_matrix(q, order)
        numerators = [[sum(f[i][t] * aq[t][j] for t in range(3))
                       for j in range(3)] for i in range(3)]
        if all(v % de == 0 for row in numerators for v in row):
            return True
    return False


def direct_short_covectors(points, r):
    """Independent primal coordinate-box search using pairwise inequalities.

    |u_j| <= r sum_i |adj(E)_ij|/|det E| follows from |uE|<=r.
    Loop over two coordinates; intersect exact intervals for the third
    imposed by -r <= u dot (p_i-p_j) <= r for all six vertex pairs.
    No edge-value tuples or divisibility filtering are used here.
    """
    e = enum.edge_matrix(points)
    det = abs(enum.determinant(e))
    adj = enum.adjugate(e)
    bounds = tuple(r * sum(abs(adj[i][j]) for i in range(3)) // det
                   for j in range(3))
    order = sorted(range(3), key=lambda j: bounds[j])
    i, j, dependent = order
    differences = tuple(tuple(p[t] - q[t] for t in range(3))
                        for p, q in itertools.combinations(points, 2))
    for a in range(-bounds[i], bounds[i] + 1):
        for b in range(-bounds[j], bounds[j] + 1):
            lo, hi = -bounds[dependent], bounds[dependent]
            valid = True
            for diff in differences:
                constant = a * diff[i] + b * diff[j]
                coefficient = diff[dependent]
                if coefficient == 0:
                    if abs(constant) > r:
                        valid = False
                        break
                    continue
                if coefficient < 0:
                    coefficient, constant = -coefficient, -constant
                lo = max(lo, -((r + constant) // coefficient))
                hi = min(hi, (r - constant) // coefficient)
                if lo > hi:
                    valid = False
                    break
            if valid:
                for c in range(lo, hi + 1):
                    u = [0, 0, 0]
                    u[i], u[j], u[dependent] = a, b, c
                    if next((v for v in u if v), 0) > 0:
                        yield tuple(u)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--finite', type=Path, required=True)
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    certificate = json.loads(args.finite.read_text())
    projections = tuple(tuple(tuple(p) for p in projection)
                        for projection in certificate['projections'])
    summary = dict(mechanism='direct integer affine maps over 24 vertex bijections',
                   pairwise_comparisons=0, accepted_lift_comparisons=0,
                   dataset_comparisons=0, dataset_width_checks=0,
                   independent_width_representatives=0, branches=[])
    source_path = args.dataset
    source = source_path.read_text()
    records = re.findall(r'ID:\s*(\d+)\s*Vertices:\s*(\[\[.*?\]\])', source)
    dataset = {}
    for identifier, vertices in records:
        k = int(identifier[4:6])
        if identifier[:4] == '0202':
            dataset.setdefault(k, []).append(tuple(tuple(p) for p in ast.literal_eval(vertices)))
    for item in certificate['results']:
        k = item['k']
        reps = tuple(tuple(tuple(p) for p in record['box_vertices'])
                     for record in item['representatives'])
        by_key = {enum.canonical_key(p): p for p in reps}
        for p in reps:
            # Compare complete short-direction sets (up to sign) at every
            # relevant width threshold, rather than merely compare ranks.
            for r in sorted({1, 2, k - 1, k}):
                edge_directions = {u if next(v for v in u if v) > 0 else tuple(-v for v in u)
                                   for u in enum.short_covectors(p, r)}
                direct_directions = set(direct_short_covectors(p, r))
                if edge_directions != direct_directions:
                    raise AssertionError(('independent width mismatch', k, r, p,
                                          edge_directions ^ direct_directions))
            summary['independent_width_representatives'] += 1
        for p, q in itertools.combinations(reps, 2):
            summary['pairwise_comparisons'] += 1
            if direct_equivalent(p, q):
                raise AssertionError(('HNF failed to merge equivalent reps', k, p, q))
        bound = item['lift_bound']
        for projection in projections:
            for h in itertools.product(range(bound + 1), repeat=4):
                if min(h) != 0 or max(h) < k:
                    continue
                points = tuple(p + (z,) for p, z in zip(projection, h))
                if enum.determinant(enum.edge_matrix(points)) == 0:
                    continue
                if not enum.has_multi_width_22k(points, k):
                    continue
                summary['accepted_lift_comparisons'] += 1
                if not direct_equivalent(points, by_key[enum.canonical_key(points)]):
                    raise AssertionError(('false HNF collision', k, points))
        source_keys = set()
        for p in dataset[k]:
            key = enum.canonical_key(p)
            source_keys.add(key)
            summary['dataset_width_checks'] += 1
            if enum.multi_width(p, k) != (2, 2, k):
                raise AssertionError(('dataset incorrect multiwidth', k, p))
            summary['dataset_comparisons'] += 1
            if key not in by_key or not direct_equivalent(p, by_key[key]):
                raise AssertionError(('dataset mismatch', k, p))
        if source_keys != set(by_key):
            raise AssertionError(('unequal class sets', k))
        summary['branches'].append(dict(k=k, enumerated_count=len(reps),
                                       dataset_count=len(dataset[k]),
                                       dataset_unique_count=len(source_keys),
                                       canonical_sets_equal=True))
    summary['elapsed_seconds'] = time.monotonic() - started
    summary['status'] = 'COMPUTATIONALLY VERIFIED IN SCOPE'
    summary['hashes'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (Path(__file__), args.finite, args.dataset,
                                      Path(__file__).with_name('enumerate.py'))}
    summary['replay_command'] = ' '.join(sys.argv)
    args.output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
