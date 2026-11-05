#!/usr/bin/env python3
"""Exact bounded enumeration for multi-width (2,2,k), standard library only.

No theorem-wide or formal computation gate is claimed by this script.
Canonical keys are row-HNF minima over all 24 orders of the four vertices.
All short integral covectors are enumerated via bounded edge evaluations;
directions differing by sign are represented once.
"""
import argparse
import itertools
import json
import time
from pathlib import Path

PERMS = tuple(itertools.permutations(range(4)))


def determinant(a):
    if len(a) == 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def adjugate(a):
    if len(a) == 2:
        return ((a[1][1], -a[0][1]), (-a[1][0], a[0][0]))
    return tuple(tuple((-1) ** (i + j) * determinant([
        [a[r][c] for c in range(3) if c != i]
        for r in range(3) if r != j])
        for j in range(3)) for i in range(3))


def edge_matrix(points, order=(0, 1, 2, 3)):
    base = points[order[0]]
    return tuple(tuple(points[j][i] - base[i] for j in order[1:])
                 for i in range(len(base)))


def row_hnf(matrix):
    """Unique upper row Hermite form of the row lattice, full row rank.

    Pivot columns increase; pivots are positive; entries above a pivot lie
    in [0,pivot). Only invertible integral row operations are performed.
    Rectangular 2x3 matrices retain the full projected four-point multiset.
    """
    a = [list(row) for row in matrix]
    m, n, pivot = len(a), len(a[0]), 0
    for col in range(n):
        nonzero = next((r for r in range(pivot, m) if a[r][col]), None)
        if nonzero is None:
            continue
        a[pivot], a[nonzero] = a[nonzero], a[pivot]
        for r in range(pivot + 1, m):
            while a[r][col]:
                q = a[pivot][col] // a[r][col]
                a[pivot] = [x - q * y for x, y in zip(a[pivot], a[r])]
                a[pivot], a[r] = a[r], a[pivot]
        if a[pivot][col] < 0:
            a[pivot] = [-x for x in a[pivot]]
        for r in range(pivot):
            q = a[r][col] // a[pivot][col]
            a[r] = [x - q * y for x, y in zip(a[r], a[pivot])]
        pivot += 1
        if pivot == m:
            break
    if pivot != m:
        raise ValueError('row_hnf requires full row rank')
    return tuple(tuple(row) for row in a)


def canonical_key(points):
    return min(row_hnf(edge_matrix(points, order)) for order in PERMS)


def vertices_from_key(key):
    d = len(key)
    return ((0,) * d,) + tuple(tuple(key[i][j] for i in range(d))
                              for j in range(3))


def independent_edge_data(points):
    e = edge_matrix(points)
    d = len(e)
    for columns in itertools.combinations(range(3), d):
        square = tuple(tuple(e[i][j] for j in columns) for i in range(d))
        det = determinant(square)
        if det:
            return e, square, det, adjugate(square)
    raise ValueError('points are not full dimensional')


def bounded_edge_values(d, r):
    """All nonzero tuples up to sign with osc(0,a_1,...,a_d)<=r."""
    for values in itertools.product(range(-r, r + 1), repeat=d):
        first = next((v for v in values if v), 0)
        if first <= 0:
            continue
        if max(0, *values) - min(0, *values) <= r:
            yield values


_VALUE_CACHE = {}


def short_covectors(points, r, data=None):
    """Exhaustive finite integral width<=r directions, one per sign.

    Choose d independent edge columns E. u E = a has |a_i|<=r.
    u = a adj(E)/det(E) is integral exactly when every numerator is divisible
    by det(E). For planar multisets the unused edge is checked afterwards.
    """
    e, square, det, adj = data or independent_edge_data(points)
    d = len(e)
    if (d, r) not in _VALUE_CACHE:
        _VALUE_CACHE[(d, r)] = tuple(bounded_edge_values(d, r))
    values_list = _VALUE_CACHE[(d, r)]
    for values in values_list:
        numerators = tuple(sum(values[i] * adj[i][j] for i in range(d))
                           for j in range(d))
        if any(v % det for v in numerators):
            continue
        u = tuple(v // det for v in numerators)
        if d != 3:
            all_values = tuple(sum(u[i] * e[i][j] for i in range(d))
                               for j in range(3))
            if max(0, *all_values) - min(0, *all_values) > r:
                continue
        yield u


def span_rank(vectors, d):
    basis = []
    for v in vectors:
        if not any(v):
            continue
        if not basis:
            basis.append(v)
        elif len(basis) == 1:
            if any(basis[0][i] * v[j] != basis[0][j] * v[i]
                   for i in range(d) for j in range(i + 1, d)):
                basis.append(v)
        elif d == 3 and determinant((basis[0], basis[1], v)):
            return 3
        if len(basis) == d:
            return d
    return len(basis)


def multi_width(points, upper=None):
    """Successive ranks of all short integral directions, exact."""
    d = len(points[0])
    if upper is None:
        upper = max(max(p[i] for p in points) - min(p[i] for p in points)
                    for i in range(d))
    data = independent_edge_data(points)
    result = []
    for r in range(1, upper + 1):
        rank = span_rank(short_covectors(points, r, data), d)
        result.extend([r] * (rank - len(result)))
        if len(result) == d:
            return tuple(result)
    raise ValueError('upper bound does not supply d independent directions')


def planar_representatives():
    """All rank-two four-point multisets in the 2x2 box with mwidth=(2,2)."""
    grid = tuple(itertools.product(range(3), repeat=2))
    representatives = {}
    examined = 0
    for points in itertools.combinations_with_replacement(grid, 4):
        try:
            data = independent_edge_data(points)
        except ValueError:
            continue
        examined += 1
        # The coordinate axes already supply rank two at width two.
        if next(short_covectors(points, 1, data), None) is not None:
            continue
        key = canonical_key(points)
        representatives.setdefault(key, points)
    return tuple(representatives[key] for key in sorted(representatives)), examined


def has_multi_width_22k(points, k):
    data = independent_edge_data(points)
    if next(short_covectors(points, 1, data), None) is not None:
        return False
    # Projected x,y axes supply rank two by width 2. For k>=3, a short
    # direction is outside their span precisely when its z coordinate is nonzero.
    if k >= 3:
        if any(u[2] for u in short_covectors(points, k - 1, data)):
            return False
        if max(p[2] for p in points) - min(p[2] for p in points) <= k:
            return True
        return any(u[2] for u in short_covectors(points, k, data))
    return span_rank(short_covectors(points, 2, data), 3) == 3


def enumerate_22k(k, projections):
    bound = min(k + 1, max(4, k))
    heights = tuple(h for h in itertools.product(range(bound + 1), repeat=4)
                    if min(h) == 0 and max(h) >= k)
    representatives = {}
    candidates = degenerate = accepted = 0
    for projection in projections:
        for h in heights:
            points = tuple(p + (z,) for p, z in zip(projection, h))
            candidates += 1
            if determinant(edge_matrix(points)) == 0:
                degenerate += 1
                continue
            if not has_multi_width_22k(points, k):
                continue
            accepted += 1
            key = canonical_key(points)
            representatives.setdefault(key, points)
    records = [dict(key=key, vertices=vertices_from_key(key), box_vertices=representatives[key],
                    normalized_volume=abs(determinant(key)))
               for key in sorted(representatives)]
    return dict(k=k, lift_bound=bound, candidates=candidates,
                degenerate=degenerate, accepted_lifts=accepted,
                count=len(records), representatives=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ks', type=int, nargs='+', default=[2, 3, 4, 5, 6])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    projections, examined = planar_representatives()
    print(f'planar (2,2): {len(projections)} classes; rank-two box multisets {examined}', flush=True)
    result = dict(evidence='exact finite enumeration, conditional on lift finite-reduction proof',
                  planar_count=len(projections), projections=projections,
                  planar_rank_two_multisets=examined, results=[])
    for k in args.ks:
        item = enumerate_22k(k, projections)
        result['results'].append(item)
        print(f'k={k}: count={item["count"]}, candidates={item["candidates"]}, '
              f'accepted_lifts={item["accepted_lifts"]}, elapsed={time.monotonic()-start:.3f}s', flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    result['elapsed_seconds'] = time.monotonic() - start
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
