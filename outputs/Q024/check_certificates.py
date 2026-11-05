#!/usr/bin/env python3
"""Independent exact checker for the finite ingredients of Q024 lift circuits.

No HNF routine, polygon-vertex routine, or edge-value enumeration is used.
Only Python standard-library integer/Fraction arithmetic is used. All-k
theorem correctness remains with the separately supplied mathematical proof.
"""
import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
from fractions import Fraction as F
from functools import reduce
from math import ceil, floor, gcd
from pathlib import Path

PERMS = tuple(itertools.permutations(range(4)))


def det(matrix):
    a = [[F(x) for x in row] for row in matrix]
    result = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        v = a[j][j]
        result *= v
        for i in range(j + 1, len(a)):
            q = a[i][j] / v
            a[i] = [x - q * y for x, y in zip(a[i], a[j])]
    return result


def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError('singular matrix')
        a[j], a[pivot] = a[pivot], a[j]
        v = a[j][j]
        a[j] = [x / v for x in a[j]]
        for i in range(n):
            if i != j:
                q = a[i][j]
                a[i] = [x - q * y for x, y in zip(a[i], a[j])]
    return tuple(tuple(row[n:]) for row in a)


def edges(points):
    return tuple(tuple(points[j][i] - points[0][i] for j in range(1, 4))
                 for i in range(len(points[0])))


def independent_columns(matrix):
    d = len(matrix)
    for columns in itertools.combinations(range(3), d):
        square = tuple(tuple(row[j] for j in columns) for row in matrix)
        if det(square):
            return columns, square
    raise ValueError('not full dimensional')


def affine_map(p, q):
    """Exact U,t mapping corresponding indexed points; no canonical forms."""
    d = len(p[0])
    ep, eq = edges(p), edges(q)
    columns, square = independent_columns(ep)
    inv = inverse(square)
    target = tuple(tuple(row[j] for j in columns) for row in eq)
    u = tuple(tuple(sum(target[i][t] * inv[t][j] for t in range(d))
                    for j in range(d)) for i in range(d))
    translation = tuple(F(q[0][i]) - sum(u[i][j] * p[0][j] for j in range(d))
                        for i in range(d))
    if any(x.denominator != 1 for row in u for x in row):
        return None
    if any(x.denominator != 1 for x in translation) or abs(det(u)) != 1:
        return None
    if any(tuple(sum(u[i][j] * point[j] for j in range(d)) + translation[i]
                 for i in range(d)) != tuple(image)
           for point, image in zip(p, q)):
        return None
    return u, translation


def equivalent(p, q):
    if len(p[0]) == 3 and abs(det(edges(p))) != abs(det(edges(q))):
        return False
    return any(affine_map(p, tuple(q[i] for i in pi)) is not None for pi in PERMS)


def pair_bound_box(inequalities):
    """Bounds a,b from every independent pair of bounded linear forms.

    Input (dx,dy,lo,hi) means lo<=dx*a+dy*b<=hi. Mapping the four
    interval corners through an inverse 2x2 matrix gives valid bounds.
    Intersecting the resulting coordinate bounds is safe and exact.
    """
    lower, upper = [None, None], [None, None]
    for f, g in itertools.combinations(inequalities, 2):
        matrix = (f[:2], g[:2])
        if not det(matrix):
            continue
        inv = inverse(matrix)
        corners = tuple(tuple(sum(inv[i][j] * value[j] for j in range(2))
                              for i in range(2))
                        for value in itertools.product(f[2:], g[2:]))
        for i in range(2):
            lo, hi = min(v[i] for v in corners), max(v[i] for v in corners)
            lower[i] = lo if lower[i] is None else max(lower[i], lo)
            upper[i] = hi if upper[i] is None else min(upper[i], hi)
    assert all(v is not None for v in lower + upper), 'unbounded coefficient search'
    return ceil(lower[0]), floor(upper[0]), ceil(lower[1]), floor(upper[1])


def all_width_one_directions(points):
    constraints = [(p[0] - q[0], p[1] - q[1], F(-1), F(1))
                   for p, q in itertools.combinations(points, 2)]
    lo_a, hi_a, lo_b, hi_b = pair_bound_box(constraints)
    return {(a, b) for a in range(lo_a, hi_a + 1) for b in range(lo_b, hi_b + 1)
            if (a or b) and all(lo <= a * x + b * y <= hi
                               for x, y, lo, hi in constraints)}


def classify_planar(expected):
    grid = tuple(itertools.product(range(3), repeat=2))
    classes, count, rejected = [], 0, 0
    for points in itertools.combinations_with_replacement(grid, 4):
        try:
            independent_columns(edges(points))
        except ValueError:
            continue
        count += 1
        if all_width_one_directions(points):
            rejected += 1
            continue
        matches = [i for i, q in enumerate(classes) if equivalent(points, q)]
        assert len(matches) <= 1
        if not matches:
            classes.append(points)
    assert count == 354 and len(classes) == 13
    mapping = []
    for p in classes:
        matches = [i for i, q in enumerate(expected) if equivalent(p, q)]
        assert len(matches) == 1
        mapping.append(matches[0])
    assert sorted(mapping) == list(range(13))
    return dict(rank_two_box_multisets=count, rejected_width_one=rejected,
                direct_affine_classes=len(classes), mapping_to_enumerator=mapping)


def circuit(points):
    # Solve the affine relation using the first three columns; all full
    # circuits have every three planar points affinely independent.
    a = tuple(tuple((1, *points[j])[i] for j in range(3)) for i in range(3))
    inv = inverse(a)
    rhs = (1, *points[3])
    rational = tuple(sum(inv[i][j] * rhs[j] for j in range(3)) for i in range(3)) + (F(-1),)
    denominator = 1
    for v in rational:
        denominator = denominator * v.denominator // gcd(denominator, v.denominator)
    integer = tuple(int(v * denominator) for v in rational)
    common = abs(reduce(gcd, integer))
    integer = tuple(v // common for v in integer)
    if integer[0] < 0:
        integer = tuple(-v for v in integer)
    return integer


def shear_membership(points, difference):
    """difference_i = t + a*x_i+b*y_i, using Fraction Gaussian solve."""
    square = tuple((1, *p) for p in points[:3])
    inv = inverse(square)
    coefficients = tuple(sum(inv[i][j] * difference[j] for j in range(3)) for i in range(3))
    return all(v.denominator == 1 for v in coefficients) and all(
        coefficients[0] + coefficients[1] * p[0] + coefficients[2] * p[1] == value
        for p, value in zip(points, difference))


def quotient_membership(points, v, difference):
    # Removing translations gives the 3x3 lattice generated by x,y,v.
    basis = tuple(tuple((p[0] - points[0][0], p[1] - points[0][1], v[i] - v[0]))
                  for i, p in enumerate(points) if i)
    inv = inverse(basis)
    anchored = tuple(difference[i] - difference[0] for i in range(1, 4))
    coefficients = tuple(sum(inv[i][j] * anchored[j] for j in range(3)) for i in range(3))
    return all(x.denominator == 1 for x in coefficients)


def offsets(points, c, r):
    pos = [i for i, value in enumerate(c) if value > 0]
    neg = [i for i, value in enumerate(c) if value < 0]
    s = sum(c[i] for i in pos)
    t = sum(a * b for a, b in zip(c, r))
    upper = max(r[i] for i in pos) - min(r[j] for j in neg)
    constraints = []
    for i in pos:
        for j in neg:
            weight = F(c[i] * (-c[j]), s * s)
            lower = (F(t, s) - (1 - weight) * upper) / weight
            delta = r[i] - r[j]
            constraints.append((points[i][0] - points[j][0],
                                points[i][1] - points[j][1],
                                lower - delta, F(upper - delta)))
    # Because the weighted mean of cross-pair heights is t/s, the
    # upper constraints imply the lower constraints displayed above.
    bounds = pair_bound_box(constraints)
    lo_a, hi_a, lo_b, hi_b = bounds
    scored = []
    for a in range(lo_a, hi_a + 1):
        for b in range(lo_b, hi_b + 1):
            h = tuple(value + a * p[0] + b * p[1] for value, p in zip(r, points))
            score = max(h[i] for i in pos) - min(h[j] for j in neg)
            if score <= upper:
                threshold = max(0, max(h[j] for j in neg) - max(h[i] for i in pos),
                                min(h[j] for j in neg) - min(h[i] for i in pos))
                scored.append((score, threshold, a, b))
    assert scored
    chosen = min(scored)
    return chosen, bounds, len(scored)


def strip_width(points, heights):
    upper = max(heights) - min(heights)
    constraints = [(p[0] - q[0], p[1] - q[1],
                    F(-upper - h + k), F(upper - h + k))
                   for (p, h), (q, k) in itertools.combinations(zip(points, heights), 2)]
    lo_a, hi_a, lo_b, hi_b = pair_bound_box(constraints)
    return min(max(values) - min(values)
               for a in range(lo_a, hi_a + 1) for b in range(lo_b, hi_b + 1)
               for values in [tuple(h + a * x + b * y for h, (x, y) in zip(heights, points))])


def check_projection(cert):
    points = tuple(tuple(p) for p in cert['points'])
    c = circuit(points)
    assert tuple(cert['relation']) == c
    assert all(c) and sum(c) == 0
    assert all(sum(c[i] * points[i][j] for i in range(4)) == 0 for j in range(2))
    v = tuple(int(value > 0) for value in c)
    s = sum(value for value in c if value > 0)
    minors = [abs(int(det(tuple((1, *points[j]) for j in indices))))
              for indices in itertools.combinations(range(4), 3)]
    D = reduce(gcd, minors)
    basis = edges(points) + (tuple(v[i] - v[0] for i in range(1, 4)),)
    index = abs(int(det(basis)))
    assert index == D * s == len(cert['residuals'])
    assert D == cert['torsion_order'] and s == cert['positive_sum']
    residuals = tuple(tuple(r) for r in cert['residuals'])
    for r in residuals:
        assert 0 <= sum(a * b for a, b in zip(c, r)) < s
    for r, t in itertools.combinations(residuals, 2):
        assert not quotient_membership(points, v, tuple(a - b for a, b in zip(r, t)))
    # Index-many distinct residual classes prove completeness of this finite
    # quotient, using the elementary determinant-index theorem.
    automorphisms = [pi for pi in PERMS
                     if affine_map(points, tuple(points[i] for i in pi)) is not None]
    assert set(automorphisms) == {tuple(a['permutation']) for a in cert['automorphisms']}
    actions = []
    for pi in automorphisms:
        sign = c[pi[0]] // c[0]
        assert sign in (-1, 1) and all(c[pi[i]] == sign * c[i] for i in range(4))
        root = next(a for a in cert['automorphisms'] if tuple(a['permutation']) == pi)
        assert root['orientation'] == sign
        action = []
        for r in residuals:
            image = tuple(sign * r[i] for i in pi)
            matches = [j for j, t in enumerate(residuals)
                       if quotient_membership(points, v, tuple(a - b for a, b in zip(image, t)))]
            assert len(matches) == 1
            target = residuals[matches[0]]
            assert sum(a * b for a, b in zip(c, image)) == sum(a * b for a, b in zip(c, target))
            assert shear_membership(points, tuple(a - b for a, b in zip(image, target)))
            action.append(matches[0])
        assert action == root['residual_action']
        actions.append(action)
    remaining, orbits = set(range(index)), []
    while remaining:
        todo = [min(remaining)]
        orbit = set(todo)
        while todo:
            i = todo.pop()
            for action in actions:
                j = action[i]
                if j not in orbit:
                    orbit.add(j)
                    todo.append(j)
        assert orbit <= remaining
        remaining -= orbit
        orbits.append(sorted(orbit))
    assert orbits == cert['residual_orbits']
    widths = []
    for r, width in zip(residuals, cert['width_certificates']):
        (C, q0, a, b), bounds, candidate_count = offsets(points, c, r)
        assert (C, q0, a, b) == (width['C'], width['q0'], *width['minimizer'])
        assert candidate_count == width['candidates']
        assert q0 + C <= 5
        small = []
        for q in range(q0):
            h = tuple(q * bit + value for bit, value in zip(v, r))
            if sum(a * b for a, b in zip(c, h)) == 0:
                small.append(dict(q=q, degenerate=True))
            else:
                w = strip_width(points, h)
                assert w < 5
                small.append(dict(q=q, width=w))
        assert small == width['small_exceptions']
        widths.append(dict(C=C, q0=q0, minimizer=[a, b],
                           independent_coefficient_box=bounds,
                           coefficient_candidates=candidate_count,
                           small_exceptions=small))
    assert len(orbits) == cert['contribution_for_every_k_ge_5']
    return dict(projection_index=cert['projection_index'], primitive_relation=c,
                index=index, torsion_order=D, residual_count=index,
                direct_automorphism_count=len(automorphisms), residual_orbits=orbits,
                width_offsets=widths, orbit_count=len(orbits))


def degenerate_recipes(projections, k):
    """The six supplied theoretical recipes, checked here only at finite k."""
    heights = []
    heights.append([(0, k, 0, t) for t in range(3) if t <= (2 * k - t) % 3])
    heights.append([(0, k, (-t) % 4, 0) for t in range(4) if t <= (k - t) % 4])
    heights.append([(0, t, 0, k) for t in range(4) if t <= (2 * k - t) % 4])
    remaining, representatives = set(itertools.product(range(2), repeat=2)), []
    while remaining:
        start = min(remaining)
        orbit, todo = {start}, [start]
        while todo:
            a, b = todo.pop()
            for image in ((b, a), ((a + k) % 2, (b + k) % 2)):
                if image not in orbit:
                    orbit.add(image)
                    todo.append(image)
        remaining -= orbit
        representatives.append(start)
    heights.append([(0, k, a, b) for a, b in representatives])
    heights.append([(0, 0, d, t) for d in (2 * k - 1, 2 * k) for t in (0, 1)])
    heights.append([(0, 0, 2 * k - 1, 0), (0, 0, 2 * k, 0), (0, 0, 2 * k, 1)])
    generated = []
    for projection, options in zip(projections[:6], heights):
        for h in options:
            generated.append(tuple(tuple(p) + (z,) for p, z in zip(projection, h)))
    return generated, [len(options) for options in heights]


def recipe_matching(circuits, finite, ks=(5, 6)):
    results = []
    for k in ks:
        record = next(item for item in finite['results'] if item['k'] == k)
        finite_full = []
        for r in record['representatives']:
            points = tuple(tuple(p) for p in r['box_vertices'])
            projected = tuple(p[:2] for p in points)
            if any(equivalent(projected, tuple(tuple(p) for p in c['points'])) for c in circuits):
                finite_full.append(points)
        generated = []
        for c in circuits:
            points = tuple(tuple(p) for p in c['points'])
            relation = tuple(c['relation'])
            v = tuple(int(value > 0) for value in relation)
            for orbit in c['residual_orbits']:
                i = orbit[0]
                r = tuple(c['residuals'][i])
                C = c['width_certificates'][i]['C']
                h = tuple((k - C) * bit + value for bit, value in zip(v, r))
                generated.append(tuple(p + (z,) for p, z in zip(points, h)))
        assert len(generated) == len(finite_full) == 29
        hits = []
        for p in generated:
            matches = [i for i, q in enumerate(finite_full) if equivalent(p, q)]
            assert len(matches) == 1
            hits.append(matches[0])
        assert sorted(hits) == list(range(29))
        projections = finite['projections']
        degenerate, contributions = degenerate_recipes(projections, k)
        assert len(degenerate) == (18 if k % 2 == 0 else 16)
        complete_recipes = generated + degenerate
        all_finite = [tuple(tuple(p) for p in r['box_vertices'])
                      for r in record['representatives']]
        all_hits = []
        for p in complete_recipes:
            matches = [i for i, q in enumerate(all_finite) if equivalent(p, q)]
            assert len(matches) == 1
            all_hits.append(matches[0])
        assert len(complete_recipes) == len(all_finite)
        assert sorted(all_hits) == list(range(len(all_finite)))
        results.append(dict(k=k, recipe_count=29, finite_full_circuit_count=29,
                            direct_affine_bijection=True,
                            degenerate_recipe_contributions=contributions,
                            complete_recipe_count=len(complete_recipes),
                            finite_tetrahedra_count=len(all_finite),
                            complete_direct_affine_bijection=True))
    return results


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--circuit', type=Path, required=True)
    parser.add_argument('--finite', type=Path, required=True)
    parser.add_argument('--proof', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    circuits = json.loads(args.circuit.read_text())
    finite = json.loads(args.finite.read_text())
    expected = tuple(tuple(tuple(p) for p in projection) for projection in finite['projections'])
    planar = classify_planar(expected)
    print('independent planar classes', planar['direct_affine_classes'], flush=True)
    checked = [check_projection(c) for c in circuits['projections']]
    contributions = [item['orbit_count'] for item in checked]
    assert contributions == [5, 4, 4, 6, 3, 3, 4]
    assert sum(contributions) == circuits['contribution'] == 29
    print('independent residual orbits', contributions, flush=True)
    recipes = recipe_matching(circuits['projections'], finite)
    hashes = {str(path): sha(path) for path in
              (Path(__file__), args.circuit, args.finite)}
    proof_identity = None
    if args.proof is not None:
        proof_text = args.proof.read_text()
        hashes[str(args.proof)] = sha(args.proof)
        proof_identity = dict(path=str(args.proof), characters_read=len(proof_text),
                              role='supplied theoretical bridge; not independently theorem-audited')
    result = dict(GATE_NAME='COMPUTATION_GATE',
                  GATE_SCOPE='finite planar classification, 34 full-circuit residuals, '
                             'seven exact orbit/strip-offset certificates, all thirteen recipes at k=5,6',
                  GATE_STATUS='PASS-WITH-LIMITATIONS',
                  COMPUTATION_VERDICT='EXACT_EXHAUSTIVE',
                  EVIDENCE='independent direct affine maps, Fraction Gaussian solves, '
                           'determinant indices, weighted interval coefficient bounds',
                  LIMITATIONS=['No theorem-wide all-k MATH_GATE is issued.',
                               'The supplied theoretical finite-to-all bridge is outside this computation gate.',
                               'The six degenerate recipes are compared only at k=5,6; '
                               'their all-k equivalence/width proofs remain outside this gate.'],
                  HANDOFF='initiating mathematical-research workflow',
                  software=dict(python=sys.version, executable=sys.executable,
                                platform=platform.platform(), packages='standard library only',
                                arithmetic='integers and fractions.Fraction', processes=1),
                  replay_command=' '.join(sys.argv), hashes=hashes,
                  proof_input=proof_identity, planar_classification=planar,
                  full_circuit_certificates=checked,
                  contribution=29, general_recipe_finite_matches=recipes,
                  warnings=[], failures=[], elapsed_seconds=time.monotonic() - started)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print('finite computation gate PASS-WITH-LIMITATIONS', flush=True)
    print('elapsed', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
