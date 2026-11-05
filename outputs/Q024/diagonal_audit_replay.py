#!/usr/bin/env python3
"""Fresh independent audit code; Python stdlib only; no project imports.

Finite diagnostics, not a proof for unbounded b.  Row-HNF enumeration is complete
within the elementary determinant bounds stated in the audit report.
"""
from fractions import Fraction
from itertools import permutations, product
from math import gcd
import json
from pathlib import Path
import sys


def egcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


def hnf(rows, chars=(0, 0)):
    (a, b), (c, d) = rows
    g, x, y = egcd(a, c)
    det = a * d - b * c
    h = abs(det) // g
    height = x * b + y * d
    t = height % h
    cv = ((-c // g) * chars[0] + (a // g) * chars[1]) % 2
    ch = (x * chars[0] + y * chars[1] - ((height - t) // h) * cv) % 2
    return g, h, t, ch, cv


def key_b(rows, chars):
    keys = []
    for sa, sb, sw in product((-1, 1), (-1, 1), (0, 1)):
        rr, cc = [], []
        for (a, b), chi in zip(rows, chars):
            u, v = sa * a, sb * b
            if sw:
                u, v = v, u
            rr.append((u, v))
            cc.append((chi + (sa == -1) * a + (sb == -1) * b) % 2)
        keys.append(hnf(rr, cc))
    return min(keys)


def key_c(rows):
    return min(hnf(rows)[:3], hnf(tuple((-a, b + a) for a, b in rows))[:3])


def norm_b(a, b, chi):
    return max(abs(a), abs(b)) + (
        (chi + min(0, b) - min(0, a)) % 2 if abs(a) == abs(b) else 0
    )


def norm_c(d, e):
    return max(abs(d), (abs(e) + 1) // 2, (abs(d + e) + 1) // 2)


def valid_hnf(g, h, t, b, family, chars=(0, 0)):
    first = None
    bound = b if family == 'B' else 2 * b
    for i in range(b // g + 1):
        d = i * g
        lo, hi = -bound, bound
        if family == 'C':
            lo, hi = max(lo, -2 * b - d), min(hi, 2 * b - d)
        start = -((i * t - lo) // h)
        end = (hi - i * t) // h
        for j in range(start, end + 1):
            e = i * t + j * h
            if d == 0 and e <= 0:
                continue
            n = norm_c(d, e) if family == 'C' else norm_b(d, e, (i * chars[0] + j * chars[1]) % 2)
            if n < b:
                return False
            if n == b:
                if first is None:
                    first = (d, e)
                elif first[0] * e != first[1] * d:
                    # Keep looking: a later point can still violate the first minimum.
                    first = (0, 0)
    return first == (0, 0)


def b_candidates(b):
    groups = {k: set() for k in ('I', 'S', 'P', 'E')}
    bases = []
    for a in range(1, b):
        for c in range(a, b):
            bases.append(('I', ((b, a), (-c, b))))
    for c in range(b // 2 + 1):
        bases.append(('S', ((b, 0), (c, b))))
    for d in range(b + 1, 2 * b + 1):
        for a in range(d - b, d // 2 + 1):
            bases.append(('P', ((b, a), (0, d))))
    d = b - 1
    bases.append(('E', ((d, d), (-d, d))))
    for D in range(d + 2, 2 * d + 3):
        ps = {d - (D - d) + 1, d - 1} if D <= 2 * d else ({0, d - 1} if D == 2 * d + 1 else {d - 1})
        for p in ps:
            bases.append(('E', ((d, d), (p, p + D))))
    for name, rows in bases:
        for chars in product(range(2), repeat=2):
            g, h, t, ch, cv = hnf(rows, chars)
            if valid_hnf(g, h, t, b, 'B', (ch, cv)):
                groups[name].add(key_b(rows, chars))
    return groups


def c_candidates(b):
    groups = {k: set() for k in ('G', 'I', 'II', 'III')}
    for h in (2 * b - 1, 2 * b):
        for t in range(h):
            groups['G'].add(key_c(((b, t), (0, h))))
    for h in range(2 * b + 1, 3 * b + 1):
        for a in range(-2 * b, b - h + 1):
            groups['G'].add(key_c(((b, a), (b, a + h))))
    for d in range(1, b):
        for delta in range(2):
            E = 2 * b - d - delta
            for e in range(-2 * b, -d - delta):
                groups['I'].add(key_c(((b, e), (d, E))))
    for d, f, delta, eps in product(range(1, b), range(1, b), range(2), range(2)):
        rows = ((d, 2 * b - d - delta), (f, -2 * b + eps))
        if d + f == b:
            groups['II'].add(key_c(rows))
        if d + f >= b + 1:
            groups['III'].add(key_c(rows))
    return groups


def complete_hnf(b, family):
    cap = (2 if family == 'B' else 4) * b * b
    keys = set()
    for g in range(1, b + 1):
        for h in range(1, cap // g + 1):
            for t in range(h):
                for chars in (product(range(2), repeat=2) if family == 'B' else [(0, 0)]):
                    if valid_hnf(g, h, t, b, family, chars):
                        rows = ((g, t), (0, h))
                        keys.add(key_b(rows, chars) if family == 'B' else key_c(rows))
    return keys


def triangle_normals(a, c):
    if a == c:
        for y in range(a // 2 + 1):
            yield ((0, 0), (a, y), (0, a))
        for x in range(1, a // 2 + 1):
            for y in range(x, a - x + 1):
                yield ((0, 0), (a, y), (x, a))
    else:
        for y in range(a):
            if y <= (c - y) % a:
                yield ((0, 0), (a, y), (0, c))
        for x in range(1, a // 2 + 1):
            for y in range(a - x + 1):
                yield ((0, 0), (a, y), (x, c))
        for x in range(2, (a + 1) // 2):
            for y in range(1, x):
                yield ((0, y), (a, 0), (x, c))


def mv(A, v):
    return tuple(sum(A[i][j] * v[j] for j in range(2)) for i in range(2))


def automorphisms(V):
    E = ((V[1][0] - V[0][0], V[2][0] - V[0][0]), (V[1][1] - V[0][1], V[2][1] - V[0][1]))
    det = E[0][0] * E[1][1] - E[0][1] * E[1][0]
    Ei = ((Fraction(E[1][1], det), Fraction(-E[0][1], det)), (Fraction(-E[1][0], det), Fraction(E[0][0], det)))
    out = []
    for p in permutations(range(3)):
        W = [V[j] for j in p]
        C = ((W[1][0] - W[0][0], W[2][0] - W[0][0]), (W[1][1] - W[0][1], W[2][1] - W[0][1]))
        U = tuple(tuple(sum(C[i][k] * Ei[k][j] for k in range(2)) for j in range(2)) for i in range(2))
        if any(z.denominator != 1 for row in U for z in row):
            continue
        U = tuple(tuple(int(z) for z in row) for row in U)
        t = tuple(W[0][i] - mv(U, V[0])[i] for i in range(2))
        out.append((U, t, p))
    return out


def check_ad(b):
    A = D0 = 0
    for V in triangle_normals(b, b):
        aut = automorphisms(V)
        residues = set(product(range(2), repeat=2))
        while residues:
            r = min(residues)
            residues -= {tuple((mv(U, r)[i] + t[i]) % 2 for i in range(2)) for U, t, p in aut}
            A += 1
        D0 += len({min(p[i] for U, t, p in aut) for i in range(3)})
    D1 = sum(len({tuple(x % 2 for x in v) for v in V}) == 2
             for ac in ((2*b-1, 2*b-1), (2*b-1, 2*b), (2*b, 2*b))
             for V in triangle_normals(*ac))
    assert A == (b*b + b//2 + 3 if b % 2 == 0 else b*b + (b+1)//2)
    assert D0 + D1 == ((6*b*b-3*b+4)//2 if b % 2 == 0 else (6*b*b-3*b+3)//2)
    return {'b': b, 'A': A, 'D_evenq': D0, 'D_oddq': D1, 'D': D0 + D1}


def check_scalar_formulas():
    for d, e in product(range(-30, 31), repeat=2):
        for lab, calc in [('C', norm_c(d, e)), ('D', max(abs(d), (abs(e)+1)//2, (abs(e-2*d)+1)//2))]:
            ans = min(max(V)-min(V) for a in range(-100, 101)
                      for V in [(0, d, a, 2*a+e) if lab == 'C' else (0, a, a+d, 2*a+e)])
            assert ans == calc, (lab, d, e, ans, calc)
    for A, B, chi in product(range(-12, 13), range(-12, 13), range(2)):
        ans = min(max(0, A, 2*a+chi, 2*a+chi+B)-min(0, A, 2*a+chi, 2*a+chi+B) for a in range(-40, 41))
        assert ans == norm_b(A, B, chi)
    return 'C/D d,e=-30..30; B A,B=-12..12, both parities; exact direct integer minimization'


def main():
    maximum = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    evidence = {'python': sys.version, 'scalar_minimization': check_scalar_formulas(),
                'AD': [check_ad(b) for b in range(3, 25)], 'BC': []}
    for b in range(3, maximum + 1):
        row = {'b': b}
        for fam, generate in [('B', b_candidates), ('C', c_candidates)]:
            groups = generate(b)
            union = set().union(*groups.values())
            assert len(union) == sum(map(len, groups.values())), (fam, b, 'overlap')
            fresh = complete_hnf(b, fam)
            assert union == fresh, (fam, b, 'missing', len(fresh-union), 'extra', len(union-fresh))
            expected = ((15*b*b+2*b+24)//8 if b % 2 == 0 else (15*b*b+4*b+5)//8) if fam == 'B' else ((17*b*b-8*b+8)//4 if b % 2 == 0 else (17*b*b-8*b+3)//4)
            assert len(fresh) == expected, (fam, b, len(fresh), expected)
            row[fam] = {'total': len(fresh), 'parts': {k: len(v) for k, v in groups.items()}}
        ad = evidence['AD'][b-3]
        row['total'] = ad['A'] + row['B']['total'] + row['C']['total'] + ad['D']
        assert row['total'] == ((81*b*b-22*b+80)//8 if b % 2 == 0 else (81*b*b-20*b+27)//8)
        evidence['BC'].append(row)
        print(json.dumps(row), flush=True)
    destination = Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(evidence, indent=2) + '\n')
    print('PASS; finite evidence only; output:', destination)


if __name__ == '__main__':
    main()
