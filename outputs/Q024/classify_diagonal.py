#!/usr/bin/env python3
"""Canonical tetrahedron representatives of multi-width (2,b,b), b>2.

Exact standard-library implementation of the A/B/C/D diagonal constructions.
Hamm triangle forms supply A/D. B uses twisted D4 character orbits; C uses
four explicit row-lattice classes modulo J. All final representatives are
compared by the complete 24-order affine-unimodular row-HNF key. The count
assertions are consistency checks, not evidence for an infinite theorem.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import enumerate as E


def triangles(a, c):
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

def automorphisms(p):
    a, c = (p[1][i] - p[0][i] for i in range(2))
    b, d = (p[2][i] - p[0][i] for i in range(2))
    det = a * d - b * c
    for perm in itertools.permutations(range(3)):
        f = tuple(p[perm[1]][i] - p[perm[0]][i] for i in range(2))
        g = tuple(p[perm[2]][i] - p[perm[0]][i] for i in range(2))
        nums = ((f[0]*d-g[0]*c, -f[0]*b+g[0]*a),
                (f[1]*d-g[1]*c, -f[1]*b+g[1]*a))
        if any(x % det for row in nums for x in row):
            continue
        u = tuple(tuple(x // det for x in row) for row in nums)
        if abs(u[0][0]*u[1][1]-u[0][1]*u[1][0]) != 1:
            continue
        t = tuple(p[perm[0]][i] - sum(u[i][j]*p[0][j] for j in range(2))
                  for i in range(2))
        yield perm, u, t

def paired_parity(p):
    return len({tuple(x % 2 for x in v) for v in p}) == 2

def orbit_reps(elements, maps):
    elements = tuple(elements)  # Materialize a generator before repeated membership checks.
    remaining = set(elements)
    while remaining:
        x = min(remaining)
        orbit = {f(x) for f in maps}
        assert orbit <= set(elements)
        yield min(orbit)
        remaining -= orbit

def a_reps(b):
    for p in triangles(b, b):
        autos = list(automorphisms(p))
        maps = [lambda x, u=u, t=t: tuple(
            (sum(u[i][j]*x[j] for j in range(2))+t[i]) % 2
            for i in range(2)) for _, u, t in autos]
        for t in orbit_reps(itertools.product(range(2), repeat=2), maps):
            yield tuple((0, *v) for v in p) + ((2, *t),)

def d_reps(b):
    for p in triangles(b, b):
        maps = [lambda x, perm=perm: perm[x] for perm, _, _ in automorphisms(p)]
        for mark in orbit_reps(range(3), maps):
            i, j = (x for x in range(3) if x != mark)
            v = tuple(p[j][c]-p[i][c] for c in range(2))
            q = tuple(2*(p[mark][c]-p[i][c]) for c in range(2))
            yield ((0,0,0),(1,0,0),(1,*v),(2,*q))
    for a, c in ((2*b-1,2*b-1),(2*b-1,2*b),(2*b,2*b)):
        for p in triangles(a, c):
            if not paired_parity(p):
                continue
            residues = [tuple(x % 2 for x in z) for z in p]
            mark = next(i for i in range(3) if residues.count(residues[i]) == 1)
            i, j = (x for x in range(3) if x != mark)
            assert all((p[j][c]-p[i][c]) % 2 == 0 for c in range(2))
            v = tuple((p[j][c]-p[i][c])//2 for c in range(2))
            q = tuple(p[mark][c]-p[i][c] for c in range(2))
            yield ((0,0,0),(1,0,0),(1,*v),(2,*q))

def objects(b):
    """All displayed bases and all four characters, with explicit paper filter."""
    for a in range(1, b):
        for c in range(a, b):
            for s, t in itertools.product(range(2), repeat=2):
                yield dict(family='I', parameters=dict(a=a, c=c),
                           basis=((b, a), (-c, b)), chi=(s, t),
                           allowed=True, reason='all I characters')
    for c in range(b // 2 + 1):
        for s, t in itertools.product(range(2), repeat=2):
            yield dict(family='S', parameters=dict(c=c),
                       basis=((b, 0), (c, b)), chi=(s, t),
                       allowed=True, reason='all S characters')
    for d in range(b + 1, 2 * b + 1):
        for a in range(d - b, d // 2 + 1):
            for s, t in itertools.product(range(2), repeat=2):
                if d == 2 * b:
                    allowed = s == 0 and t == b % 2
                    reason = 'P diamond: both corners good'
                elif a == d - b:
                    allowed = (s + t - b) % 2 == 0
                    reason = 'P boundary: anti-corner good'
                else:
                    allowed = True
                    reason = 'P has two non-corner contacts'
                yield dict(family='P', parameters=dict(d=d, a=a),
                           basis=((b, a), (0, d)), chi=(s, t),
                           allowed=allowed, reason=reason)
    d = b - 1
    for s, t in itertools.product(range(2), repeat=2):
        yield dict(family='E_double', parameters=dict(d=d),
                   basis=((d, d), (-d, d)), chi=(s, t),
                   allowed=s == 1 and t == (1 + d) % 2,
                   reason='both length b-1 corners bad')
    for D in range(d + 2, 2 * d + 3):
        if D <= 2 * d:
            ps = {2 * d - D + 1, d - 1}
        elif D == 2 * d + 1:
            ps = {0, d - 1}
        else:
            ps = {d - 1}
        for p in sorted(ps):
            for s, t in itertools.product(range(2), repeat=2):
                allowed = s == 1
                reason = 'unique length b-1 corner bad'
                if D == 2 * d + 2:
                    allowed = allowed and t == b % 2
                    reason += '; independent length b anti-corner good'
                yield dict(family='E_single', parameters=dict(d=d, D=D, p=p),
                           basis=((d, d), (p, p + D)), chi=(s, t),
                           allowed=allowed, reason=reason)

def vertices(basis, chi):
    p, q = basis
    v = (p[0], q[0])
    w = (p[1], q[1])
    r = tuple(int(x % 2) for x in chi)
    return ((0, 0, 0), (0, *v), (2, *r),
            (2, r[0] + w[0], r[1] + w[1]))

def d4_images(basis, chi):
    """The exact affine character action of all eight endpoint/layer symmetries."""
    rows = tuple((a, b, c) for (a, b), c in zip(basis, chi))
    for sa in (-1, 1):
        for sb in (-1, 1):
            for swap in (False, True):
                result = []
                for a, b, c in rows:
                    aa, bb = sa * a, sb * b
                    if swap:
                        aa, bb = bb, aa
                    cc = (c + (sa == -1) * a + (sb == -1) * b) % 2
                    result.append((aa, bb, cc))
                yield (sa, sb, swap), tuple(result)

def d4_key(basis, chi):
    return min(E.row_hnf(rows + ((0, 0, 2),))
               for _, rows in d4_images(basis, chi))

def b_reps(b):
    """One actual tetrahedron per displayed twisted-D4 orbit, for any b>2."""
    if b <= 2:
        raise ValueError('the B diagonal classification requires b>2')
    seen = set()
    for item in objects(b):
        if not item['allowed']:
            continue
        key = d4_key(item['basis'], item['chi'])
        if key not in seen:
            seen.add(key)
            yield vertices(item['basis'], item['chi'])

def c_reps(b):
    seen=set()
    def add(v,w):
        ordinary=E.row_hnf((v,w))
        reflected=E.row_hnf(tuple((-d,e+d) for d,e in (v,w)))
        seen.add(min(ordinary,reflected))
    for h in (2*b-1,2*b):
        for t in range(h):
            add((b,t),(0,h))
    for h in range(2*b+1,3*b+1):
        for a in range(-2*b,b-h+1):
            add((b,a),(b,a+h))
    for d in range(1,b):
        for delta in (0,1):
            top=(d,2*b-d-delta)
            for e in range(-2*b,-d-delta):
                add((b,e),top)
            for f in range(1,b):
                for epsilon in (0,1):
                    if d+f>=b:
                        add(top,(f,-2*b+epsilon))
    for (d1,e1),(d2,e2) in sorted(seen):
        yield ((0,0,0),(0,d1,d2),(1,0,0),(2,e1,e2))


def expected(b):
    return (81*b*b-(22 if b%2==0 else 20)*b+(80 if b%2==0 else 27))//8


def classify(b,verify=False):
    if b<=2:
        raise ValueError('Use enumerate.py for the separately established b=2 branch.')
    seen=set()
    reps=[]
    counts={}
    for family,build in (('A',a_reps),('B',b_reps),('C',c_reps),('D',d_reps)):
        before=len(reps)
        for p in build(b):
            if verify:
                assert E.multi_width(p)==(2,b,b),(b,family,p,E.multi_width(p))
            key=E.canonical_key(p)
            assert key not in seen,('duplicate full affine class',b,family,key)
            seen.add(key)
            reps.append(dict(family=family,vertices=p,canonical_key=key))
        counts[family]=len(reps)-before
    assert len(reps)==expected(b),(b,counts,len(reps),expected(b))
    return dict(b=b,k=b,count=len(reps),by_family=counts,
                finite_run_scope=dict(b=b),
                all_3D_widths_independently_checked=verify,
                representatives=reps)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b',type=int,required=True)
    parser.add_argument('--verify-directions',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=classify(args.b,args.verify_directions)
    folder=Path(__file__).resolve().parent
    result['code_sha256']={n:hashlib.sha256((folder/n).read_bytes()).hexdigest()
                           for n in ('classify_diagonal.py','enumerate.py')}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('b','k','count','by_family','all_3D_widths_independently_checked')}))


if __name__=='__main__':
    main()

