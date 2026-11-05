#!/usr/bin/env python3
"""Construct all affine-unimodular classes with multi-width (2,b,k), k>b>2.

Standard-library exact arithmetic only. Mathematical bridge: Q024-S in
proofs.md, frozen SHA256 8e0de53dd27c185973e0757c0624e508199097d57844555a2c842b07dd5a6ab1.
This script implements the finite height quotient; it does not infer an
infinite theorem from agreement of finite counts. All returned tetrahedra
use the existing full 24-order row-HNF affine class key.

Example: python3 classify_strict.py --b 3 --k 4 --verify-directions --output out.json
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as F
from math import ceil, floor
from pathlib import Path

import enumerate as E
import lift_circuits as C


def normal_forms(b):
    groups = {name: [] for name in 'ABCD'}
    for m in range(b//2+1):
        for t in (0,1):
            groups['A'].append(((0,0),(0,m),(0,b),(2,t)))
    for w in range(b):
        for t in (0,1):
            if t <= (b-w-t) % 2:
                groups['B'].append(((0,0),(0,b),(2,t),(2,t+w)))
    groups['B'] += [((0,0),(0,b),(2,0),(2,b)),
                    ((0,0),(0,b-1),(2,1),(2,b))]
    for d in range(b+1):
        for e in (range(-(b//2),b+1) if d == b
                  else (2*b-1-d,2*b-d)):
            groups['C'].append(((0,0),(0,d),(1,0),(2,e)))
    for d in range(b+1):
        for e in (range(b,2*b+1) if d == b else (2*b-1,2*b)):
            groups['D'].append(((0,0),(1,0),(1,d),(2,e)))
    return groups


def strip_leq(p, h, bound):
    """Exist integers a,c with range(h_i+a*x_i+c*y_i)<=bound?

    Pairwise differences give necessary and sufficient linear inequalities.
    Projecting intersections of any two independent edge strips bounds c;
    vertical edges can provide the bound directly. All improved integer
    solutions are exhaustively included; no floating approximation is used.
    """
    edges = [(p[i][0]-p[j][0],p[i][1]-p[j][1],h[i]-h[j])
             for i,j in itertools.combinations(range(4),2)]
    vertical = [(min(F(-bound-z,y),F(bound-z,y)),max(F(-bound-z,y),F(bound-z,y)))
                for x,y,z in edges if x == 0 and y]
    if vertical:
        lo,hi = max(x for x,y in vertical),min(y for x,y in vertical)
    else:
        ranges = []
        for (x,y,z),(xx,yy,zz) in itertools.combinations(edges,2):
            det = x*yy-xx*y
            if not det:
                continue
            vals = [F(-xx*a+x*c,det) for a in (-bound-z,bound-z)
                    for c in (-bound-zz,bound-zz)]
            ranges.append((min(vals),max(vals)))
        lo,hi = max(x for x,y in ranges),min(y for x,y in ranges)
    for c in range(ceil(lo),floor(hi)+1):
        lower,upper = None,None
        feasible = True
        for x,y,z in edges:
            if x == 0:
                if abs(z+c*y) > bound:
                    feasible = False
                    break
                continue
            if x < 0:
                x,y,z = -x,-y,-z
            left,right = ceil(F(-bound-z-c*y,x)),floor(F(bound-z-c*y,x))
            lower = left if lower is None else max(lower,left)
            upper = right if upper is None else min(upper,right)
        if feasible and lower <= upper:
            return True
    return False


def orbit_reps(actions, size):
    remaining = set(range(size))
    while remaining:
        start = min(remaining)
        orbit = {a[start] for a in actions}
        assert orbit <= remaining
        remaining -= orbit
        yield start


def ordinary_lifts(p, k):
    lam,_ = C.relation(p)
    s = sum(x for x in lam if x > 0)
    delta = tuple(int(x > 0) for x in lam)
    edge = E.edge_matrix(p)
    hnf = E.row_hnf(edge+(tuple(delta[i]-delta[0] for i in range(1,4)),))
    keys = tuple(itertools.product(*(range(hnf[i][i]) for i in range(3))))
    lookup = {key:i for i,key in enumerate(keys)}
    residues = []
    for key in keys:
        r = (0,)+key
        q = sum(l*x for l,x in zip(lam,r))//s
        residues.append(tuple(x-q*d for x,d in zip(r,delta)))
    actions = []
    for perm in E.PERMS:
        if E.row_hnf(E.edge_matrix(p,perm)) != E.row_hnf(edge):
            continue
        j = next(i for i,l in enumerate(lam) if l)
        sign = lam[perm[j]]//lam[j]
        image = []
        for r in residues:
            target = tuple(sign*r[j] for j in perm)
            key = C.reduce_coset(tuple(target[i]-target[0] for i in range(1,4)),hnf)
            image.append(lookup[key])
        actions.append(image)
    for index in orbit_reps(actions,len(residues)):
        r = residues[index]
        t = sum(l*x for l,x in zip(lam,r))
        for q in range(k+1):
            if not 0 < s*q+t <= s*k:
                continue
            h = tuple(q*d+x for d,x in zip(delta,r))
            if strip_leq(p,h,k) and not strip_leq(p,h,k-1):
                yield h


def duplicate_lifts(p,k):
    pair = next((i,j) for i,j in itertools.combinations(range(4),2) if p[i] == p[j])
    others = [i for i in range(4) if i not in pair]
    pt = p[pair[0]]
    lattice = tuple(tuple(p[j][axis]-pt[axis] for j in others) for axis in range(2))
    hnf = E.row_hnf(lattice)
    reps = tuple(itertools.product(range(hnf[0][0]),range(hnf[1][1])))
    lookup = {r:i for i,r in enumerate(reps)}
    swap = E.row_hnf(tuple(tuple(reversed(row)) for row in lattice)) == hnf
    for d in range(1,k+1):
        actions = []
        for flip in (False,True):
            for exchange in ((False,True) if swap else (False,)):
                image = []
                for a,b in reps:
                    if exchange:
                        a,b = b,a
                    if flip:
                        a,b = d-a,d-b
                    image.append(lookup[C.reduce_coset((a,b),hnf)])
                actions.append(image)
        for index in orbit_reps(actions,len(reps)):
            a,c = reps[index]
            h = [0]*4
            h[pair[1]],h[others[0]],h[others[1]] = d,a,c
            if strip_leq(p,h,k) and not strip_leq(p,h,k-1):
                yield tuple(h)


def expected(b,k):
    const = (76 if k%2 == 0 else 56) if b%2 == 0 else (37 if k%2 == 0 else 25)
    return (81*b*b-18*b+const)//4


def classify(b,k,verify=False):
    if not k > b > 2:
        raise ValueError('This constructor requires k>b>2; diagonal is a separate problem.')
    seen = set()
    result = []
    counts = {}
    for family,projections in normal_forms(b).items():
        before = len(result)
        for p in projections:
            lifts = duplicate_lifts(p,k) if len(set(p)) < 4 else ordinary_lifts(p,k)
            for h in lifts:
                vertices = tuple(point+(z,) for point,z in zip(p,h))
                if verify:
                    assert not any(u[2] for u in E.short_covectors(vertices,k-1))
                key = E.canonical_key(vertices)
                assert key not in seen, ('duplicate full affine class',b,k,family,key)
                seen.add(key)
                result.append(dict(family=family, vertices=vertices, canonical_key=key))
        counts[family] = len(result)-before
    assert len(result) == expected(b,k), (b,k,counts,len(result),expected(b,k))
    return dict(b=b,k=k,count=len(result),by_family=counts,
                classification_scope='Q024-S; all k>b>2 by the cited proof',
                finite_run_scope=dict(b=b,k=k),
                all_short_outer_directions_checked=verify,
                representatives=result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b',type=int,required=True)
    parser.add_argument('--k',type=int,required=True)
    parser.add_argument('--verify-directions',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = classify(args.b,args.k,args.verify_directions)
    folder = Path(__file__).resolve().parent
    result['code_sha256'] = {name:hashlib.sha256((folder/name).read_bytes()).hexdigest()
                             for name in ('classify_strict.py','enumerate.py','lift_circuits.py')}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('b','k','count','by_family','all_short_outer_directions_checked')}))


if __name__ == '__main__':
    main()
