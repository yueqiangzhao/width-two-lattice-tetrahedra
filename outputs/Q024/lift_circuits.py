#!/usr/bin/env python3
"""Exact finite certificates for the seven full-circuit (2,2) projections.

The mathematics justifying finite-to-all lifting is in paper/sections/boundary.tex.
No sampling, fitting, floating-point arithmetic, or external modules are used.
"""
import argparse
import itertools
import json
from fractions import Fraction
from functools import reduce
from math import gcd, ceil, floor
from pathlib import Path
import enumerate as E


def relation(points):
    a = ((1,)*4, tuple(p[0] for p in points), tuple(p[1] for p in points))
    c = [(-1)**i * E.determinant(tuple(tuple(row[j] for j in range(4) if j != i)
                                          for row in a)) for i in range(4)]
    d = abs(reduce(gcd, c))
    c = tuple(v//d for v in c)
    if next(v for v in c if v) < 0:
        c = tuple(-v for v in c)
    return c, d


def reduce_coset(v, hnf):
    v = list(v)
    for i, row in enumerate(hnf):
        q = v[i] // row[i]
        v = [x-q*y for x,y in zip(v,row)]
    return tuple(v)


def polygon_vertices(inequalities):
    """Vertices of bounded {a*x+b*y<=c}; coefficients/rhs are exact."""
    vertices = set()
    for (a,b,c),(d,e,f) in itertools.combinations(inequalities,2):
        determinant = a*e-b*d
        if not determinant:
            continue
        x, y = Fraction(c*e-b*f,determinant), Fraction(a*f-c*d,determinant)
        if all(A*x+B*y <= C for A,B,C in inequalities):
            vertices.add((x,y))
    if not vertices:
        raise AssertionError('no polygon vertices')
    return tuple(sorted(vertices))


def integer_polygon_points(inequalities):
    vertices = polygon_vertices(inequalities)
    limits = (ceil(min(x for x,y in vertices)), floor(max(x for x,y in vertices)),
              ceil(min(y for x,y in vertices)), floor(max(y for x,y in vertices)))
    lo_x,hi_x,lo_y,hi_y = limits
    points = tuple((a,b) for a in range(lo_x,hi_x+1) for b in range(lo_y,hi_y+1)
                   if all(A*a+B*b <= C for A,B,C in inequalities))
    return points, limits, vertices


def oscillation(heights):
    return max(heights)-min(heights)


def strip_width(points, heights):
    """All integer shears potentially improving the zero-shear upper bound."""
    upper = oscillation(heights)
    inequalities = [(points[i][0]-points[j][0],points[i][1]-points[j][1],
                     upper-heights[i]+heights[j])
                    for i in range(4) for j in range(4) if i!=j]
    candidates,_,_ = integer_polygon_points(inequalities)
    return min(oscillation(tuple(h+a*x+b*y for h,(x,y) in zip(heights,points)))
               for a,b in candidates)


def eventual_width(points, c, residual):
    positive = [i for i,v in enumerate(c) if v>0]
    negative = [i for i,v in enumerate(c) if v<0]
    upper = max(residual[i] for i in positive)-min(residual[j] for j in negative)
    inequalities = [(points[i][0]-points[j][0],points[i][1]-points[j][1],
                     upper-residual[i]+residual[j])
                    for i in positive for j in negative]
    candidates,limits,vertices = integer_polygon_points(inequalities)
    scored = []
    for a,b in candidates:
        h = tuple(r+a*x+b*y for r,(x,y) in zip(residual,points))
        g = max(h[i] for i in positive)-min(h[j] for j in negative)
        q0 = max(0,max(h[j] for j in negative)-max(h[i] for i in positive),
                 min(h[j] for j in negative)-min(h[i] for i in positive))
        scored.append((g,q0,a,b))
    C = min(row[0] for row in scored)
    chosen = min(row for row in scored if row[0]==C)
    _,q0,a,b = chosen
    v = tuple(int(x>0) for x in c)
    small = []
    for q in range(q0):
        h = tuple(q*x+r for x,r in zip(v,residual))
        if sum(x*y for x,y in zip(c,h))==0:
            small.append(dict(q=q,degenerate=True))
        else:
            small.append(dict(q=q,width=strip_width(points,h)))
    # For every k>=5 there is exactly one positive-volume q with width k.
    assert q0+C <= 5, (points,residual,C,q0)
    assert all(s.get('degenerate') or s['width']<5 for s in small)
    return dict(C=C,q0=q0,minimizer=(a,b),shear_box=limits,
                polygon_vertices=[[str(x),str(y)] for x,y in vertices],
                candidates=len(candidates),small_exceptions=small)


def projection_certificate(index,points):
    c,D = relation(points)
    assert all(c)
    s = sum(x for x in c if x>0)
    v = tuple(int(x>0) for x in c)
    edge = E.edge_matrix(points)
    basis = edge+(tuple(v[i]-v[0] for i in range(1,4)),)
    hnf = E.row_hnf(basis)
    assert abs(E.determinant(basis))==D*s
    residues = []
    keys = tuple(itertools.product(*(range(hnf[i][i]) for i in range(3))))
    for key in keys:
        r = (0,)+key
        t = sum(x*y for x,y in zip(c,r))
        shift = t//s
        r = tuple(x-shift*y for x,y in zip(r,v))
        assert 0<=sum(x*y for x,y in zip(c,r))<s
        residues.append(r)
    lookup = {key:i for i,key in enumerate(keys)}
    automorphisms = [pi for pi in E.PERMS
                     if E.row_hnf(E.edge_matrix(points,pi))==E.row_hnf(edge)]
    actions = []
    for pi in automorphisms:
        i = next(i for i in range(4) if c[i])
        sign = c[pi[i]]//c[i]
        assert sign in (-1,1) and all(c[pi[j]]==sign*c[j] for j in range(4))
        action = []
        for r in residues:
            image = tuple(sign*r[j] for j in pi)
            image = tuple(image[j]-image[0] for j in range(1,4))
            action.append(lookup[reduce_coset(image,hnf)])
        actions.append(dict(permutation=pi,orientation=sign,residual_action=action))
    remaining = set(range(len(residues)))
    orbits = []
    while remaining:
        i = min(remaining)
        orbit = {a['residual_action'][i] for a in actions}
        assert orbit<=remaining
        remaining-=orbit
        orbits.append(sorted(orbit))
    widths = [eventual_width(points,c,r) for r in residues]
    for orbit in orbits:
        assert len({widths[i]['C'] for i in orbit})==1
    return dict(projection_index=index,points=points,relation=c,torsion_order=D,
                positive_sum=s,residue_hnf=hnf,residuals=residues,
                automorphisms=actions,residual_orbits=orbits,width_certificates=widths,
                contribution_for_every_k_ge_5=len(orbits))


def duplicate_lifts(points,k):
    pair=next((i,j) for i,j in itertools.combinations(range(4),2) if points[i]==points[j])
    other=[i for i in range(4) if i not in pair]
    p=points[pair[0]]
    lattice=tuple(tuple(points[j][axis]-p[axis] for j in other) for axis in range(2))
    hnf=E.row_hnf(lattice)
    cosets=tuple(itertools.product(range(hnf[0][0]),range(hnf[1][1])))
    # The only possible nontrivial projected automorphism fixes p and swaps q,r.
    swapped=tuple(tuple(reversed(row)) for row in lattice)
    swap=E.row_hnf(swapped)==hnf
    remaining=set(cosets)
    lifts=[]
    while remaining:
        s,t=min(remaining)
        orbit={reduce_coset((s,t),hnf),reduce_coset((k-s,k-t),hnf)}
        if swap:
            orbit|={reduce_coset((t,s),hnf),reduce_coset((k-t,k-s),hnf)}
        assert orbit<=remaining
        remaining-=orbit
        h=[0]*4
        h[pair[1]]=k
        h[other[0]],h[other[1]]=s,t
        lifts.append(tuple(h))
    return lifts


def representative_recipes(k):
    """One representative per class for all k>=5, using the proved reduction."""
    if k<5:
        raise ValueError('the structural recipe is for k>=5')
    projections,_=E.planar_representatives()
    records=[]
    for i,p in enumerate(projections,1):
        if i<=4:
            heights=duplicate_lifts(p,k)
        elif i<=6:
            heights=[(0,0,d,t) for d in (2*k-1,2*k)
                     for t in ((0,1) if i==5 or d%2==0 else (0,))]
        else:
            cert=projection_certificate(i,p)
            v=tuple(int(x>0) for x in cert['relation'])
            heights=[]
            for orbit in cert['residual_orbits']:
                j=min(orbit)
                q=k-cert['width_certificates'][j]['C']
                heights.append(tuple(q*x+r for x,r in zip(v,cert['residuals'][j])))
        for h in heights:
            records.append(dict(projection_index=i,heights=h,
                                vertices=tuple((x,y,z) for (x,y),z in zip(p,h))))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    projections,_=E.planar_representatives()
    certs=[projection_certificate(i+1,p) for i,p in enumerate(projections) if all(relation(p)[0])]
    result=dict(evidence='exact finite certificates with theoretical finite-to-all bridge in paper/sections/boundary.tex',
                scope='seven full-circuit projection types; k>=5',projections=certs,
                contribution=sum(c['contribution_for_every_k_ge_5'] for c in certs))
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print([(c['projection_index'],len(c['residuals']),len(c['automorphisms']),
            c['contribution_for_every_k_ge_5']) for c in certs])
    print('full-circuit contribution',result['contribution'])


if __name__=='__main__':
    main()
