#!/usr/bin/env python3
"""Exact finite source-set and all-direction check; no all-parameter inference."""
import ast
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'outputs/Q024'))
import classify_strict as S
import enumerate as E

cases=tuple((b,k) for b in range(3,7) for k in (b+1,b+2))+((3,10),(4,10))
source=ROOT/'work/q024/sources/width_2_tetrahedra_v2.txt'
expected={case:set() for case in cases}
current=None
for line in source.read_text().splitlines():
    if line.startswith('ID: '):
        identity=line[4:]
        current=(int(identity[:2]),int(identity[2:4]),int(identity[4:6]))
    elif line.startswith('Vertices: ') and current is not None:
        first,b,k=current
        if first==2 and (b,k) in expected:
            p=tuple(tuple(v) for v in ast.literal_eval(line[10:]))
            expected[b,k].add(E.canonical_key(p))
rows=[]
for b,k in cases:
    result=S.classify(b,k,verify=True)
    keys={tuple(tuple(row) for row in p['canonical_key']) for p in result['representatives']}
    assert keys==expected[b,k],(b,k,len(keys),len(expected[b,k]))
    row={key:result[key] for key in ('b','k','count','by_family')}
    row.update(full_affine_class_sets_equal=True,all_short_outer_directions_checked=True)
    rows.append(row)
    print(json.dumps(row),flush=True)
files=[ROOT/'outputs/Q024'/n for n in ('classify_strict.py','enumerate.py','lift_circuits.py')]
files += [Path(__file__),source]
Path(__file__).with_suffix('.json').write_text(json.dumps(dict(
    scope='exact finite cases displayed; source agreement is not an infinite proof',
    sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    results=rows),indent=2)+'\n')
