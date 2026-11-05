#!/usr/bin/env python3
"""Finite full-class replay and constructor computation certificate for Q024.

New replay: production classify_diagonal.py for b=3..10, full GL3 class sets
against the complete corresponding local Zenodo records, and full short-
direction multi-width checks for b=3..5.  Hash-valid prior strict ten-case
evidence is reused without claiming to rerun it.  No count formula accepts
or filters a comparison result.
"""
import argparse
import ast
import collections
import datetime
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUTS = ROOT / 'outputs/Q024'
sys.path.insert(0, str(OUTPUTS))
import classify_diagonal as D
import classify_strict as S
import enumerate as E

STRICT_CASES = tuple((b, k) for b in range(3, 7) for k in (b + 1, b + 2)) + ((3, 10), (4, 10))
DIAGONAL_CASES = tuple(range(3, 11))
SOURCE = ROOT / 'work/q024/sources/width_2_tetrahedra_v2.txt'
PRIOR_SCRIPT = ROOT / 'work/q024/check_strict_constructor.py'
PRIOR_OUTPUT = PRIOR_SCRIPT.with_suffix('.json')
PRODUCTION = ('classify_strict.py', 'classify_diagonal.py', 'enumerate.py', 'lift_circuits.py')
PATTERNS = {(0, 0, 0, 2): 'A', (0, 0, 2, 2): 'B',
            (0, 0, 1, 2): 'C', (0, 1, 1, 2): 'D'}
PINNED = {
    'outputs/Q024/classify_strict.py': '2596e13f84d03234bb015feb93ef701e40f52ec3691da1f047bf508919b4d795',
    'outputs/Q024/classify_diagonal.py': '593285d4ca7a0eb74096e8ab2b2af7208fd6124fb88147a3f871fc62a772c1db',
    'outputs/Q024/enumerate.py': '66ea4c0cb1e674d80ce34913039c1708295fd6ea67ed9dc71dd83a6125188adb',
    'outputs/Q024/lift_circuits.py': '64dc5ad9f44a20c82e0443e92efa84120eea67ea0af026e2b18d424afad9915c',
    'work/q024/check_strict_constructor.py': 'da6088f588dfb0bb1add8642c8dccf2cd5a9c96e436287798ae8f94412359026',
    'work/q024/check_strict_constructor.json': '5d97b9b14c40ce6bf632a6d1ce2fc46b00efaab4b974b62cf53f0bc3fbc45ca0',
    'work/q024/sources/width_2_tetrahedra_v2.txt': '3a12fb3f19a3ecdc2bb8500dbf67ab84a615af6f54ce2bd906eb19285023fd0a',
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def key_tuple(value):
    return tuple(tuple(row) for row in value)


def hash_manifest(paths):
    return {relative(path): sha256(path) for path in paths}


def inspect_runtime_dependencies(failures):
    names = set(Path(name).stem for name in PRODUCTION)
    graph = {}
    for name in PRODUCTION:
        path = OUTPUTS / name
        tree = ast.parse(path.read_text())
        imports = set()
        dynamic = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split('.')[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    failures.append(dict(kind='relative_production_import', file=name))
                imports.add((node.module or '').split('.')[0])
            elif isinstance(node, ast.Call):
                func = node.func
                if isinstance(func, ast.Name) and func.id in ('__import__', 'eval', 'exec'):
                    dynamic.append(func.id)
                if isinstance(func, ast.Attribute) and func.attr in ('import_module', 'exec_module'):
                    dynamic.append(func.attr)
        nonstandard = imports - sys.stdlib_module_names
        if nonstandard - names or dynamic:
            failures.append(dict(kind='production_import_closure', file=name,
                                 nonstandard=sorted(nonstandard), dynamic=dynamic))
        graph[name] = dict(imports=sorted(imports),
                           local_source_dependencies=sorted(x + '.py' for x in nonstandard),
                           uninspected_dynamic_loaders=dynamic)
    modules = {'classify_strict': S, 'classify_diagonal': D, 'enumerate': E,
               'lift_circuits': S.C}
    actual = {name: relative(module.__file__) for name, module in modules.items()}
    for name, path in actual.items():
        if path != 'outputs/Q024/' + name + '.py':
            failures.append(dict(kind='runtime_module_location', module=name, actual=path))
    if S.E is not E or D.E is not E or S.C.E is not E:
        failures.append(dict(kind='runtime_dependency_identity'))
    return dict(source_import_graph=graph, actual_project_module_files=actual,
                runtime_requires_work_sources=False,
                expected_standalone_source_files=['outputs/Q024/' + n for n in PRODUCTION])


def function_asts(path):
    return {node.name: ast.dump(node, include_attributes=False)
            for node in ast.parse(Path(path).read_text()).body
            if isinstance(node, ast.FunctionDef)}


def development_provenance(failures):
    """Verify exact transplanted functions; provenance is not a runtime import."""
    target = function_asts(OUTPUTS / 'classify_diagonal.py')
    specifications = (
        ('work/q024/check_ad_diagonal.py', ('triangles', 'automorphisms', 'paired_parity')),
        ('work/q024/check_ad_representatives.py', ('orbit_reps', 'a_reps', 'd_reps')),
        ('work/q024/check_b_diagonal_representatives.py',
         ('objects', 'vertices', 'd4_images', 'd4_key', 'b_reps')),
    )
    comparisons = []
    for source, functions in specifications:
        original = function_asts(ROOT / source)
        for name in functions:
            equal = target[name] == original[name]
            comparisons.append(dict(function=name, development_source=source,
                                    normalized_AST_equal=equal))
            if not equal:
                failures.append(dict(kind='transplanted_function_changed',
                                     function=name, source=source))
    origins = [ROOT / source for source, _ in specifications]
    origins += [ROOT / 'work/q024/c_diagonal_check.py',
                ROOT / 'work/q024/planar_normal_forms.py',
                ROOT / 'work/q024/near_diagonal.py']
    return dict(
        exact_function_transplants=comparisons,
        generation_source_hashes=hash_manifest(origins),
        C_transcription='classify_diagonal.c_reps transcribes the four row-lattice branches of c_diagonal_check.construction; II and III are united by d+f>=b, and full row-HNF modulo J replaces its two-coordinate HNF implementation. This transcription is checked here by full tetrahedron class sets in the stated eight cases, not by an AST identity claim.',
        strict_transcription='classify_strict transcribes the four planar normal forms and finite height quotient/shell construction. The prior ten exact source-set checks are reused only at the pinned unchanged hashes.',
        distinction='These development files document source generation; production constructors do not import or execute them.')


def reuse_strict_evidence(failures):
    prior = json.loads(PRIOR_OUTPUT.read_text())
    current_hashes = {}
    for name, recorded in prior['sha256'].items():
        actual = sha256(ROOT / name)
        current_hashes[name] = actual
        if actual != recorded:
            failures.append(dict(kind='prior_evidence_hash_mismatch', file=name,
                                 recorded=recorded, current=actual))
    rows = prior['results']
    cases = tuple((row['b'], row['k']) for row in rows)
    if cases != STRICT_CASES:
        failures.append(dict(kind='prior_evidence_case_scope', expected=STRICT_CASES, actual=cases))
    for row in rows:
        if not row.get('full_affine_class_sets_equal') or not row.get('all_short_outer_directions_checked'):
            failures.append(dict(kind='prior_evidence_missing_check', case=row))
        if row['count'] != sum(row['by_family'].values()):
            failures.append(dict(kind='prior_evidence_count_integrity', case=row))
    return dict(
        evidence_state='REUSED_HASH_MATCHING_PRIOR_EXECUTED_OUTPUT',
        rerun_in_this_certificate=False,
        reason_for_reuse='All production, checker, dependency and dataset hashes match the already executed ten-case evidence; no change or unresolved failure requires a repeated full run.',
        checker=relative(PRIOR_SCRIPT), checker_sha256=sha256(PRIOR_SCRIPT),
        output=relative(PRIOR_OUTPUT), output_sha256=sha256(PRIOR_OUTPUT),
        recorded_and_current_input_hashes=current_hashes,
        exact_cases=[dict(b=b, k=k) for b, k in STRICT_CASES],
        results=rows,
        checked_properties=['complete full-affine class sets equal the corresponding supplied source sets',
                            'all outer directions of width<k exhaustively checked by exact edge-value divisibility'],
        replay_command=[sys.executable, relative(PRIOR_SCRIPT)],
        original_runtime_metadata='The prior output records hashes and exact checks but does not record its Python version, platform, elapsed time or timestamp. Current runtime metadata below is not relabelled as its historical runtime.',
        bridge='Actual tetrahedra have vertices (p_i,h_i). The checker calls the unchanged production constructor, compares all unordered 24-order affine-GL3 HNF keys with the full corresponding source sets, and enables the exhaustive exact outer-direction check.')


def layer(points):
    directions = list(E.short_covectors(points, 2))
    if len(directions) != 1:
        raise ValueError(f'width<=2 directions not unique up to sign: {directions}')
    values = [sum(a * z for a, z in zip(directions[0], p)) for p in points]
    lo, hi = min(values), max(values)
    if hi - lo != 2:
        raise ValueError(f'first direction has width {hi-lo}')
    pattern = min(tuple(sorted(x - lo for x in values)),
                  tuple(sorted(hi - x for x in values)))
    if pattern not in PATTERNS:
        raise ValueError(f'unexpected first-width layer pattern {pattern}')
    return PATTERNS[pattern]


def source_classes(failures):
    expected = {b: collections.defaultdict(list) for b in DIAGONAL_CASES}
    families = {b: {} for b in DIAGONAL_CASES}
    current = identity = None
    for line in SOURCE.read_text().splitlines():
        if line.startswith('ID: '):
            identity = line[4:].strip()
            current = tuple(int(identity[i:i+2]) for i in (0, 2, 4))
        elif line.startswith('Vertices: ') and current:
            first, b, k = current
            if first != 2 or b != k or b not in expected:
                continue
            points = tuple(tuple(v) for v in ast.literal_eval(line[10:]))
            key = E.canonical_key(points)
            expected[b][key].append(identity)
            try:
                family = layer(points)
            except ValueError as exc:
                failures.append(dict(kind='source_first_layer', source_id=identity, error=str(exc)))
                continue
            if key in families[b] and families[b][key] != family:
                failures.append(dict(kind='source_family_invariance', b=b, key=key))
            families[b][key] = family
    return expected, families


def replay_diagonal(failures):
    expected, source_families = source_classes(failures)
    rows = []
    for b in DIAGONAL_CASES:
        start = time.monotonic()
        try:
            built = D.classify(b, verify=False)
        except Exception as exc:
            failures.append(dict(kind='production_diagonal_exception', b=b,
                                 exception_type=type(exc).__name__, error=str(exc)))
            rows.append(dict(b=b, k=b, completed=False))
            continue
        constructed = {}
        per_family = collections.Counter()
        width_checks = 0
        actual_layer_checks = 0
        matches = []
        for record in built['representatives']:
            points = tuple(tuple(p) for p in record['vertices'])
            key = E.canonical_key(points)
            if key != key_tuple(record['canonical_key']):
                failures.append(dict(kind='stored_production_key', b=b, record=record))
            if key in constructed:
                failures.append(dict(kind='duplicate_full_GL3_class', b=b, key=key))
            constructed[key] = record
            family = record['family']
            per_family[family] += 1
            try:
                actual_family = layer(points)
                actual_layer_checks += 1
                if actual_family != family:
                    failures.append(dict(kind='constructed_layer_bridge', b=b,
                                         declared=family, actual=actual_family, vertices=points))
            except ValueError as exc:
                failures.append(dict(kind='constructed_first_layer', b=b,
                                     family=family, error=str(exc), vertices=points))
            if key in source_families[b] and source_families[b][key] != family:
                failures.append(dict(kind='source_constructed_family', b=b,
                                     key=key, declared=family,
                                     source_family=source_families[b][key]))
            width = None
            if b <= 5:
                width_checks += 1
                try:
                    width = E.multi_width(points, upper=b)
                except ValueError as exc:
                    failures.append(dict(kind='full_3D_width_exception', b=b,
                                         family=family, error=str(exc), vertices=points))
                if width != (2, b, b):
                    failures.append(dict(kind='full_3D_width', b=b, family=family,
                                         width=width, vertices=points))
            matches.append(dict(full_affine_key=key, family=family, vertices=points,
                                source_ids=expected[b].get(key, []),
                                independent_full_multi_width=width))
        built_keys, source_keys = set(constructed), set(expected[b])
        missing, extra = sorted(source_keys - built_keys), sorted(built_keys - source_keys)
        if missing or extra:
            failures.append(dict(kind='diagonal_full_source_set_difference', b=b,
                                 missing=missing, extra=extra))
        if built['count'] != len(constructed) or dict(per_family) != built['by_family']:
            failures.append(dict(kind='production_result_integrity', b=b))
        source_counts = collections.Counter(source_families[b].values())
        per_family_equal = {}
        for family in 'ABCD':
            a = {key for key, rec in constructed.items() if rec['family'] == family}
            z = {key for key, name in source_families[b].items() if name == family}
            per_family_equal[family] = a == z
            if a != z:
                failures.append(dict(kind='diagonal_family_set_difference', b=b,
                                     family=family, missing=sorted(z-a), extra=sorted(a-z)))
        row = dict(b=b, k=b, completed=True, count=len(constructed),
                   production_by_family=dict(per_family), source_by_family=dict(source_counts),
                   source_records=sum(map(len, expected[b].values())),
                   source_unique_classes=len(source_keys),
                   full_affine_class_sets_equal=built_keys == source_keys,
                   family_class_sets_equal=per_family_equal, missing_keys=missing, extra_keys=extra,
                   independently_recomputed_canonical_keys=len(matches),
                   exact_first_layer_checks=actual_layer_checks,
                   full_3D_multi_width_checks=width_checks,
                   matches=sorted(matches, key=lambda item: item['full_affine_key']),
                   elapsed_seconds=time.monotonic()-start)
        rows.append(row)
        print(json.dumps({key:row[key] for key in (
            'b','count','production_by_family','full_affine_class_sets_equal',
            'family_class_sets_equal','full_3D_multi_width_checks')}), flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description="Portable replay of the historical diagonal mathematical comparison only")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    failures = []
    rows = replay_diagonal(failures)
    result = dict(scope="Historical mathematical replay function only; not the original host-dependent provenance gate",
                  results=rows, failures=failures,
                  source_sha256={"work/q024/check_diagonal_constructor.py": sha256(__file__),
                                 "work/q024/sources/width_2_tetrahedra_v2.txt": sha256(SOURCE)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    if failures:
        raise SystemExit(1)
    print("Diagonal mathematical replay: PASS")

if __name__ == "__main__":
    main()
