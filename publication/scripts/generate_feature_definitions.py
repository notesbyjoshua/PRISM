#!/usr/bin/env python3
import ast
import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FACEKIT = ROOT / "FaceKit"
COMMIT = "e474628f3cb17a83edad59ac6e1d9e5438753032"
SOURCE_PATH = "src/facekit/core/geometric/extractor.py"
SCHEMA_PATH = ROOT / "revision_results/final_20260914/exact_feature_names.csv"
MAPPING_PATH = ROOT / "local_review/bibm2026_evidence/feature_region_family_mapping.csv"
OUTPUT = ROOT / "publication/feature_definitions"
SUPPLEMENT = ROOT / "publication/supplementary/Supplementary_Table_S6_exact_feature_definitions.csv"
METHODS = [
    "_feat_eyebrow", "_feat_eye", "_feat_nose", "_feat_philtrum", "_feat_mouth_lip", "_feat_chin_jaw",
    "_feat_forehead", "_feat_face_global", "_feat_midface_cheek",
]
COLUMNS = [
    "feature", "anatomical_region", "feature_family", "historical_source_file", "historical_source_function",
    "historical_source_line_or_code_location", "measurement_type", "landmark_indices_or_named_points",
    "numerator_definition", "denominator_definition", "full_formula_description", "units",
    "normalization_reference", "left_right_rule", "bilateral_mean_rule", "asymmetry_rule", "area_rule",
    "angle_rule", "notes", "extraction_status",
]


def git_show(path):
    return subprocess.check_output(["git", "-C", str(FACEKIT), "show", f"{COMMIT}:{path}"], text=True)


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as handle: return list(csv.DictReader(handle))


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Recorder(ast.NodeTransformer):
    def __init__(self, source):
        self.source = source
        self.function = ""

    def visit_FunctionDef(self, node):
        previous, self.function = self.function, node.name
        node = self.generic_visit(node)
        self.function = previous
        return node

    def visit_Assign(self, node):
        node = self.generic_visit(node)
        if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Subscript): return node
        target = node.targets[0]
        if not isinstance(target.value, ast.Name) or target.value.id != "f": return node
        expression = ast.unparse(node.value)
        node.value = ast.Call(
            func=ast.Name(id="_record", ctx=ast.Load()),
            args=[target.slice, node.value, ast.Constant(self.function), ast.Constant(node.lineno),
                  ast.Constant(node.end_lineno), ast.Constant(expression), ast.Call(func=ast.Name(id="locals", ctx=ast.Load()), args=[], keywords=[])],
            keywords=[],
        )
        return ast.copy_location(node, node)


def safe_locals(values):
    safe = {}
    for key, value in values.items():
        if isinstance(value, (str, int, bool)): safe[key] = value
        elif isinstance(value, (list, tuple)) and all(isinstance(item, (str, int)) for item in value): safe[key] = list(value)
    return safe


def feature_references(expression, local_values):
    refs = set()
    tree = ast.parse(expression, mode="eval")
    for node in ast.walk(tree):
        if not isinstance(node, ast.Subscript) or not isinstance(node.value, ast.Name) or node.value.id != "f": continue
        try: refs.add(str(eval(compile(ast.Expression(node.slice), "<feature-ref>", "eval"), {}, local_values)))
        except (NameError, TypeError, ValueError): pass
    return sorted(refs)


class Instantiate(ast.NodeTransformer):
    def __init__(self, values): self.values = values

    def visit_Name(self, node):
        aliases = {"bz": "bizygomatic_width", "fh": "face_height", "mfh": "midface_height", "mw": "mouth_width"}
        if node.id in aliases: return ast.copy_location(ast.Name(id=aliases[node.id], ctx=node.ctx), node)
        value = self.values.get(node.id)
        if isinstance(value, (str, int, bool)): return ast.copy_location(ast.Constant(value=value), node)
        if isinstance(value, list) and all(isinstance(item, (str, int)) for item in value):
            return ast.copy_location(ast.List(elts=[ast.Constant(item) for item in value], ctx=ast.Load()), node)
        return node

    def visit_Subscript(self, node):
        node = self.generic_visit(node)
        if isinstance(node.value, ast.Name) and node.value.id == "s" and isinstance(node.slice, ast.Constant):
            aliases = {"bizyg": "bizygomatic_width", "face_h": "face_height", "midface_h": "midface_height", "mouth_w": "mouth_width"}
            if node.slice.value in aliases: return ast.copy_location(ast.Name(id=aliases[node.slice.value], ctx=ast.Load()), node)
        return node


def instantiate(expression, values):
    tree = Instantiate(values).visit(ast.parse(expression, mode="eval"))
    ast.fix_missing_locations(tree)
    return ast.unparse(tree.body)


def assignment_definitions(tree):
    definitions = defaultdict(lambda: defaultdict(list))
    for node in tree.body:
        if not isinstance(node, ast.ClassDef) or node.name != "GeometricFeatureExtractor": continue
        for function in node.body:
            if not isinstance(function, ast.FunctionDef): continue
            for item in ast.walk(function):
                if not isinstance(item, ast.Assign) or len(item.targets) != 1 or not isinstance(item.targets[0], ast.Name): continue
                definitions[function.name][item.targets[0].id].append((item.lineno, ast.unparse(item.value)))
    return definitions


def supporting_definitions(function, line, expression, definitions):
    excluded = {"f", "lm", "s", "np", "float", "abs", "len", "min", "max", "range"}
    pending = [node.id for node in ast.walk(ast.parse(expression, mode="eval")) if isinstance(node, ast.Name) and node.id not in excluded]
    found, seen = [], set()
    while pending:
        name = pending.pop(0)
        if name in seen or name in {"bz", "fh", "mfh", "mw", "bizygomatic_width", "face_height", "midface_height", "mouth_width"}: continue
        seen.add(name)
        candidates = [item for item in definitions[function].get(name, []) if item[0] <= line]
        if not candidates: continue
        def_line, value = candidates[-1]
        found.append((name, def_line, value))
        pending.extend(node.id for node in ast.walk(ast.parse(value, mode="eval")) if isinstance(node, ast.Name) and node.id not in excluded)
    return found


def load_historical_source():
    source = git_show(SOURCE_PATH)
    records = {}

    def record(key, value, function, line, end_line, expression, values):
        local_values = safe_locals(values)
        records.setdefault(key, {
            "function": function, "line": line, "end_line": end_line, "expression": expression,
            "instantiated_expression": instantiate(expression, local_values), "locals": local_values,
            "direct_features": feature_references(expression, local_values),
        })
        return value

    tree = ast.parse(source)
    transformed = Recorder(source).visit(ast.parse(source))
    ast.fix_missing_locations(transformed)
    namespace = {"_record": record, "__name__": "historical_facekit_extractor"}
    exec(compile(transformed, SOURCE_PATH, "exec"), namespace)
    cls = namespace["GeometricFeatureExtractor"]
    rng = np.random.default_rng(7)
    landmarks = rng.uniform(0, 1000, size=(478, 2))
    instance = cls.__new__(cls)
    canonical, scales = instance._canonicalize(landmarks)
    for method in METHODS: getattr(instance, method)(canonical, scales)
    return source, tree, namespace, records


def compute_features(namespace, landmarks):
    cls = namespace["GeometricFeatureExtractor"]
    instance = cls.__new__(cls)
    canonical, scales = instance._canonicalize(landmarks)
    features = {}
    for method in METHODS: features.update(getattr(instance, method)(canonical, scales))
    return features


def landmark_dependencies(namespace, features, seeds=12):
    dependencies = {feature: set() for feature in features}
    for seed in range(seeds):
        rng = np.random.default_rng(1000 + seed)
        landmarks = rng.uniform(100, 900, size=(478, 2))
        baseline = compute_features(namespace, landmarks)
        for index in range(478):
            changed = set()
            for axis in (0, 1):
                perturbed = landmarks.copy()
                perturbed[index, axis] += .01
                result = compute_features(namespace, perturbed)
                for feature in features:
                    a, b = baseline[feature], result[feature]
                    if not np.isclose(a, b, rtol=1e-9, atol=1e-11, equal_nan=True): changed.add(feature)
            for feature in changed: dependencies[feature].add(index)
    return dependencies


def constant_names(namespace):
    names = defaultdict(list)
    for name, value in namespace.items():
        if name.isupper() and isinstance(value, int): names[value].append(name)
    return names


def denominator_parts(expression):
    tree = ast.parse(expression, mode="eval").body
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Div): return ast.unparse(tree.left), ast.unparse(tree.right)
    return ast.unparse(tree), "none beyond any already normalized component features"


def measurement_type(feature, expression):
    if feature.endswith("_mean"): return "bilateral mean"
    if feature.endswith("_asym") or feature.startswith("gaze_asym") or feature == "face_asymmetry": return "asymmetry metric"
    if "angle" in feature or "slant" in feature: return "angle"
    if "area" in feature or "roundness" in feature or "thickness" in feature: return "area ratio/index"
    if any(term in feature for term in ("ratio", "aspect", "fill", "offset", "flare", "curvature", "tenting", "elevation", "drop", "uniformity", "triangularity", "evers", "std", "bulge")):
        return "local or normalized ratio/index"
    if "_dist(" in expression: return "normalized distance"
    return "derived scalar"


def units(feature):
    return "degrees" if "angle" in feature or "slant" in feature else "unitless"


def normalizer(feature, expression, direct_features):
    checks = [
        (r"bizygomatic_width\s*\*\*\s*2|bizygomatic_width \*\* 2", "bizygomatic width squared"),
        (r"\bbizygomatic_width\b", "bizygomatic width"), (r"\bmidface_height\b", "midface height"),
        (r"\bface_height\b", "face height"), (r"\bmouth_width\b", "mouth width"),
    ]
    for pattern, label in checks:
        if re.search(pattern, expression): return label
    if feature == "nose_to_face_area_ratio": return "face-oval polygon area"
    if direct_features: return "inherited from directly referenced component features"
    if any(term in feature for term in ("angle", "slant")): return "none; angular result in degrees"
    return "local denominator or no additional normalization; see exact expression"


def side_rule(feature):
    if feature.endswith("_r"): return "Subject's right side under the historical MediaPipe convention."
    if feature.endswith("_l"): return "Subject's left side under the historical MediaPipe convention."
    return "Not applicable."


def formula_rows(source, tree, namespace, records, dependencies):
    schema = [row["feature"] for row in read_csv(SCHEMA_PATH)]
    mapping = {row["feature"]: row for row in read_csv(MAPPING_PATH)}
    if len(schema) != 125 or len(set(schema)) != 125: raise ValueError("Historical schema must contain 125 unique features")
    if set(records) != set(schema):
        missing, extra = sorted(set(schema) - set(records)), sorted(set(records) - set(schema))
        raise ValueError(f"Historical source/schema mismatch; missing={missing}, extra={extra}")
    names, definitions = constant_names(namespace), assignment_definitions(tree)
    rows = []
    for feature in schema:
        record = records[feature]
        expression = record["instantiated_expression"]
        numerator, denominator = denominator_parts(expression)
        indices = sorted(dependencies[feature])
        named = [f"{index} ({'/'.join(sorted(names[index]))})" if names[index] else str(index) for index in indices]
        support = supporting_definitions(record["function"], record["line"], record["expression"], definitions)
        support_text = "; ".join(f"L{line} {name} = {value}" for name, line, value in support)
        direct = record["direct_features"]
        description = f"Historical expression: {expression}."
        if support_text: description += f" Supporting definitions: {support_text}."
        mean_rule = expression if feature.endswith("_mean") else "Not applicable."
        asym_rule = expression if feature.endswith("_asym") or feature.startswith("gaze_asym") or feature == "face_asymmetry" else "Not applicable."
        area_rule = "Convex-hull or ordered shoelace area exactly as shown in the historical expression and supporting definitions." if any(
            term in description for term in ("_hull_area", "_polygon_area_ordered", "area", "poly_pts")) else "Not applicable."
        angle_rule = "Historical helper `_angle_at` or `degrees(arctan2(...))`, exactly as shown." if any(
            term in description for term in ("_angle_at", "degrees", "arctan2")) else "Not applicable."
        note_parts = ["Landmarks are first centered at the inner-canthal midpoint and roll-corrected by the inner-canthal axis."]
        if direct: note_parts.append("Direct component features: " + ", ".join(direct) + ".")
        rows.append({
            "feature": feature, "anatomical_region": mapping[feature]["anatomical_region"],
            "feature_family": mapping[feature]["feature_family"], "historical_source_file": SOURCE_PATH,
            "historical_source_function": record["function"],
            "historical_source_line_or_code_location": f"{COMMIT}:{SOURCE_PATH}:L{record['line']}-L{record['end_line']}",
            "measurement_type": measurement_type(feature, expression), "landmark_indices_or_named_points": "; ".join(named),
            "numerator_definition": numerator, "denominator_definition": denominator,
            "full_formula_description": description, "units": units(feature),
            "normalization_reference": normalizer(feature, expression, direct), "left_right_rule": side_rule(feature),
            "bilateral_mean_rule": mean_rule, "asymmetry_rule": asym_rule, "area_rule": area_rule,
            "angle_rule": angle_rule, "notes": " ".join(note_parts), "extraction_status": "EXACT_FROM_SOURCE",
        })
    return rows, records


def dependency_outputs(rows, records, dependencies):
    by_feature = {row["feature"]: row for row in rows}
    landmark_rows = []
    for row in rows:
        feature = row["feature"]
        landmark_rows.append({
            "feature": feature, "landmark_count": len(dependencies[feature]),
            "landmark_indices": "; ".join(map(str, sorted(dependencies[feature]))),
            "landmark_indices_or_named_points": row["landmark_indices_or_named_points"],
            "normalization_reference": row["normalization_reference"],
            "direct_component_features": "; ".join(records[feature]["direct_features"]),
        })
    edges = []
    features = [row["feature"] for row in rows]
    for i, source in enumerate(features):
        for target in features[i + 1:]:
            shared = sorted(dependencies[source] & dependencies[target])
            shared_normalizer = by_feature[source]["normalization_reference"] == by_feature[target]["normalization_reference"]
            direct = source in records[target]["direct_features"] or target in records[source]["direct_features"]
            if not shared and not shared_normalizer and not direct: continue
            relationships = []
            if shared: relationships.append("shared_landmarks")
            if shared_normalizer: relationships.append("shared_normalizer")
            if direct: relationships.append("direct_derivation")
            edges.append({
                "source_feature": source, "target_feature": target, "relationship_types": "; ".join(relationships),
                "shared_landmark_count": len(shared), "shared_landmark_indices": "; ".join(map(str, shared)),
                "shared_normalization_reference": by_feature[source]["normalization_reference"] if shared_normalizer else "",
                "direct_derivation": "yes" if direct else "no",
            })
    return landmark_rows, edges


def write_formula_report(rows):
    lines = [
        "# Historical FaceKit feature formulas", "",
        f"Source revision: `{COMMIT}`", "", f"Source file: `{SOURCE_PATH}`", "",
        "All 125 rows were extracted from the exact historical source revision. The table records the historical Python",
        "expression, supporting local definitions, landmarks detected by deterministic perturbation of the exact code,",
        "normalization, units, bilateral rules, and source line location. It does not use the current 120-feature checkout.", "",
        "The canonicalization shared by every feature centers landmarks at the midpoint of inner-canthal landmarks 133 and",
        "362, rotates the coordinates to make that axis horizontal, and then applies the feature-specific denominator shown",
        "in the CSV. A small `1e-12` term in the source protects divisions from a zero denominator and is preserved in the",
        "recorded expression.", "", "## Extraction summary", "",
        f"- Features in corrected schema: {len(rows)}", f"- Unique definitions: {len({row['feature'] for row in rows})}",
        f"- Exact from source: {sum(row['extraction_status'] == 'EXACT_FROM_SOURCE' for row in rows)}",
        "- Newer-only features: 0", "", "## Definitions", "",
        "| Feature | Type | Source | Normalization | Landmarks |", "|---|---|---|---|---:|",
    ]
    for row in rows:
        count = len(row["landmark_indices_or_named_points"].split("; ")) if row["landmark_indices_or_named_points"] else 0
        lines.append(f"| `{row['feature']}` | {row['measurement_type']} | `{row['historical_source_function']}` "
                     f"{row['historical_source_line_or_code_location'].rsplit(':', 1)[-1]} | {row['normalization_reference']} | {count} |")
    lines += ["", "The complete exact expressions and supporting definitions are in `historical_feature_formulas.csv` and",
              "`Supplementary_Table_S6_exact_feature_definitions.csv`."]
    (OUTPUT / "HISTORICAL_FEATURE_FORMULAS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_dependency_report(rows, landmark_rows, edges):
    counts = sorted(((row["feature"], int(row["landmark_count"])) for row in landmark_rows), key=lambda item: (-item[1], item[0]))
    normalizers = Counter(row["normalization_reference"] for row in rows)
    direct = [edge for edge in edges if edge["direct_derivation"] == "yes"]
    overlap = sorted((edge for edge in edges if edge["shared_landmark_count"]), key=lambda edge: (-edge["shared_landmark_count"], edge["source_feature"], edge["target_feature"]))
    lines = [
        "# Historical feature dependency report", "",
        f"This descriptive map uses the exact FaceKit source revision `{COMMIT}`. It complements the empirical reference-",
        "cohort correlation analysis; it does not alter the historical 6,125-test Benjamini–Hochberg family.", "",
        "## Landmark use", "", f"Across 125 features, landmark counts range from {counts[-1][1]} to {counts[0][1]}.",
        "The largest dependency sets are:", "",
    ]
    lines += [f"- `{feature}`: {count} landmarks" for feature, count in counts[:10]]
    lines += ["", "## Shared normalizers", ""]
    lines += [f"- {normalizer}: {count} features" for normalizer, count in normalizers.most_common()]
    lines += ["", "## Direct derivations", "", f"The source contains {len(direct)} pairwise edges in which one feature directly references another.", ""]
    for edge in direct[:30]: lines.append(f"- `{edge['source_feature']}` ↔ `{edge['target_feature']}`")
    lines += ["", "## Highest landmark overlap", "", "| Feature 1 | Feature 2 | Shared landmarks |", "|---|---|---:|"]
    for edge in overlap[:25]: lines.append(f"| `{edge['source_feature']}` | `{edge['target_feature']}` | {edge['shared_landmark_count']} |")
    lines += ["", "Complete pairwise edges are in `feature_dependency_edges.csv`; per-feature landmark sets are in",
              "`feature_landmark_dependencies.csv`."]
    (OUTPUT / "FEATURE_DEPENDENCY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    source, tree, namespace, records = load_historical_source()
    schema = [row["feature"] for row in read_csv(SCHEMA_PATH)]
    dependencies = landmark_dependencies(namespace, schema)
    rows, records = formula_rows(source, tree, namespace, records, dependencies)
    landmark_rows, edges = dependency_outputs(rows, records, dependencies)
    write_csv(OUTPUT / "historical_feature_formulas.csv", rows, COLUMNS)
    write_csv(SUPPLEMENT, rows, COLUMNS)
    write_csv(OUTPUT / "feature_landmark_dependencies.csv", landmark_rows, list(landmark_rows[0]))
    write_csv(OUTPUT / "feature_dependency_edges.csv", edges, list(edges[0]))
    audit = [{
        "schema_position": position, "feature": row["feature"], "present_in_corrected_schema": "yes",
        "captured_from_exact_source": "yes", "source_commit": COMMIT,
        "source_location": row["historical_source_line_or_code_location"], "landmark_dependency_count": len(dependencies[row["feature"]]),
        "newer_only_feature": "no", "extraction_status": row["extraction_status"], "audit_result": "PASS",
    } for position, row in enumerate(rows, 1)]
    write_csv(OUTPUT / "formula_extraction_audit.csv", audit, list(audit[0]))
    write_formula_report(rows)
    write_dependency_report(rows, landmark_rows, edges)
    print(f"Generated exact definitions for {len(rows)} historical features and {len(edges)} dependency edges")


if __name__ == "__main__": main()
