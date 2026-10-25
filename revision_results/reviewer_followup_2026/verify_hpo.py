#!/usr/bin/env python3
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import urllib.request

import pandas as pd

OUT = Path(__file__).resolve().parent
URL = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ontology-file', type=Path, help='Verify a retained public ontology snapshot without network access')
    args = parser.parse_args()
    verified_at = datetime.now(timezone.utc)
    if args.ontology_file:
        payload = args.ontology_file.read_bytes()
    else:
        with urllib.request.urlopen(URL, timeout=60) as response:
            payload = response.read()
    (OUT/'local_only').mkdir(exist_ok=True)
    (OUT/'local_only/hp.obo').write_bytes(payload)
    text = payload.decode()
    terms = {}
    for stanza in text.split('[Term]')[1:]:
        id_match = re.search(r'^id: (HP:\d+)$', stanza, re.M)
        name_match = re.search(r'^name: (.+)$', stanza, re.M)
        if id_match and name_match:
            definition = re.search(r'^def: "(.*?)"', stanza, re.M)
            terms[id_match.group(1)] = {'label': name_match.group(1), 'definition': definition.group(1) if definition else '',
                                      'obsolete': 'is_obsolete: true' in stanza}
    proposals = [
        ('philtrum_length', 'higher', 'HP:0000343', 'Long philtrum', 'normalized length uses midface height; no absolute or age-specific clinical threshold'),
        ('philtrum_length', 'lower', 'HP:0000322', 'Short philtrum', 'same ratio can change through denominator; requires clinical review'),
        ('inter_pupillary_distance', 'higher', 'HP:0000316', 'Hypertelorism', 'IPD/bizygomatic ratio does not establish bony orbital spacing or age-standardized hypertelorism'),
        ('inter_pupillary_distance', 'lower', 'HP:0000601', 'Hypotelorism', 'lower normalized IPD is not sufficient for clinical hypotelorism'),
        ('eb_thickness_mean', 'higher', 'HP:0000574', 'Thick eyebrow', 'convex-hull area/bizygomatic width squared; does not measure hair density or diameter'),
        ('eb_thickness_r;eb_thickness_l', 'higher', 'HP:0000574', 'Thick eyebrow', 'bilateral landmark polygon area differs from clinical hair-density definition'),
        ('malar_bulge_mean', 'lower', 'HP:0000272', 'Malar flattening', '2D contour residual, not measured bone projection; no one-to-one equivalence'),
        ('malar_bulge_mean', 'higher', 'HP:0010620', 'Malar prominence', '2D contour residual cannot establish clinical malar prominence')]
    rows = []
    for feature, direction, identifier, expected, limitation in proposals:
        term = terms.get(identifier)
        verified = term is not None and not term['obsolete'] and term['label'].lower() == expected.lower()
        rows.append({'facekit_feature': feature, 'candidate_direction': direction, 'hpo_id': identifier,
                     'expected_label': expected, 'verified_label': term['label'] if term else '',
                     'id_label_verified': verified, 'ontology_definition': term['definition'] if term else '',
                     'evidence_url': f'https://hpo.jax.org/browse/term/{identifier}', 'ontology_source': URL,
                     'metric_limitation': limitation, 'mapping_status': 'candidate_requires_clinician_approval',
                     'relationship': 'uncertain; many-to-many; not a validated phenotype mapping'})
    pd.DataFrame(rows).to_csv(OUT/'candidate_hpo_mapping.csv', index=False)
    version = re.search(r'^data-version: (.+)$', text, re.M)
    (OUT/'hpo_source_provenance.json').write_text(json.dumps({'url': URL, 'sha256': hashlib.sha256(payload).hexdigest(),
        'data_version': version.group(1) if version else None, 'verified_candidates': sum(r['id_label_verified'] for r in rows),
        'candidate_rows': len(rows), 'verified_utc': verified_at.isoformat(),
        'retrieved_date': verified_at.date().isoformat() if not args.ontology_file else None,
        'acquisition': 'local_snapshot' if args.ontology_file else 'public_download', 'retained_snapshot': 'local_only/hp.obo'}, indent=2))
    print(f'HPO identifier/label checks: {sum(r["id_label_verified"] for r in rows)}/{len(rows)}; no patient metadata transmitted')


if __name__ == '__main__':
    main()
