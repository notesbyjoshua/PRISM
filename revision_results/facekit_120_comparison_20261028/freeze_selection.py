import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'stats_testing'))
from prepare_cohorts import deduplicate_reference, select_patients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    out = parser.parse_args().output.resolve()
    if HERE/'local_only' not in out.parents or not out.is_dir():
        raise ValueError('An existing local-only audit directory is required')
    target = out/'private_selected_image_manifest.csv'
    if target.exists():
        raise FileExistsError(target)
    disease = pd.read_csv(ROOT/'data/master_dataset_disease.csv', dtype={'patient_id': 'string', 'image_id': 'string'}, low_memory=False)
    reference = pd.read_csv(ROOT/'data/master_dataset_healthy.csv', low_memory=False)
    selected, _ = select_patients(disease, '3-9', 'pose')
    ref, _ = deduplicate_reference(reference)
    rows = []
    for source, frame, column in [('GMDB', selected, 'facekit_relative_path'), ('FairFace', ref, 'file')]:
        for item in frame.to_dict('records'):
            locator = str(item[column])
            path = Path(locator).expanduser()
            candidates = [path] if path.is_absolute() else [ROOT/path, ROOT/'data'/path]
            resolved = next((candidate for candidate in candidates if candidate.is_file()), None)
            rows.append({'source': source, 'patient_id': item.get('patient_id', ''), 'image_id': item['image_id'],
                         'recorded_image_locator': locator, 'original_file_available': resolved is not None,
                         'original_sha256': hashlib.sha256(resolved.read_bytes()).hexdigest() if resolved else '',
                         'input_use_authorized': 'unconfirmed; no extraction performed'})
    pd.DataFrame(rows).to_csv(target, index=False)
    print('Private selected-input manifest frozen locally; no image extraction performed.')


if __name__ == '__main__':
    main()
