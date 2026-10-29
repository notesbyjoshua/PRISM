import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[2]


def current_source(schema_path=None,source_path=None):
    config=json.loads((ROOT/'ACTIVE_ANALYSIS.json').read_text())
    schema=pd.read_csv(schema_path or ROOT/'publication/feature_definitions/current_120_feature_formulas.csv')
    if len(schema)!=120 or schema.feature.nunique()!=120 or not schema.source_commit.eq(config['facekit_commit']).all():
        raise ValueError('Exactly 120 current pinned FaceKit definitions are required')
    source=Path(source_path or ROOT/'FaceKit/src/facekit/core/geometric/extractor.py').read_bytes()
    if not schema.extractor_sha256.eq(hashlib.sha256(source).hexdigest()).all():
        raise ValueError('Current source differs from the pinned feature definitions')
    return source.decode()


if __name__=='__main__':
    from current_120 import export
    current_source()
    export('schema')
