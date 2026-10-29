import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def active_configuration(require_inputs=True):
    path=ROOT/'ACTIVE_ANALYSIS.json'
    if not path.exists():
        raise FileNotFoundError('No validated active analysis manifest; run the private 120-feature rebuild first')
    config=json.loads(path.read_text())
    if config.get('status')!='VALIDATED_ACTIVE' or config.get('feature_count')!=120 or config.get('synthetic_rows_in_primary')!=0:
        raise ValueError('The active analysis is not a validated real-source 120-feature configuration')
    for field in ['disease_csv','reference_csv','disease_master_csv','reference_master_csv','feature_schema','private_run','results_directory']:
        config[field]=str(ROOT/config[field])
    if require_inputs:
        for field in ['disease_csv','reference_csv','disease_master_csv','reference_master_csv','feature_schema']:
            if not Path(config[field]).is_file():
                raise FileNotFoundError('Authorized private current inputs are unavailable; historical inputs will not be substituted')
    return config
