import argparse
from pathlib import Path

import pandas as pd

from analysis_config import active_configuration


def main():
    parser = argparse.ArgumentParser(description='Export the validated current 120-feature primary score')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    config = active_configuration(require_inputs=False)
    source = Path(config['results_directory'])/'stats_each_disease_vs_healthy.csv'
    frame = pd.read_csv(source)
    if frame.feature.nunique() != 120 or len(frame) != config['primary_hypotheses']:
        raise ValueError('Current validated 120-feature results required')
    if args.output is None:
        print(f'Current score: {source}; historical standalone PDS v2 is retired in the private archive')
        return
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame[['disease', 'feature', 'hedges_g', 'mannwhitney_q_fdr_global', 'phenotype_difference_score_pdf']].to_csv(args.output, index=False)
    print('Current aggregate score export written')


if __name__ == '__main__':
    main()
