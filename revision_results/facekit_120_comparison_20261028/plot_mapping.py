import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    out = parser.parse_args().output
    crosswalk = pd.read_csv(out/'feature_crosswalk.csv')
    summary = crosswalk.groupby('mapping_status').size().rename('old_feature_count').reset_index()
    summary.to_csv(out/'mapping_summary_data.csv', index=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, ax = plt.subplots(figsize=(8.5, 2.7), constrained_layout=True)
    labels = ['Changed full measurement definition', 'Removed without direct match']
    values = [int(crosswalk.mapping_status.eq(status).sum()) for status in ['changed_definition', 'removed']]
    bars = ax.barh(labels, values, color=['#247c87', '#b1bdc7'], height=.55)
    ax.invert_yaxis()
    ax.set_xlim(0, 130)
    ax.set_xlabel('Historical features (n = 125; includes one hairline-to-forehead merge)')
    ax.set_title('Feature-definition audit, not patient-level agreement', loc='left', fontweight='bold', pad=16)
    for bar, value in zip(bars, values):
        ax.text(value+1, bar.get_y()+bar.get_height()/2, str(value), va='center', fontweight='bold')
    for side in ['top', 'right', 'left']:
        ax.spines[side].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.savefig(out/'feature_mapping_summary.png', dpi=220, bbox_inches='tight')
    fig.savefig(out/'feature_mapping_summary.svg', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
