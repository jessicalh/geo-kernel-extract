#!/usr/bin/env python3
"""Render Stage 1.1 per-atom-type results as JPG tables for email."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

out = Path(__file__).parent / 'tables'
out.mkdir(exist_ok=True)

# ── Table 1: Non-normalised (raw kernels) ──

cols1 = ['Atom type', 'Atoms', 'R²']
data1 = [
    ['H',              '230,135', '0.921'],
    ['C  (all)',       '133,488', '0.514'],
    ['    CA',          '29,944', '0.541'],
    ['    C=O',         '29,944', '0.361'],
    ['    CB',          '27,429', '0.578'],
    ['    C side',      '46,171', '0.660'],
    ['N  (all)',        '39,954', '0.210'],
    ['    N bb',        '29,944', '0.212'],
    ['    N side',      '10,010', '0.588'],
    ['O  (all)',        '42,429', '0.231'],
    ['    O bb',        '29,944', '0.253'],
    ['    O side',      '12,485', '0.241'],
]

# ── Table 2: Normalised + progressive scalars ──

cols2 = ['Atom type', 'Atoms', 'Base', '+Scales', '+Mut', 'Fair']
data2 = [
    ['H',              '230,135', '0.848', '0.924', '0.861', '0.928'],
    ['C  (all)',       '133,488', '0.471', '0.529', '0.512', '0.562'],
    ['    CA',          '29,944', '0.539', '0.577', '0.597', '0.627'],
    ['    C=O',         '29,944', '0.372', '0.411', '0.430', '0.463'],
    ['    CB',          '27,429', '0.544', '0.597', '0.604', '0.647'],
    ['    C side',      '46,171', '0.646', '0.690', '0.694', '0.729'],
    ['N  (all)',        '39,954', '0.245', '0.292', '0.345', '0.380'],
    ['    N bb',        '29,944', '0.254', '0.301', '0.351', '0.387'],
    ['    N side',      '10,010', '0.720', '0.762', '0.870', '0.887'],
    ['O  (all)',        '42,429', '0.274', '0.303', '0.358', '0.382'],
    ['    O bb',        '29,944', '0.303', '0.338', '0.395', '0.422'],
    ['    O side',      '12,485', '0.334', '0.373', '0.543', '0.566'],
]

def render_table(cols, data, path, col_widths=None):
    n_cols = len(cols)
    n_rows = len(data)

    if col_widths is None:
        col_widths = [1.0 / n_cols] * n_cols

    fig_w = 5.5
    row_h = 0.30
    fig_h = (n_rows + 1) * row_h + 0.15

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis('off')

    tbl = ax.table(
        cellText=data,
        colLabels=cols,
        colWidths=col_widths,
        loc='center',
        cellLoc='center',
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.35)

    for j in range(n_cols):
        cell = tbl[0, j]
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold', fontsize=8)
        cell.set_edgecolor('white')
        cell.set_linewidth(1.5)

    summary_rows = {1, 2, 7, 10}
    subtype_rows = {3, 4, 5, 6, 8, 9, 11, 12}

    for i in range(1, n_rows + 1):
        for j in range(n_cols):
            cell = tbl[i, j]
            cell.set_edgecolor('#bdc3c7')
            cell.set_linewidth(0.8)
            if i in summary_rows:
                cell.set_facecolor('#eaf2f8')
                cell.set_text_props(fontweight='bold')
            elif i in subtype_rows:
                cell.set_facecolor('#f8f9fa')
            else:
                cell.set_facecolor('white')
            if j == 0:
                cell.set_text_props(ha='left')
                cell._loc = 'left'

    fig.tight_layout(pad=0.1)
    fig.savefig(str(path), dpi=200, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f'Wrote {path}')

render_table(
    cols1, data1,
    out / 'raw_kernels.jpg',
    col_widths=[0.35, 0.30, 0.35],
)

render_table(
    cols2, data2,
    out / 'normalised_fair.jpg',
    col_widths=[0.24, 0.16, 0.15, 0.15, 0.15, 0.15],
)

print('Done.')
