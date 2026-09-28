"""Select the maximum unrounded AutoXZ Informedness for k >= 3."""
import math
from pathlib import Path


def select_best(rows):
    candidates = [row for row in rows
                  if row['chromosome_class'] == 'AutoXZ'
                  and int(row['Kmer']) >= 3
                  and math.isfinite(float(row['Informedness']))]
    if not candidates:
        raise ValueError('No finite AutoXZ Informedness for k >= 3')
    # Stable min retains original height-column order when score and k both tie.
    return min(candidates, key=lambda row: (-float(row['Informedness']), int(row['Kmer'])))


def write_selection(rows, destination):
    best = select_best(rows)
    Path(destination).write_text(f"{int(best['Kmer'])}\t{best['height_class']}\n")
