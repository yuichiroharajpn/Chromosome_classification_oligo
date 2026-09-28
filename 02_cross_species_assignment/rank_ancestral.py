#!/usr/bin/env python3
'Score dendrogram height classes against ancestral linkage groups.'

import glob
import os
import re
import sys
import pandas as pd
import numpy as np
from select_best import write_selection


def cal_dif(a, b):
    # The first retained chromosome defines the macrochromosome cluster.
    nn = a.iloc[0]
    df_cmp = pd.DataFrame({'class': a})
    df_cmp['Mm'] = b
    df_cmp['rMm'] = np.where(df_cmp['class'] == nn, 'Larger', 'Smaller')
    a1_true = len(df_cmp[(df_cmp['Mm'] == 'Larger') & (df_cmp['rMm'] == 'Larger')])
    a2_true = len(df_cmp[(df_cmp['Mm'] == 'Smaller') & (df_cmp['rMm'] == 'Smaller')])
    a1_false = len(df_cmp[(df_cmp['Mm'] == 'Larger') & (df_cmp['rMm'] == 'Smaller')])
    a2_false = len(df_cmp[(df_cmp['Mm'] == 'Smaller') & (df_cmp['rMm'] == 'Larger')])
    if a1_true + a2_true + a1_false + a2_false > 0:
        accuracy = (a1_true + a2_true) / (a1_true + a2_true + a1_false + a2_false)
        if a1_true + a1_false > 0:
            accuracy2 = a1_true / (a1_true + a1_false)
        else:
            accuracy2 = 0
        if a2_true + a2_false > 0:
            accuracy3 = a2_true / (a2_true + a2_false)
        else:
            accuracy3 = 0
        accuracy4 = accuracy2 + accuracy3 - 1
    else:
        accuracy = 0
        accuracy2 = 0
        accuracy3 = 0
        accuracy4 = 0
    return pd.Series([a1_true + a2_false, a2_true + a1_false, accuracy, accuracy2, accuracy3, accuracy4], index=['Num_Macro', 'Num_Micro', 'Accuracy', 'Sensitivity', 'Specificity', 'Informedness'])

def natural_keys(path):
    name = os.path.basename(path)
    return [int(s) if s.isdigit() else s for s in re.split('(\\d+)', name)]

def main():
    selection_file = None
    if "--selection-file" in sys.argv:
        position = sys.argv.index("--selection-file")
        if position + 1 == len(sys.argv):
            sys.exit("--selection-file requires a path")
        selection_file = sys.argv[position + 1]
        del sys.argv[position:position + 2]
    if len(sys.argv) != 4:
        sys.exit("Usage: rank_ancestral.py CLUSTER_DIR ANCESTRAL_SYNTENY.tsv ANCESTRAL_CLASSES.tsv [--selection-file FILE]")
    args = sys.argv
    df_ref = pd.read_csv(args[2], sep='\t')
    df_ref = df_ref.rename(columns={'sp1.Chr': 'rChr', 'sp2.Chr': 'Chr'})
    df_ref = df_ref.loc[df_ref.groupby('rChr')['orthologs'].idxmax()]
    df_rMm = pd.read_csv(args[3], sep='\t', names=('rChr', 'rMm'))
    df_ref = pd.merge(df_ref, df_rMm, on='rChr', how='left')
    pattern = os.path.join(args[1], '*mer_mid_relfreq.dendro.grp.h.txt')
    files = glob.glob(pattern)
    l = sorted(files, key=natural_keys)
    if not l:
        raise ValueError('No *mer_mid_relfreq.dendro.grp.h.txt files found')
    selection_rows = []
    for f in l:
        df = pd.read_csv(f, sep='\t')
        df = df.drop(columns='Length')
        df = df[df['Chr'].isin(df_ref['Chr'])]
        df = df[~df['Chr'].isin(['chrY', 'chrW'])]
        df1 = pd.merge(df, df_ref.loc[:, ['Chr', 'rChr', 'rMm']], on='Chr', how='inner')
        df1 = df1.set_index('Chr')
        x = df1.iloc[:, 0:-2].apply(cal_dif, b=df1['rMm']).T
        x['chromosome_class'] = 'AutoXZ'
        x.reset_index(inplace=True)
        x = x.rename(columns={'index': 'height_class'})
        m = re.search('(\\d+)mer', os.path.basename(f))
        x['Kmer'] = m.group(1)
        x['ratio_Macro'] = x['Num_Macro'] / len(df)
        x = x.astype({'Num_Macro': int, 'Num_Micro': int})
        selection_rows.extend(x.to_dict('records'))
        # Keep the legacy column order and four-decimal serialization.
        if l.index(f) == 0:
            x.iloc[:, [8, 0, 1, 2, 7, 9, 3, 4, 5, 6]].to_csv(sys.stdout, header=True, index=False, sep='\t', float_format='%.4f')
        else:
            x.iloc[:, [8, 0, 1, 2, 7, 9, 3, 4, 5, 6]].to_csv(sys.stdout, header=False, index=False, sep='\t', float_format='%.4f')
    if selection_file is not None:
        write_selection(selection_rows, selection_file)


if __name__ == "__main__":
    main()
