#!/usr/bin/env python3
'Compare chromosome assignments with ancestral linkage groups.'

import sys, re, glob, math
import pandas as pd
import numpy as np




def main():
    if len(sys.argv) != 6:
        sys.exit("Usage: match_anc_chrom.py CLUSTERS.tsv ANCESTRAL_SYNTENY.tsv ANCESTRAL_CLASSES.tsv HEIGHT_CLASS HUANG_ANCESTRAL.tsv")
    args = sys.argv
    clstr = args[1]
    df_ref = pd.read_csv(args[2], sep='\t')
    df_ref = df_ref.rename(columns={'sp1.Chr': 'rChr', 'sp2.Chr': 'Chr'})
    df_rMm = pd.read_csv(args[3], sep='\t', names=('rChr', 'rMm'))
    df_ref = pd.merge(df_ref, df_rMm, on='rChr', how='left')
    layer = args[4]
    df = pd.read_csv(clstr, sep='\t')
    df = df.drop(columns='Length')
    nn = df.at[0, layer]
    df = df[df['Chr'].isin(df_ref['Chr'])]
    df = df[~df['Chr'].isin(['chrY', 'chrW'])]
    df['Mm'] = np.where(df[layer] == nn, 'Larger', 'Smaller')
    df_ref = pd.merge(df_ref, df.loc[:, ['Chr', 'Mm']], on='Chr', how='left')
    df20 = pd.read_csv(args[5], sep='\t')
    df2 = df20.dropna(subset=['Chicken.chromosome'])['anc_jv'].value_counts().reset_index()
    df2.columns = ['rChr', 'ngene']
    df_ref = pd.merge(df_ref, df2, on='rChr', how='left')
    df_ref['rortho'] = df_ref['orthologs'] / df_ref['ngene']
    df_ref = df_ref.assign(flg=df_ref['rMm'] == df_ref['Mm'])
    df_ref = df_ref.astype({'flg': int})
    df_ref['qflg'] = df_ref['flg'] * df_ref['rortho']
    df_ref.loc[:, ['rChr', 'Chr', 'orthologs', 'rMm', 'Mm', 'ngene', 'rortho', 'flg', 'qflg']].to_csv(sys.stdout, header=True, index=False, sep='\t', float_format='%.4f')


if __name__ == "__main__":
    main()
