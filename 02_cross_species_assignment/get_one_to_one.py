#!/usr/bin/env python3
'Filter one-to-one orthologs without changing their input order.'

import sys, re, glob
import pandas as pd




def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: get_one_to_one.py PAIRS.tsv")
    args = sys.argv
    df = pd.read_csv(args[1], sep='\t', header=None, names=['Sp_1', 'Sp_2'])
    s1_uniq = df['Sp_1'].value_counts().where(lambda d: d == 1).dropna()
    s2_uniq = df['Sp_2'].value_counts().where(lambda d: d == 1).dropna()
    df_o = df[df['Sp_1'].isin(s1_uniq.index) & df['Sp_2'].isin(s2_uniq.index)]
    df_o.replace('_[A-Za-z]{6}$', '', regex=True).to_csv(sys.stdout, header=False, index=False, sep='\t')


if __name__ == "__main__":
    main()
