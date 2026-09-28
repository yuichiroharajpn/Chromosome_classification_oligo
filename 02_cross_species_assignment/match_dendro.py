#!/usr/bin/env python3
'Assign Larger/Smaller using the first chromosome at the selected dendrogram height.'

import sys, re, glob, math
import pandas as pd
import numpy as np




def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: match_dendro.py CLUSTERS.tsv HEIGHT_CLASS")
    args = sys.argv
    clstr = args[1]
    layer = args[2]
    df = pd.read_csv(clstr, sep='\t')
    df = df.drop(columns='Length')
    nn = df.at[0, layer]
    df['Mm'] = np.where(df[layer] == nn, 'Larger', 'Smaller')
    df.loc[:, ['Chr', 'Mm']].to_csv(sys.stdout, header=True, index=False, sep='\t', float_format='%.4f')


if __name__ == "__main__":
    main()
