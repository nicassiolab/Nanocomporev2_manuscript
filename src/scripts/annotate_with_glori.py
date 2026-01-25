#!python

import sys

import pandas as pd


reference_file = sys.argv[1]
glori_file = sys.argv[2]


reference = pd.read_csv(reference_file, sep='\t')
glori = pd.read_csv(glori_file, sep='\t', usecols=['chrom', 'strand', 'start'])
glori['modified'] = 1
# The high confidence set contains the start-end pairs indicating a 5mer.
# We want to annotate tho central position as the modified one.
glori['start'] += 2


joined = reference.join(glori.set_index(['chrom', 'strand', 'start']),
                        on=['chr', 'strand', 'genomicPos'],
                        how='left')
joined['modified'] = joined['modified'].fillna(0).astype(int)

joined.to_csv(sys.stdout, index=False, sep='\t')

