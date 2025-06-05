"""v0.1b
Processing of duplicated DSs that fell into inners.
According to the list obtained during the creation of the summary of BREX regions.
Manual check of duplicates.

Only 4 cases of duplications. According to e-val, in three cases the PDC definitely remains.
In one case: it is also the same, although the e-val is slightly higher, but the coverage is much better.
In all cases, DSs consist of duplicated proteins entirely.
"""

import json
import pandas as pd

from pathlib import Path

pd.options.mode.copy_on_write = True

# Annotations
input_annot = Path('/home/holydiver/Main/2024_BREX/Data/20250324_dupl_brex_processing/'
                   'protein_annotations_dupl_brex_proc.tsv')
output_annot = Path('/home/holydiver/Main/2024_BREX/Data/20250406_brex_regions_summary/'
                    'protein_annotations_dupl_proc.tsv')

df_ann = pd.read_csv(input_annot, 
                     sep='\t',
                     dtype={'Protein': str, 'Padloc_ann': str, 'System': str, 'Nucleotide': str, 
                            'Start': int, 'End': int, 'Strand': str, 'Localisation': str})

idxs_prots_to_drop = [1332145, 56265, 56267, 1120637, 1120639, 1095724, 1095726]

# Drop
df_ann = df_ann.drop(idxs_prots_to_drop)
# Save
df_ann.to_csv(output_annot, sep='\t', index=False)

# Defsys regions summary
input_summary_path = Path('/home/holydiver/Main/2024_BREX/Data/20250511_brex_regions_summary/inner_summary_brex.json')
output_summary_path = Path('/home/holydiver/Main/2024_BREX/Data/20250511_brex_regions_summary/'
                           'inner_summary_brex_dupl_proc.json')

# Target defsys regions summary
with open(input_summary_path) as f:
    regions_summary = json.load(f)

brex = ['brex_type_I%21%CP000482.1', 'brex_type_I%53%NZ_CP113517.1',
        'brex_type_I%28%NZ_CP077135.1', 'brex_type_I%26%NZ_CP084075.1']
defsys_to_drop = ['SoFic%20%CP000482.1', 'HEC-01%9%NZ_CP113517.1',
                  'AbiL%2%NZ_CP077135.1', 'AbiL%2%NZ_CP084075.1']

# Drop
for k, ds in zip(brex, defsys_to_drop):
    del regions_summary[k]['Inner_DS_full'][ds]
# Save
with open(output_summary_path, mode='w') as f:
    json.dump(regions_summary, f, indent=4)
