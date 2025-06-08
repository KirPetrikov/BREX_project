"""
Add assembly accession IDs to defsys summary
"""
import json
import pandas as pd

input_summary_path = '/home/holydiver/Main/2024_BREX/Data/20250323_defsys_summary/defsys_summary_all.json'
input_accessions = '/home/holydiver/Main/2024_BREX/Data/Accessions_summary.tsv'

with open(input_summary_path) as f:
    data: dict = json.load(f)

df_acc = pd.read_csv(input_accessions, sep='\t')
n_to_a = dict(zip(df_acc.Nucleotide, df_acc.Accession))

for k, v in data.items():
    acc = n_to_a[v['Nucleotide']]
    data[k]['Accession'] = acc

with open(input_summary_path, mode='w') as f:
    json.dump(data, f, indent=4)
