"""Creates protein annotation table with:
Padloc and DF annotations;
systems;
unique ids;
nucleotides and accessios.

NB: Contains all duplicates and Padloc/DF conflicts
"""
import pandas as pd
import json

pd.options.mode.copy_on_write = True

# Padloc protein annotations source table
input_padloc_ann = ('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/'
                    '20250323_defsys_summary/protein_annotations_all.tsv')

# Padloc source summary
input_padloc_summary = ('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/'
                        '20250323_defsys_summary/defsys_summary_all.json')

# DF genes/protein annotations updated table
input_dfnfnr_genes_path = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/dfnfnr_genes_acc.tsv'

input_accessions = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/Accessions_summary.tsv'

df_ann_padloc = (
    pd.read_csv(
        input_padloc_ann,
        sep='\t'
    ).rename({'System': 'Padloc_System'}, axis=1)
     .drop(['Localisation', 'Start', 'End', 'Strand'], axis=1)
)

df_ann_dfnfnr = pd.read_csv(
    input_dfnfnr_genes_path,
    sep='\t',
    dtype={'Nucleotide': str, 'Protein': str, 'DF_ann': str, 'DF_System': str,
           'DF_System_sub': str, 'Accession': str, 'DF_DS_ID': str},
).drop('Accession', axis=1)
# Drop Accession as not all nucleotides are in DF results, and adding Accessions still neccesary after merge

df_ann = pd.merge(
    df_ann_dfnfnr, df_ann_padloc, on=['Protein', 'Nucleotide'], how='outer'
).fillna('miss')

# Check all nucl
a = set(df_ann_padloc.Nucleotide.tolist())
b = set(df_ann_dfnfnr.Nucleotide.tolist())
c = a.union(b)
d = set(df_ann.Nucleotide.tolist())

assert not c.difference(d), 'Some nucleotides missed in combined dataframe'

# Check all prots
a = set(df_ann_padloc.Protein.tolist())
b = set(df_ann_dfnfnr.Protein.tolist())
c = a.union(b)
d = set(df_ann.Protein.tolist())

assert not c.difference(d), 'Some proteins missed in combined dataframe'

# --- Map unique Padloc DS_IDs to proteins
# Make Padloc tmp-ids
df_ann['tmpid_Padloc'] = df_ann.apply(lambda x: f'{x.Protein}_{x.Padloc_System}', axis=1)

with open(input_padloc_summary) as f:
    data_defsys = json.load(f)

proteins_padloc = {}

for defsys_id, curr_reg in data_defsys.items():
    for prot_id in curr_reg['DS_Prots']:
        tmp_id = f'{curr_reg["Nucleotide"]}_{prot_id}_{curr_reg["System"]}'
        proteins_padloc[tmp_id] = defsys_id

df_ann['Padloc_DS_ID'] = df_ann.tmpid_Padloc.map(proteins_padloc)

# --- Finalaze df
df_ann = df_ann.fillna('miss')
df_ann = df_ann.where(df_ann != '', 'miss')

df_ann = df_ann.drop(['tmpid_Padloc'], axis=1)

# --- Add Accessions
df_acc = pd.read_csv(input_accessions, sep='\t')[['Nucleotide', 'Accession']]

df_ann = df_ann.merge(df_acc, on='Nucleotide', how='left')

print(df_ann.head(2))

df_ann.to_csv('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/'
              'protein_annotations_20250821.tsv',
              sep='\t',
              index=False)
