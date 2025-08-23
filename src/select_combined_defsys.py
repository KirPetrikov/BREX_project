"""v0.2
Takes summary of Padloc and DF results, and table with they combined output.
Create combined summary by choosing unique DS id.
In case of ambiguity priority choice is Padloc variant.
Considers DFs, which are joined in combined output table ('DS1/DS2'), just selects first.
Checks duplicated proteins.
Also takes list of DS ids from genomes, processed only by Padloc, and adds them to combined summary.
"""
import json
import pandas as pd

from pathlib import Path
from defsys import find_dupl_defsys


data_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data')

input_ann_path = data_path / 'protein_annotations_20250822.tsv'
input_accessions_path = data_path / 'Accessions_summary.tsv'

input_pdf_comb_path = data_path / '20250822_PDF_Comb/PDF_comb_results'

input_padloc_summary_path = data_path / '20250323_defsys_summary/defsys_summary_all.json'
input_padloc_dupl_path = data_path / '20250323_defsys_summary/duplicated_defsys.tsv'
input_dfnfnr_summary_path = data_path / '20250820_dfnfnr_summary/defsys_summary_all.json'

input_padloc_uniq_path = data_path / 'padloc_uniq_ds_ids.txt'

output_results_path = data_path / '20250822_combined_summary'
output_results_path.mkdir(parents=True, exist_ok=True)

# --- Add unique DS_ID for each DF/row in ann-table
df_ann = (
    pd.read_csv(
        input_ann_path,
        sep='\t',
        dtype={
            'Protein': str, 'Nucleotide': str,
            'Padloc_ann': str, 'Padloc_System': str,
            'DF_ann': str, 'DF_System': str, 'DF_System_sub': str,
            'Padloc_DS_ID': str, 'DF_DS_ID': str
        }
    )
)

dfnfnr_mask = (df_ann.DF_DS_ID != 'miss').values
padloc_mask = (df_ann.Padloc_System != 'miss').values

# "Tmp-ids - DS_IDs" correspondence for DF.
print('---Tmp-ids for DF---')
df_ann['df_tmp_id'] = df_ann.loc[dfnfnr_mask].apply(
    lambda x: f'{x.Protein}%{x.DF_System_sub}',
    axis=1
)
df_ids_code = dict(zip(
    df_ann.loc[dfnfnr_mask].df_tmp_id,
    df_ann.loc[dfnfnr_mask].DF_DS_ID
))

# "Tmp-ids - DS_IDs" correspondence for Padloc
print('---Tmp-ids for Padloc---')
df_ann['p_tmp_id'] = df_ann.loc[padloc_mask].apply(
    lambda x: f'{x.Protein}%{x.Padloc_System}',
    axis=1
)
padloc_ids_code = dict(zip(
    df_ann.loc[padloc_mask].p_tmp_id,
    df_ann.loc[padloc_mask].Padloc_DS_ID
))

padloc_sel_ids = []
dfnfnr_sel_ids = []
duplicated_ds_ids = {'Accession': [], 'DS_ID': []}
splitted_ds_ids = {'Accession': [], 'DS_ID': []}

# --- Select unique DS_ID from all combined tables
print('---Read combined tables---')
for pdf_comb_table in input_pdf_comb_path.iterdir():
    print(f'---Processing {pdf_comb_table.name}---')

    df_pdfc = pd.read_csv(
        pdf_comb_table,
        sep=',',
        dtype={'system family DefenseFinder': str, 'system family PADLOC': str,
               'subtype DefenseFinder': str, 'subtype PADLOC': str,
               'genes': str, ' genes DefenseFinder': str, ' genes PADLOC': str},
        names=['DF_System', 'Padloc_System', 'DF_System_sub',
               'Padloc_System_sub', 'Proteins', 'DF_anns', 'Padloc_anns'],
        skiprows=1
    )
    defsys_number = df_pdfc.shape[0]

    # Process merged defsys
    idxs_to_split = []
    for tool_col in ('DF_System_sub', 'Padloc_System_sub'):
        s_mask = df_pdfc[tool_col].apply(lambda x: '/' in x).values

        if s_mask.any():
            curr_idxs_split = df_pdfc.loc[s_mask].index
            idxs_to_split.extend(curr_idxs_split)

            for i in curr_idxs_split:
                df_pdfc.loc[i, tool_col] = df_pdfc.loc[i, tool_col].split('/')[0]

    # --- Mapping of DS_ID
    # Tmp-ids for Padloc
    pdfc_padloc_mask = (df_pdfc.Padloc_System != 'N.A.').values

    df_pdfc['tmpid'] = ''
    df_pdfc.loc[pdfc_padloc_mask, 'tmpid'] = df_pdfc.apply(
        lambda x: f'{x.Proteins.split(";")[0]}%{x.Padloc_System_sub}',
        axis=1
    )

    # Tmp-ids for DF
    df_pdfc.loc[~pdfc_padloc_mask, 'tmpid'] = df_pdfc.apply(
        lambda x: f'{x.Proteins.split(";")[0]}%{x.DF_System_sub}',
        axis=1
    )

    df_pdfc['DS_ID'] = ''
    curr_padloc_sel_ids = [padloc_ids_code[i] for i in df_pdfc.loc[pdfc_padloc_mask, 'tmpid'].values]
    df_pdfc.loc[pdfc_padloc_mask, 'DS_ID'] = curr_padloc_sel_ids
    curr_df_sel_ids = [df_ids_code[i] for i in df_pdfc.loc[~pdfc_padloc_mask, 'tmpid'].values]
    df_pdfc.loc[~pdfc_padloc_mask, 'DS_ID'] = curr_df_sel_ids

    # Update final lists of DS_ID
    padloc_sel_ids.extend(curr_padloc_sel_ids)
    dfnfnr_sel_ids.extend(curr_df_sel_ids)

    # --- Find duplicated proteins
    prots_ids = {'Protein': [], 'DS_ID': []}

    tmp_prots = dict(zip(df_pdfc.DS_ID, df_pdfc.Proteins))
    for ds_id, prots in tmp_prots.items():
        curr_prots = prots.split(';')
        prots_ids['Protein'].extend(curr_prots)
        prots_ids['DS_ID'].extend([ds_id] * len(curr_prots))

    curr_dupl_ds_ids = (
        find_dupl_defsys(
            pd.DataFrame(prots_ids)
        ).DS_ID
         .unique()
         .tolist()
    )

    duplicated_ds_ids['DS_ID'].extend(curr_dupl_ds_ids)
    duplicated_ds_ids['Accession'].extend([pdf_comb_table.stem] * len(curr_dupl_ds_ids))

    if idxs_to_split:
        curr_splitted_ds_ids = df_pdfc.loc[idxs_to_split, 'DS_ID']
        splitted_ds_ids['DS_ID'].extend(curr_splitted_ds_ids)
        splitted_ds_ids['Accession'].extend([pdf_comb_table.stem] * len(curr_splitted_ds_ids))

    df_pdfc = df_pdfc.dropna()
    if defsys_number != df_pdfc.shape[0]:
        print(f'Warning: DS missed in {pdf_comb_table.stem}')

# --- Save intermediate results
with open(output_results_path / 'tmp_padloc_sel_ids.txt', mode='w') as f:
    f.write('\n'.join(padloc_sel_ids))
with open(output_results_path / 'tmp_dfnfnr_sel_ids.txt', mode='w') as f:
    f.write('\n'.join(dfnfnr_sel_ids))

pd.DataFrame(duplicated_ds_ids).to_csv(
    output_results_path / 'duplicated_ds_ids.tsv',
    sep='\t',
    index=False
)

pd.DataFrame(splitted_ds_ids).to_csv(
    output_results_path / 'splitted_ds_ids.tsv',
    sep='\t',
    index=False
)

print('---Combined tables processing completed---')

# --- Make final combined summary
with open(input_padloc_summary_path) as f:
    padloc_data = json.load(f)
with open(input_dfnfnr_summary_path) as f:
    dfnfnr_data = json.load(f)

# List of only Padloc DS ids (from genomes absent in DF results)
padloc_uniq_ds_ids = []
with open(input_padloc_uniq_path) as f:
    for line in f:
        padloc_uniq_ds_ids.append(line.strip())

combined_summary = {}
for ds_id in padloc_sel_ids:
    combined_summary[ds_id] = padloc_data[ds_id]
for ds_id in dfnfnr_sel_ids:
    combined_summary[ds_id] = dfnfnr_data[ds_id]
for ds_id in padloc_uniq_ds_ids:
    combined_summary[ds_id] = padloc_data[ds_id]

with open(output_results_path / 'combined_summary.json', mode='w') as f:
    json.dump(combined_summary, f, indent=4)

print('---Combined DS summary selection completed---')

# --- Add duplicates from only Padloc DS ids
df_padloc_dupl = pd.read_csv(input_padloc_dupl_path, sep='\t')
df_padloc_dupl = df_padloc_dupl.loc[df_padloc_dupl.DS_ID.isin(padloc_uniq_ds_ids)]

# Add Accessions
df_acc = pd.read_csv(input_accessions_path, sep='\t')[['Nucleotide', 'Accession']]
df_padloc_dupl = df_padloc_dupl.merge(df_acc, on='Nucleotide', how='left')

duplicated_ds_ids['DS_ID'].extend(df_padloc_dupl.DS_ID.to_list())
duplicated_ds_ids['Accession'].extend(df_padloc_dupl.Accession.to_list())

pd.DataFrame(duplicated_ds_ids).to_csv(
    output_results_path / 'duplicated_ds_ids.tsv',
    sep='\t',
    index=False
)

print('---Unique duplicated DS added---')
print('--- Finished ---')
