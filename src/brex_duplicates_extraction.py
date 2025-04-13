"""v0.1
   Duplicate DS summary
"""
import json
import pandas as pd

from pathlib import Path
from collections import defaultdict


def create_defsys_dupl_summary(path_to_csv: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Input: Padlock csv
    Output: duplicated DS summary table
    """

    df_all = pd.read_csv(path_to_csv,
                         names=['SysNo', 'Nucleotide', 'System', 'Protein',
                                'Padloc_ann', 'Start', 'End', 'Strand'],
                         dtype={'SysNo': str, 'Nucleotide': str, 'System': str, 'Protein': str,
                                'Padloc_ann': str, 'Start': int, 'End': int, 'Strand': str},
                         usecols=[0, 1, 2, 3, 6, 11, 12, 13],
                         skiprows=1
                         )

    # Create unique DefSystems IDs
    df_all['DS_ID'] = df_all['System'] + '%' + df_all['SysNo'] + '%' + df_all['Nucleotide']

    # --- Log duplicate DS assignment ---
    duplicated_prots = df_all['Protein'].value_counts()[df_all['Protein'].value_counts() > 1].index
    dupl_defsys_ids = df_all[df_all['Protein'].isin(duplicated_prots)].DS_ID.unique()
    df_dupl_full = df_all[df_all.DS_ID.isin(dupl_defsys_ids)]

    df_duplicated_proteins = (
        df_dupl_full[df_dupl_full['Protein'].isin(duplicated_prots)
                     ].groupby(['Protein', 'Nucleotide'])[['DS_ID', 'Padloc_ann']]
                      .agg(lambda x: x.unique().tolist())
                      .reset_index()[['Nucleotide', 'Protein', 'DS_ID', 'Padloc_ann']]
                              )

    return df_dupl_full, df_duplicated_proteins


def check_brex(frame):
    defsys_l = frame.unique().tolist()
    result = any([i.startswith('brex') for i in defsys_l])
    return result


def find_brex_dupl_from_df(df_all: pd.DataFrame):
    duplicated_prots = df_all['Protein'].value_counts()[df_all['Protein'].value_counts() > 1].index

    m_brex = (df_all[
                  df_all['Protein'].isin(duplicated_prots)
              ].groupby('Protein')[['DS_ID']]
              .agg(check_brex).values
              )

    duplicated_prots_brex = (df_all[
                                 df_all['Protein'].isin(duplicated_prots)
                             ].groupby('Protein')[['DS_ID']]
                             .agg(check_brex)
                             )[m_brex].index

    dupl_defsys_ids_brex = df_all[df_all['Protein'].isin(duplicated_prots_brex)].DS_ID.unique()

    df_dupl_brex = df_all[df_all.DS_ID.isin(dupl_defsys_ids_brex)]

    return df_dupl_brex


input_folder_path = Path('/home/niagara/Storage/MetaRus/Common_dir/collections/complete_bacteria_collection/padloc')
output_path_results = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250311_duplicates_new')

duplicate_defsys_full = []
duplicate_proteins = []

for folder in input_folder_path.iterdir():
    print(f'---Processing {folder.name}---')
    dupl_curr_full, dupl_curr_prots = create_defsys_dupl_summary(folder / (folder.name + '_padloc.csv'))
    duplicate_defsys_full.append(dupl_curr_full)
    duplicate_proteins.append(dupl_curr_prots)

# Check existence at least one duplicated DS case
if duplicate_defsys_full:
    df_dupl_all = pd.concat(duplicate_defsys_full).reset_index(drop=True)
    df_dupl_prots = pd.concat(duplicate_proteins).reset_index(drop=True)

    df_dupl_all.to_csv(output_path_results / 'duplicated_defsys_extracted.tsv', sep='\t', index=False)
    df_dupl_prots.to_json(output_path_results / 'duplicated_proteins_summary.json', orient='index')

    brex_dupl_df = find_brex_dupl_from_df(df_dupl_all)
    brex_dupl_df.to_csv(output_path_results / 'brex_duplicated_extracted.tsv', sep='\t',
                        index=False)

    # --- Count duplicateted annotations for prots ---
    dupl_prots = df_dupl_prots.copy().to_dict()
    # Rename systems
    for k, v in dupl_prots['DS_ID'].items():
        dupl_prots['DS_ID'][k] = [i.split('%')[0] for i in v]

    ann_count_prots_brex = defaultdict(int)
    for v in dupl_prots['DS_ID'].values():
        check_brex = any([i.startswith('brex') for i in v])
        if check_brex:
            k = ', '.join(v)
            ann_count_prots_brex[k] += 1

    with open(output_path_results / 'brex_duplicated_proteins_counts.json', mode='w') as f:
        json.dump(ann_count_prots_brex, f, indent=4)

    with open(output_path_results / 'brex_report_duplicated_proteins_counts.txt', mode='w') as f:
        for k, v in ann_count_prots_brex.items():
            f.write(f'{k}:\t{v}\n')

