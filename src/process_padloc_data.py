"""v0.5
Parse Padloc results (csv-files) to create json-defense systems summary,
table of Padloc proteins annotations, table of DS where there is protein duplication,
and table with DS where there art non-unidirectional genes.
"""
import argparse
import json
import pandas as pd

from pathlib import Path
from collections import defaultdict
from defsys import parse_padloc_csv, find_dupl_defsys, create_defsys_summary

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Parse Padloc results (csv-files) to create json-defense systems summary, '
                    'table of Padloc proteins annotations, table of DS where there is protein duplication, '
                    'and table with DS where there art non-unidirectional genes.'
                                     )
    parser.add_argument('input_folder_path', type=str,
                        help='Path to the data folder')
    parser.add_argument('output_results_filder_path', type=str,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def process_padloc_data(path_to_csv) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """
    Reads Padloc csv, creates proteins annotations table, finds DS where there is proteins duplication,
    chooses uniform strandness for DS with non-unidirectional genes, create summary for all DS.
    Does not take into account possible proteins duplications.

    Params:
    - path_to_csv

    Return:
    - defsys_summary: Summary table with general DS data;
    - df_proteins: Table with proteins annotations;
    - dupl_curr_full: Table with DS where there is protein duplications
    - nonunidir_defsys_names: DS_IDs joined by ","
    """

    df_all = parse_padloc_csv(path_to_csv)

    df_proteins = df_all.iloc[:, [3, 4, 2, 1, 5, 6, 7]]
    df_proteins['Localisation'] = 'no'

    dupl_curr_full = find_dupl_defsys(df_all)

    # Systems with non-unidirectional/antiparallel genes ---
    nonunidir_defsys_ids = (df_all.loc[:, ('Strand', 'DS_ID')]
                                  .groupby('DS_ID')['Strand']
                                  .agg(pd.Series.nunique)
                                  .loc[lambda x: x > 1]
                                  .index)
    # Choose uniform strandness/direction by voting; '+' in case of equality
    if nonunidir_defsys_ids.empty:
        nonunidir_defsys = ()
    else:
        strand = (df_all.loc[df_all.DS_ID.isin(nonunidir_defsys_ids)]
                        .groupby(['DS_ID'])['Strand']
                        .agg(lambda x: x.value_counts().sort_index().idxmax()))

        df_all.loc[df_all.DS_ID.isin(nonunidir_defsys_ids), 'Strand'] = (
            df_all.loc[df_all.DS_ID.isin(nonunidir_defsys_ids), 'DS_ID'].map(strand)
                                                                         )

        nonunidir_defsys = ','.join(nonunidir_defsys_ids)

    defsys_summary = create_defsys_summary(df_all)

    return defsys_summary, df_proteins, dupl_curr_full, nonunidir_defsys


args = parse_arguments()
input_folder_path = Path(args.input_folder_path)
output_results_filder_path = Path(args.output_results_filder_path)

output_results_filder_path.mkdir(parents=True, exist_ok=True)

summary_defsys_all = {}
duplicate_defsys_all = []
nonunidir_defsys_all = defaultdict(str)  # List of DS with non-uniderectional/antiparallel genes
protein_annotations_all = []

for folder in input_folder_path.iterdir():
    print(f'---Processing {folder.name}---')
    summary_curr, prot_curr, dupl_defsys_curr, nonunidir_curr = process_padloc_data(
                                                                    folder / (folder.name + '_padloc.csv')
                                                                                    )

    curr_accession_result_path = output_results_filder_path / f'By_Accessions/{folder.name}'
    curr_accession_result_path.mkdir(parents=True, exist_ok=True)

    protein_annotations_all.append(prot_curr)

    if not dupl_defsys_curr.empty:
        duplicate_defsys_all.append(dupl_defsys_curr)
    if nonunidir_curr:
        nonunidir_defsys_all[folder.name] = nonunidir_curr

    summary_curr.to_json(curr_accession_result_path / (folder.name + '_summary.json'), orient='index')
    summary_defsys_all.update(summary_curr.to_dict(orient='index'))

# --- Write results ---
with open(output_results_filder_path / 'non_unidirectional_defsys.tsv', mode='w') as f:
    for k, v in nonunidir_defsys_all.items():
        f.write(f'{k}\t{v}\n')

if duplicate_defsys_all:
    (pd.concat(duplicate_defsys_all)
       .reset_index(drop=True)
       .to_csv(output_results_filder_path / 'duplicated_defsys.tsv', sep='\t', index=False))

(pd.concat(protein_annotations_all)
   .reset_index(drop=True)
   .to_csv(output_results_filder_path / 'protein_annotations_all.tsv', sep='\t', index=False))

with open(output_results_filder_path / 'defsys_summary_all.json', mode='w') as f:
    json.dump(summary_defsys_all, f, indent=4)
