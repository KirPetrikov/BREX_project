"""v0.1
Parse DefenseFinder results ('systems.tsv' files) to create
- json defense systems summary;
- tsv-table of proteins annotations;
- tsv-table of DS where there is protein duplication.
"""
import argparse
import json
import pandas as pd

from pathlib import Path
from defsys import (parse_defense_finder_tsv,
                    make_unidir_genes_defsys,
                    parse_prodigal_gff,
                    find_dupl_defsys,
                    create_defsys_summary)

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
                    'table of Padloc proteins annotations, table of DS where there is protein duplication, '
                    'and table with DS where there art non-unidirectional genes.'
                                     )
    parser.add_argument('input_dfnfnr_path', type=str,
                        help='Path to the DefenseFinder data folder')
    parser.add_argument('input_padloc_path', type=str,
                        help='Path to the Padloc data folder, for gff-files')
    parser.add_argument('input_accessions_path', type=str,
                        help='Path to the Nucleotide-Accession table, for gff-files')
    parser.add_argument('output_results_path', type=str,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def process_dfnfnr_data(path_to_systems,
                        path_to_gff,
                        path_to_accession) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reads DefenseFinder systems-tsv,
    chooses uniform strandness for DS with non-unidirectional genes,
    create summary for all DS.
    Does not take into account possible proteins duplications, but
    finds DS where there is proteins duplication.

    Params:
    - path_to_tsv

    Return:
    - defsys_summary: Summary table with general DS data;
    - dupl_curr_full: Table with DS where there is protein duplications
    """

    defsys_data = (
        parse_defense_finder_tsv(path_to_systems, add_accession=path_to_accession)[
            ['DS_ID', 'Nucleotide', 'Accession', 'DF_System_sub', 'All_Proteins']
        ].to_dict(orient='index')
    )

    new_columns = ('DS_ID', 'Accession', 'Nucleotide', 'DF_System_sub')
    list_sys_new = {'DS_ID': [],
                    'Accession': [],
                    'Nucleotide': [],
                    'DF_System_sub': [],
                    'Protein': []}

    for defsys in defsys_data.values():
        if ',' in defsys['All_Proteins']:
            all_prots: list = defsys['All_Proteins'].split(',')
        else:
            all_prots = [defsys['All_Proteins']]
        for curr_prot in all_prots:
            for col in new_columns:
                list_sys_new[col].append(defsys[col])
            list_sys_new['Protein'].append(curr_prot)

    df_new = pd.DataFrame.from_dict(list_sys_new)

    df_gff = parse_prodigal_gff(path_to_gff).rename(
        {'Gene_ID': 'Protein', 'Chrom': 'Nucleotide'},
        axis=1
    )

    # Create new DS-Proteins table
    df_new = (df_new.merge(df_gff[['Protein', 'Start', 'End', 'Strand']])
                    .rename({'DF_System_sub': 'System'}, axis=1))

    make_unidir_genes_defsys(df_new)

    dupl_defsys = find_dupl_defsys(df_new)

    df_summary = create_defsys_summary(df_new)

    return df_summary, dupl_defsys


args = parse_arguments()
input_dfnfnr_path = Path(args.input_dfnfnr_path)
input_padloc_path = Path(args.input_padloc_path)
input_accessions_path = args.input_accessions_path
output_results_path = Path(args.output_results_path)

output_results_path.mkdir(parents=True, exist_ok=True)

summary_defsys_all = {}
duplicate_defsys_all = []

for folder in input_dfnfnr_path.iterdir():
    print(f'---Processing {folder.name}---')

    curr_dfnfnr_path = folder / f'{folder.name}_systems.tsv'
    curr_padloc_path = input_padloc_path / f'{folder.name}/{folder.name}_prodigal.gff'

    summary_curr, dupl_defsys_curr = process_dfnfnr_data(
        curr_dfnfnr_path, curr_padloc_path, input_accessions_path
    )

    if not dupl_defsys_curr.empty:
        duplicate_defsys_all.append(dupl_defsys_curr)

    curr_accession_result_path = output_results_path / f'By_Accessions/{folder.name}'
    curr_accession_result_path.mkdir(parents=True, exist_ok=True)

    summary_curr.to_json(curr_accession_result_path / (folder.name + '_summary.json'), orient='index')
    summary_defsys_all.update(summary_curr.to_dict(orient='index'))

# --- Write results ---
if duplicate_defsys_all:
    (pd.concat(duplicate_defsys_all)
       .reset_index(drop=True)
       .to_csv(output_results_path / 'duplicated_defsys.tsv', sep='\t', index=False))

with open(output_results_path / 'defsys_summary_all.json', mode='w') as f:
    json.dump(summary_defsys_all, f, indent=4)
