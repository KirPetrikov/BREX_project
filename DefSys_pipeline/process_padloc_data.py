"""v0.6p
Parse all given Padloc results (csv-files) to create:
- defense systems json-summary
- table of Padloc proteins annotations
- table of DS where there are protein redundancy

Folder structure required:
input_data_path/
    Accession1/
        Accession1_padloc.csv
    Accession2/
        Accession2_padloc.csv
    etc.
"""
import argparse
import json
import pandas as pd

from pathlib import Path
from commons import make_unidir_genes_defsys, find_redundancy_defsys, create_defsys_summary

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        'Parse Padloc results (csv-files) to create DS-summary json,'
        'table of Padloc proteins annotations, table of DS where there is protein redundancy.'
    )
    parser.add_argument('-i', '--input_data_path', type=str,
                        help='Path to the Padloc folder')
    parser.add_argument('-o', '--output_path', type=str,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def parse_padloc_csv(path_to_csv: str | Path, unifyed: bool = True) -> pd.DataFrame:
    """
    Reads Padloc csv and create table in convinient form.
    Adds columns:
    - 'DS_ID' with unique ID to prevent merging of DS with the same name
    - 'Accession' with assembly GenBank ID, gets from folder name

    Params:
    - path_to_csv: Padloc systems csv-table
    - unifyed: If 'True' returns table in unifyed form with only selected columns

    Return:
    - DS-table, unifyed by default
    """

    df = pd.read_csv(path_to_csv,
                     names=['SysNo', 'Nucleotide', 'System', 'Protein',
                            'HMM_profile', 'HMM_prot_name', 'Annotation', 'Eval_full',
                            'Eval_i', 'Cov_target', 'Cov_hmm', 'Start', 'End', 'Strand',
                            'Description', 'Protein_ID', 'Contig_end', 'All'],
                     index_col=None,
                     dtype={'SysNo': str, 'Nucleotide': str, 'System': str, 'Protein': str,
                            'HMM_profile': str, 'HMM_prot_name': str, 'Annotation': str,
                            'Eval_full': float, 'Eval_i': float, 'Cov_target': float, 'Cov_hmm': float,
                            'Start': int, 'End': int, 'Strand': str, 'Description': str,
                            'Protein_ID': int, 'Contig_end': int, 'All': str},
                     usecols=[_ for _ in range(18)],
                     skiprows=1
                     )

    # Create unique DefSystems IDs
    df['DS_ID'] = df['System'] + '%' + df['SysNo'] + '%' + df['Nucleotide']
    df['Accession'] = path_to_csv.name[:-11]

    if unifyed:
        return df[
            ['DS_ID', 'Accession', 'Nucleotide', 'Protein', 'Annotation', 'System', 'Start', 'End', 'Strand']
        ]
    else:
        return df


def process_single_padloc_table(path_to_csv: str | Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    From Padloc csv creates summary for all DS, proteins annotations table,
    finds DS where there is proteins redundancy (but does not take it into account).

    Params:
    - path_to_csv: Padloc systems csv-table

    Return:
    - defsys_summary: Summary table with general data for each DS
    - df_proteins: Shorten table with annotations of DSs' proteins
    - df_rdn: Table with DSs where there is protein redundansy
    """

    df_all = parse_padloc_csv(path_to_csv)

    df_rdn = find_redundancy_defsys(df_all)

    df_proteins = df_all[
        ['DS_ID', 'Accession', 'Nucleotide', 'Protein', 'Annotation', 'System']
    ]

    make_unidir_genes_defsys(df_all)

    defsys_summary = create_defsys_summary(df_all)

    return defsys_summary, df_proteins, df_rdn


def process_padloc_data(
        input_data_path: Path,
        results_path: Path
) -> None:
    results_path.mkdir(parents=True, exist_ok=True)

    summary_defsys_all = {}
    redundancy_defsys_all = []
    protein_annotations_all = []

    for folder in input_data_path.iterdir():
        print(f'---Processing {folder.name}---')

        summary_curr, prot_curr, rdn_curr = process_single_padloc_table(folder / f'{folder.name}_padloc.csv')

        protein_annotations_all.append(prot_curr)

        redundancy_defsys_all.append(rdn_curr)

        summary_defsys_all.update(summary_curr.to_dict(orient='index'))

        # --- Write current results ---
        curr_accession_result_path = results_path / f'By_Accessions/{folder.name}'
        curr_accession_result_path.mkdir(parents=True, exist_ok=True)
        summary_curr.to_json(curr_accession_result_path / f'{folder.name}_summary.json', orient='index')

    # --- Write results ---
    if redundancy_defsys_all:
        (
            pd.concat(redundancy_defsys_all)
              .reset_index(drop=True)
              .to_csv(results_path / 'redundant_defsys.tsv', sep='\t', index=False)
        )

    (
        pd.concat(protein_annotations_all)
          .reset_index(drop=True)
          .to_csv(results_path / 'protein_annotations.tsv', sep='\t', index=False)
    )

    with open(results_path / 'defsys_summary.json', mode='w') as f:
        json.dump(summary_defsys_all, f, indent=4)


if __name__ == '__main__':
    args = parse_arguments()

    process_padloc_data(
        Path(args.input_data_path),
        Path(args.output_path)
    )
