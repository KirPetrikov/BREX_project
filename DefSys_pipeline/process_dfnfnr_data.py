"""v0.2p
Parse all given DefenseFinder results (genes tsv-files) to create:
- defense & anti-defense systems json-summary
- table with anti-defense systems, accessions and DS ids
- table with DefenseFinder proteins annotations
- table with DS where there are protein redundancy

Folder structure required:
input_data_path/
    Accession1/
        Accession1_genes.tsv
    Accession2/
        Accession1_genes.tsv
    etc.

For gff-files:
input_gff_path/
    Accession1/
        Accession1_prodigal.gff
    Accession2/
        Accession1_prodigal.gff
    etc.
"""
import argparse
import json
import pandas as pd
import re

from pathlib import Path
from commons import (make_unidir_genes_defsys,
                     find_redundancy_defsys,
                     create_defsys_summary,
                     parse_gff)

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        'Parse DefenseFinder results (genes files) to create DS-summary json,'
        'table of DefenseFinder proteins annotations, table of DS where there is protein redundancy,'
        'and separate table for antidefense systems.'
    )
    parser.add_argument('-i', '--input_data_path', type=str,
                        help='Path to the DefenseFinder folder')
    parser.add_argument('-g', '--input_gff_path', type=str,
                        help='Path to the folder with Prodigal gff-files')
    parser.add_argument('-o', '--output_path', type=str,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def parse_dfnfnr_genes(
        path_to_tsv: str | Path,
        path_to_gff: str | Path = '',
        pattern: str = '',
        rough: bool = False
) -> pd.DataFrame:
    """
    Reads DefenseFinder genes tsv and create table in convinient form.
    Separate poteins with antidefense activity.

    Adds columns:
    - 'DS_ID' with unique ID to prevent merging of DS with the same name
    - 'Accession' with assembly GenBank ID, gets from folder name
    If specified gff-file adds genes coordinates:
    - 'Start', 'End', 'Strand'

    Params:
    - path_to_tsv
    - path_to_gff: path to Prodigal gff; if isn't specified - just skip

    Return:
    - unifyed table for defense systems genes/proteins
    """
    df = pd.read_csv(path_to_tsv, sep='\t',
                     usecols=[0, 1, 2, 5, 22, 23, 24],
                     dtype={'replicon': str, 'hit_id': str, 'gene_name': str,
                            'sys_id': str, 'type': str, 'subtype': str, 'activity': str},
                     names=['replicon', 'Protein', 'Annotation', 'tmp_ID',
                            'System', 'System_sub', 'Activity'],
                     skiprows=1)

    df = df.drop_duplicates()

    # Create unique DS_IDs
    df['DS_ID'] = df.apply(
        lambda x: f'{x.System_sub}%{x.tmp_ID.split("_")[-1]}%{x.replicon}',
        axis=1
    )

    # Add coords for protein genes
    if path_to_gff:
        assert pattern, 'Pattern for protein id search in gff does not provided!'
        df_gff = parse_gff(path_to_gff, pattern).rename({'ID': 'Protein'}, axis=1)
        df_gff = df_gff.loc[df_gff.Protein != 0]

    if rough:
        df = df.rename({'replicon': 'Accession'}, axis=1)
        if path_to_gff:
            df_gff['tmp_Nucleotide'] = df_gff.Comment.apply(lambda x: re.search(r'ID=(\w+)_\d+', x).group(1))
            df_gff['tmp_num'] = df_gff.Nucleotide.apply(lambda x: x.split('_')[-1])
            df_gff.loc[:, 'Protein'] = df_gff.tmp_Nucleotide + '_' + df_gff.astype({'Protein': str}).Protein
            df_gff.loc[:, 'Nucleotide'] = df_gff.tmp_Nucleotide + df_gff.tmp_num

            df = df.merge(df_gff[['Nucleotide', 'Protein', 'Start', 'End', 'Strand']], on='Protein', how='left')

    else:
        df = df.rename({'replicon': 'Nucleotide'}, axis=1)
        df['Accession'] = Path(path_to_tsv).parent.name
        if path_to_gff:
            df_gff.loc[:, 'Protein'] = df_gff.Nucleotide + '_' + df_gff.astype({'Protein': str}).Protein

            df = df.merge(df_gff[['Protein', 'Start', 'End', 'Strand']], on='Protein', how='left')

    return df.drop(['tmp_ID'], axis=1)


# def parse_dfnfnr_genes_rough(
#         path_to_tsv: str | Path,
#         path_to_gff: str | Path,
#         pattern: str
# ) -> pd.DataFrame:
#     """
#     Reads DefenseFinder genes tsv and create table in convinient form.
#     Separate poteins with antidefense activity.
#
#     Adds columns:
#     - 'DS_ID' with unique ID to prevent merging of DS with the same name
#     - 'Accession' with assembly GenBank ID, gets from folder name
#     If 'path_to_gff' specified adds genes coordinates:
#     - 'Start', 'End', 'Strand'
#     If 'rough_sample' specified rewrire proteins IDs in spesified manner
#
#     Params:
#     - path_to_tsv
#     - path_to_gff: path to Prodigal gff; if isn't specified - just skip
#
#     Return:
#     - unifyed table for defense systems genes/proteins
#     """
#     df = pd.read_csv(path_to_tsv, sep='\t',
#                      usecols=[0, 1, 2, 5, 22, 23, 24],
#                      dtype={'replicon': str, 'hit_id': str, 'gene_name': str,
#                             'sys_id': str, 'type': str, 'subtype': str, 'activity': str},
#                      names=['Accession', 'Protein', 'Annotation', 'tmp_ID',
#                             'System', 'System_sub', 'Activity'],
#                      skiprows=1)
#
#     df = df.drop_duplicates()
#
#     # Create unique DS_IDs
#     df['DS_ID'] = df.apply(
#         lambda x: f'{x.System_sub}%{x.tmp_ID.split("_")[-1]}%{x.Accession}',
#         axis=1
#     )
#
#     # Add coords for protein genes
#     df_gff = parse_gff(path_to_gff, pattern).rename({'ID': 'Protein'}, axis=1)
#
#     df_gff = df_gff.loc[df_gff.Protein != 0]
#     df_gff = df_gff.astype({'Protein': str})
#     df_gff['tmp_Nucleotide'] = df_gff.Comment.apply(lambda x: re.search(r'ID=(\w+)_\d+', x).group(1))
#     df_gff['tmp_num'] = df_gff.Nucleotide.apply(lambda x: x.split('_')[-1])
#     df_gff.loc[:, 'Protein'] = df_gff.tmp_Nucleotide + '_' + df_gff.Protein
#     df_gff.loc[:, 'Nucleotide'] = df_gff.tmp_Nucleotide + df_gff.tmp_num
#
#     df = df.merge(df_gff[['Protein', 'Start', 'End', 'Strand', 'Nucleotide']], on='Protein', how='left')
#
#     return df.drop(['tmp_ID'], axis=1)


def process_single_dfnfnr_table(
        path_to_tsv: str | Path,
        path_to_gff: str | Path,
        pattern: str = '',
        rough: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    From DefenseFinder genes tsv creates summary for all DS, proteins annotations table,
    finds DS where there is proteins redundancy (but does not take it into account).
    Save poteins with antidefense activity to separate summary.
    Adds genes coordinates using corresponding gff-file.

    Params:
    - path_to_tsv: DefenseFinder genes table
    - path_to_gff: path to Prodigal gff

    Return:
    - defsys_summary: Summary table with general data for each DS
    - df_proteins: Shorten table with annotations of DSs' proteins
    - df_rdn: Table with DSs where there is protein redundansy
    """

    df_defsys = parse_dfnfnr_genes(path_to_tsv,
                                   path_to_gff,
                                   pattern,
                                   rough)

    df_anti = df_defsys.loc[df_defsys.Activity == 'Antidefense', ['Accession', 'DS_ID']]

    df_rdn = find_redundancy_defsys(df_defsys)

    df_proteins = df_defsys[
        ['DS_ID', 'Accession', 'Nucleotide', 'Protein', 'Annotation', 'System', 'System_sub']
    ]

    make_unidir_genes_defsys(df_defsys)

    defsys_summary = create_defsys_summary(df_defsys)

    cols_order = ['Accession', 'Nucleotide', 'Protein', 'Annotation', 'System', 'System_sub',
                  'Start', 'End', 'Strand', 'DS_Prots', 'Have_inner']
    return defsys_summary[cols_order], df_proteins, df_rdn, df_anti


def process_dfnfnr_data(
        input_data_path,
        input_gff_path,
        pattern,
        results_path
) -> None:
    Path(results_path).mkdir(parents=True, exist_ok=True)

    summary_defsys_all = {}
    redundancy_defsys_all = []
    protein_annotations_all = []
    antidefense_all = []

    for folder in Path(input_data_path).iterdir():
        sample_id = folder.name
        print(f'---Processing {sample_id}---')

        summary_curr, prot_curr, rdn_curr, anti_curr = process_single_dfnfnr_table(
            folder / f'{sample_id}_genes.tsv',
            Path(input_gff_path) / f'{sample_id}/{sample_id}_prodigal.gff',
            pattern
        )

        protein_annotations_all.append(prot_curr)

        redundancy_defsys_all.append(rdn_curr)

        summary_defsys_all.update(summary_curr.to_dict(orient='index'))

        antidefense_all.append(anti_curr)

        # --- Write current results ---
        curr_accession_result_path = results_path / f'By_Accessions/{sample_id}'
        curr_accession_result_path.mkdir(parents=True, exist_ok=True)
        summary_curr.to_json(curr_accession_result_path / f'{sample_id}_summary.json', orient='index')

    # --- Write results ---
    if redundancy_defsys_all:
        (
            pd.concat(redundancy_defsys_all)
              .reset_index(drop=True)
              .to_csv(results_path / 'redundant_defsys.tsv', sep='\t', index=False)
        )

    if antidefense_all:
        (
            pd.concat(antidefense_all).to_csv(
                results_path / 'antidefense_systems.tsv',
                sep='\t',
                index=False
            )
        )

    (
        pd.concat(protein_annotations_all)
          .reset_index(drop=True)
          .to_csv(results_path / 'protein_annotations.tsv', sep='\t', index=False)
    )

    with open(results_path / 'defsys_summary.json', mode='w') as f:
        json.dump(summary_defsys_all, f, indent=4)


def process_dfnfnr_data_rough(
        input_samples: list[tuple],
        pattern: str,
        results_path: Path
) -> None:
    results_path.mkdir(parents=True, exist_ok=True)

    summary_defsys_all = {}
    redundancy_defsys_all = []
    protein_annotations_all = []
    antidefense_all = []

    for sample_id, dfnfnr_dir, gff_file in input_samples:
        for file_name in Path(dfnfnr_dir).iterdir():
            if str(file_name).endswith('genes.tsv'):
                print(f'---Processing {sample_id}---')

                summary_curr, prot_curr, rdn_curr, anti_curr = process_single_dfnfnr_table(
                    file_name,
                    gff_file,
                    pattern,
                    True
                )

                protein_annotations_all.append(prot_curr)

                redundancy_defsys_all.append(rdn_curr)

                summary_defsys_all.update(summary_curr.to_dict(orient='index'))

                antidefense_all.append(anti_curr)

                # --- Write current results ---
                curr_accession_result_path = Path(results_path) / 'By_Accessions'
                curr_accession_result_path.mkdir(parents=True, exist_ok=True)
                summary_curr.to_json(curr_accession_result_path / f'{sample_id}_summary.json', orient='index')

    # --- Write results ---
    if redundancy_defsys_all:
        (
            pd.concat(redundancy_defsys_all)
              .reset_index(drop=True)
              .to_csv(results_path / 'redundant_defsys.tsv', sep='\t', index=False)
        )

    if antidefense_all:
        (
            pd.concat(antidefense_all).to_csv(
                results_path / 'antidefense_systems.tsv',
                sep='\t',
                index=False
            )
        )

    (
        pd.concat(protein_annotations_all)
          .reset_index(drop=True)
          .to_csv(Path(results_path) / 'protein_annotations.tsv', sep='\t', index=False)
    )

    with open(Path(results_path) / 'defsys_summary.json', mode='w') as f:
        json.dump(summary_defsys_all, f, indent=4)


if __name__ == '__main__':
    args = parse_arguments()

    process_dfnfnr_data(
        Path(args.input_data_path),
        Path(args.input_gff_path),
        Path(args.output_path)
    )
