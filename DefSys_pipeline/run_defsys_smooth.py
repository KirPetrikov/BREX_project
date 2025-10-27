"""v0.1p
Processed Padloc and DefenseFinder data to make summaries
for each tool and summary which combine both results
"""
import argparse
import pandas as pd

from pathlib import Path
from process_padloc_data import process_padloc_data_smooth
from process_dfnfnr_data import process_dfnfnr_data_smooth
from merge_padloc_and_defensefinder import merge_padloc_and_defensefinder
from make_combined_summary import make_combined_summary

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        'Processed Padloc and DefenseFinder data to make summaries '
        'for each tool and summary which combine both results.'
    )
    parser.add_argument('-p', '--padloc', type=Path,
                        help='Path to the Padloc folder')
    parser.add_argument('-d', '--defensefinder', type=Path,
                        help='Path to the DefenseFinder folder')
    parser.add_argument('-g', '--gff', type=Path,
                        help='Path to the folder with Prodigal gff-files for DefenseFinder')
    parser.add_argument('-x', '--regex', type=str,
                        help='Regex pattern to searg genes IDs in ggf-file Comment')
    parser.add_argument('-o', '--output', type=Path,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def merge_annotations(
        padloc_ann,
        dfnfnr_ann,
        results_dir
) -> set:
    df_ann_padloc = (
        pd.read_csv(padloc_ann, sep='\t')
        .rename(
            {'DS_ID': 'Padloc_DS_ID',
             'System': 'Padloc_System',
             'Annotation': 'Padloc_Annotation'},
            axis=1
        )
    )

    df_ann_dfnfnr = (
        pd.read_csv(dfnfnr_ann, sep='\t')
        .rename(
            {'DS_ID': 'DF_DS_ID',
             'System': 'DF_System',
             'System_sub': 'DF_System_sub',
             'Annotation': 'DF_Annotation'},
            axis=1
        )
    )

    # --- Check Padloc/DF ID uniqueness
    non_unq_ds_ids = set(df_ann_padloc.Padloc_DS_ID).intersection(set(df_ann_dfnfnr.DF_DS_ID))
    if non_unq_ds_ids:
        print(
            f'### Non-unique IDs:\n'
            f'{non_unq_ds_ids}\n'
            f'---'
        )
    else:
        print('### All DS_ID are unique')

    common_accessions = set(df_ann_padloc.Accession).intersection(set(df_ann_dfnfnr.Accession))
    padloc_unq_acc = set(df_ann_padloc.Accession).difference(common_accessions)
    dfnfnr_unq_acc = set(df_ann_dfnfnr.Accession).difference(common_accessions)

    # --- Write results
    df_ann = pd.merge(
        df_ann_dfnfnr,
        df_ann_padloc,
        on=['Protein', 'Accession', 'Nucleotide'], how='outer'
    ).fillna('miss')

    df_ann.to_csv(f'{results_dir}/Merged_annotations.tsv', sep='\t', index=False)

    with open(f'{results_dir}/Unique_accessions_Padloc.txt', mode='w') as f:
        for acc in padloc_unq_acc:
            f.write(f'{acc}\n')

    with open(f'{results_dir}/Unique_accessions_DefenseFinder.txt', mode='w') as f:
        for acc in dfnfnr_unq_acc:
            f.write(f'{acc}\n')

    with open(f'{results_dir}/Common_accessions.txt', mode='w') as f:
        for acc in common_accessions:
            f.write(f'{acc}\n')

    return common_accessions


def run_defsys_smooth(padloc,
                      dfnfnr,
                      gff,
                      pattern: str,
                      results_path):

    print('\n>>> Run smooth DefSys pipeline <<<')
    results_path = Path(results_path)
    results_path.mkdir(parents=True, exist_ok=True)

    print('\n>>> Padloc processing')

    padloc_results = results_path / 'Summary_Padloc'
    padloc_results.mkdir(parents=True, exist_ok=True)
    process_padloc_data_smooth(
        Path(padloc),
        Path(padloc_results)
    )

    print('\n>>> DefenseFinder processing')

    dfnfnr_results = results_path / 'Summary_DefenseFinder'
    dfnfnr_results.mkdir(parents=True, exist_ok=True)
    process_dfnfnr_data_smooth(
        dfnfnr,
        gff,
        pattern,
        dfnfnr_results
    )

    print('\n>>> Merge annotations')
    accessions = merge_annotations(
        padloc_results / 'protein_annotations.tsv',
        dfnfnr_results / 'protein_annotations.tsv',
        results_path
    )

    print('\n>>> Make samples list')
    samples_list = []
    for acc in accessions:
        samples_list.append(
            f'{acc},{dfnfnr}/{acc},{padloc}/{acc}\n'
        )
    samples = ''.join(samples_list)
    with open(results_path / 'samples.csv', mode='w') as f:
        f.write(samples)

    print('\n>>> Merge Padloc and DefenseFinder defense systems')
    merge_results = results_path / 'PDF_merge'
    merge_results.mkdir(parents=True, exist_ok=True)
    merge_padloc_and_defensefinder(
        results_path / 'samples.csv',
        merge_results
    )

    print('\n>>> Make combined defense system summary')
    make_combined_summary(
        results_path / 'Merged_annotations.tsv',
        results_path / 'PDF_merge/By_Accessions',
        results_path / 'Summary_Padloc/defsys_summary.json',
        results_path / 'Summary_Padloc/redundant_defsys.tsv',
        results_path / 'Unique_accessions_Padloc.txt',
        results_path / 'Summary_DefenseFinder/defsys_summary.json',
        results_path / 'Summary_DefenseFinder/redundant_defsys.tsv',
        results_path / 'Unique_accessions_DefenseFinder.txt',
        results_path
    )

    print("That's all")


if __name__ == '__main__':
    args = parse_arguments()

    run_defsys_smooth(args.padloc,
                      args.defensefinder,
                      args.gff,
                      args.regex,
                      args.output)
