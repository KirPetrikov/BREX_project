"""v0.1p
WARNING: docstring for rough pipeline do not updadted!

# TODO Update docstrings

Processed Padloc and DefenseFinder data to make summaries
for each tool and summary which combine both results
"""
import argparse
import csv
import pandas as pd

from pathlib import Path
from process_padloc_data import process_padloc_data_rough
from process_dfnfnr_data import process_dfnfnr_data_rough
from merge_padloc_and_defensefinder import merge_padloc_and_defensefinder
from make_combined_summary import make_combined_summary

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        'Processed Padloc and DefenseFinder data to make summaries '
        'for each tool and summary which combine both results.'
    )
    # TODO Add description
    parser.add_argument('-s', '--samples_list', type=Path,
                        help='___')
    parser.add_argument('-x', '--regex', type=str,
                        help='Regex pattern to searg genes IDs in ggf-file Comment')
    parser.add_argument('-o', '--output_path', type=Path,
                        help='Path to the results folder. Will be created if it does not exist')
    return parser.parse_args()


def test_csv_input(input_samples_file):
    print('\n >>> Check files <<<')

    with open(input_samples_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            flag = False
            for file_name in Path(row[2]).iterdir():
                if str(file_name).endswith('genes.tsv'):
                    flag = True
            assert flag, f'Missed DefenseFinder "genes.tsv" file in {row[2]}'

            flag = False
            for file_name in Path(row[1]).iterdir():
                if str(file_name).endswith('padloc.csv'):
                    flag = True
            assert flag, f'Missed Padloc "padloc.csv" file in {row[1]}'

            assert Path(row[3]).exists(), f'Missed gff-file in {row[1]}'

    print('---Checking is complete')


def parse_csv_input(input_samples_file, output):
    padloc = []
    dfnfnr = []
    merging_list = []
    with open(input_samples_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            padloc.append((row[0], row[1]))
            dfnfnr.append((row[0], row[2], row[3]))
            merging_list.append((row[0],
                                row[2],
                                row[1]))

    with open(output / 'samples_for_merging.csv', mode='w') as f:
        for line in merging_list:
            f.write(f'{",".join(line)}\n')

    return padloc, dfnfnr


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

    df_ann[
        ['Accession', 'Nucleotide', 'Protein',
         'Padloc_Annotation', 'Padloc_System', 'Padloc_DS_ID',
         'DF_Annotation', 'DF_System', 'DF_System_sub', 'DF_DS_ID']
    ].to_csv(f'{results_dir}/Merged_annotations.tsv', sep='\t', index=False)

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


def run_defsys_rough(samples_list, pattern, output_path):

    print('\n >>> Run rough DefSys pipeline <<<')

    output_path.mkdir(parents=True, exist_ok=True)

    padloc_samples, dfnfnr_samples = parse_csv_input(samples_list, output_path)

    print('\n>>> Padloc processing')
    padloc_results = Path(output_path) / 'Summary_Padloc'
    padloc_results.mkdir(parents=True, exist_ok=True)
    process_padloc_data_rough(padloc_samples, padloc_results)

    print('\n>>> DefenseFinder processing')
    dfnfnr_results = Path(output_path) / 'Summary_DefenseFinder'
    dfnfnr_results.mkdir(parents=True, exist_ok=True)
    process_dfnfnr_data_rough(dfnfnr_samples, pattern, dfnfnr_results)

    print('\n>>> Merge annotations')
    _ = merge_annotations(
        padloc_results / 'protein_annotations.tsv',
        dfnfnr_results / 'protein_annotations.tsv',
        output_path
    )

    print('\n>>> Merge Padloc and DefenseFinder defense systems')
    merge_results = output_path / 'PDF_merge'
    merge_results.mkdir(parents=True, exist_ok=True)
    merge_padloc_and_defensefinder(
        output_path / 'samples_for_merging.csv',
        merge_results
    )

    print('\n>>> Make combined defense system summary')
    make_combined_summary(
        output_path / 'Merged_annotations.tsv',
        output_path / 'PDF_merge/By_Accessions',
        output_path / 'Summary_Padloc/defsys_summary.json',
        output_path / 'Summary_Padloc/redundant_defsys.tsv',
        output_path / 'Unique_accessions_Padloc.txt',
        output_path / 'Summary_DefenseFinder/defsys_summary.json',
        output_path / 'Summary_DefenseFinder/redundant_defsys.tsv',
        output_path / 'Unique_accessions_DefenseFinder.txt',
        output_path
    )

    print("That's all")


if __name__ == '__main__':
    args = parse_arguments()

    test_csv_input(args.samples_list)

    run_defsys_rough(args.samples_list,
                     args.regex,
                     args.output_path)
