"""v0.4
Create combined summary by choosing unique DS id.
In case of ambiguity priority choice is DefenseFinder variant.
Considers DSs, which are joined in combined output table ('DS1/DS2'), and just selects first.
Checks protein redundancy.
Also takes lists of accessios uniquely processed by each tool, and adds them to combined summary.

Takes , mergeded protein annotations table, lists of Padloc and DF unique accessions,
lists of Padloc and DF redundant defsystems, 
unified json-summaries of Padloc and DF, tables with their merged results.

Output:
- "Combined_summary.json"
- "Combined_redundant_defsys.tsv"
- "Joined_defsys_from_merged_results.tsv"
"""
import argparse
import json
import pandas as pd

from pathlib import Path
from commons import find_redundancy_defsys


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Makes combined summary from Padloc and DefenseFinder merged tables.'
    )

    parser.add_argument('-a', '--input_ann_path', type=Path, required=True,
                        help='Path to merged annotation file')
    parser.add_argument('-t', '--input_merged_table_path', type=Path, required=True,
                        help='Path to combined Padloc-DefenceFinder tables')

    parser.add_argument('-p', '--input_padloc_summary_path', type=Path, required=True,
                        help='Path to Padloc summary JSON file')
    parser.add_argument('-l', '--input_padloc_redund_path', type=Path, required=True,
                        help='Path to Padloc redundant systems file')
    parser.add_argument('-c', '--input_padloc_uniq_path', type=Path, required=True,
                        help='Path to unique Padloc accessions file')

    parser.add_argument('-d', '--input_dfnfnr_summary_path', type=Path, required=True,
                        help='Path to DefenseFinder summary JSON file')
    parser.add_argument('-f', '--input_dfnfnr_redund_path', type=Path, required=True,
                        help='Path to DefenseFinder redundant systems file')
    parser.add_argument('-n', '--input_dfnfnr_uniq_path', type=str, required=True,
                        help='Path to unique DefenseFinder accessions file')

    parser.add_argument('-o', '--output_results_path', type=Path, required=True,
                        help='Output results directory')

    return parser.parse_args()


def make_ids_dicts_and_unq(
    ann_path,
    padloc_uniq_path,
    dfnfnr_uniq_path
) -> tuple[dict, dict, dict, dict]:

    print('---Tmp-ID dictionary making---')

    # --- Add unique DS_ID for each DS/row in ann-table
    df_ann = (
        pd.read_csv(
            ann_path,
            sep='\t',
            dtype={
                'Accession': str, 'Protein': str, 'Nucleotide': str,
                'Padloc_ann': str, 'Padloc_System': str,
                'DF_ann': str, 'DF_System': str, 'DF_System_sub': str,
                'Padloc_DS_ID': str, 'DF_DS_ID': str
            }
        )
    )

    padloc_uniq_ds_accs = []
    with open(padloc_uniq_path) as file:
        for line in file:
            padloc_uniq_ds_accs.append(line.strip())
    padloc_uniq = df_ann.loc[df_ann.Accession.isin(padloc_uniq_ds_accs), 'Padloc_DS_ID']

    dfnfnr_uniq_ds_accs = []
    with open(dfnfnr_uniq_path) as file:
        for line in file:
            dfnfnr_uniq_ds_accs.append(line.strip())
    dfnfnr_uniq = df_ann.loc[df_ann.Accession.isin(dfnfnr_uniq_ds_accs), 'DF_DS_ID']

    dfnfnr_mask = (df_ann.DF_DS_ID != 'miss').values
    padloc_mask = (df_ann.Padloc_System != 'miss').values

    # "Tmp-ids - DS_IDs" correspondence for DF
    df_ann['df_tmp_id'] = df_ann.loc[dfnfnr_mask].apply(
        lambda x: f'{x.Protein}%{x.DF_System_sub}',
        axis=1
    )
    df_dict = dict(zip(
        df_ann.loc[dfnfnr_mask].df_tmp_id,
        df_ann.loc[dfnfnr_mask].DF_DS_ID
    ))

    # "Tmp-ids - DS_IDs" correspondence for Padloc
    df_ann['p_tmp_id'] = df_ann.loc[padloc_mask].apply(
        lambda x: f'{x.Protein}%{x.Padloc_System}',
        axis=1
    )
    padloc_dict = dict(zip(
        df_ann.loc[padloc_mask].p_tmp_id,
        df_ann.loc[padloc_mask].Padloc_DS_ID
    ))

    del df_ann

    return padloc_dict, df_dict, padloc_uniq, dfnfnr_uniq


def process_merged_table(
    table_path,
    padloc_dict,
    df_dict
):
    df_pdfc = pd.read_csv(
        table_path,
        sep=',',
        dtype={'DF_System': str, 'Padloc_System': str, 'DF_System_sub': str,
               'Padloc_System_sub': str, 'Proteins': str, 'DF_anns': str, 'Padloc_anns': str},
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
    pdfc_dfnfnr_mask = (df_pdfc.DF_System != 'N.A.').values
    # Tmp-ids for Padloc
    df_pdfc['tmpid'] = ''
    df_pdfc.loc[~pdfc_dfnfnr_mask, 'tmpid'] = df_pdfc.loc[~pdfc_dfnfnr_mask].apply(
        lambda x: f'{x.Proteins.split(";")[0]}%{x.Padloc_System_sub}',
        axis=1
    )

    # Tmp-ids for DF
    df_pdfc.loc[pdfc_dfnfnr_mask, 'tmpid'] = df_pdfc.loc[pdfc_dfnfnr_mask].apply(
        lambda x: f'{x.Proteins.split(";")[0]}%{x.DF_System_sub}',
        axis=1
    )

    df_pdfc['DS_ID'] = ''
    padloc_sel_ids = [padloc_dict[i] for i in df_pdfc.loc[~pdfc_dfnfnr_mask, 'tmpid'].values]
    df_pdfc.loc[~pdfc_dfnfnr_mask, 'DS_ID'] = padloc_sel_ids
    dfnfnr_sel_ids = [df_dict[i] for i in df_pdfc.loc[pdfc_dfnfnr_mask, 'tmpid'].values]
    df_pdfc.loc[pdfc_dfnfnr_mask, 'DS_ID'] = dfnfnr_sel_ids

    # --- Find redundant proteins
    prots_ids = {'Protein': [], 'DS_ID': []}

    tmp_prots = dict(zip(df_pdfc.DS_ID, df_pdfc.Proteins))
    for ds_id, prots in tmp_prots.items():
        curr_prots = prots.split(';')
        prots_ids['Protein'].extend(curr_prots)
        prots_ids['DS_ID'].extend([ds_id] * len(curr_prots))

    redund_ds_ids = (
        find_redundancy_defsys(
            pd.DataFrame(prots_ids)
        ).DS_ID
         .unique()
         .tolist()
    )

    if idxs_to_split:
        splitted_ds_ids = df_pdfc.loc[idxs_to_split, 'DS_ID'].to_list()
    else:
        splitted_ds_ids = []

    df_pdfc = df_pdfc.dropna()
    if defsys_number != df_pdfc.shape[0]:
        print(f'### Warning: DS missed in {table_path.stem}')

    return padloc_sel_ids, dfnfnr_sel_ids, redund_ds_ids, splitted_ds_ids


def make_combined_summary(
    ann_path,
    merged_table_path: Path,
    padloc_summary_path,
    padloc_redund_path,
    padloc_uniq_path,
    dfnfnr_summary_path,
    dfnfnr_redund_path,
    dfnfnr_uniq_path,
    results_path: Path
) -> None:
    padloc_sel_ids = []
    dfnfnr_sel_ids = []
    redundant_ds = {'Accession': [], 'DS_ID': []}
    splitted_ds = {'Accession': [], 'DS_ID': []}

    results_path.mkdir(parents=True, exist_ok=True)

    padloc_ids_code, dfnfnr_ids_code, padloc_uniq_ds_ids, dfnfnr_uniq_ds_ids = make_ids_dicts_and_unq(
        ann_path,
        padloc_uniq_path,
        dfnfnr_uniq_path
    )

    # --- Select unique DS_ID from all combined tables
    print('---Processing merged tables---')
    for pdf_merged_table in merged_table_path.iterdir():
        print(f'---Processing {pdf_merged_table.name}---')

        curr_padloc_sel_ids, curr_df_sel_ids, curr_redund_ds_ids, curr_splitted_ds_ids = process_merged_table(
            pdf_merged_table,
            padloc_ids_code,
            dfnfnr_ids_code
        )

        # Update final lists of DS_ID
        padloc_sel_ids.extend(curr_padloc_sel_ids)
        dfnfnr_sel_ids.extend(curr_df_sel_ids)

        splitted_ds['DS_ID'].extend(curr_splitted_ds_ids)
        splitted_ds['Accession'].extend([pdf_merged_table.stem] * len(curr_splitted_ds_ids))

        redundant_ds['DS_ID'].extend(curr_redund_ds_ids)
        redundant_ds['Accession'].extend([pdf_merged_table.stem] * len(curr_redund_ds_ids))

    print('---Merged tables processing completed---')

    results_path.mkdir(parents=True, exist_ok=True)

    # --- --- Save joined/splitted defsys
    pd.DataFrame(splitted_ds).to_csv(
        results_path / 'Joined_defsys_from_merged_results.tsv',
        sep='\t',
        index=False
    )

    # --- --- Add unique redundant defsys and save
    df_padloc_redund = pd.read_csv(padloc_redund_path, sep='\t')
    df_padloc_redund = df_padloc_redund.loc[df_padloc_redund.DS_ID.isin(padloc_uniq_ds_ids)]
    redundant_ds['DS_ID'].extend(df_padloc_redund.DS_ID.to_list())
    redundant_ds['Accession'].extend(df_padloc_redund.Accession.to_list())

    df_dfnfnr_redund = pd.read_csv(dfnfnr_redund_path, sep='\t')
    df_dfnfnr_redund = df_dfnfnr_redund.loc[df_dfnfnr_redund.DS_ID.isin(padloc_uniq_ds_ids)]
    redundant_ds['DS_ID'].extend(df_dfnfnr_redund.DS_ID.to_list())
    redundant_ds['Accession'].extend(df_dfnfnr_redund.Accession.to_list())

    pd.DataFrame(redundant_ds).to_csv(
        results_path / 'Combined_redundant_defsys.tsv',
        sep='\t',
        index=False
    )

    combined_summary = {}

    # --- --- Add Padloc to combined summary
    with open(padloc_summary_path) as f:
        padloc_data = json.load(f)
    for ds_id in padloc_sel_ids:
        combined_summary[ds_id] = padloc_data[ds_id]
        combined_summary[ds_id]['System_sub'] = combined_summary[ds_id]['System']
        combined_summary[ds_id]['System'] = combined_summary[ds_id]['System_sub'].split('_')[0]

    # --- --- Add DF to combined summary
    with open(dfnfnr_summary_path) as f:
        dfnfnr_data = json.load(f)
    for ds_id in dfnfnr_sel_ids:
        combined_summary[ds_id] = dfnfnr_data[ds_id]

    # --- --- Add unique DS ids
    for ds_id in padloc_uniq_ds_ids:
        combined_summary[ds_id] = padloc_data[ds_id]
        combined_summary[ds_id]['System_sub'] = combined_summary[ds_id]['System']
        combined_summary[ds_id]['System'] = combined_summary[ds_id]['System_sub'].split('_')[0]

    for ds_id in dfnfnr_uniq_ds_ids:
        combined_summary[ds_id] = dfnfnr_data[ds_id]

    with open(results_path / 'Combined_summary.json', mode='w') as f:
        json.dump(combined_summary, f, indent=4)

    print('---Making combined summary is complet---')


if __name__ == '__main__':
    args = parse_arguments()

    make_combined_summary(args.input_ann_path, args.input_merged_table_path, args.input_padloc_summary_path,
                          args.input_padloc_redund_path, args.input_padloc_uniq_path, args.input_dfnfnr_summary_path,
                          args.input_dfnfnr_redund_path, args.input_dfnfnr_uniq_path, args.output_results_path)
