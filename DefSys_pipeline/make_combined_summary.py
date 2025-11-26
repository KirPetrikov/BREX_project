import json
import os
import sys

import pandas as pd

from pathlib import Path
from commons import find_redundancy_defsys


def make_ids_dicts_and_unq(
    ann_path,
    padloc_uniq_path,
    dfnfnr_uniq_path
) -> tuple[dict, dict, dict, dict]:
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

    # Select final DS_ID from merged tables to add in final combined summary
    print('>>> Processing merged tables')
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

    results_path.mkdir(parents=True, exist_ok=True)

    print('>>> Add DS to final combined summary')
    combined_summary = {}

    ref_file = Path(os.path.dirname(sys.argv[0]), 'immune_system_list_reference.json')
    with open(ref_file) as f:
        ref_dict = json.load(f)

    # Make final summary
    # --- Add Padloc DSs to combined summary
    with open(padloc_summary_path) as f:
        padloc_data = json.load(f)
    for ds_id in padloc_sel_ids:
        combined_summary[ds_id] = padloc_data[ds_id]
        # Unify DS names
        # Padloc DS names need splittitg
        old_sys_name: str = combined_summary[ds_id]['System'].split('_')[0]
        combined_summary[ds_id]['System_sub'] = combined_summary[ds_id]['System'].replace(
            old_sys_name, ref_dict[old_sys_name]
        )
        combined_summary[ds_id]['System'] = ref_dict[old_sys_name]
        combined_summary[ds_id]['Activity'] = 'miss'

    # --- Add DF to combined summary

    print('\n\n')

    with open(dfnfnr_summary_path) as f:
        dfnfnr_data = json.load(f)
    for ds_id in dfnfnr_sel_ids:
        combined_summary[ds_id] = dfnfnr_data[ds_id]
        # Unify DS names
        # DFs DS names do not need splitting
        old_sys_name: str = combined_summary[ds_id]['System']
        combined_summary[ds_id]['System_sub'] = combined_summary[ds_id]['System_sub'].replace(
            old_sys_name, ref_dict[old_sys_name]
        )
        combined_summary[ds_id]['System'] = ref_dict[old_sys_name]

    # --- --- Add unique DS ids
    for ds_id in padloc_uniq_ds_ids:
        combined_summary[ds_id] = padloc_data[ds_id]
        combined_summary[ds_id]['System_sub'] = combined_summary[ds_id]['System']
        combined_summary[ds_id]['System'] = combined_summary[ds_id]['System_sub'].split('_')[0]
        combined_summary[ds_id]['Activity'] = 'miss'

    for ds_id in dfnfnr_uniq_ds_ids:
        combined_summary[ds_id] = dfnfnr_data[ds_id]

    # Save results
    # --- Save combined summary
    with open(results_path / 'Combined_summary.json', mode='w') as f:
        json.dump(combined_summary, f, indent=4)

    # --- Save joined/splitted defsys
    pd.DataFrame(splitted_ds).to_csv(
        results_path / 'Joined_defsys_from_merged_results.tsv',
        sep='\t',
        index=False
    )

    # --- Add unique redundant defsys and save
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

    print('--- Making combined summary is complet')
