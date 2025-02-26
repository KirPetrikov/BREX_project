"""v0.4a WIP
   No extraction functionality
"""
import argparse
import json
import pandas as pd

from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Parse Padloc results, csv-files, to create defense systems summary.'
                                     )

    parser.add_argument('input_folder_path', type=str,
                        help='Path to the data folder')
    parser.add_argument('output_path_results', type=str,
                        help='Path to the results folder. Will be created if it does not exist')
    parser.add_argument('-ds', type=str, default='brex',
                        help='Name or prefix of the Padloc protection system.\nDefault: brex')

    return parser.parse_args()


def coords_selector(frame):
    """
    Auxiliary function for groupby-agg
    Selects the coordinates and strand of the DS's region
    """
    if frame.name == 'start':
        return min(frame)
    elif frame.name == 'end':
        return max(frame)
    elif frame.name == 'strand':
        return frame.unique()[0]


# def get_proteins_numbers(frame):
#     """
#     Auxiliary function for groupby-agg
#     Return only protein's number w/o nucleotide
#     """
#     return [int(i.split('_')[-1]) for i in frame]

def get_defsys_with_proteins(cell: list, df: pd.DataFrame) -> dict:
    """
    Auxiliary function for groupby-agg
    Return only protein's number w/o nucleotide
    """
    tmp_d = {}
    for ds in cell:
        tmp_d[ds] = df.loc[ds, 'DS_Prots']
    return tmp_d


def create_defsys_summary(path_to_csv: str | Path) -> tuple[pd.DataFrame, tuple, list, pd.DataFrame]:
    """
    Input: Padlock csv
    Output: DefSys summary table

    Ignores duplicate-systems in summary table
    The region direction is selected by the brxC direction or by simple voting

    Field 'DS_ID':
        x.split('%')[-1] contig/nucleotide name;
        x.split('%')[-2] Padlock number of the system
    """

    df_all = pd.read_csv(path_to_csv, dtype={'system.number': str},
                         usecols=[0, 1, 2, 3, 6, 11, 12, 13])

    # Create unique DefSystems IDs
    df_all['DS_ID'] = df_all['system'] + '%' + df_all['system.number'] + '%' + df_all['seqid']

    # Collection of systems with antiparallel genes if needed
    anti_defsys_ids = (df_all.loc[:, ('strand', 'DS_ID')]
                       .groupby('DS_ID')
                       .agg(lambda x: x.nunique())
                       .query('strand > 1')
                       .index)
    anti_defsys_names = [path_to_csv.stem[:15] + '_' + i for i in anti_defsys_ids]

    # Choose uniform strandness/direction for non-unidirectional DS
    if anti_defsys_names:
        for curr_defsys in df_all.DS_ID.unique():
            # Checking for BrxC
            have_brxc = (df_all['DS_ID'] == curr_defsys) & (df_all['protein.name'] == 'BrxC')
            # If present, select its direction
            if not df_all.loc[have_brxc].empty:
                df_all.loc[(df_all.loc[:, 'DS_ID'] == curr_defsys), ['strand']
                           ] = df_all.loc[have_brxc, 'strand'].values[0]
            else:
                # If not, by voting
                df_all.loc[(df_all.DS_ID == curr_defsys), ['strand']
                           ] = df_all.loc[(df_all.DS_ID == curr_defsys), ['strand']
                                          ].value_counts().sort_index().index[0][0]

    # Log duplicate DS assignment
    duplicated_prots = df_all['target.name'].value_counts()[df_all['target.name'].value_counts() > 1].index
    dupl_defsys_ids = df_all[df_all['target.name'].isin(duplicated_prots)].DS_ID.unique()
    if dupl_defsys_ids.size != 0:
        dupl_defsys_names = path_to_csv.stem[:15] + '\t' + '\t'.join([i for i in dupl_defsys_ids]),
    else:
        dupl_defsys_names = ()

    # --- Table of protein annotations ---
    df_proteins = df_all.iloc[:, [3, 4, 2, 1, 5, 6, 7]].copy()
    df_proteins.columns = ('Protein', 'Padloc_ann', 'System', 'Nucleotide', 'Start', 'End', 'Strand')
    df_proteins['Localisation'] = 'no'
    # Merge duplicate DS assignment
    df_dupl_merged = (df_proteins.loc[df_proteins.Protein.isin(duplicated_prots)]
                      .groupby(
        ['Protein', 'Padloc_ann', 'Nucleotide', 'Start', 'End', 'Strand', 'Localisation']
    )
                      .agg({'System': lambda x: f'DUPL_{",".join(x)}'})
                      .reset_index()
                      )
    df_proteins_fin = pd.concat([df_proteins.loc[~df_proteins.Protein.isin(duplicated_prots)],
                                 df_dupl_merged])

    # Ignore duplicate DS assignment
    df_all = df_all[~df_all.DS_ID.isin(dupl_defsys_ids)]
    if df_all.empty:
        return df_all, dupl_defsys_names, anti_defsys_names, df_proteins_fin

    # --- Summary table ---
    df_result: pd.DataFrame = df_all.loc[:, ['DS_ID', 'start', 'end', 'strand']].groupby('DS_ID').agg(coords_selector)
    df_result = df_result.astype({'start': 'int32',
                                  'end': 'int32'})
    # Get only proteins numbers of all DS
    df_result['DS_Prots'] = (df_all.loc[:, ['DS_ID', 'target.name']]
                                   .groupby('DS_ID')
                                   .agg(lambda x: [int(i.split('_')[-1]) for i in x]))
    df_result['Nucleotide'] = (df_all.loc[:, ['DS_ID', 'seqid']]
                                     .groupby('DS_ID')
                                     .agg(lambda x: x.unique()[0]))

    # Check presence of any inner genes in every DS
    df_result['Have_inner'] = df_result.DS_Prots.apply(lambda x: set(x) != set(range(min(x), max(x) + 1)))

    # Find intersections of DS having inners with DS without inners
    df_result.reset_index(inplace=True)

    # For different nucleotides there may be a false intersection of coordinates
    # Process them separately
    df_by_nucl = []
    for nucl in df_result.Nucleotide.unique():
        df_curr_nucl = df_result.loc[df_result.Nucleotide == nucl].copy()

        merged = df_curr_nucl.loc[df_curr_nucl['Have_inner']].merge(df_curr_nucl.loc[~df_curr_nucl['Have_inner']],
                                                                    how='cross')
        merged_inners = merged.query('(start_x <= start_y <= end_x) or (start_x <= end_y <= end_x)')

        df_curr_nucl.set_index('DS_ID', inplace=True)

        df_curr_nucl['Inner_DS'] = merged_inners.loc[:, ['DS_ID_x', 'DS_ID_y']].groupby('DS_ID_x').agg(
            lambda x: x.unique().tolist())

        df_by_nucl.append(df_curr_nucl.copy())

    df_final = pd.concat(df_by_nucl)
    df_final.fillna({'Inner_DS': ''}, inplace=True)

    # For DS having inners write proteins numbers of inner DS
    df_final['Inner_DS_Prots'] = ''
    df_final.loc[df_final['Inner_DS'] != '', 'Inner_DS'] = (
            df_final.loc[df_final['Inner_DS'] != '', 'Inner_DS']
                    .apply(lambda x: get_defsys_with_proteins(x, df_final))
                                                            )

    df_final.rename({'start': 'Start', 'end': 'End', 'strand': 'Strand'}, axis=1, inplace=True)

    return df_final, dupl_defsys_names, anti_defsys_names, df_proteins_fin


args = parse_arguments()
input_folder_path = Path(args.input_folder_path)
assert input_folder_path.exists(), 'Input folder does not exist!'
target_defsys = args.ds
output_path_results = Path(args.output_path_results)

summary_all_defsys = {}
summary_target_defsys = {}
duplicate_defsys = []
antiparallel_defsys = []  # List of DS with antiparallel genes
no_defsys_genomes = []  # List of accessions w/o DS (i.e. all DS were with duplicated assignment)
no_target_defsys = [] # List of accessions w/o target DS
protein_annotations = pd.DataFrame(columns=('Protein', 'Padloc_ann',
                                            'System', 'Nucleotide',
                                            'Start', 'End',
                                            'Strand', 'Localisation'))


for folder in input_folder_path.iterdir():
    print(f'---Processing {folder.name}---')
    df_summary, dupl_curr, anti_curr, prot_curr = create_defsys_summary(folder / (folder.name + '_padloc.csv'))

    curr_accession_folder = output_path_results / f'By_Accessions/{folder.name}'
    curr_accession_folder.mkdir(parents=True, exist_ok=True)

    protein_annotations = pd.concat([protein_annotations, prot_curr])
    if dupl_curr:
        duplicate_defsys.extend(dupl_curr)
    if anti_curr:
        antiparallel_defsys.extend(anti_curr)

    if df_summary.empty:
        no_defsys_genomes.append(folder.name)
    else:
        df_summary.to_json(curr_accession_folder / (folder.name + '_summary.json'), orient='index')
        summary_all_defsys.update(df_summary.to_dict(orient='index'))

        # Select only DefSys of interest
        df_selected = df_summary.loc[df_summary.index.str.startswith(target_defsys)]
        if df_selected.empty:
            no_target_defsys.append(folder.name)
        else:
            summary_target_defsys.update(df_selected.to_dict(orient='index'))
            df_selected.to_csv(curr_accession_folder / (folder.name + f'_{target_defsys}.tsv'), sep='\t')



logs = (duplicate_defsys, antiparallel_defsys, no_defsys_genomes, no_target_defsys)
f_names = ('duplicate_defsys.tsv',
           'antiparallel_defsys.txt',
           'no_defsys_genomes.txt',
           'no_target_defsys.txt')
for log, f_name in zip(logs, f_names):
    with open(output_path_results / f_name, mode='w') as f:
        for i in log:
            f.write(i + '\n')

protein_annotations.to_csv(output_path_results / 'protein_annotations.tsv', sep='\t', index=False)

with open(output_path_results / 'All_DefSys_summary.json', mode='w') as f:
    json.dump(summary_all_defsys, f, indent=4)

with open(output_path_results / f'{target_defsys}_summary.json', mode='w') as f:
    json.dump(summary_target_defsys, f, indent=4)
