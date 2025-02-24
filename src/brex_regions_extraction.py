"""v0.4 WIP
   No extraction functionality!
"""
import argparse
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


def get_proteins_numbers(frame):
    """
    Auxiliary function for groupby-agg
    Return only protein's number w/o nucleotide
    """
    return [int(_.split('_')[-1]) for _ in frame]


def create_defsys_summary(path_to_csv: str | Path) -> tuple:
    """
    Input: Padlock csv
    Output: DefSys summary table

    Ignores duplicate-systems in summary table
    The region direction is selected by the brxC direction or by simple voting

    Field 'DS_ID':
        x.split('%')[-1] contig/nucleotide name;
        x.split('%')[-2] Padlock number of the system
    """

    df = pd.read_csv(path_to_csv, dtype={'system.number': str},
                     usecols=[0, 1, 2, 3, 6, 11, 12, 13])

    # Create unique DefSystems IDs
    df['DS_ID'] = df['system'] + '%' + df['system.number'] + '%' + df['seqid']

    # Collection of systems with antiparallel genes if needed
    anti_defsys_ids = (df.loc[:, ('strand', 'DS_ID')]
                         .groupby('DS_ID')
                         .agg(lambda x: x.nunique())
                         .query('strand > 1')
                         .index)
    anti_defsys_names = [path_to_csv.stem[:15] + '_' + i for i in anti_defsys_ids]

    # Choose uniform strandness/direction for non-unidirectional DS
    if anti_defsys_names:
        for curr_defsys in df.DS_ID.unique():
            # Checking for BrxC
            have_brxc = (df['DS_ID'] == curr_defsys) & (df['protein.name'] == 'BrxC')
            # If present, select its direction
            if not df.loc[have_brxc].empty:
                df.loc[(df.loc[:, 'DS_ID'] == curr_defsys), ['strand']
                       ] = df.loc[have_brxc, 'strand'].values[0]
            else:
                # If not, by voting
                df.loc[(df.DS_ID == curr_defsys), ['strand']
                       ] = df.loc[(df.DS_ID == curr_defsys), ['strand']
                                  ].value_counts().index[0][0]

    # Log and ignore DS with proteins with duplicated DS assignment
    duplicated_prots = df['target.name'].value_counts()[df['target.name'].value_counts() > 1].index
    dupl_defsys_ids = df[df['target.name'].isin(duplicated_prots)].DS_ID.unique()
    dupl_defsys_names = [path_to_csv.stem[:15] + '_' + i for i in dupl_defsys_ids]

    # Add protein annotation
    df_proteins = df.iloc[:, [3, 4, 2, 1, 5, 6, 7]].copy()
    df_proteins.columns = ('Protein', 'Padloc_ann', 'System', 'Nucleotide', 'Start', 'End', 'Strand')
    df_proteins['Localisation'] = 'no'
    # --- Merge duplicated DS assignment
    df_dupl_merged = (df_proteins.loc[df_proteins.Protein.isin(duplicated_prots)]
                                 .groupby(
        ['Protein', 'Padloc_ann', 'Nucleotide', 'Start', 'End', 'Strand', 'Localisation']
                                         )
                                 .agg({'System': lambda x: f'DUPL_{",".join(x)}'})
                                 .reset_index()
                      )
    df_proteins_fin = pd.concat([df_proteins.loc[~df_proteins.Protein.isin(duplicated_prots)],
                                 df_dupl_merged])

    df = df[~df.DS_ID.isin(dupl_defsys_ids)]

    # Summary table
    df_result = df.loc[:, ['DS_ID', 'start', 'end', 'strand']].groupby('DS_ID').agg(coords_selector)
    df_result = df_result.astype({'start': 'int32',
                                  'end': 'int32'})
    # --- Get only proteins numbers of all DS
    df_result['Target_DS_Prots'] = df.loc[:, ['DS_ID', 'target.name']].groupby('DS_ID').agg(get_proteins_numbers)

    # --- Check presence of any inner genes in every DS
    df_result['Have_inner'] = df_result.Target_DS_Prots.apply(lambda x: set(x) != set(range(min(x), max(x) + 1)))

    # --- Find intersections of DS having inners with DS without inners
    df_result.reset_index(inplace=True)

    merged = df_result.loc[df_result['Have_inner']].merge(df_result.loc[~df_result['Have_inner']], how='cross')
    merged_inners = merged.query('(start_x <= start_y <= end_x) or (start_x <= end_y <= end_x)')

    df_result.set_index('DS_ID', inplace=True)

    df_result['Inner_DS'] = merged_inners.loc[:, ['DS_ID_x', 'DS_ID_y']].groupby('DS_ID_x').agg(
        lambda x: x.unique().tolist())

    df_result.fillna({'Inner_DS': ''}, inplace=True)

    # For DS having inners write proteins numbers of inner DS
    df_result['Inner_DS_Prots'] = None
    for sys_id, defsys in df_result.loc[df_result['Inner_DS'] != '', 'Inner_DS'].items():
        inner_defsys_prots = [x for xs in df_result.loc[defsys, 'Target_DS_Prots'] for x in xs]
        df_result.at[sys_id, 'Inner_DS_Prots'] = inner_defsys_prots

    df_result.rename({'start': 'Start', 'end': 'End', 'strand': 'Strand'}, axis=1, inplace=True)

    return df_result, dupl_defsys_names, anti_defsys_names, df_proteins_fin

"""
input_folder_path = Path('/home/holydiver/Main/2024_BREX/Data/TMP/Exmpls/Padloc_big')
assert input_folder_path.exists(), 'Input folder does not exist!'
target_defsys = 'brex'
output_path_results = Path('/home/holydiver/Main/2024_BREX/Data/New_data')
"""

args = parse_arguments()
input_folder_path = Path(args.input_folder_path)
assert input_folder_path.exists(), 'Input folder does not exist!'
target_defsys = args.ds
output_path_results = Path(args.output_path_results)

duplicate_defsys = []
antiparallel_defsys = []  # List of nucleotides w/o DS of interest
no_target_defsys_genomes = []  # List of nucleotides w/o DS of interest
protein_annotations = pd.DataFrame(columns=('Protein', 'Padloc_ann',
                                            'System', 'Nucleotide',
                                            'Start', 'End',
                                            'Strand', 'Localisation'))

for folder in input_folder_path.iterdir():
    print(f'---Processing {folder.name}---')
    df_all, dupl_curr, anti_curr, prot_curr = create_defsys_summary(folder / (folder.name + '_padloc.csv'))

    curr_accession_folder = output_path_results / folder.name
    curr_accession_folder.mkdir(parents=True, exist_ok=True)

    df_all.to_csv(curr_accession_folder / (folder.name + '_DS_summary.tsv'), sep='\t')

    # Select only DefSys of interest
    df_selected = df_all.loc[df_all.index.str.startswith(target_defsys)]

    if df_selected.empty:
        no_target_defsys_genomes.append(folder.name)
    else:
        df_selected.to_csv(curr_accession_folder / (folder.name + f'_{target_defsys}.tsv'), sep='\t')

    if dupl_curr:
        duplicate_defsys.extend(dupl_curr)
    if anti_curr:
        antiparallel_defsys.extend(anti_curr)

    protein_annotations = pd.concat([protein_annotations, prot_curr])

logs = (duplicate_defsys, antiparallel_defsys, no_target_defsys_genomes)

f_names = ('duplicate_defsys.txt',
           'antiparallel_defsys.txt',
           'no_target_defsys_genomes.txt')

for log, f_name in zip(logs, f_names):
    with open(output_path_results / f_name, mode='w') as f:
        for i in log:
            f.write(i + '\n')

protein_annotations.to_csv(output_path_results / 'protein_annotations.tsv', sep='\t', index=False)
