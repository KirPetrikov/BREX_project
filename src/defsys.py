"""v0.2
Scripts collection for import
"""
import numpy as np
import pandas as pd
import re

from collections import defaultdict
from itertools import combinations

pd.options.mode.copy_on_write = True


def co_occurence_matrix(data_items: list | tuple, to_dataframe=True) -> np.ndarray | pd.DataFrame:
    """
    Create co-occurence matrix for collection of items
    Return numpy array or pandas dataframe

    :param data_items: List of lists with items
    :param to_dataframe: Set True (default) to return pd.DataFrame, else - np.array
    :return: Co-occurence array or dataframe
    """
    all_items = sorted(set(x for xs in data_items for x in xs))

    co_occurrence = defaultdict(lambda: defaultdict(int))

    for curr_item in data_items:
        for item_1, item_2 in combinations(sorted(set(curr_item)), 2):
            co_occurrence[item_1][item_2] += 1
            co_occurrence[item_2][item_1] += 1

        # Self co-occurrence
        for ds in curr_item:
            co_occurrence[ds][ds] += 1

    all_items = list(all_items)
    co_occurrence_mtx = np.zeros((len(all_items), len(all_items)))

    codes_idxs = {}
    for num, val in enumerate(all_items):
        codes_idxs[val] = num

    for row_idx, in_val in co_occurrence.items():
        for col_idx, c_val in in_val.items():
            co_occurrence_mtx[
                codes_idxs[row_idx], codes_idxs[col_idx]
                              ] = c_val

    if to_dataframe:
        df = (pd.DataFrame(co_occurrence_mtx.astype(int))
              .set_axis(all_items, axis=1)
              .set_axis(all_items, axis=0))

        return df

    return co_occurrence_mtx


def top_co_occurred(co_occurrence_mtx: np.ndarray | pd.DataFrame,
                    n_top: int = 10) -> dict[tuple: int] | dict[str: int]:
    """
    Create dictionary with top n co-occured pairs of items

    :param co_occurrence_mtx: Co-occurence matrix
    :param n_top: number of top entries to get
                  If set to -1 - get all pairs with non-zero co-occurrence
    :return: Pairs of ID's with them co-occurrence value
    """

    dataframe = False

    if isinstance(co_occurrence_mtx, pd.DataFrame):
        items = co_occurrence_mtx.index
        co_occurrence_mtx = co_occurrence_mtx.values
        dataframe = True

    triu_indices = np.triu_indices_from(co_occurrence_mtx, k=1)
    upper_values = co_occurrence_mtx[triu_indices]

    if n_top == -1:
        n_top = np.count_nonzero(upper_values)

    top_n_idxs = np.argpartition(-upper_values, n_top)[:n_top]
    top_n_sorted = top_n_idxs[np.argsort(-upper_values[top_n_idxs])]

    top_n_pairs = [(triu_indices[0][i], triu_indices[1][i]) for i in top_n_sorted]

    top_n_values = [co_occurrence_mtx[i, j] for i, j in top_n_pairs]

    top_n_fin = {}

    if dataframe:
        for i, j in zip(top_n_pairs, top_n_values):
            curr_pair = tuple(items[list(i)])
            top_n_fin[curr_pair] = int(j)
    else:
        for i, j in zip(top_n_pairs, top_n_values):
            top_n_fin[i] = int(j)

    return top_n_fin


def find_dupl_defsys(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select rows corresponding to DS where there is protein duplications

    Params:
    - dataframe: Padloc csv

    Return:
    - dataframe: only duplicated DS
    """

    dupl_prots = df.Protein.value_counts()[df.Protein.value_counts() > 1].index
    dupl_defsys_ids = df[df.Protein.isin(dupl_prots)].DS_ID.unique()
    df_dupl = df[df.DS_ID.isin(dupl_defsys_ids)]

    return df_dupl


def coords_selector(frame: pd.DataFrame):
    """
    Auxiliary function for groupby-agg in defsys_summary()
    Selects the coordinates of the DS's region
    """
    if frame.name == 'Start':
        return min(frame)
    elif frame.name == 'End':
        return max(frame)
    else:
        return frame.unique()[0]


def create_defsys_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary for defense systems

    Params:
    - dataframe: table of DFs' proteins.
                 Must contain cols:
                 ('DS_ID', 'Accession', 'Nucleotide', 'Protein', 'System', 'Start', 'End', 'Strand')
                 The strandness in each system must be uniform.

     Return:
    - dataframe: Every row - unique DS; index - 'DS_ID'.
                 Columns:
                 - 'Accession', 'Nucleotide', 'DS_ID', 'Strand': corr. values;
                 - 'Start', 'End': DS's region boundary coordinates;
                 - 'DS_Prots': list, IDs of DS's proteins
                 - 'Have_inner': True/False
    """

    df_result = (
        df.loc[:, ['DS_ID', 'Accession', 'Nucleotide', 'System', 'Start', 'End', 'Strand']]
          .groupby('DS_ID')
          .agg(coords_selector)
    )

    # Get only proteins numbers of all DS
    df_result['DS_Prots'] = (
        df.loc[:, ['DS_ID', 'Protein']]
          .groupby('DS_ID')
          .agg(
            lambda x: [int(i.split('_')[-1]) for i in x]
          )
    )

    # Check presence of any inner genes in every DS
    df_result['Have_inner'] = df_result['DS_Prots'].apply(
        lambda x: set(x) != set(
            range(min(x), max(x) + 1)
        )
    )

    return df_result


def select_target_dupl_defsys(df: pd.DataFrame, target_defsys: str, dupl_ready: bool = False) -> pd.DataFrame:
    """
    Select those DS which intersect with target DS by duplicatedt proteins

    - dupl_ready: Set to True if input df contains not only duplicates
    """
    if not dupl_ready:
        df = find_dupl_defsys(df)

    # Mask for proteins from target DS
    dupl_prots_from_target_defsys = (df.groupby('Protein')['DS_ID']
                                       .agg(
        lambda x: any([i.startswith(target_defsys) for i in x.unique()])
                                            )
                                     )
    # Select proteins
    duplicated_prots_target = dupl_prots_from_target_defsys[dupl_prots_from_target_defsys].index

    # Select IDs of all DS which contain duplicated proteins
    dupl_defsys_ids_target = df[df['Protein'].isin(duplicated_prots_target)].DS_ID.unique()

    # Select full DS
    df_dupl_target = df[df.DS_ID.isin(dupl_defsys_ids_target)]

    return df_dupl_target


def create_dupl_proteins_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary for duplicated proteins only

    Params:
    - df: Combined table of defense systems' proteins
          (as parsed Padloc-csv table with added assembly accessions)

    Return:
    - dataframe: Every row - unique duplicated protein.
                 Columns:
                 - 'Protein', 'Nucleotide', 'Accession': corresponding values;
                 - 'Padloc_ann', 'System': collected as lists in corresponding columns;
                 - 'DS_ID': all DS_IDs joined by ","
    """

    prots_counts = df.Protein.value_counts()
    dupl_prots = prots_counts[prots_counts > 1].index
    df = (df[
             df.Protein.isin(dupl_prots)
             ].groupby(['Protein', 'Nucleotide', 'Accession'])[['DS_ID', 'Padloc_ann', 'System']]
              .agg(lambda x: x.tolist())
              .reset_index()
          )
    df.loc[:, 'DS_ID'] = df.DS_ID.apply(lambda x: ','.join(x))
    return df


def parse_defense_finder_tsv(path_to_tsv, /,
                             short: bool = True,
                             add_accession='') -> pd.DataFrame:
    """
    Reads DefenseFinder systems tsv in a convenient form.
    Adds column 'DS_ID' with unique ID to prevent merging of DS with the same name.
    Add column 'Accession' with assembly accessions if
    path to table 'Nucleotide-Accession' provided

    Params:
    - path_to_tsv
    - short: If 'True' returns only selected columns.
    - add_accession: If provided path add 'Accession'

    Return:
    - dataframe
    """

    df = pd.read_csv(path_to_tsv,
                     sep='\t',
                     dtype={'sys_id': str, 'type': str, 'subtype': str,
                            'activity': str, 'sys_beg': str, 'sys_end': str,
                            'protein_in_syst': str, 'genes_count': int,
                            'name_of_profiles_in_sys': str, 'replicon': str},
                     names=['tmp_ID', 'DF_System', 'DF_System_sub',
                            'activity', 'sys_beg', 'sys_end',
                            'All_Proteins', 'genes_count', 'All_Annot', 'Nucleotide'],
                     skiprows=1
                     )

    # Create unique DS_IDs
    df['DS_ID'] = df.apply(
        lambda x: f'{x.DF_System_sub}%{x.tmp_ID.split("_")[-1]}%{x.Nucleotide}',
        axis=1
    )

    cols_for_short = ['DS_ID', 'DF_System', 'DF_System_sub',
                      'All_Proteins', 'All_Annot', 'Nucleotide']

    if add_accession:
        df_acc = pd.read_csv(add_accession, sep='\t')
        df = df.merge(df_acc[['Nucleotide', 'Accession']], on='Nucleotide', how='left')
        cols_for_short.append('Accession')

    if short:
        df = df[cols_for_short]

    return df


def parse_padloc_csv(path_to_csv, short: bool = True) -> pd.DataFrame:
    """
    Reads Padloc csv in a convenient form.
    Adds columns 'DS_ID' with unique ID to prevent merging of DS with the same name.

    Params:
    - path_to_csv
    - short: If 'True' returns only selected columns.

    Return:
    - dataframe
    """

    df = pd.read_csv(path_to_csv,
                     names=['SysNo', 'Nucleotide', 'System', 'Protein',
                            'HMM_profile', 'HMM_prot_name', 'Padloc_ann', 'Eval_full',
                            'Eval_i', 'Cov_target', 'Cov_hmm', 'Start', 'End', 'Strand',
                            'Description', 'Protein_ID', 'Contig_end', 'All'],
                     index_col=None,
                     dtype={'SysNo': str, 'Nucleotide': str, 'System': str, 'Protein': str,
                            'HMM_profile': str, 'HMM_prot_name': str, 'Padloc_ann': str,
                            'Eval_full': float, 'Eval_i': float, 'Cov_target': float, 'Cov_hmm': float,
                            'Start': int, 'End': int, 'Strand': str, 'Description': str,
                            'Protein_ID': int, 'Contig_end': int, 'All': str},
                     usecols=[_ for _ in range(18)],
                     skiprows=1
                     )

    # Create unique DefSystems IDs
    df['DS_ID'] = df['System'] + '%' + df['SysNo'] + '%' + df['Nucleotide']

    if short:
        return df.iloc[:, [0, 1, 2, 3, 6, 11, 12, 13, 18]]
    else:
        return df


def make_unidir_genes_defsys(df: pd.DataFrame):
    """
    Modifies table by choosing uniform strandness/direction of genes within every DS.
    By voting, or '+' in case of equality
    """

    nonunidir_defsys_ids = (df.loc[:, ('Strand', 'DS_ID')]
                            .groupby('DS_ID')['Strand']
                            .agg(pd.Series.nunique)
                            .loc[lambda x: x > 1]
                            .index)

    if not nonunidir_defsys_ids.empty:
        strand = (
            df.loc[df.DS_ID.isin(nonunidir_defsys_ids)]
              .groupby(['DS_ID'])['Strand']
              .agg(
                lambda x: x.value_counts().sort_index().idxmax()
              )
        )

        df.loc[df.DS_ID.isin(nonunidir_defsys_ids), 'Strand'] = (
            df.loc[df.DS_ID.isin(nonunidir_defsys_ids), 'DS_ID'].map(strand)
        )


def parse_hmmer_domtblout(input_dom_path) -> pd.DataFrame:
    """
    Read HMMER-format domtblout-file

    Params
    - input_dom_path

    Return
    - dataframe
    """
    colnames = ['Protein', 'Target_accession', 'Target_len', 'HMM_prot_name', 'HMM_profile',
                'HMM_len', 'Eval_full', 'Score_full', 'Bias_full', 'N_dom', 'Total_dom',
                'Eval_c', 'Eval_i', 'Score_dom', 'Bias_dom', 'HMM_start', 'HMM_end',
                'Ali_start', 'Ali_end', 'Envelope_start', 'Envelope_end', 'Prob_acc', 'Description']
    
    df = pd.read_csv(input_dom_path,
                     sep=r'[ ]{1,}',
                     names=colnames,
                     comment='#',
                     dtype={0: str, 1: str, 2: int, 3: str, 4: str, 5: int, 6: float,
                            7: float, 8: float, 9: int, 10: int, 11: float, 12: float,
                            13: float, 14: float, 15: int, 16: int, 17: int, 18: int,
                            19: int, 20: int, 21: float, 22: str},
                     engine='python')

    return df


def parse_prodigal_gff(path_to_gff,
                       add_id: bool = True,
                       with_nucl: bool = True) -> pd.DataFrame:
    """
    Parse Prodigal gff-file to pandas DataFrame.
    Can add gene id as just number or as number with nucleotide preffix
    """
    gff_cols_names = ('Chrom', 'Sourse', 'Feature', 'Start', 'End',
                      'Score', 'Strand', 'Frame', 'Comment')
    df = pd.read_csv(path_to_gff,
                     sep='\t',
                     names=gff_cols_names,
                     dtype={'Chrom': str, 'Sourse': str, 'Feature': str,
                            'Start': int, 'End': int, 'Score': float,
                            'Strand': str, 'Frame': str, 'Comment': str}
                     )
    if add_id:
        pattern = re.compile(r'ID=\d+_(\d+)')
        df['Gene_ID'] = df.Comment.apply(lambda x: int(re.search(pattern, x).group(1)))
        if with_nucl:
            df = df.astype({'Gene_ID': str})
            df.loc[:, 'Gene_ID'] = df.Chrom + '_' + df.Gene_ID

    return df
