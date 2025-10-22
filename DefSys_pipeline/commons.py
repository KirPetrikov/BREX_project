"""v0.3p
Common scripts for pipeline
"""
import pandas as pd
import re

pd.options.mode.copy_on_write = True


def find_redundancy_defsys(df: pd.DataFrame) -> pd.DataFrame:
    """
    Finds redundancy of proteins annotation to DS:
    when the same protein was annotated in different DS

    Params:
    - dataframe: unifyed DS-table

    Return:
    - dataframe: only rows corresponding to DS where there are redundant proteins
    """

    dupl_prots = df.Protein.value_counts()[df.Protein.value_counts() > 1].index
    dupl_defsys_ids = df[df.Protein.isin(dupl_prots)].DS_ID.unique()
    df_redund = df[df.DS_ID.isin(dupl_defsys_ids)]

    return df_redund


def make_unidir_genes_defsys(df: pd.DataFrame) -> None:
    """
    Inplace modifies unifyed DS-table by choosing
    uniform strandness/direction of genes within every DS.
    By simple voting, or '+' in case of equality.
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


def create_defsys_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary for defense systems from unifyed DS-table

    Params:
    - dataframe: unifyed DS-table.
                 Must contain cols:
                 ('DS_ID', 'Accession', 'Nucleotide', 'Protein', 'System', 'Start', 'End', 'Strand')
                 The strandness in each system must be uniform.

     Return:
    - dataframe: Every row - unique DS; index - 'DS_ID'.
                 Columns:
                 - 'Accession', 'Nucleotide', 'DS_ID', 'Strand': corr. values;
                 - 'Start', 'End': DS's region boundary coordinates;
                 - 'DS_Prots': list, only numbers of proteins' IDs
                 - 'Have_inner': True/False
    """

    def coords_selector(frame: pd.DataFrame) -> str:
        """
        Auxiliary function for groupby-agg
        Selects boundary coordinates of the DS's region
        """
        if frame.name == 'Start':
            return min(frame)
        elif frame.name == 'End':
            return max(frame)
        else:
            return frame.unique()[0]

    df_result = (
        df.groupby('DS_ID')
          .agg(coords_selector)
    )

    # Get only proteins numbers of all DS
    df_result['DS_Prots'] = (
        df.loc[:, ['DS_ID', 'Protein']]
          .groupby('DS_ID')
          .agg(
            lambda x: [
                int(i.split('_')[-1]) for i in x
            ]
          )
    )

    # Check presence of any inner genes in every DS
    df_result['Have_inner'] = df_result['DS_Prots'].apply(
        lambda x: set(x) != set(
            range(min(x), max(x) + 1)
        )
    )

    return df_result


def parse_prodigal_gff(path_to_gff,
                       add_id: bool = True,
                       with_nucl: bool = True) -> pd.DataFrame:
    """
    Parse Prodigal gff-file to pandas DataFrame.
    Can add gene id as just number or as number with nucleotide prefix
    """
    gff_cols_names = ('Nucleotide', 'Sourse', 'Feature', 'Start', 'End',
                      'Score', 'Strand', 'Frame', 'Comment')
    df = pd.read_csv(path_to_gff,
                     sep='\t',
                     names=gff_cols_names,
                     dtype={'Nucleotide': str, 'Sourse': str, 'Feature': str,
                            'Start': int, 'End': int, 'Score': float,
                            'Strand': str, 'Frame': str, 'Comment': str})
    if add_id:
        pattern = re.compile(r'ID=\d+_(\d+)')

        df['Protein'] = df.Comment.apply(
            lambda x: int(re.search(pattern, x).group(1))
        )

        if with_nucl:
            df = df.astype({'Protein': str})
            df.loc[:, 'Protein'] = df.Nucleotide + '_' + df.Protein

    return df
