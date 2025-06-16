"""v 0.1
Add to table with Padloc annotation the ones by DefenseFinder"""
import pandas as pd

input_dfnfnr_genes_path = ('/home/niagara/Storage/MetaRus/Common_dir/collections/complete_bacteria_collection/'
                           'defensfinder_out/concat_defense_finder_genes.tsv')
df_ann_dfnfnr = pd.read_csv(input_dfnfnr_genes_path, sep='\t',
                            usecols=[0, 1, 2, 22, 23],
                            dtype={'replicon': str, 'hit_id': str, 'gene_name': str, 'type': str, 'subtype': str},
                            names=['Nucleotide', 'Protein', 'DF_ann', 'DF_System', 'DF_System_sub'],
                            skiprows=1
                            )

input_annotpadloc_path = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/curr_protein_annotations_prev.tsv'
df_ann_padloc = (
    pd.read_csv(
        input_annotpadloc_path,
        sep='\t',
        dtype={'Protein': str, 'Padloc_ann': str,
               'System': str, 'Nucleotide': str,
               'Start': int, 'End': int,
               'Strand': str, 'Localisation': str}
    )
    [['Protein', 'Padloc_ann', 'System', 'Nucleotide']]
    .rename({'System': 'P_System'}, axis=1)
)

df_combine = pd.merge(
    df_ann_dfnfnr, df_ann_padloc, on=['Protein', 'Nucleotide'], how='outer'
).fillna('miss')

# Check all nucl
a = set(df_ann_padloc.Nucleotide.tolist())
b = set(df_ann_dfnfnr.Nucleotide.tolist())
c = a.union(b)
d = set(df_combine.Nucleotide.tolist())

assert not c.difference(d), 'Some nucleotides missed in combined dataframe'

# Check all prots
a = set(df_ann_padloc.Protein.tolist())
b = set(df_ann_dfnfnr.Protein.tolist())
c = a.union(b)
d = set(df_combine.Protein.tolist())

assert not c.difference(d), 'Some proteins missed in combined dataframe'

df_combine.to_csv('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/curr_protein_annotations.tsv',
                  sep='\t',
                  index=False)

# Check if there are redundant protein annotations
print(df_ann_dfnfnr.groupby('Protein').count().Nucleotide.value_counts())
