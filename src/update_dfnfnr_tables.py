"""
Add Accesions and unique defsystem IDs to 'genes' and 'systems' table.
Add Nucleotide to 'systems' table.
"""
import pandas as pd

from pathlib import Path

pd.options.mode.copy_on_write = True

input_dfnfnr_sys_path = ('/home/niagara/Storage/MetaRus/Common_dir/collections/complete_bacteria_collection/'
                         'defensfinder_out/concat_defense_finder_systems.tsv')

input_dfnfnr_genes_path = ('/home/niagara/Storage/MetaRus/Common_dir/collections/complete_bacteria_collection/'
                           'defensfinder_out/concat_defense_finder_genes.tsv')

output_genes = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/dfnfnr_genes_acc.tsv'
output_systems = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/dfnfnr_sys_acc.tsv'

input_accessions = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/Accessions_summary.tsv'
df_acc = pd.read_csv(input_accessions, sep='\t')

df_g = pd.read_csv(input_dfnfnr_genes_path, sep='\t',
                   usecols=[0, 1, 2, 22, 23],
                   dtype={'replicon': str, 'hit_id': str, 'gene_name': str, 'type': str, 'subtype': str},
                   names=['Nucleotide', 'Protein', 'DF_ann', 'DF_System', 'DF_System_sub'],
                   skiprows=1
                   )

# Show number of duplicates and drop
print(f'Genes duplicates: {df_g.duplicated().sum()}')
df_g = df_g.drop_duplicates()

# Add Accession
df_g = df_g.merge(df_acc[['Nucleotide', 'Accession']], on='Nucleotide', how='left')

# Tmp ID
df_g['tmpid'] = df_g.apply(lambda x: f'{x.Protein},{x.DF_System_sub},{x.Nucleotide}',
                           axis=1
                           )

# Systems table
df_sys = pd.read_csv(input_dfnfnr_sys_path, sep='\t',
                     usecols=[0, 1, 2, 6, 8],
                     dtype={'sys_id': str,
                               'type': str,
                               'subtype': str,
                               'protein_in_syst': str,
                               'name_of_profiles_in_sys': str},
                     names=['DS_ID', 'DF_System', 'DF_System_sub', 'All_Proteins', 'All_Annot'],
                     skiprows=1
                     )

# Show number of duplicates and drop
print(f'Systems duplicates: {df_sys.duplicated().sum()}')
df_sys = df_sys.drop_duplicates()

# Add Nucleotide
df_sys['Nucleotide'] = ''

df_sys.loc[
    (df_sys.DS_ID.str.startswith('NZ_'))
    |
    (df_sys.DS_ID.str.startswith('NC_')),
    'Nucleotide'
] = df_sys.loc[
        (df_sys.DS_ID.str.startswith('NZ_'))
        |
        (df_sys.DS_ID.str.startswith('NC_')), 'DS_ID'
].apply(lambda x: '_'.join(x.split('_')[:2]))

df_sys.loc[
    ~(
        (df_sys.DS_ID.str.startswith('NZ_'))
        |
        (df_sys.DS_ID.str.startswith('NC_'))
      ),
    'Nucleotide'
] = df_sys.loc[
        ~(
            (df_sys.DS_ID.str.startswith('NZ_'))
            |
            (df_sys.DS_ID.str.startswith('NC_'))
          ),
        'DS_ID'
].apply(lambda x: x.split('_')[0])

# Add Accession
df_sys = df_sys.merge(df_acc[['Nucleotide', 'Accession']], on='Nucleotide', how='left')

# Create unique defsys IDs as for Padloc
df_sys['number'] = df_sys.DS_ID.apply(lambda x: x.split('_')[-1])  # Добавить временный ИД - номер ЗС
df_sys['DF_DS_ID'] = df_sys.apply(lambda x: f'{x.DF_System_sub}%{x.number}%{x.Nucleotide}', axis=1)

# Assign unique defsys IDs to each protein
dftmp_data = df_sys[
    ['All_Proteins', 'DF_System_sub', 'number', 'Nucleotide', 'DF_DS_ID']
].to_dict(orient='index')

prots_from_sys = {}

for val in dftmp_data.values():
    curr_prots = val['All_Proteins'].split(',')
    for prot in curr_prots:
        tmp_id = f'{prot},{val["DF_System_sub"]},{val["Nucleotide"]}'
        prots_from_sys[tmp_id] = val['DF_DS_ID']

df_g['DF_DS_ID'] = df_g.tmpid.map(prots_from_sys)

# Save to files
df_g.drop('tmpid', axis=1).to_csv(output_genes, sep='\t', index=False)
df_sys.drop('number', axis=1).to_csv(output_systems, sep='\t', index=False)
