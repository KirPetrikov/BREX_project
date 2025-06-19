"""
Parse DefenseFinder table with complete defsys to get protein annotations
"""
import json
import pandas as pd

pd.options.mode.copy_on_write = True

input_dfnfnr_sys_path = '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/dfnfnr_sys_acc.tsv'
df_dr_sys = pd.read_csv(input_dfnfnr_sys_path, sep='\t',
                        dtype={'DS_ID': str, 'DF_System': str, 'DF_System_sub': str, 'All_Proteins': str,
                               'All_Annot': str, 'Nucleotide': str, 'Accession': str, 'My_ID': str},
                        )

dftmp_data = df_dr_sys.to_dict(orient='index')

prots_from_sys = {'Protein': [], 'DF_ann': [], 'DF_System': [], 'DF_System_sub': [], 'Nucleotide': [], 'Accession': []}

for val in dftmp_data.values():
    curr_prots = val['All_Proteins'].split(',')
    curr_anns = val['All_Annot'].split(',')
    for c_prot, c_ann in zip(curr_prots, curr_anns):
        prots_from_sys['Protein'].append(c_prot)
        prots_from_sys['DF_ann'].append(c_ann)
        prots_from_sys['DF_System'].append(val['DF_System'])
        prots_from_sys['DF_System_sub'].append(val['DF_System_sub'])
        prots_from_sys['Nucleotide'].append(val['Nucleotide'])
        prots_from_sys['Accession'].append(val['Accession'])
        prots_from_sys['DS_ID'].append(val['DS_ID'])
    
pd.DataFrame(prots_from_sys).to_csv('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/dfnfnr_prots_from_sys.tsv', sep='\t', index=False)

