"""v0.1

Обработка дуплицирующихся БРЕКС-систем

Входные данные - на основе ручной проверки таблицы дупликатов.

**Сложные случаи**
Просто дропать нуклеотиды целиком, т.к. мы не знаем, где и как могут пересечься ЗС с двойным учётом.
Всего таких: 30

**Простые случаи**
Всего: 104
Где для случаев аннотации БРЕКС/RM есть единственный дуплицирующийся PglX/REase_MTase_IIG.
PglX - обязательный во всех БРЕКСсистемах, так что там, где он единственный,
его нельзя засчитать, как систему "RM_type_IIG".
Сл-но, просто дропаем эту альтернативную аннотацию, оставляя PglX для БРЕКС-сситем.
В аннотациях белков - дропаем соответствующие системы "RM_type_IIG", т.е. этот белок для них единственный.
"""

import json
import pandas as pd

from pathlib import Path

output_path_results = Path('/home/holydiver/Main/2024_BREX/Data/20250324_dupl_brex_processing')

duplicates_processing_json = Path('dupl_brex_processing_data.json')

input_summary_path = Path('/home/holydiver/Main/2024_BREX/Data/20250323_defsys_summary/defsys_summary_all.json')
input_prot_annot_path = Path('/home/holydiver/Main/2024_BREX/Data/20250323_defsys_summary/protein_annotations_all.tsv')

with open(duplicates_processing_json) as f:
    duplicates_processing_data = json.load(f)

df_sum = pd.read_json(input_summary_path, orient='index')

# Add accession
input_accessions = Path('/home/holydiver/Main/2024_BREX/Data/Accessions_summary.tsv')
df_acc = pd.read_csv(input_accessions, sep='\t')
df_sum = (df_sum.reset_index()
                .rename({'index': 'DS_ID'}, axis=1)
                .merge(df_acc[['Nucleotide', 'Accession']], on='Nucleotide', how='left')
                .set_index('DS_ID'))

df_ann = pd.read_csv(input_prot_annot_path, sep='\t')

# Easy cases
df_sum = df_sum.drop(duplicates_processing_data['dsid_to_drop_rm'])
df_ann = df_ann.loc[~(
    (df_ann.Padloc_ann == 'REase_MTase_IIG')
    &
    (df_ann.Protein.isin(duplicates_processing_data['proteins_to_drop_rm']))
                      )]

# CP001337.1
df_sum.at['RM_type_II%12%CP001337.1', 'DS_Prots'] = [1153, 1154]
idx_to_drop = df_ann.loc[(
    (df_ann.Padloc_ann == 'MTase_II')
    &
    (df_ann.Protein == 'CP001337.1_1151')
                          )].index
df_ann = df_ann.drop(idx_to_drop)

# Hard cases
df_sum = df_sum.loc[~df_sum.Nucleotide.isin(duplicates_processing_data['nucleotides_to_drop_hard'])]
df_ann = df_ann.loc[~df_ann.Nucleotide.isin(duplicates_processing_data['nucleotides_to_drop_hard'])]

# Write full results
df_sum.to_json(output_path_results / 'defsys_summary_dupl_brex_proc.json', orient='index')
df_ann.to_csv(output_path_results / 'protein_annotations_dupl_brex_proc.json', sep='\t', index=False)

# Select only nucleotides with target DS
target_defsys = 'brex'
nucl_with_target_defsys = df_sum.loc[df_sum.System.str.startswith(target_defsys)].Nucleotide.values

df_sum = df_sum.loc[df_sum.Nucleotide.isin(nucl_with_target_defsys)]
df_ann = df_ann.loc[df_ann.Nucleotide.isin(nucl_with_target_defsys)]

# Write selected results
df_sum.to_json(output_path_results / 'defsys_summary_dupl_brex_proc_select.json', orient='index')
df_ann.to_csv(output_path_results / 'protein_annotations_dupl_brex_proc_select.tsv', sep='\t', index=False)
