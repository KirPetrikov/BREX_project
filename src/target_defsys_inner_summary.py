"""v0.1
"""

import pandas as pd

from pathlib import Path

pd.options.mode.copy_on_write = True

input_summary_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250324_dupl_brex_processing/'
                          'defsys_summary_dupl_brex_proc.json')
input_dupl_defsys_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250323_defsys_summary/'
                              'duplicated_defsys.tsv')
output_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250331_brex_inner_summary')

target_defsys = 'brex'

df_sum = pd.read_json(input_summary_path,
                      orient='index',
                      dtype={'Nucleotide': str, 'System': str, 'Start': int, 'End': int,
                             'Strand': str, 'DS_Prots': list, 'Have_inner': bool})
nucl_with_target_defsys = df_sum.loc[df_sum.System.str.startswith(target_defsys)].Nucleotide.values
df_target_nucl = df_sum.loc[df_sum.Nucleotide.isin(nucl_with_target_defsys)]

df_target_nucl['Inner_DS'] = ''
df_target_nucl['Inner_DS_Names'] = ''
df_target_nucl['Inner_DS_Prots'] = ''
df_target_nucl['All_DS_Prots'] = ''
df_target_nucl['Non_DS'] = ''

duplicates = pd.read_csv(input_dupl_defsys_path, sep='\t')
duplicates = duplicates.loc[duplicates.Nucleotide.isin(nucl_with_target_defsys)]
defsys_dupl = duplicates.loc[duplicates.Nucleotide.isin(nucl_with_target_defsys)].DS_ID.values

target_defsys_for_summary = df_target_nucl.loc[
    df_target_nucl.System.str.startswith(target_defsys), ['Nucleotide', 'Start', 'End', 'Have_inner']
                                               ].to_dict(orient='index')

target_defsys_collect = []
inner_defsys_prots_outbonds = []
inner_defsys_coords_outbonds = []
inner_defsys_duplicates = []

# Для каждой целевой ЗС собираем саммари отдельно
for ds_id, values in target_defsys_for_summary.items():
    curr_target_defsys_df = df_target_nucl.loc[[ds_id]]

    # Сразу создаём список белков всех ЗС, пока что - только целевой ЗС
    curr_all_defsys_prots = curr_target_defsys_df.DS_Prots.iat[0]

    if values['Have_inner']:
        extracted_inner_defsys = df_target_nucl.loc[
            # Проверяем, что целевой нуклеотид
            (df_target_nucl.Nucleotide == values['Nucleotide'])
            &
            # Смотрим пересечения
            (
                ((values['Start'] <= df_target_nucl.Start) & (df_target_nucl.Start <= values['End']))
                |
                ((values['Start'] <= df_target_nucl.End) & (df_target_nucl.End <= values['End']))
             )
            &
            # Саму целевую ЗС не включаем
            (~df_target_nucl.System.str.startswith(target_defsys))
            ]
        if not extracted_inner_defsys.empty:
            curr_inner_defsys = {}
            curr_inner_prots = []
            curr_inner_names = []
            # Добавляем список внутренних ЗС
            for i, frame in extracted_inner_defsys[['System', 'DS_Prots']].iterrows():
                curr_inner_defsys[frame.System] = frame.DS_Prots
                curr_inner_prots.extend(frame['DS_Prots'])
                curr_inner_names.append(frame.System)
            curr_target_defsys_df.at[ds_id, 'Inner_DS'] = curr_inner_defsys
            curr_target_defsys_df.at[ds_id, 'Inner_DS_Names'] = curr_inner_names
            curr_target_defsys_df.at[ds_id, 'Inner_DS_Prots'] = curr_inner_prots

            # Обновляем список белков всех ЗС
            curr_all_defsys_prots = curr_target_defsys_df.DS_Prots.iat[0] + curr_inner_prots

            # Проверка выхода ИД белков внутренних ЗС за границы
            max_target_prot = max(curr_target_defsys_df.DS_Prots.iloc[0])
            min_target_prot = min(curr_target_defsys_df.DS_Prots.iloc[0])
            inner_prots = [i for ii in extracted_inner_defsys.DS_Prots.to_list() for i in ii]
            max_inner_prot = max(inner_prots)
            min_inner_prot = min(inner_prots)
            cond1 = (min_inner_prot < min_target_prot) or (max_inner_prot > max_target_prot)
            if cond1:
                inner_defsys_prots_outbonds.append(ds_id)

            # Проверка выхода координат внутренних ЗС за границы
            max_inner_coord = extracted_inner_defsys.loc[:, ['Start', 'End']].max().max()
            min_inner_coord = extracted_inner_defsys.loc[:, ['Start', 'End']].min().min()
            cond2 = ((min_inner_coord < curr_target_defsys_df.Start.iloc[0])
                     or
                     (max_inner_coord > curr_target_defsys_df.End.iloc[0]))
            if cond2:
                if ds_id not in set(inner_defsys_prots_outbonds):
                    inner_defsys_coords_outbonds.append(ds_id)

            # Проверка наличия ЗС с дупликацией
            curr_dupl_defsys = extracted_inner_defsys.index.intersection(defsys_dupl)
            cond3 = curr_dupl_defsys.size != 0
            if cond3:
                inner_defsys_duplicates.extend(curr_dupl_defsys.tolist())

        # Список белков не-ЗС
        curr_non_defsys_prots = set(
            i for i in range(
                min(curr_target_defsys_df.DS_Prots.iat[0]),
                max(curr_target_defsys_df.DS_Prots.iat[0]) + 1
                             )
                                    ) - set(curr_all_defsys_prots)
        if curr_non_defsys_prots:
            curr_target_defsys_df.at[ds_id, 'Non_DS'] = list(curr_non_defsys_prots)

    # Добавляем список белков всех ЗС
    curr_target_defsys_df.at[ds_id, 'All_DS_Prots'] = curr_all_defsys_prots

    target_defsys_collect.append(curr_target_defsys_df)

df_target = pd.concat(target_defsys_collect).astype({'Have_inner': bool})

df_target.to_json(output_path / f'target_{target_defsys}_summary.json', orient='index')

logs = (
    inner_defsys_prots_outbonds,
    inner_defsys_coords_outbonds,
    inner_defsys_duplicates
)

logs_names = (
    f'target_{target_defsys}_inner_prots_outbonds.txt',
    f'target_{target_defsys}_inner_coords_outbonds.txt',
    f'target_{target_defsys}inner_duplicates.txt'
)

for data, name in zip(logs, logs_names):
    with open(output_path / name, mode='w') as f:
        for line in data:
            f.write(line + '\n')
