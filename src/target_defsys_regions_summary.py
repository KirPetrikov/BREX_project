"""v0.2b
Parse defense systems summary and gff-file to create summary for target DS region inner genes
and get ID for up/down genes.
Genes that intersect with the outermost genes of the target system are considered external (up- or downstream).
Logs cases when among inner DS occur intersections and duplicates.
Can consider "DMS_other" genes as genes that are not related to defense systems.
"""
import argparse
import json
import pandas as pd

from pathlib import Path
from defsys import parse_gff, add_gene_id_to_gff

pd.options.mode.copy_on_write = True


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Parse defense systems summary and gff-file to create for target DS region '
                    'inner genes summary and get up/down genes. Logs cases when among inner DS occur duplicates. '
                    'Can consider "DMS_other" genes as genes that are not related to defense systems.'
                                     )
    parser.add_argument('-p', '--padloc', type=str,
                        help='Path to the Padloc data (for gff-file).')
    parser.add_argument('-s', '--summary', type=str,
                        help='Path to the DS summary.')
    parser.add_argument('-d', '--duplicates', type=str,
                        help='Path to the DS duplicates table.')
    parser.add_argument('-o', '--out', type=str,
                        help='Path to the results folder. Will be created if it does not exist.')
    parser.add_argument('-t', '--target', type=str,
                        help='Name or prefix of target defense system(s).')
    parser.add_argument('-D', '--drop_dms', action='store_true',
                        help='Set up flag if needed to consider DMS_other genes as non DS.')
    return parser.parse_args()


args = parse_arguments()
input_padloc_path = Path(args.padloc)
input_summary_path = Path(args.summary)
input_dupl_defsys_path = Path(args.duplicates)
output_path_results = Path(args.out)
output_path_results.mkdir(parents=True, exist_ok=True)
target_defsys = args.target
drop_dms = args.drop_dms

with open(input_summary_path) as f:
    data_summary = json.load(f)

df_summary = pd.read_json(input_summary_path,
                          orient='index',
                          dtype={'Nucleotide': str, 'System': str, 'Start': int, 'End': int,
                                 'Strand': str, 'DS_Prots': list, 'Have_inner': bool})

duplicates = pd.read_csv(input_dupl_defsys_path, sep='\t').DS_ID.values

boundary_proteins = ['DS_ID\tAccession\tUp_Prot\tUp_Int\tDw_Prot\tDw_Int\n']
data_target = {}
inner_defsys_duplicates = ['DS_ID\tAccession\n']

for ds_id, values in data_summary.items():
    gff_path = input_padloc_path / f'{values["Accession"]}/{values["Accession"]}_prodigal.gff'
    if gff_path.exists() and ds_id.startswith(target_defsys):
        print(f'---Processing {ds_id}---')

        data_target[ds_id] = values
        data_target[ds_id]['Inner_DS'] = ''
        data_target[ds_id]['Inner_DS_Prots'] = ''

        curr_gff = parse_gff(gff_path)
        curr_gff = add_gene_id_to_gff(curr_gff, add_nucl=False)

        # Пересечения с границами региона целевой ЗС считаем внешними генами
        # Делаем проверку на такие гены
        left_intersect = (curr_gff.loc[
            (curr_gff.Chrom == values['Nucleotide'])
            &
            ((curr_gff.Start <= values['Start']) & (curr_gff.End >= values['Start']))
            &
            (~curr_gff.Gene_ID.isin(values['DS_Prots']))
                                        ]).Gene_ID

        right_intersect = (curr_gff.loc[
            (curr_gff.Chrom == values['Nucleotide'])
            &
            ((curr_gff.Start <= values['End']) & (curr_gff.End >= values['End']))
            &
            (~curr_gff.Gene_ID.isin(values['DS_Prots']))
                                        ]).Gene_ID

        max_id = curr_gff.loc[curr_gff.Chrom == values['Nucleotide']].Gene_ID.max()

        # В зависимости от направления целевой ЗС назначаем ап/даунстрим белки
        # intersect_upstr, intersect_dwstr: метки о наличии пересечения
        # Проверяем ID на выход за границы всего нуклеотида (левая: 1, правая: наибольший ИД)
        if values['Strand'] == '+':
            if not left_intersect.empty:
                prot_upstr = left_intersect.iat[0]
                intersect_upstr = True
            else:
                prot_upstr = min(values['DS_Prots']) - 1
                intersect_upstr = False

            if not right_intersect.empty:
                prot_downstr = right_intersect.iat[0]
                intersect_dwstr = True
            else:
                prot_downstr = max(values['DS_Prots']) + 1
                # Если правый белок - граница региона, то даунстрим будет 0, т.е. несуществующим белком
                if prot_downstr > max_id:
                    prot_downstr = 0
                intersect_dwstr = False

        elif values['Strand'] == '-':
            if not left_intersect.empty:
                prot_downstr = left_intersect.iat[0]
                intersect_dwstr = True
            else:
                prot_downstr = min(values['DS_Prots']) - 1
                intersect_dwstr = False

            if not right_intersect.empty:
                prot_upstr = right_intersect.iat[0]
                intersect_upstr = True
            else:
                prot_upstr = max(values['DS_Prots']) + 1
                # Если правый белок - граница региона, то апстрим будет 0, т.е. несуществующим белком
                if prot_upstr > max_id:
                    prot_upstr = 0
                intersect_upstr = False

        else:
            raise ValueError(f'In {ds_id} "Strand" value is not + or -')

        p_u = values['Nucleotide'] + '_' + str(prot_upstr)
        p_d = values['Nucleotide'] + '_' + str(prot_downstr)
        boundary_proteins.append(
            f'{ds_id}\t{values["Accession"]}\t{p_u}\t{intersect_upstr}\t{p_d}\t{intersect_dwstr}\n'
                                 )

        # Смотрим внутренние ЗС
        curr_target_defsys = df_summary.loc[[ds_id]]

        if values['Have_inner']:
            extracted_inner_defsys = df_summary.loc[
                # Проверяем, что целевой нуклеотид
                (df_summary.Nucleotide == values['Nucleotide'])
                &
                # Смотрим полностью вложенные системы
                ((values['Start'] <= df_summary.Start) & (df_summary.End <= values['End']))
                &
                # Саму целевую ЗС не включаем
                (~(df_summary.index == ds_id))
                ]

            if drop_dms:
                extracted_inner_defsys = extracted_inner_defsys.loc[
                                           ~(extracted_inner_defsys.System == 'DMS_other')
                                                                    ]
            if not extracted_inner_defsys.empty:
                data_target[ds_id]['Inner_DS'] = extracted_inner_defsys.System.tolist()
                data_target[ds_id]['Inner_DS_Prots'] = [x for xs in extracted_inner_defsys.DS_Prots for x in xs]

                # Проверка наличия ЗС с дупликацией
                curr_dupl_defsys = extracted_inner_defsys.index.intersection(duplicates)
                if curr_dupl_defsys.size != 0:
                    l = [f'{i}\t{values["Accession"]}\n' for i in curr_dupl_defsys]
                    inner_defsys_duplicates.extend(l)

# Write results
with open(output_path_results / f'inner_summary_{target_defsys}.json', 'w') as f:
    json.dump(data_target, f, indent=4)

with open(output_path_results / f'boundary_proteins_{target_defsys}.tsv', mode='w') as f:
    for i in boundary_proteins:
        f.write(i)

with open(output_path_results / f'inner_summary_{target_defsys}_dupl.tsv', mode='w') as f:
    for i in inner_defsys_duplicates:
        f.write(i)

