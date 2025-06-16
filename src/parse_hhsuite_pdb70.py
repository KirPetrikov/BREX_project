"""
v0.1c
"""
import json
import re
import sys
import pandas as pd

from pathlib import Path


input_hh_suite_results = Path(sys.argv[1])
input_n_hits = int(sys.argv[2])
output_path = Path(sys.argv[3])
prefix_for_save = input_hh_suite_results.name.lower()

output_path.mkdir(parents=True, exist_ok=True)

output_tsv_name = Path(f'{prefix_for_save}_pdb_top_annot.tsv')
output_json_name = Path(f'{prefix_for_save}_pdb_annot.json')

descriptions_catalog = {}  # Collections of all top annotations with descriptions
only_first_annotations = {'Cluster': [],  # Top-1 annotations for table
                          'PDB_ID': [],
                          'PDB_ann': [],
                          'EVal': [],
                          'Prob': []}

pattern_table = re.compile(r'\d+ +(\w+) ')
pattern_descr = re.compile(r'>(\w+) +(.+);.+.$')

for file_name in input_hh_suite_results.iterdir():
    if file_name.suffix == '.hhr':
        with open(file_name) as file:
            cluster_id = file_name.stem
            print(f'Parse {cluster_id}')
            annotations = []  # List of top-hits annotations
            descriptions_catalog[cluster_id] = {}

            # Skip headers line
            count = 0
            while count < 9:
                file.readline()
                count += 1

            curr_line = file.readline()

            # Only top hit
            only_first_annotations['Cluster'].append(cluster_id)
            only_first_annotations['PDB_ID'].append(re.search(pattern_table, curr_line[:35]).group(1))
            only_first_annotations['EVal'].append(float(curr_line[41:49].strip()))
            only_first_annotations['Prob'].append(float(curr_line[35:41].strip()))

            # Top N hits with additional description
            count_n_hits = 0
            for line in file:
                if count_n_hits == input_n_hits:
                    break
                elif line.startswith('>'):
                    count_n_hits += 1

                    curr_search = re.search(pattern_descr, line)
                    curr_pdb_id = curr_search.group(1)
                    curr_descr = curr_search.group(2)

                    if count_n_hits == 1:
                        only_first_annotations['PDB_ann'].append(
                            curr_search.group(2).split(';')[0]
                        )

                    stat_line = file.readline()

                    curr_description = {
                        f'Hit_No_{count_n_hits}': {
                            'PDB_ID': curr_pdb_id,
                            'Description': curr_descr,
                            'EVal': float(stat_line.split('  ')[1][8:]),
                            'Prob': float(stat_line.split('  ')[0][7:])
                        }
                    }
                    descriptions_catalog[cluster_id].update(curr_description)

                else:
                    continue

(pd.DataFrame.from_dict(only_first_annotations, orient='columns')
             .to_csv(output_path / output_tsv_name, sep='\t', index=False))

with open(output_path / output_json_name, mode='w') as file:
    json.dump(descriptions_catalog, file)
