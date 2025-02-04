"""
v0.1
"""
import json
import sys
import pandas as pd

from pathlib import Path

input_hh_suite_results = Path(sys.argv[1])
input_n_hits = int(sys.argv[2])
output_path = Path(sys.argv[3])

output_tsv_name = Path('pdb_top_annot.tsv')
output_json_name = Path('pdb_annot.json')

descriptions_catalog = {}  # Collections of all top annotations with descriptions
only_first_annotations = {'Cluster': [],  # Top-1 annotations for table
                          'pdb70_id': [],
                          'Annotation': []}

for file_name in input_hh_suite_results.iterdir():
    if file_name.suffix == '.hhr':
        with open(file_name) as file:
            cluster_id = file_name.stem
            top_pdb_ids = []  # List of read top-hits ids
            annotations = []  # List of top-hits annotations

            # Skip headers line
            count = 0
            while count < 9:
                file.readline()
                count += 1

            curr_line = file.readline()

            # Only top hit
            curr_first_pdb_id = curr_line[4:11].strip(' ')
            curr_first_annot = curr_line[11:36].strip(' ')
            only_first_annotations['Cluster'].append(cluster_id)
            only_first_annotations['pdb70_id'].append(curr_first_pdb_id)
            only_first_annotations['Annotation'].append(curr_first_annot)

            # Top N hits with additional description
            count_n_hits = 0
            for line in file:
                if count_n_hits == input_n_hits:
                    break
                elif line.startswith('>'):
                    count_n_hits += 1
                    curr_pdb_id = line[1:8].strip()
                    curr_descr = '; '.join(line[8:].strip().split(';')[:-2])
                    curr_e_val = file.readline().strip().split('  ')[1][8:]
                    annotations.append((curr_pdb_id, curr_descr, curr_e_val))
                else:
                    continue

            all_top_annotations = {
                f"Hit_No_{num + 1}": {'ID': ann[0], 'Description': ann[1], "E-val": ann[2]} for
                num, ann in enumerate(annotations)}

            descriptions_catalog[cluster_id] = all_top_annotations

(pd.DataFrame.from_dict(only_first_annotations, orient='columns')
   .to_csv(output_path / output_tsv_name, sep='\t', index=False))

with open(output_path / output_json_name, mode='w') as file:
    json.dump(descriptions_catalog, file)
