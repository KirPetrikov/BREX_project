"""v0.3
Drops inner DS by name to further ignore them.
Calculates different compositions of the region, how many of what elements are in the region.

`intersect_log` - a list where internal DS are logged,
for which the list of protein IDs goes beyond the boundaries of the protein IDs of the target DS
"""

import json
import pandas as pd

from pathlib import Path

pd.options.mode.copy_on_write = True

input_accessions = Path('/home/holydiver/Main/2024_BREX/Data/Accessions_summary.tsv')
input_json_path = Path('/home/holydiver/Main/2024_BREX/Data/20250511_brex_regions_summary/'
                       'inner_summary_brex_dupl_proc.json')

output_path = Path('/home/holydiver/Main/2024_BREX/Data/20250511_brex_regions_summary')

defsys_to_drop = 'DMS'

with open(input_json_path) as f:
    data_defsys = json.load(f)

intersect_log = []

for brex_id, curr_reg in data_defsys.items():
    curr_reg['Inner_DS'] = ''
    curr_reg['Inner_DS_Prots'] = ''
    curr_reg['Others_Prots'] = ''

    if curr_reg['Have_inner']:
        all_prots_in_reg = set(
            [i for i in range(min(curr_reg['DS_Prots']), max(curr_reg['DS_Prots']) + 1)]
        )
        inner_defsys_prots = []

        if curr_reg['Inner_DS_full']:

            # Check for cases in which some proteins
            # included in the inner DS are outside the region.
            # Take the IDs of all the inner DS proteins
            # and keep only those that intersect with the set of all proteins in the region.
            for ds in curr_reg['Inner_DS_full'].copy():
                if defsys_to_drop and ds.startswith(defsys_to_drop):
                    del curr_reg['Inner_DS_full'][ds]
                    continue
                    
                curr_ds_prots = set(curr_reg['Inner_DS_full'][ds])

                # Logging intersected DS
                check_outbonds = curr_ds_prots.difference(all_prots_in_reg)
                if check_outbonds:
                    # print(f'{brex_id}\t{curr_reg["Nucleotide"]}\t{ds}')
                    intersect_log.append(brex_id)

                # Check how many IDs of the current inner DS belong to the region
                # If none - drop DS
                # Those proteins that fell into the region - save
                check_innner = all_prots_in_reg.intersection(curr_ds_prots)
                if not check_innner:
                    del curr_reg['Inner_DS_full'][ds]
                else:
                    inner_defsys_prots.extend(list(check_innner))

        if inner_defsys_prots:
            curr_reg['Inner_DS_Prots'] = inner_defsys_prots

        if curr_reg['Inner_DS_full']:
            curr_reg['Inner_DS'] = [i.split('%')[0] for i in curr_reg['Inner_DS_full']]
        else:
            curr_reg['Inner_DS_full'] = ''

        # Other proteins
        all_ds_prots = curr_reg['DS_Prots'] + inner_defsys_prots

        other_prots = list(
            all_prots_in_reg.difference(set(all_ds_prots))
        )

        if other_prots:
            curr_reg['Others_Prots'] = other_prots

    curr_reg['Others_Count'] = len(curr_reg['Others_Prots'])
    curr_reg['Inner_DS_Count'] = len(curr_reg['Inner_DS'])
    curr_reg['Prots_Count'] = curr_reg['Others_Count'] + curr_reg['Inner_DS_Count']

with open(output_path / 'inner_summary_brex_recalc.json', mode='w') as f:
    json.dump(data_defsys, f, indent=4)
