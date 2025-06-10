#!/usr/bin/env python3
"""v0.1
Add last common ancestor to summary. Rewrite source files!
"""
import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Add last common ancestor to summary.'
                                                 'Rewrite source files!')
    parser.add_argument('-j', '--json_data', type=Path, required=True,
                        help='Description summary, .json')
    parser.add_argument('-t', '--table_data', type=Path, required=True,
                        help='Description summary table, .tsv')
    return parser.parse_args()


args = parse_args()

with open(args.json_data) as f:
    data_clusters = json.load(f)

results = {}

pipeline = (
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit name2taxid | cut -f 2 | paste -sd ',' | "
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit lca -s ',' | "
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit lineage -Ln -i 2 | cut -f 3"
)

for num, (key, val) in enumerate(data_clusters.items()):
    print(f'---Process {key}---')
    tax_names = '\n'.join(val['Tax'])

    proc = subprocess.Popen(
        pipeline,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    stdout, stderr = proc.communicate(input=tax_names)

    if proc.returncode != 0:
        print(f'Error for cluster {key}: {stderr}')
    else:
        results[key] = stdout.strip()

    data_clusters[key]['Tax_LCA'] = stdout.strip()

pd.DataFrame.from_dict(data_clusters, orient='index').to_csv(args.table_data, sep='\t')

with open(args.json_data, mode='w') as f:
    json.dump(data_clusters, f, indent=4)
