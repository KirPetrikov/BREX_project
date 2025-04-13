"""v0.1
Creates a "Cluster - Protein - Cluster Size" table and an summary histogram.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Creates a "Cluster - Protein - Cluster Size" table '
                    'and an summary histogram.'
                                     )
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Path to the mmseqs clusters table.')
    parser.add_argument('-o', '--output', type=str, required=False,
                        help='Optional. Path to save result files. '
                             'If not specified, the parent directory of the fasta-file is used. '
                             'Will be created if it does not exist.')
    return parser.parse_args()


args = parse_arguments()
input_mmseqs_tsv = Path(args.input)
assert input_mmseqs_tsv.is_file(), 'The path to the input file is incorrect or the file does not exist.'
if args.output:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    output_clusters_tsv = output / 'clusters_table.tsv'
    output_clusters_hist = output / 'clusters_hist.png'
else:
    output_clusters_tsv = input_mmseqs_tsv.parents[0] / 'clusters_table.tsv'
    output_clusters_hist = input_mmseqs_tsv.parents[0] / 'clusters_hist.png'

df = pd.read_csv(input_mmseqs_tsv, sep='\t', names=['Cluster_p', 'Protein'])

cluster_names = {clust: 'CLUS_' + clust for clust in df.Cluster_p.unique()}

df['Cluster'] = df.Cluster_p.map(cluster_names)

df['Cluster_size'] = df.Cluster.map(df.Cluster.value_counts())

df[['Protein', 'Cluster', 'Cluster_size']].to_csv(output_clusters_tsv, sep='\t', index=False)

clss_size = df.Cluster.value_counts().value_counts().sort_index()

if len(clss_size) > 10:
    clss_size = pd.concat([clss_size.iloc[:5], clss_size[-5:]])

clss_size.plot.bar()

for i in range(len(clss_size)):
    plt.text(i - 0.3,
             clss_size.iloc[i] + 400,
             clss_size.iloc[i],
             fontsize=12)

plt.text(len(clss_size) - 3,
         clss_size.iloc[0] * 0.9,
         f'Total clusters:\n{df.Cluster.nunique()}',
         fontsize=12)

plt.tick_params(left=False, labelleft=False)
plt.xticks(size=12,
           rotation=45,
           ha='center'
           )
           
plt.ylabel('Number of clusters')
plt.xlabel('Cluster size')
plt.title('Number of five smallest and five largest clusters', loc='right')

plt.savefig(output_clusters_hist, bbox_inches='tight', dpi=100)

