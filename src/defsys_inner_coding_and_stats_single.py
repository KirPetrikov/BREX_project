"""v0.1
Summary for most abundant clusters of inner other proteins.

Coding of the composition of a region: protein to corresponding cluster.
Updates DS summary by removing from 'Others_' fields proteins from annotated clusters
Find last common ancestor taxon.

Input clusters table must be with all annotations added (DefenseFinder, PDB70, PfamA)
"""

import argparse
import pandas as pd
import json

from collections import defaultdict
from pathlib import Path

pd.options.mode.copy_on_write = True


def parse_args():
    parser = argparse.ArgumentParser(description='Clusters coding of other genes'
                                                 ' and description summary')
    parser.add_argument('-c', '--clusters_table', type=Path, required=True,
                        help='Clusters table with all annotations added, .tsv')
    parser.add_argument('-s', '--defsys_summary', type=Path, required=True,
                        help='Defense systems summary, .json')
    parser.add_argument('-a', '--accessions', type=Path, required=True,
                        help='Table Nucleotide - Assembly accession - Taxon, .tsv')
    parser.add_argument('-b', '--pdb70_parent', type=Path, required=True,
                        help='PDB70 annotations directory')
    parser.add_argument('-f', '--pfama_parent', type=Path, required=True,
                        help='Pfam-A annotations directory')
    parser.add_argument('-p', '--padloc_data', type=Path, required=True,
                        help='Padloc results directory')
    parser.add_argument('-t', '--cluster_size_threshold', type=int, default=2,
                        help='Minimal size of clusters to be considered (default: 2)')
    parser.add_argument('-m', '--top_clusters_number', type=int, default=10,
                        help='Number of most abundant clusters to select for descriptoin (default: 10)')
    parser.add_argument('-o', '--output_folder', type=Path, required=True,
                        help='Directory to save result. Will be created')

    return parser.parse_args()


def clusters_coding(data, clusters, threshold=2):
    """
    Update defense systems summary by
    adding field 'Other_Cluster' with list of corresponding clusters.
    Only for Other proteins.

    :param data: Defense systems summary
    :param clusters: Clusters dataframe
    :param threshold: Minimal size of clusters to be considered (default: 2)
    :return: DataFrame: Cluster | List of DS_IDs with it | Number of regions
    """

    codes = dict(
        zip(clusters.loc[clusters.Cluster_size >= threshold].Protein,
            clusters.loc[clusters.Cluster_size >= threshold].Cluster)
    )

    regs_composition = defaultdict(list)

    for ds_id, val in data.items():
        data[ds_id]['Other_Cluster'] = ''
        curr_other_prots = [f'{val["Nucleotide"]}_{i}' for i in val['Others_Prots']]
        curr_clusters = []
        for prot in curr_other_prots:
            if prot in codes:
                curr_clusters.append(codes[prot])
            else:
                pass
        if curr_clusters:
            data[ds_id]['Other_Cluster'] = curr_clusters

            for cluster in data[ds_id]['Other_Cluster']:
                regs_composition[cluster].append(ds_id)

    df = pd.DataFrame(
        pd.Series(data=regs_composition,
                  index=[i for i in regs_composition],
                  name='DS_IDs').rename_axis('Cluster')
    )

    df['Regs_Count'] = df.DS_IDs.apply(len)

    return df


def get_accession(nucl, df_acc):
    accession = df_acc.loc[df_acc.Nucleotide == nucl].iat[0, 0]
    return accession


def read_fasta(fasta_file_path) -> tuple[str, str]:
    """
    Read fasta-file
    Returns (>fasta_header, sequence) as in source
    (keep one/multi-lines format)
    """
    name, seq = None, []
    for line in fasta_file_path:
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


def extract_sequence(fasta_file_path,
                     prot_ids: list | tuple) -> list:
    """
    Extract sequences with specific IDs/headers from fasta-file
    """

    count, stop_num, selected_seqs = 0, len(prot_ids), []

    with open(fasta_file_path, mode='r') as file:
        for name, seq in read_fasta(file):
            header = name.strip()[1:]
            if header in prot_ids:
                seq = seq.replace("\n", "")
                selected_seqs.append(f'{header}:{seq}')

                count += 1
                if count == stop_num:
                    break

    return selected_seqs


args = parse_args()

input_clusters_table_path = args.clusters_table
input_defsys_summary_path = args.defsys_summary
input_accessions_path = args.accessions
input_pdb70_parent_path = args.pdb70_parent
input_pfama_parent_path = args.pfama_parent
input_padloc_data_path = args.padloc_data
output_folder = args.output_folder
cluster_size_threshold = args.cluster_size_threshold
top_clusters_number = args.top_clusters_number

output_folder.mkdir(parents=True, exist_ok=True)

pdb70_path = input_pdb70_parent_path / 'region_pdb_top_annot.tsv'
pfama_path = input_pfama_parent_path / 'region_pfama_top_annot.tsv'

# Read data
with open(input_defsys_summary_path) as f:
    data_defsys = json.load(f)

df_clusters = pd.read_csv(input_clusters_table_path, sep='\t')
df_accessions = pd.read_csv(input_accessions_path, sep='\t')

df_pdb70 = pd.read_csv(pdb70_path, sep='\t',
                       dtype={'Cluster': str, 'PDB_ID': str, 'PDB_ann': str, 'EVal': float, 'Prob': float})

df_pfama = pd.read_csv(pfama_path, sep='\t',
                       dtype={'Cluster': str, 'Pfam_ID': str, 'Pfam_ann': str, 'EVal': float, 'Prob': float})

# Select clusters w/o any defsys annotation
df_to_sel_clust_no_ann = df_clusters.groupby('Cluster')[['DF_ann', 'Padloc_ann']].agg(set)

clust_no_ann = df_to_sel_clust_no_ann.loc[
    (df_to_sel_clust_no_ann.DF_ann == {'miss'})
    &
    (df_to_sel_clust_no_ann.Padloc_ann == {'miss'})
].index

# Proteins from clusters w/o defsys annotation
prots_without_defsys_ann = df_clusters.loc[df_clusters.Cluster.isin(clust_no_ann)].Protein.tolist()

# Update 'Others_' fields
for curr_reg in data_defsys.values():
    new_other_prots = []

    for prot_id in curr_reg['Others_Prots']:
        curr_prot = f'{curr_reg["Nucleotide"]}_{prot_id}'
        if curr_prot in prots_without_defsys_ann:
            new_other_prots.append(prot_id)

    curr_reg['Others_Prots'] = new_other_prots
    curr_reg['Others_Count'] = len(new_other_prots)

# Clusters coding
df_prots_regs = clusters_coding(data_defsys, df_clusters, threshold=cluster_size_threshold)

selected_clusters = df_prots_regs.nlargest(top_clusters_number, 'Regs_Count').index

top_clusters_data = {}

for curr_cluster in selected_clusters:
    print(f'---Process {curr_cluster}---')

    ann_1_pdb = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].PDB_ann.iat[0]
    id_1_pdb = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].PDB_ID.iat[0]
    eval_1_pdb = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].EVal.iat[0]
    ann_1_pfama = df_pfama.loc[df_pfama.Cluster == curr_cluster].Pfam_ann.iat[0]
    id_1_pfama = df_pfama.loc[df_pfama.Cluster == curr_cluster].Pfam_ID.iat[0]
    eval_1_pfama = df_pfama.loc[df_pfama.Cluster == curr_cluster].EVal.iat[0]

    # Representatives sequences
    fasta_reprs_acc = get_accession(df_clusters.loc[df_clusters.Cluster == curr_cluster].Nucleotide.iat[0],
                                    df_accessions)
    fasta_reprs_id = df_clusters.loc[df_clusters.Cluster == curr_cluster].Protein.iat[0]

    path_to_fasta_reprs = input_padloc_data_path / f'{fasta_reprs_acc}/{fasta_reprs_acc}_prodigal.faa'

    fasta_reprs = extract_sequence(path_to_fasta_reprs, (fasta_reprs_id,))

    top_clusters_data[curr_cluster] = {}
    top_clusters_data[curr_cluster]['Size'] = int(df_prots_regs.loc[curr_cluster, 'Regs_Count'])
    top_clusters_data[curr_cluster]['PDB_ID'] = id_1_pdb
    top_clusters_data[curr_cluster]['PDB_Descr'] = ann_1_pdb
    top_clusters_data[curr_cluster]['PDB_EVal'] = eval_1_pdb
    top_clusters_data[curr_cluster]['Pfam_ID'] = id_1_pfama
    top_clusters_data[curr_cluster]['Pfam_Descr'] = ann_1_pfama
    top_clusters_data[curr_cluster]['Pfam_EVal'] = eval_1_pfama
    top_clusters_data[curr_cluster]['Repr'] = fasta_reprs[0]

# Save data
(pd.DataFrame.from_dict(top_clusters_data, orient='index')
 .rename_axis('Cluster')
 .to_csv(output_folder / 'top_clusters_data.tsv', sep='\t'))

with open(output_folder / 'top_clusters_data.json', mode='w') as f:
    json.dump(top_clusters_data, f, indent=4)
