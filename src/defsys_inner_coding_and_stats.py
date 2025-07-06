"""v0.3
Summary for inner other proteins:
- for pairs of most co-represented clusters
- for most abundant clusters

Coding of the composition of a region: protein to corresponding cluster.
Updates DS summary by removing from 'Others_' fields proteins from annotated clusters
Find last common ancestor taxon.

Input clusters table must be with all annotations added (DefenseFinder, PDB70, PfamA)
"""

import argparse
import subprocess
import pandas as pd
import json

from collections import defaultdict
from pathlib import Path
from defsys import co_occurence_matrix, top_co_occurred

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
    parser.add_argument('-m', '--top_coocc_number', type=int, default=10,
                        help='Number of top cooccured pairs to select for descriptoin (default: 10)')
    parser.add_argument('-n', '--top_clusters_number', type=int, default=10,
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
        other_prots = [f'{val["Nucleotide"]}_{i}' for i in val['Others_Prots']]
        curr_clusters = []
        for prot in other_prots:
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


def jaccard_index(id_1, id_2, df) -> tuple[float, set]:
    """
    Calculate Jaccard index, find set of common regions.

    :param id_1: Cluster_1 ID
    :param id_2: Cluster_2 ID
    :param df: Result of 'clusters_coding()'
    :return: Jaccard index, set of common regions (DS_IDs)
    """
    set_1 = set(df.loc[id_1].to_list()[0])
    set_2 = set(df.loc[id_2].to_list()[0])

    intrs = len(set_1.intersection(set_2))
    un = len(set_1.union(set_2))

    return round(intrs / un, 3), set_1.intersection(set_2)


def get_accession(nucl, df_acc):
    accession = df_acc.loc[df_acc.Nucleotide == nucl].iat[0, 0]
    return accession


def get_dist_between_proteins(frame):
    """
    Returns number of genes between two other genes
    based on sequential IDs difference

    For use in .apply() on pd.Series:
        |    index   |       'Protein'      |
        |------------|----------------------|
        | nucleotide | [Prot1_ID, Prot2_ID] |
    """
    prot_ids = [int(i.split('_')[-1]) for i in frame]

    return abs(prot_ids[0] - prot_ids[1]) - 1


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


PIPELINE_LCA = (
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit name2taxid | cut -f 2 | paste -sd ',' | "
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit lca -s ',' | "
    "/home/niagara/Storage/MetaRus/k_petrikov/bin/taxonkit lineage -Ln -i 2 | cut -f 3"
)

args = parse_args()

input_clusters_table_path = args.clusters_table
input_defsys_summary_path = args.defsys_summary
input_accessions_path = args.accessions
input_pdb70_parent_path = args.pdb70_parent
input_pfama_parent_path = args.pfama_parent
input_padloc_data_path = args.padloc_data
output_folder = args.output_folder
cluster_size_threshold = args.cluster_size_threshold
top_coocc_number = args.top_coocc_number
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

# Proteins w/o defsys annotation
prots_with_dfnfnr_ann = df_clusters.loc[(df_clusters.DF_ann != 'miss')].Protein.tolist()

# Update 'Others_' fields
for brex_id, curr_reg in data_defsys.items():
    curr_other_prots = []

    for prot_id in curr_reg['Others_Prots']:
        curr_prot = f'{curr_reg["Nucleotide"]}_{prot_id}'
        if curr_prot not in prots_with_dfnfnr_ann:
            curr_other_prots.append(prot_id)

    data_defsys[brex_id]['Others_Prots'] = curr_other_prots
    data_defsys[brex_id]['Others_Count'] = len(curr_other_prots)

# Clusters coding
df_prots_regs = clusters_coding(data_defsys, df_clusters, threshold=cluster_size_threshold)
df_defsys = pd.DataFrame.from_dict(data_defsys, orient='index').drop(
    ['Start', 'End', 'Strand', 'Inner_DS_full', 'Accession'], axis=1)

# Pairs of most co-represented clusters
# Co-occurrence matrix
clusters_by_regions = df_defsys.loc[df_defsys.Others_Count > 0].Other_Cluster.to_list()
df_co_occurrence = co_occurence_matrix(clusters_by_regions)
coocc_data = top_co_occurred(df_co_occurrence, top_coocc_number)

coocc_prots_data = {}

for pair_clusters in coocc_data:
# for pair_clusters in [('CLUS_CP028720.1_1266', 'CLUS_NZ_OZ061338.1_1399')]:  # Test
    print(f'---Process {pair_clusters}---')
    j_idx, regs = jaccard_index(pair_clusters[0], pair_clusters[1], df_prots_regs)

    # Nucleotides and protein pairs
    pair_nucleotides = [i.split('%')[-1] for i in regs]

    prot_to_defsys = defaultdict(dict)

    for curr_defsys in regs:
        for curr_prot_id in data_defsys[curr_defsys]['Others_Prots']:
            curr_prot = f"{data_defsys[curr_defsys]['Nucleotide']}_{curr_prot_id}"
            prot_to_defsys[curr_prot]['DS_ID'] = curr_defsys

    pair_proteins = (pd.DataFrame.from_dict(prot_to_defsys, orient='index')
                       .rename_axis('Protein', axis=0)
                       .reset_index())
    pair_proteins = (pair_proteins.loc[
                         pair_proteins.Protein.isin(
                             df_clusters.loc[
                                 df_clusters.Cluster.isin(pair_clusters)
                             ].Protein
                         )
                     ].groupby('DS_ID').agg(list)
                     )

    # Check protein pairs lenght
    if not (pair_proteins.Protein.apply(len) == 2).all():
        print(f'Warning! For {pair_clusters} cluster appears more than once in single region!')

    # Representatives sequences
    fasta_reprs_acc = get_accession(pair_proteins.index[0].split('%')[-1], df_accessions)
    fasta_reprs_ids = pair_proteins.iat[0, 0]

    path_to_fasta_reprs = input_padloc_data_path / f'{fasta_reprs_acc}/{fasta_reprs_acc}_prodigal.faa'

    fasta_reprs = extract_sequence(path_to_fasta_reprs, fasta_reprs_ids)

    # Intergenic distances
    distances_all = pair_proteins.Protein.apply(get_dist_between_proteins)
    if distances_all.nunique() > 1:
        dist_stat = ','.join([str(i) for i in distances_all.unique()])
    else:
        dist_stat = str(distances_all.unique()[0])

    # Taxa, all genera
    pair_taxa = (df_accessions.loc[df_accessions.Nucleotide.isin(pair_nucleotides)]
                 .Species
                 .apply(lambda x: x.split(' ')[0])
                 .unique()
                 .tolist())

    # Find LCA
    tax_names = '\n'.join(pair_taxa)

    proc = subprocess.Popen(
        PIPELINE_LCA,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    stdout, stderr = proc.communicate(input=tax_names)

    if proc.returncode != 0:
        print(f'Error for cluster {pair_clusters}: {stderr}')
        tax_lca = 'ERROR'
    else:
        tax_lca = stdout.strip()

    # Write dict for JSON
    json_valid_key = '%'.join(pair_clusters)

    coocc_prots_data[json_valid_key] = {}
    coocc_prots_data[json_valid_key]['J_index'] = j_idx
    coocc_prots_data[json_valid_key]['Regs_number'] = len(regs)
    coocc_prots_data[json_valid_key]['Distance'] = dist_stat
    coocc_prots_data[json_valid_key]['Tax'] = pair_taxa
    coocc_prots_data[json_valid_key]['Tax_LCA'] = tax_lca
    coocc_prots_data[json_valid_key]['Size_1'] = int(
        df_clusters.loc[df_clusters.Cluster == pair_clusters[0]].iat[0, 2]
    )
    coocc_prots_data[json_valid_key]['PDB_ID_1'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[0]
    ].PDB_ID.iat[0]
    coocc_prots_data[json_valid_key]['PDB_Descr_1'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[0]
    ].PDB_ann.iat[0]
    coocc_prots_data[json_valid_key]['PDB_EVal_1'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[0]
    ].EVal.iat[0]
    coocc_prots_data[json_valid_key]['Size_2'] = int(
        df_clusters.loc[df_clusters.Cluster == pair_clusters[1]].iat[0, 2]
    )
    coocc_prots_data[json_valid_key]['PDB_ID_2'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[1]
    ].PDB_ID.iat[0]
    coocc_prots_data[json_valid_key]['PDB_Descr_2'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[1]
    ].PDB_ann.iat[0]
    coocc_prots_data[json_valid_key]['PDB_EVal_2'] = df_pdb70.loc[
        df_pdb70.Cluster == pair_clusters[1]
    ].EVal.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_ID_1'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[0]
    ].Pfam_ID.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_Descr_1'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[0]
    ].Pfam_ann.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_EVal_1'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[0]
    ].EVal.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_ID_2'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[1]
    ].Pfam_ID.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_Descr_2'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[1]
    ].Pfam_ann.iat[0]
    coocc_prots_data[json_valid_key]['Pfam_EVal_2'] = df_pfama.loc[
        df_pfama.Cluster == pair_clusters[1]
    ].EVal.iat[0]
    coocc_prots_data[json_valid_key]['Repr_1'] = fasta_reprs[0]
    coocc_prots_data[json_valid_key]['Repr_2'] = fasta_reprs[1]

# Most abundant clusters
selected_clusters = df_prots_regs.nlargest(top_clusters_number, 'Regs_Count').index

top_clusters_data = {}

for curr_cluster in selected_clusters:
    print(f'---Process {curr_cluster}---')

    # Representatives sequences
    fasta_reprs_acc = get_accession(df_clusters.loc[df_clusters.Cluster == curr_cluster].Nucleotide.iat[0],
                                    df_accessions)
    fasta_reprs_id = df_clusters.loc[df_clusters.Cluster == curr_cluster].Protein.iat[0]

    path_to_fasta_reprs = input_padloc_data_path / f'{fasta_reprs_acc}/{fasta_reprs_acc}_prodigal.faa'

    fasta_reprs = extract_sequence(path_to_fasta_reprs, (fasta_reprs_id,))

    top_clusters_data[curr_cluster] = {}
    top_clusters_data[curr_cluster]['Size'] = int(df_prots_regs.loc[curr_cluster, 'Regs_Count'])
    top_clusters_data[curr_cluster]['PDB_ID'] = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].PDB_ID.iat[0]
    top_clusters_data[curr_cluster]['PDB_Descr'] = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].PDB_ann.iat[0]
    top_clusters_data[curr_cluster]['PDB_EVal'] = df_pdb70.loc[df_pdb70.Cluster == curr_cluster].EVal.iat[0]
    top_clusters_data[curr_cluster]['Pfam_ID'] = df_pfama.loc[df_pfama.Cluster == curr_cluster].Pfam_ID.iat[0]
    top_clusters_data[curr_cluster]['Pfam_Descr'] = df_pfama.loc[df_pfama.Cluster == curr_cluster].Pfam_ann.iat[0]
    top_clusters_data[curr_cluster]['Pfam_EVal'] = df_pfama.loc[df_pfama.Cluster == curr_cluster].EVal.iat[0]
    top_clusters_data[curr_cluster]['Repr'] = fasta_reprs[0]

# Save data
(pd.DataFrame.from_dict(top_clusters_data, orient='index')
 .rename_axis('Cluster')
 .to_csv(output_folder / 'top_clusters_stat.tsv', sep='\t'))

with open(output_folder / 'top_clusters_stat.json', mode='w') as f:
    json.dump(top_clusters_data, f, indent=4)

(pd.DataFrame.from_dict(coocc_prots_data, orient='index')
 .rename_axis('Pairs')
 .to_csv(output_folder / 'cluster_pairs_stat.tsv', sep='\t'))

with open(output_folder / 'cluster_pairs_stat.json', mode='w') as f:
    json.dump(coocc_prots_data, f, indent=4)
