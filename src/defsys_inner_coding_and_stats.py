"""v0.1
Coding of the composition of a region from inner genes into corresponding clusters.
Statistics on clusters and proteins.
"""

import pandas as pd
import json

from collections import defaultdict
from pathlib import Path
from defsys import co_occurence_matrix, top_co_occurred, parse_gff, add_gene_id_to_gff

pd.options.mode.copy_on_write = True


def clusters_coding(data, clusters, threshold=2):
    """
    Add to all DS_ID in defense systems summary
    field 'Other_Cluster' with list of regs_composition.
    Only for other_proteins.

    Create dataframe:
    Cluster - List of DS_IDs with it - Number of regions

    :param data: Defense systems summary
    :param clusters: Clusters dataframe
    :param threshold: Minimal size of clusters to be considered
    :return: DataFrame
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


def jaccard_index(id_1, id_2, df) -> tuple[float, set]:
    """

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


def get_dist_between_proteins(frame, path_to_gff_parent, df_acc):
    """
    Look in gff-file for two proteins and returns
    distance between them (in b.p.)

    For use in .apply() on DataFrame:
        |    index   |       'Protein'      |
        |------------|----------------------|
        | nucleotide | [Prot1_ID, Prot2_ID] |
    """
    prot_pair = frame.Protein
    # nucl = frame.name
    # accession = df_acc.loc[df_acc.Nucleotide == nucl].iat[0, 0]

    accession = get_accession(frame.name, df_acc)

    path_to_gff = path_to_gff_parent / f'{accession}/{accession}_prodigal.gff'
    gff = parse_gff(path_to_gff)
    gff = add_gene_id_to_gff(gff)

    coords = gff.loc[gff.Gene_ID.isin(sorted([prot_pair[0], prot_pair[1]]))
                     ][['Start', 'End']].values.flatten()

    return coords[2] - coords[1]


def read_fasta(fasta_file_path):
    """
    Read fasta-file
    Returns (>fasta_header, sequence) as in source
    (does not create one-line-fasta)
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

    with open(fasta_file_path, mode='r') as f:
        for name, seq in read_fasta(f):
            header = name.strip()[1:]
            if header in prot_ids:
                seq = seq.replace("\n", "")
                selected_seqs.append(f'{header}:{seq}')

                count += 1
                if count == stop_num:
                    break

    return selected_seqs


input_clusters_table_path = Path(
    '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250413_brex_regions_clusters/Region/clusters_table.tsv')

input_defsys_summary_path = Path(
    '/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250511_brex_regions_summary/inner_summary_brex_recalc.json')

input_accessions_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/Accessions_summary.tsv')

# Only Padloc annotations
input_annot_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/Curr_protein_annotations.tsv')

input_reg_type = 'region'
input_pdb70_parent_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250503_brex_annot_pdb70')
input_pfama_parent_path = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250503_brex_annot_pfama')

input_padloc_data_path = Path('/home/niagara/Storage/MetaRus/Common_dir/collections/complete_bacteria_collection/padloc')

output_folder = Path('/home/niagara/Storage/MetaRus/k_petrikov/2024_BREX/Data/20250511_brex_regions_summary/'
                     '20250605_inner_clus_stat')
output_folder.mkdir(parents=True, exist_ok=True)

pdb70_path = input_pdb70_parent_path / f'{input_reg_type}_pdb_top_annot.tsv'
pfama_path = input_pfama_parent_path / f'{input_reg_type}_pfama_top_annot.tsv'

# Read data
with open(input_defsys_summary_path) as f:
    data_defsys = json.load(f)

df_clusters = pd.read_csv(input_clusters_table_path, sep='\t')
df_accessions = pd.read_csv(input_accessions_path, sep='\t')

df_annotations = pd.read_csv(input_annot_path,
                             sep='\t',
                             dtype={'Protein': str, 'Padloc_ann': str, 'System': str, 'Nucleotide': str,
                                    'Start': int, 'End': int, 'Strand': str, 'Localisation': str})

df_pdb70 = pd.read_csv(pdb70_path, sep='\t',
                       dtype={'Cluster': str, 'PDB_ID': str, 'PDB_ann': str, 'EVal': float, 'Prob': float})

df_pfama = pd.read_csv(pfama_path, sep='\t',
                       dtype={'Cluster': str, 'Pfam_ID': str, 'Pfam_ann': str, 'EVal': float, 'Prob': float})

# Add annotations, Padloc, PDB70, PfamA
df_clusters = df_clusters.merge(df_annotations,
                                how='left',
                                on=['Protein']
                                ).fillna('miss')

m = df_pdb70[['Cluster', 'PDB_ann']].set_index('Cluster').to_dict()['PDB_ann']
df_clusters['PDB_ann'] = df_clusters['Cluster'].map(m)

m = df_pfama[['Cluster', 'Pfam_ann']].set_index('Cluster').to_dict()['Pfam_ann']
df_clusters['Pfam_ann'] = df_clusters['Cluster'].map(m)

# Add nucleotides
df_clusters.loc[df_clusters.Nucleotide == 'miss', 'Nucleotide'] = df_clusters.loc[
    df_clusters.Nucleotide == 'miss'].Protein.apply(lambda x: '_'.join(x.split('_')[:-1]))

# Drop unnecessary columns
# df_clusters = df_clusters.drop(['Start', 'End', 'Strand', 'Localisation'], axis=1)

# Clusters coding
df_prots_regs = clusters_coding(data_defsys, df_clusters, threshold=4)

df_defsys = pd.DataFrame.from_dict(data_defsys, orient='index').drop(
    ['Start', 'End', 'Strand', 'Inner_DS_full', 'Accession'], axis=1)

# Co-occurrence matrix
clusters_by_regions = df_defsys.loc[df_defsys.Others_Count > 0].Other_Cluster.to_list()

df_co_occurrence = co_occurence_matrix(clusters_by_regions)

coocc_data = top_co_occurred(df_co_occurrence, 20)

# Skip Paris-proteins
# del coocc_data[('CLUS_NZ_CP113504.1_1637', 'CLUS_NZ_CP149360.1_435')]

coocc_prots_data = {}

# for pair in [('CLUS_CP028720.1_1266', 'CLUS_NZ_OZ061338.1_1399')]:  # Test
for pair in coocc_data:
    j_idx, regs = jaccard_index(pair[0], pair[1], df_prots_regs)

    ann_1_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[0]].PDB_ann.iat[0]
    ann_2_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[1]].PDB_ann.iat[0]
    id_1_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[0]].PDB_ID.iat[0]
    id_2_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[1]].PDB_ID.iat[0]
    eval_1_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[0]].EVal.iat[0]
    eval_2_pdb = df_pdb70.loc[df_pdb70.Cluster == pair[1]].EVal.iat[0]
    ann_1_pfama = df_pfama.loc[df_pfama.Cluster == pair[0]].Pfam_ann.iat[0]
    ann_2_pfama = df_pfama.loc[df_pfama.Cluster == pair[1]].Pfam_ann.iat[0]
    id_1_pfama = df_pfama.loc[df_pfama.Cluster == pair[0]].Pfam_ID.iat[0]
    id_2_pfama = df_pfama.loc[df_pfama.Cluster == pair[1]].Pfam_ID.iat[0]
    eval_1_pfama = df_pfama.loc[df_pfama.Cluster == pair[0]].EVal.iat[0]
    eval_2_pfama = df_pfama.loc[df_pfama.Cluster == pair[1]].EVal.iat[0]

    clust1_size, clust2_size = (
        df_clusters.loc[df_clusters.Cluster == pair[0]].iat[0, 2],
        df_clusters.loc[df_clusters.Cluster == pair[1]].iat[0, 2]
    )

    # Список нуклеотидов и пар белков
    pair_nucleotides = [i.split('%')[-1] for i in regs]

    pair_proteins = (
        df_clusters.iloc[:, [0, 1, 5]]
                   .loc[
            df_clusters.Cluster.isin(pair)
            &
            df_clusters.Nucleotide.isin(pair_nucleotides)
                        ]
                   .groupby('Nucleotide')
                   .agg(list)
                   .iloc[:, [0]]
    )

    # Все белки должы отобраться парами
    try:
        assert (pair_proteins.Protein.apply(len) == 2).all()
    except AssertionError:
        print(f'Fail for {pair}')

    # Representatives sequences
    fasta_reprs_acc = get_accession(pair_proteins.index[0], df_accessions)
    fasta_reprs_ids = pair_proteins.iat[0, 0]

    path_to_fasta_reprs = input_padloc_data_path / f'{fasta_reprs_acc}/{fasta_reprs_acc}_prodigal.faa'

    fasta_reprs = extract_sequence(path_to_fasta_reprs, fasta_reprs_ids)

    # Intergenic distances
    distances_all = pair_proteins.apply(
        get_dist_between_proteins, args=(input_padloc_data_path, df_accessions), axis=1
    )
    dist_stat = f'{distances_all.mean():.0f}\u00B1{distances_all.sem():.0f}'

    # Taxa, only genus
    pair_taxa = (df_accessions.loc[df_accessions.Nucleotide.isin(pair_nucleotides)]
                 .Species
                 .apply(lambda x: x.split(' ')[0])
                 .unique()
                 .tolist())

    # Full taxa
    # df_acc.loc[df_acc.Nucleotide.isin(curr_nucls)].Species.unique().tolist()

    json_valid_key = '%'.join(pair)

    coocc_prots_data[json_valid_key] = {}
    coocc_prots_data[json_valid_key]['J_index'] = j_idx
    coocc_prots_data[json_valid_key]['Regs_number'] = len(regs)
    coocc_prots_data[json_valid_key]['Distance'] = dist_stat
    coocc_prots_data[json_valid_key]['Tax'] = pair_taxa
    coocc_prots_data[json_valid_key]['Size_1'] = int(clust1_size)
    coocc_prots_data[json_valid_key]['PDB_ID_1'] = id_1_pdb
    coocc_prots_data[json_valid_key]['PDB_Descr_1'] = ann_1_pdb
    coocc_prots_data[json_valid_key]['PDB_EVal_1'] = eval_1_pdb
    coocc_prots_data[json_valid_key]['Size_2'] = int(clust2_size)
    coocc_prots_data[json_valid_key]['PDB_ID_2'] = id_2_pdb
    coocc_prots_data[json_valid_key]['PDB_Descr_2'] = ann_2_pdb
    coocc_prots_data[json_valid_key]['PDB_EVal_2'] = eval_2_pdb
    coocc_prots_data[json_valid_key]['Pfam_ID_1'] = id_1_pfama
    coocc_prots_data[json_valid_key]['Pfam_Descr_1'] = ann_1_pfama
    coocc_prots_data[json_valid_key]['Pfam_EVal_1'] = eval_1_pfama
    coocc_prots_data[json_valid_key]['Pfam_ID_2'] = id_2_pfama
    coocc_prots_data[json_valid_key]['Pfam_Descr_2'] = ann_2_pfama
    coocc_prots_data[json_valid_key]['Pfam_EVal_2'] = eval_2_pfama
    coocc_prots_data[json_valid_key]['Repr_1'] = fasta_reprs[0].split(':')[1]
    coocc_prots_data[json_valid_key]['Repr_2'] = fasta_reprs[1].split(':')[1]

# Save data
(pd.DataFrame.from_dict(coocc_prots_data, orient='index')
   .rename_axis('Pairs')
   .to_csv(output_folder / 'cluster_pairs_stat.tsv', sep='\t'))

with open(output_folder / 'cluster_pairs_stat.json', mode='w') as file:
    json.dump(coocc_prots_data, file, indent=4)
