"""
v0.4b
"""
import json
import pandas as pd
import sys

from collections import defaultdict
from pathlib import Path

BREX_SET = {'BrxA', 'BrxB', 'BrxC',
            'BrxD', 'BrxE', 'BrxF',
            'BrxHII', 'BrxHI', 'BrxL',
            'BrxP', 'PglW', 'PglXI',
            'PglX', 'PglZ'
            }


def modify_ann_prot(df: pd.DataFrame, prot_name: str, cutoff: int = 0) -> set:
    """

    :param df:
    :param prot_name:
    :param cutoff:
    :return:
    """
    clusters_list = (
        df.query('Annotation == @prot_name & Cluster_size > @cutoff')['Cluster']
          .unique()
                     )

    annot_to_reset = {
        k: prot_name + '%CLUS_' + str(v) for k, v in zip(clusters_list,
                                                         range(clusters_list.shape[0])
                                                         )
                      }

    df.loc[df['Cluster'].isin(clusters_list), 'Annotation'] = (
        df.loc[df['Cluster'].isin(clusters_list), ['Cluster']]
          .apply(lambda x: annot_to_reset[str(x.values[0])], axis=1)
                                                                )

    return {_ for _ in annot_to_reset.values()}


def mask_singletons(df: pd.DataFrame, cutoff: int = 1) -> None:
    """

    :param df:
    :param cutoff:
    :return:
    """
    df.loc[df['Cluster_size'] < cutoff, ['Annotation']] = '%MASK' + str(cutoff)


input_annot = Path('/home/holydiver/Main/2024_BREX/Data/Prrr_Data/2024_05_07_full_annotation.tsv')
input_path = Path('/home/holydiver/Main/2024_BREX/Data/20250118_exctracted')
input_json_summary = 'regions_summary.json'

consider_upstream = True
prots_to_consider_clus = ('WYL', 'BrxC')
prot_singl_cutoff = 0


output_path = Path('/home/holydiver/Main/2024_BREX/Data/New_data')
output_tsv_name = '20250124_unique_regions_summary.tsv'
output_json_name = '20250124_unique_regions_accessions.json'


annot_df = pd.read_csv(input_annot, sep='\t').loc[:, ['Protein', 'Annotation', 'Cluster', 'Cluster_size']]

curr_brex_set = BREX_SET

# for protein in prots_to_consider_clus:
#     annot_to_add = modify_ann_prot(annot_df, protein)
#     if protein in BREX_SET:
#         curr_brex_set.update(_ for _ in annot_to_add)

# mask_singletons(annot_df)

annot_dict = dict(zip(annot_df.Protein, annot_df.Annotation))

with open(input_path / input_json_summary) as f:
    regions_summary = json.load(f)

unique_regions_summary = {}
unique_regions_accessions = defaultdict(list)

for curr_region in regions_summary:
    # For regions with upstream gene
    if curr_region['defsys_idxs'][0] != 0:

        if consider_upstream:
            upstream_idx = curr_region['defsys_idxs'][0] - 1
        else:
            upstream_idx = curr_region['defsys_idxs'][0]

        end_idx = curr_region['defsys_idxs'][-1] + 1

        # Unique region
        curr_region_selected_proteins = [
            curr_region['nucleotide'] + '_' + i for i in curr_region['proteins_ids'][upstream_idx:end_idx]
                                         ]
        annotated_region = ','.join(
            [annot_dict[protein] for protein in curr_region_selected_proteins]
                                    )

        if annotated_region in unique_regions_summary:
            unique_regions_summary[annotated_region]['Counts'] += 1
            unique_regions_summary[annotated_region]['Defsys_type'].add(curr_region['defsys_type'])
        else:
            unique_regions_summary[annotated_region] = {}
            unique_regions_summary[annotated_region]['Counts'] = 1
            unique_regions_summary[annotated_region]['Defsys_type'] = {curr_region['defsys_type']}

        unique_regions_accessions[annotated_region].append(curr_region['accession'] +
                                                           '%' +
                                                           curr_region['nucleotide']
                                                           )

    else:
        print('Warning:',
              curr_region['accession'] + '%' + curr_region['nucleotide'] + '%' + curr_region['defsys_type'],
              'does not has upstream gene!',
              file=sys.stderr
              )

unique_regions_df = (pd.DataFrame.from_dict(unique_regions_summary, orient='index')
                       .rename_axis('Region')
                       .reset_index()
                     )

if (unique_regions_df.Defsys_type.apply(lambda x: len(x)) > 1).any():
    print('Warning: Some regions are assigned a non-unique defense system type!', file=sys.stderr)
else:  # Defsys_type-set cast to string
    unique_regions_df.loc[:, ['Defsys_type']] = (
        unique_regions_df.loc[:, 'Defsys_type'].astype(str).apply(lambda x: x[2:-2])
                                                 )

unique_regions_df['Inner'] = (
    unique_regions_df.Region.apply(lambda x: bool(set(x.split(',')[1:]).difference(curr_brex_set)))
                              )

unique_regions_df.to_csv(output_path / output_tsv_name, sep='\t', index=False)

with open(output_path / output_json_name, mode='w') as f:
    json.dump(unique_regions_accessions, f, indent=4)
