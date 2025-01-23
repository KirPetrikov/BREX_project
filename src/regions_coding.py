import json
import pandas as pd
import sys

from collections import defaultdict
from pathlib import Path


def modify_ann_prot(df: pd.DataFrame, prot_name: str, cutoff: int = 0):
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


def mask_singletons(df: pd.DataFrame, cutoff: int = 1) -> None:
    """

    :param df:
    :param cutoff:
    :return:
    """
    df.loc[df['Cluster_size'] < cutoff, ['Annotation']] = '%MASK' + str(cutoff)

# input_annot = Path('/home/holydiver/Main/2024_BREX/Data/Prrr_Data/2024_05_07_full_annotation.tsv')
# input_path = Path('/home/holydiver/Main/2024_BREX/Data/20250118_exctracted')
#
# tmpa = annot_df.iloc[10:15]
# tmpb = annot_df.query('Annotation == @prot_to_consider_clus').iloc[:15]
# test_df = pd.concat([tmpa, tmpb])
#

#


input_annot = Path('/home/holydiver/Main/2024_BREX/Data/Prrr_Data/2024_05_07_full_annotation.tsv')
input_path = Path('/home/holydiver/Main/2024_BREX/Data/TMP/Exmpls/')
input_json_summary = Path('defsys_regions_summary.json')

consider_upstream = True
prot_to_consider_clus = 'WYL'
prot_singl_cutoff = 0
singl_cutoff_total = 0

output_path = Path('./')
output_json_name = Path('summary.json')


annot_df = pd.read_csv(input_annot, sep='\t').loc[:, ['Protein', 'Annotation', 'Cluster', 'Cluster_size']]

# modify_ann_prot(test_df, prot_to_consider_clus)
# mask_singletons(annot_df)

annot_dict = dict(zip(annot_df.Protein, annot_df.Annotation))

with open(input_path / input_json_summary) as f:
    regions_summary = json.load(f)

unique_regions_summary = defaultdict(lambda: {'counts': 0,
                                              'inner_presents': True,
                                              'defsys_type': [],
                                              'accession': []}
                                     )

for curr_region in regions_summary:
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
            [annot_dict.get(protein, protein) for protein in curr_region_selected_proteins]
                                    )

        unique_regions_summary[annotated_region]['counts'] += 1
        unique_regions_summary[annotated_region]['defsys_type'].append(curr_region['defsys_type'])
        unique_regions_summary[annotated_region]['accession'].append(curr_region['accession'] +
                                                                     '%' +
                                                                     curr_region['nucleotide']
                                                                     )
        unique_regions_summary[annotated_region]['inner_presents'] = bool(curr_region['inner_idxs'])

    else:
        print('Warning:',
              curr_region['accession'] + '%' + curr_region['nucleotide'],
              'does not has upstream gene!',
              file=sys.stderr
              )
