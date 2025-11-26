import json
import pandas as pd

from pathlib import Path
from commons import (make_unidir_genes_defsys,
                     find_redundancy_defsys,
                     create_defsys_summary)

try:
    from parse_single_tables import parse_single_dfnfnr_table
except ImportError as e:
    raise ImportError(f'{e}\nFunctions for parsing single tables must be implemmented!')

pd.options.mode.copy_on_write = True


def process_dfnfnr_data(
        input_samples: list[tuple],
        results_path: str | Path
) -> None:
    """

    Args:
        input_samples:
        results_path:

    Returns:

    """
    curr_accession_result_path = Path(results_path) / 'By_Accessions'
    curr_accession_result_path.mkdir(parents=True, exist_ok=True)

    summary_defsys_all = {}
    redundancy_defsys_all = []
    protein_annotations_all = []
    antidefense_all = []

    cols_order = ['Accession', 'Nucleotide', 'System', 'System_sub',
                  'Start', 'End', 'Strand', 'DS_Prots', 'Have_inner', 'Activity']

    for sample_id, dfnfnr_dir, gff_file in input_samples:
        for file_name in Path(dfnfnr_dir).iterdir():
            if str(file_name).endswith('genes.tsv'):
                print(f'---Processing {sample_id}---')

                df_curr = parse_single_dfnfnr_table(file_name, gff_file, sample_id)

                antidefense_all.append(
                    df_curr.loc[df_curr.Activity == 'Antidefense', ['Accession', 'DS_ID']]
                )

                redundancy_defsys_all.append(
                    find_redundancy_defsys(df_curr)
                )

                protein_annotations_all.append(
                    df_curr[['DS_ID', 'Accession', 'Nucleotide', 'Protein',
                             'Annotation', 'System', 'System_sub']]
                )

                make_unidir_genes_defsys(df_curr)
                summary_curr = create_defsys_summary(df_curr)[cols_order]
                summary_defsys_all.update(summary_curr.to_dict(orient='index'))

                # Write current results
                summary_curr.to_json(curr_accession_result_path / f'{sample_id}_summary.json', orient='index')

    # Write results
    if redundancy_defsys_all:
        (
            pd.concat(redundancy_defsys_all)
              .reset_index(drop=True)
              .to_csv(results_path / 'redundant_defsys.tsv', sep='\t', index=False)
        )

    if antidefense_all:
        (
            pd.concat(antidefense_all)
              .to_csv(results_path / 'antidefense_systems.tsv',
                      sep='\t',
                      index=False)
        )

    (
        pd.concat(protein_annotations_all)
          .reset_index(drop=True)
          .to_csv(Path(results_path) / 'protein_annotations.tsv', sep='\t', index=False)
    )

    with open(Path(results_path) / 'defsys_summary.json', mode='w') as f:
        json.dump(summary_defsys_all, f, indent=4)
