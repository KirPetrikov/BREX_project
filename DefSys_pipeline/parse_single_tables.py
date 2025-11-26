"""
These functions must be implemented to keep data consistency.

The formatting of output tables should be pandas DataFrame with columns:

'Accession': str - Sample identifier for a single organism (e.g., GenBank accession number or MAG number)
'Nucleotide': str - Specific nucleotide/contig identifier; each must be unique among all samples
'System': str - Defense system name
'System_sub': str [Only for DefenseFinder] - Defense system subtype name
'Protein': str - Gene/protein identifier; each must be unique, of the form '<Nucleotide identifier>_<number>'
'Start': int - gene 5'-coordinate
'End': int - gene 3'-coordinate
'Strand': str - Chain identifier; must be '+'/'-'
'Activity': str [Only for DefenseFinder] - Protection type; must be 'Defense'/'Antidefense'
"""

# def parse_single_dfnfnr_table(
#         path_to_tsv: str | Path,
#         path_to_gff: str | Path,
#         sample_id: str
# ) -> pd.DataFrame:
#
#     pass
#
#     return ...
#
# def parse_single_padloc_table(
#         path_to_csv: str | Path,
#         sample_id: str
# ) -> pd.DataFrame:
#
#     pass
#
#     return ...
