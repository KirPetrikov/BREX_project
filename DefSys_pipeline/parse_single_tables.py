"""
These functions must be implemented to keep data consistency.

Input arguments of functions must be:
For 'parse_single_dfnfnr_table':
- path to DefenseFinder genes tsv-file
- path to gff-file for getting genes coordinates
- sample id
For 'parse_single_padloc_table':
- path to Padloc systems csv-file
- sample id

NB:  'sample id' not nesserary for using in function, but must be present for consistency

The output tables should be pandas DataFrames with columns:
'Accession': str - Sample identifier for a single organism (e.g., GenBank accession number or MAG number)
'Nucleotide': str - Specific nucleotide/contig identifier; each must be unique among all samples
'System': str - Defense system name
'System_sub': str [Only for DefenseFinder] - Defense system subtype name
'Protein': str - Gene/protein identifier;
                 each must be unique, of the form '<Nucleotide/smple identifier>_<number: int>'
'Start': int - gene 5'-coordinate
'End': int - gene 3'-coordinate
'Strand': str - Chain identifier; must be '+'/'-'
'Activity': str [Only for DefenseFinder] - Protection type; must be 'Defense'/'Antidefense'

--- Pattern:

def parse_single_dfnfnr_table(
        path_to_tsv: str | Path,
        path_to_gff: str | Path,
        sample_id: str
) -> pd.DataFrame:

    pass

    return df


def parse_single_padloc_table(
        path_to_csv: str | Path,
        sample_id: str
) -> pd.DataFrame:

    pass

    return df

"""

