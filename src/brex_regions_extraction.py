"""v0.3d"""

import json
import subprocess
import pandas as pd
import re

from pathlib import Path


def coords_selector(x):
    """Auxiliary function for df.apply"""
    if x.name == 'start':
        return min(set(x))
    elif x.name == 'end':
        return max(set(x))
    elif x.name == 'strand':
        return x.unique()[0]


def count_flanks(x, add_up: int, add_down: int) -> tuple:
    """Auxiliary function for df.apply"""
    if x.strand == '+':
        return max(x.start - add_up, 1), x.end + add_down
    elif x.strand == '-':
        return max(x.start - add_down, 1), x.end + add_up


def create_defsys_summary(path_to_csv: str | Path) -> tuple:
    """
    Input: Padlock csv
    Output: DefSys table in a short, convenient format
    Ignores duplicate-systems completely
        (logged to the dupl_sys_ids)
    The region direction is selected by the brxC direction or by simple voting
        (logged to the variable antiparallel_sys)

    Field 'SysID':
        x.split('%')[-1] contig name;
        x.split('%')[-2] Padlock number of the system
    """

    sel_col_numbers = [0, 1, 2, 3, 6, 11, 12, 13]  # Just to skip unnecessary cols

    df = pd.read_csv(path_to_csv, dtype={'system.number': str}).iloc[:, sel_col_numbers]
    df['SysID'] = df['system'] + '%' + df['system.number'] + '%' + df['seqid']

    # Collection of systems with antiparallel genes if needed
    anti_sys = (df.loc[:, ('strand', 'SysID')]
                  .groupby('SysID')
                  .agg(lambda x: len(x.unique()))
                  .query('strand > 1'))
    anti_sys_names = [path_to_csv.stem[:15] + '_' + i for i in anti_sys.index.to_list()]
    
    # Choose uniform direction
    for curr_sys in df.SysID.unique():
        if not df.loc[(df.loc[:, 'SysID'] == curr_sys) & (df.loc[:, 'protein.name'] == 'BrxC')].empty:
            df.loc[(df.loc[:, 'SysID'] == curr_sys), ['strand']
                   ] = df.loc[(df.loc[:, 'SysID'] == curr_sys) & (df.loc[:, 'protein.name'] == 'BrxC'), 'strand'
                              ].values[0]
        else:
            df.loc[(df.loc[:, 'SysID'] == curr_sys), ['strand']
                   ] = (df.loc[(df.loc[:, 'SysID'] == curr_sys), ['strand', 'SysID']]
                          .groupby('SysID')
                          .agg(lambda x: x.value_counts().index[x.value_counts().argmax()])
                          .values[0][0])

    # Ignore multi systems
    prot_from_multi_sys = df['target.name'].value_counts()[df['target.name'].value_counts() > 1].index
    dupl_sys_ids = df[df['target.name'].isin(prot_from_multi_sys)].SysID.unique()
    dupl_sys_list = [path_to_csv.stem[:15] + '_' + i for i in dupl_sys_ids]
    df = df[~df.SysID.isin(dupl_sys_ids)]

    # Prot ID only number w/o nucleotide ID
    df.loc[:, ['target.name']] = df['target.name'].apply(lambda x: x.split('_')[-1])

    # Final df
    df_result = df.loc[:, ['SysID', 'start', 'end', 'strand']].groupby('SysID').agg(coords_selector)
    df_result['proteins_ids'] = (df.loc[:, ('target.name', 'SysID')]
                                 .groupby('SysID')
                                 .agg(lambda x: x.unique())
                                 )
    df_result.reset_index(inplace=True)
    
    # return df_result
    return df_result, dupl_sys_list, anti_sys_names


def extract_region_gff(gff_file: str | Path,
                       chrom_id: str,
                       start_coord: int,
                       end_coord: int,
                       extracted_gff_file: str | Path
                       ) -> None:
    """
    Extract given interval from input gff-file
    and save it to output gff-file

    :param gff_file:
    :param chrom_id:
    :param start_coord:
    :param end_coord:
    :param extracted_gff_file:
    :return:
    """

    cmd = ['bedtools', 'intersect',
           '-a', gff_file,  # Input GFF file
           '-b', 'stdin',
           '-wa'
           ]

    bed_region = f'{chrom_id}\t{start_coord - 1}\t{end_coord}\n'

    with open(extracted_gff_file, mode='w') as file:
        subprocess.run(cmd, input=bed_region, text=True, stdout=file)


def add_loc_marks_to_gff(gff_file: str, defsys_type, defsys_params):
    gff_df_columns = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Comment']
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, names=gff_df_columns)

    gff_df['GeneID'] = gff_df.Comment.apply(lambda x: re.search(r'ID=\d+_(\d+)', x).group(1))
    gff_df['Loc_mark'] = 'No'
    
    # defsys_prot_ids = [_.split('_')[-1] for _ in defsys_params['proteins_ids']]
    # gff_df.loc[gff_df['GeneID'].isin(defsys_prot_ids), ['Loc_mark']] = defsys_type.split('%')[0]
    gff_df.loc[gff_df['GeneID'].isin(defsys_params['proteins_ids']), ['Loc_mark']] = defsys_type.split('%')[0]

    if defsys_params['strand'] == '+':
        flank_types = ('upstream_', 'downstream_')
    else:  # defsys_params['strand'] == '-'
        flank_types = ('downstream_', 'upstream_')

    flank_size = gff_df.loc[gff_df['Start'] < defsys_params['start']].shape[0]
    gff_df.loc[gff_df['Start'] < defsys_params['start'], ['Loc_mark']] = [flank_types[0] + str(i) for i in
                                                                          range(flank_size, 0, -1)]

    flank_size = gff_df.loc[gff_df['End'] > defsys_params['end']].shape[0]
    gff_df.loc[gff_df['End'] > defsys_params['end'], ['Loc_mark']] = [flank_types[1] + str(i) for i in
                                                                      range(1, flank_size + 1, 1)]

    gff_df.loc[gff_df['Loc_mark'] == 'No', ['Loc_mark']] = 'inner'

    gff_df.loc[:, 'Comment'] = gff_df.Comment + 'loc_mark=' + gff_df.Loc_mark + ';'

    return gff_df


def check_bounds(df, start, end, up_flank, down_flank, assembly, defsys, strand) -> str:
    """
    Check if given flanks extend beyond the nucleotide boundaries
    """
    if strand == '+':
        bound_marks = ('upstream', 'downstream')
    else:  # strand == '-'
        bound_marks = ('downstream', 'upstream')
    
    bounds = []
    
    up_bound_interval = start - df.Start.iloc[0]
    down_bound_interval = df.End.iloc[df.shape[0] - 1] - end
    if up_bound_interval < up_flank:
        bounds.append(assembly + '_' + defsys + '\t' + bound_marks[0] + '\t' + str(up_flank - up_bound_interval))
    if down_bound_interval < down_flank:
        bounds.append(assembly + '_' + defsys + '\t' + bound_marks[1] + '\t' + str(down_flank - down_bound_interval))
   
    return '\n'.join(bounds)


def read_fasta(fasta_file):
    name, seq = None, []
    for line in fasta_file:
        if line.startswith('>'):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


def extract_region_fasta(fasta_file: Path,
                         prot_ids: list,
                         extracted_fasta_file: Path) -> None:
    """
    Extracts sequences with specific IDs from fasta-file
    """

    count, stop_num = 0, len(prot_ids)
    
    with open(fasta_file, mode='r') as f_read, open(extracted_fasta_file, mode='w') as f_write:
        for name, seq in read_fasta(f_read):
            if name.strip()[1:] in prot_ids:
                f_write.write(name + seq)
                count += 1
                if count == stop_num:
                    break


# #### Regs summary

# In[8]:


def create_region_summary(df_with_marks: pd.DataFrame, direction: str, meta_info: tuple) -> dict:
    """
    meta_info = (defsys_type, nucleotide, accession)
    """
    
    region_proteins = df_with_marks.GeneID.to_list()
    
    if direction == '-':
        region_proteins = region_proteins[::-1]
        defsys_idxs = list(df_with_marks.set_axis([_ for _ in range(df_with_marks.shape[0])][::-1], axis='index'
                                                  ).query('Loc_mark == @meta_info[0]').index)[::-1]
        inner_idxs = list(df_with_marks.set_axis([_ for _ in range(df_with_marks.shape[0])][::-1], axis='index'
                                                 ).query('Loc_mark == "inner"').index)
    else:  # direction == '+'
        defsys_idxs = list(df_with_marks.query('Loc_mark == @meta_info[0]').index)
        inner_idxs = list(df_with_marks.query('Loc_mark == "inner"').index)
        
    final_region = {'proteins_ids': region_proteins, 
                    'defsys_idxs': defsys_idxs,
                    'inner_idxs': inner_idxs,
                    'defsys_type': meta_info[0],
                    'nucleotide': meta_info[1],
                    'accession': meta_info[2]
                    }
    
    return final_region


input_def_sys_list = ['brex']
input_folder_path = Path('/home/holydiver/Main/2024_BREX/Data/Prrr_Data/01_Padloc_hmm_custom/padloc_hmm_custom')
output_path_regions = Path('/home/holydiver/Main/2024_BREX/Data/New_data/output')
output_path_summaries = Path('/home/holydiver/Main/2024_BREX/Data/New_data')
input_up_flank, input_down_flank = 10000, 10000

dupl_sys = []
antiparallel_sys = []
no_selected_defsys_genomes = []  # List of nucleotides w/o DefSys of interest
no_full_flanks = []  # List of nucleotides where the distance to the edge is less than the given flank

regions_summary_collection = []  # List of summaries for defsys regions for JSON-file

output_path_summaries.mkdir(parents=True, exist_ok=True)
output_path_regions.mkdir(parents=True, exist_ok=True)

# Logs writing
logs = (dupl_sys, antiparallel_sys, no_selected_defsys_genomes, no_full_flanks)
f_names = ('duplicates.txt',
           'antiparallel.txt',
           'no_selected_defsys.txt',
           'no_full_flanks.txt'
           )

for log, f_name in zip(logs, f_names):
    with open(output_path_summaries / f_name, mode='w') as f:
        for i in log:
            f.write(i + '\n')

with open(output_path_summaries / 'regions_summary.json', mode='w') as f:
    json.dump(regions_summary_collection, f, indent=4)
