####################
# import libraries #
####################
import argparse
import pandas as pd

from pathlib import Path
from Step1mod import step1

parser = argparse.ArgumentParser(
                    prog='df_pl_merge',
                    description='merge defensefinder and padloc results in one csv file for each sample')

parser.add_argument('-s',
                    '--samples_list_path',
                    help='Provide path to txt file with each line representing comma separated: '
                         'sample_id, path_to_defensefinder_results, path_to_padloc_results',
                    default='samples.txt')

args = parser.parse_args()
samples_list_path = args.samples_list_path

if Path(samples_list_path).exists():
    df_pl_viz = pd.read_csv(samples_list_path, names = ['sample_id', 'df_dir', 'pl_dir'])
else:
    raise Exception(f'Path {samples_list_path} not found')

step1(df_pl_viz)
