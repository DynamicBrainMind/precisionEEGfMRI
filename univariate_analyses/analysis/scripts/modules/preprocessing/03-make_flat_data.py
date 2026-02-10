import mne
import numpy as np
import pandas as pd
from pyprojroot import here
import os
from glob import glob
from pathlib import Path
import pickle
from itertools import product
import sys

'''
This script converts the nested dict data in formatted/full to a flat data
frame that can be used to compute correlations. Output is stored in 
'''



def process_subject(subject, col_names, data_path):

    file = data_path / Path(f'{subject}.pkl')
    with open(file, 'rb') as file:
        data = pickle.load(file)

    subject_d = pd.DataFrame()

    for session in data:
        for run in data[session]:

            run_data = data[session][run]
            run_df = pd.DataFrame(run_data['X'], columns=col_names)
            run_df.insert(0, 'subject', subject)
            run_df.insert(1, 'session', session)
            run_df.insert(2, 'run', run)
            trs = data[session][run]['trs']
            run_df.insert(3, 'tr', trs)
            
            networks = list(data[session][run]['y'].keys())
            for idx, network in enumerate(networks, start=4):
                run_df.insert(idx, network, data[session][run]['y'][network])

            subject_d = pd.concat([subject_d, run_df], axis=0)

    return subject_d


def make_column_names(channels, freqs, num_lags):

    lags = list(range(num_lags))
    
    col_names_tuple = list(product(channels, freqs, lags))

    col_names = []
    for name_tuple in col_names_tuple:
        name = '_'.join([str(x) for x in name_tuple])
        col_names.append(name)

    return col_names


if __name__ == '__main__':


    # Parse input args
    args = sys.argv[1:]

    # If no path to data given, use analysis/data/original
    if not args:
        data_root = Path('analysis/data/original')
        os.chdir(here())
        sys.path.append(str(here()))
    else:
        data_root = Path(args[0])

    # Get channel names
    raw = mne.io.read_raw_eeglab(f'{data_root}/sub-001/ses-001/eeg/sub-001_ses-001_bld001_eeg_Bergen_CWreg_filt_ICA_rej.set')
    channels = raw.info['ch_names']

    # Get freqs
    with open(data_root / Path('../formatted/freqs.txt'), 'r') as file:
        freqs = file.read().splitlines()

    # Get number of lags
    with open(data_root / Path('../formatted/num_lags.txt'), 'r') as file:
        num_lags = int(file.read())

    col_names = make_column_names(channels, freqs, num_lags)

    # Get subjects
    subjects = [Path(x).stem for x in glob(f'{data_root}/../formatted/full/sub[-_][0-9]*')]
    subjects = sorted(subjects, key=lambda x: int(x.split('-')[1]))

    d = pd.DataFrame()

    # Process each subject and concatenate
    for subject in subjects:
        print(f'\n---PROCESSING SUBJECT {subject} of {subjects[-1]}---\n')
        new_d = process_subject(subject, col_names,
                                data_path=Path(f'{data_root}/../formatted/full'))
        d = pd.concat([d, new_d], axis=0)

    # Write merged data and channel names to file
    print('\nWriting out full data...')
    out_path = Path(f'{data_root}/../correlation_data')
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    d.reset_index(drop=True).to_feather(out_path / Path('merged_data.feather'))
    with open(f'{data_root}/../correlation_data/ch_names.txt', 'w') as file:
        file.write('\n'.join(channels))

