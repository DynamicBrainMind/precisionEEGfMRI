import pickle
import numpy as np
from pyprojroot import here
from pathlib import Path
from glob import glob
import sys
import os
from itertools import product
import pandas as pd

# Needs commenting...


def get_col_names(data_root):

    # Import frequencies
    with open(data_root / Path('../freqs.txt'), 'r') as file:
        freqs = [float(x) for x in file.read().splitlines()]

    if all([not x % 1 for x in freqs]):
        freqs = [int(x) for x in freqs]

    freqs = [str(x) for x in freqs]

    # Import channel names
    with open(data_root / Path('../../correlation_data/ch_names.txt'), 'r') as file:
        ch_names = file.read().splitlines()

    # Import lags
    with open(data_root / Path('../num_lags.txt'), 'r') as file:
        num_lags = int(file.read())

    lags = [str(x) for x in list(range(num_lags))]

    col_names = ['_'.join(list(x)) for x in list(product(ch_names, freqs, lags))]

    return col_names


def format_gradcpt(subject, col_names, data_root):

    with open(data_root / Path(f'{subject}.pkl'), 'rb') as file:
        d = pickle.load(file)


    out = pd.DataFrame()

    for session in d:

        if 'run-001' not in d[session]:
            continue
        n = d[session]['run-001']['X'].shape[0]
        trs = d[session]['run-001']['trs']
        stem = np.column_stack([np.full(n, subject), np.full(n, session), trs])
        eeg = d[session]['run-001']['X']

        ses_ar = np.column_stack([stem, eeg])

        stem_names = ['subject', 'session', 'tr']
        full_col_names = stem_names + col_names

        ses_df = pd.DataFrame(ses_ar, columns = full_col_names)

        out = pd.concat([out, ses_df], axis=0)

    return out






if __name__ == '__main__':

    args = sys.argv[1:]

    if not args:
        os.chdir(here())
        sys.path.append(str(here()))
        data_root = Path('analysis/data/original')

    else:
        if len(args) > 2:
            raise RuntimeError('Usage: python extract_data.py path/to/data')

        data_root = Path(args[1])

    data_root = data_root / Path('../formatted/full')

    subjects = glob(str(data_root / Path('sub-*')))
    subjects = [Path(x).stem for x in subjects]
    subjects = sorted(subjects, key= lambda x: int(x.split('-')[1]))

    col_names = get_col_names(data_root)

    out = pd.DataFrame()

    for subject in subjects:

        print(f'Processing subject {subject} of {subjects[-1]}')
        sub_df = format_gradcpt(subject, col_names, data_root)
        out = pd.concat([out, sub_df], axis=0)


    dout = data_root / Path('../../correlation_data/gradcpt_eeg.feather')
    out = out.reset_index(drop=True)
    out.to_feather(dout)







