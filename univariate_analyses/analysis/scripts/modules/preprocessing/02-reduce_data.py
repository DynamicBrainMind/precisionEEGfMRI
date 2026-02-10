import mne
from pyprojroot import here
import pandas as pd
import os
from glob import glob
import numpy as np
import pickle
from pathlib import Path
from itertools import product
import sys


class Reduce:
    """
    Reduce and summarize EEG feature arrays by averaging over frequency bands and 
    limiting lags.

    This class takes fully preprocessed EEG data (from the Reformat step), reshapes 
    and aggregates the data by canonical frequency bands (Delta, Theta, Alpha, Beta, 
    Gamma), and reduces the number of lags considered. The reduced data is saved 
    per subject in pickle files.

    Parameters
    ----------
    in_path : pathlib.Path
        Directory containing fully formatted EEG/fMRI data per subject.
    out_path : pathlib.Path
        Directory where reduced EEG data will be saved.
    overwrite : str
        Whether to overwrite existing data
    n_freqs : int, default=40
        Number of frequency bins in the original EEG feature arrays.
    n_chans : int, default=31
        Number of EEG channels.
    n_lags : int, default=11
        Number of time lags in the original EEG feature arrays.

    Attributes
    ----------
    n_freqs : int
        Number of frequency bins.
    n_chans : int
        Number of EEG channels.
    n_lags : int
        Number of time lags.
    in_path : pathlib.Path
        Input directory path.
    out_path : pathlib.Path
        Output directory path.
    subjects : list of str
        Sorted list of subject identifiers to process.
    ch_names : list of str
        List of EEG channel names loaded from `ch_names.txt`.

    Methods
    -------
    run()
        Iterates over all subjects and sessions, reduces EEG data, and saves the 
        reduced arrays per subject as pickle files.
    _reduce_eeg(X)
        Converts a full EEG array into a reduced array by averaging over canonical 
        frequency bands and limiting the number of lags.
    
    Notes
    -----
    - Assumes that input EEG data has already been processed with the Reformat step 
      (e.g., `01-reformat_data.py`).
    - Reduction involves:
        1. Reshaping the array to a DataFrame with separate columns for channel, 
           frequency, and lag.
        2. Binning frequencies into canonical EEG bands.
        3. Averaging power over frequency bands for lags 0–5.
        4. Pivoting back to a wide format and returning a numpy array.
    """

    def __init__(self, in_path, out_path, overwrite, n_freqs=40, n_chans=31, n_lags=11):

        # Initialize class level variables
        self.n_freqs = n_freqs
        self.n_chans = n_chans
        self.n_lags = n_lags
        self.in_path = in_path
        self.out_path = out_path
        self.overwrite = overwrite
        completed = glob(str(out_path / Path('sub*')))
        self.completed = [Path(x).stem for x in completed]

        # Make sure full formatted data exists
        if not os.path.exists(in_path):
            raise RuntimeError('Must run 01-reformat_data.py before this script')

        subjects = [Path(x).stem for x in glob(str(in_path / Path('*')))]
        self.subjects = sorted(subjects, key=lambda x: int(x.split('-')[1]))

        # Get channel names
        ch_names_path = here() / Path('analysis/scripts/modules/preprocessing/ch_names.txt')
        if not os.path.exists(ch_names_path):
            raise RuntimeError('Need to have a ch_names.txt file containing '
            'one channel name per line in the preprocessing dir')

        with open(ch_names_path, 'r') as file:
            self.ch_names = file.read().splitlines()

    def run(self):

        # Iterate over subjects
        for subject in self.subjects:
            print('\n')
            print(f'---Processing subject {subject} of {self.subjects[-1]}---\n')
            out = {}

            # Check for overwrite
            if subject in self.completed and not self.overwrite:
                continue

            dfile = self.in_path / Path(f'{subject}.pkl')
            with open(dfile, 'rb') as file:
                d = pickle.load(file)

            # Iterate over sessions
            for session in d:
                # Iterate over runs
                for run in d[session]:
                    X = self._reduce_eeg(d[session][run]['X'])
                    if session not in out:
                        out[session] = {run: {'X': X}}
                    else:
                        out[session][run] = {'X': X}
                    out[session][run]['y'] = d[session][run]['y']

            # Write per subject
            out_fpath = self.out_path / Path(f'{subject}.pkl')

            if not os.path.exists(out_fpath.parent):
                os.makedirs(out_fpath.parent)

            with open(out_fpath, 'wb') as file:
                pickle.dump(out, file)


    def _reduce_eeg(self, X):
        # Takes as input the full EEG array X
        # Returns truncated array of only 5 lags averaged over frequency
        # bands

        # Make a list of (fully wide) column names
        freqs = [f'freq{x}' for x in range(1, self.n_freqs+1)]
        lags = [f'lag{x}' for x in range(self.n_lags)]
        cols = list(product(self.ch_names, freqs, lags))
        cols = [f'{x[0]}_{x[1]}_{x[2]}' for x in cols]

        # --- Transform array to df and reshape ---
        X_df = pd.DataFrame(X, columns=cols)
        X_df['sample'] = range(1, X_df.shape[0]+1)
        X_df = X_df.melt(id_vars='sample', var_name='col', value_name='power')
        X_df[['channel', 'frequency', 'lag']] = X_df['col'].str.split('_', expand=True)
        X_df = X_df.drop(columns='col')
        X_df['frequency'] = X_df['frequency'].str.replace('freq', '').astype(int)
        X_df['lag'] = X_df['lag'].str.replace('lag', '').astype(int)

        # Band frequency into canonical freq bands
        bins = [0, 4, 8, 13, 30, 40]
        labels = ['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma']

        X_df['eeg_band'] = pd.cut(X_df['frequency'], bins=bins,
                                  labels=labels,
                                  right=True,
                                  include_lowest=True)

        # Average over frequency bands for 5 lags only
        X_df = X_df[X_df['lag'] <= 5].groupby(['sample', 'channel', 'eeg_band',
                                        'lag'])['power']\
                                                .mean()\
                                                .reset_index()

        # Spread to wide
        X_df['col'] = X_df[['channel', 'eeg_band', 'lag']].astype(str).agg('_'.join, axis=1)
        X_df = X_df.drop(columns=['channel', 'eeg_band', 'lag'])
        X_df = X_df.pivot(index='sample', columns='col', values='power').reset_index()
        X_df = X_df.drop(columns='sample')

        # Back to array
        X = X_df.to_numpy()

        return X


if __name__ == '__main__':

    args = sys.argv[1:]

    if not args:
        in_path = Path('analysis/data/formatted/full')
        out_path = Path('analysis/data/formatted/reduced')
        os.chdir(here())
        sys.path.append(str(here()))

    else:
        in_path = Path(args[0])
        out_path = in_path.parent / Path('reduced')

    valid_response = False
    while not valid_response:
        raw = input('\nDo you want to overwrite existing data? (y/n): ')
        print('\n')
        response = raw.strip().lower()
        if response in ['y', 'n']:
            valid_response = True
            overwrite = True if response == 'y' else False


    reduce = Reduce(in_path=in_path, out_path=out_path, overwrite=overwrite)
    reduce.run()
