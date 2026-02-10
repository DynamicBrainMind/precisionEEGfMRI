import mne
from mne.time_frequency import tfr_array_morlet
from pathlib import Path
import numpy as np
import pandas as pd
from pyprojroot import here
import os
from scipy.signal import butter, filtfilt
os.chdir(here() / Path('analysis'))


subject = 'sub-006'
session = 'ses-001'
run = 2


data_root = Path('data/original')
data_leaf = Path(f'{subject}/{session}/eeg/{subject}_{session}_bld00{run}_eeg_Bergen_CWreg_filt_ICA_rej.set')

raw = mne.io.read_raw_eeglab(data_root / data_leaf)

# Decompose
freqs = np.array(range(1, 41))
njobs = int(os.cpu_count() - 1)
sfreq = raw.info['sfreq']

tf_full = tfr_array_morlet(
        raw.get_data()[np.newaxis, :, :],
        sfreq=sfreq,
        freqs=freqs,
        n_cycles=7,
        n_jobs=njobs,
        output='complex')[0]


power = np.abs(tf_full.data) ** 2
rel_power = power / power.sum(axis=2, keepdims=True)

# LPF
order = 4
cutoff = 0.25
do_filter = False
if do_filter:
    b, a  = butter(4, cutoff / (sfreq / 2), btype = 'low')
    filtered = filtfilt(b, a, rel_power)
else:
    filtered = rel_power

filtered = filtered.transpose(2, 0, 1)

time_cutoff = int(30 * sfreq)

# Get Oz channel
ch_names = raw.info['ch_names']
d = filtered[:time_cutoff, ch_names.index('Oz'), :]

out_leaf = Path('scripts/sandbox/timeseries/eeg_spectrum.csv')
pd.DataFrame(d, columns = freqs).to_csv(out_leaf, index = False)

