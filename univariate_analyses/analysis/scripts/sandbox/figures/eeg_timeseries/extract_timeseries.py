import mne
from mne.time_frequency import tfr_array_morlet
import numpy as np
import os
from pyprojroot import here
import sys
from pathlib import Path

args = sys.argv[1:]
if not args:
    os.chdir(here())
    data_root = Path('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/preprocessed')
else:
    data_root = Path(args[0])


script_root = here() / Path('analysis/scripts/sandbox/figures/eeg_timeseries')
subject = 'sub-020'
session = 'ses-001'
run = 'bld002'
fname = Path(f'{subject}_{session}_{run}_eeg_Bergen_CWreg_filt_ICA_rej.set')
data_path = data_root / Path(f'{subject}/{session}/eeg')

raw = mne.io.read_raw_eeglab(data_path / fname)
sfreq = raw.info['sfreq']
raw = raw.get_data()

time_range = 10
np.random.seed(42)
start = np.random.choice(range(raw.shape[1]))
#start = raw.shape[1] // 2
end = int(start + time_range * sfreq)


#freqs = np.logspace(np.log10(1), np.log10(40), 40)
freqs = np.linspace(1, 40, 40)
n_cycles = freqs / 2

power = tfr_array_morlet(
        raw[np.newaxis, :, :],
        freqs=freqs,
        n_cycles=n_cycles,
        sfreq=sfreq,
        n_jobs=os.cpu_count()-1,
        output='power')[0]

raw_s = raw.transpose()[np.arange(start, end)]
power = power.transpose(2, 0, 1)
power_s = power[np.arange(start, end)]

np.save(script_root / Path('timeseries.npy'), raw_s)
np.save(script_root / Path('power.npy'), power_s)
