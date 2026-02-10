import mne
import numpy as np
import pandas as pd
from pyprojroot import here
import os
from pathlib import Path

root = Path('analysis/scripts/sandbox/figures_general/cwreg_overview')
os.chdir(here() / root)

# Originally sub 2 ses 1
subject = 'sub-003'
session = 'ses-001'

dpaths = {'pre_cw': Path(f'data/{subject}_{session}_task-ExperienceSampling_run-002_eeg_Bergen.set'),
         'post_cw': Path(f'data/{subject}_{session}_task-ExperienceSampling_run-002_eeg_Bergen_CWreg.set')}

raw = mne.io.read_raw_eeglab(dpaths['pre_cw'])
channels = ['Oz', 'O1', 'O2']
channels += [x for x in raw.info['ch_names'] if x[0] == 'P']

for dtype in dpaths:
    dpath = dpaths[dtype]

    raw = mne.io.read_raw_eeglab(dpath)
    raw_picks = raw.pick_channels(channels)
    out = raw_picks.get_data().transpose(1, 0)
    out = pd.DataFrame(out, columns=channels)

    out.to_csv(f'data/{dtype}_timeseries.csv', index=False)

