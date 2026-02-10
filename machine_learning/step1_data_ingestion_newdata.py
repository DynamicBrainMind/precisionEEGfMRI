"""
This script constructs the multimodal EEG–fMRI dataset used for predictive modeling.
For each subject, session, and block, preprocessed EEG recordings are transformed into
time–frequency representations aligned to fMRI sampling markers. EEG power estimates
from a fixed pre-marker window are extracted and concatenated across blocks and sessions.
In parallel, fMRI-derived network time series are loaded for each network of interest and
aligned in time with the EEG samples. The resulting EEG feature tensors, fMRI target
variables, and subject/session identifiers are concatenated across all participants and
saved for downstream machine-learning analyses.
"""

import os
import numpy as np
from utils import compute_pre_event_tfr_segments

data_path = os.path.join(os.getcwd(), 'eegfmri_data_11212025')

# -------------------------------------------------------------------------
# Define subnetworks
# -------------------------------------------------------------------------
network_names = ['DNa', 'DNb', 'dATNa', 'dATNb', 'FPCNa', 'FPCNb', 'SAL']

all_x = []
# One list of segments per network
all_y = {net: [] for net in network_names}

subject_ids = []
session_ids = []  # store the session number for each sample

subjects = [f"sub-{i:03d}" for i in [
    1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20,
    21, 23, 24, 25
]]

sessions = ['001', '002']
blocks = ['001', '002', '003', '004']

bin_size = 1.0
n_bins = 20

for sub in subjects:
    for ses in sessions:
        for bld in blocks:
            try:
                # -----------------------------------------------------------------
                # EEG file
                # -----------------------------------------------------------------
                eegfile = os.path.join(
                    data_path, sub, f'ses-{ses}', 'eeg',
                    f'{sub}_ses-{ses}_bld{bld}_eeg_Bergen_CWreg_filt_ICA_rej.set'
                )

                # -----------------------------------------------------------------
                # fMRI subnetwork time-series files
                # Assumes filenames: <NET>_ts_MSHBM_<sub>_bld<bld>_highpass.txt
                # e.g., DNa_ts_MSHBM_sub-001_bld001_highpass.txt
                # -----------------------------------------------------------------
                fmri_files = {
                    net: os.path.join(
                        data_path, sub, f'ses-{ses}', 'func', 'general',
                        f'{net}_ts_MSHBM_{sub}_bld{bld}_highpass.txt'
                    )
                    for net in network_names
                }

                # Check that all files exist
                all_files = [eegfile] + list(fmri_files.values())
                if not all(os.path.isfile(f) for f in all_files):
                    print(f'Skipping missing files for {sub}, ses-{ses}, bld{bld}')
                    continue

                # -----------------------------------------------------------------
                # Load EEG-derived segments
                # -----------------------------------------------------------------
                x = compute_pre_event_tfr_segments(
                    eegfile,
                    bin_size=bin_size,
                    n_bins=n_bins
                )

                # -----------------------------------------------------------------
                # Load fMRI subnetwork time series
                # -----------------------------------------------------------------
                y_dict = {
                    net: np.loadtxt(fmri_files[net]).reshape(-1, 1)
                    for net in network_names
                }

                # -----------------------------------------------------------------
                # Trim to the shortest length across EEG and all subnetworks
                # -----------------------------------------------------------------
                lengths = [x.shape[0]] + [y_dict[net].shape[0] for net in network_names]
                min_len = min(lengths)

                x = x[:min_len]
                for net in network_names:
                    y_dict[net] = y_dict[net][:min_len]

                # -----------------------------------------------------------------
                # Append to master lists
                # -----------------------------------------------------------------
                all_x.append(x)
                for net in network_names:
                    all_y[net].append(y_dict[net])

                subject_ids.extend([sub] * x.shape[0])
                session_ids.extend([ses] * x.shape[0])

            except Exception as e:
                print(f"Error in {sub}, ses-{ses}, bld{bld}: {e}")

# -------------------------------------------------------------------------
# Concatenate across all subjects/sessions/blocks
# -------------------------------------------------------------------------
X = np.concatenate(all_x, axis=0)
Y = {net: np.concatenate(all_y[net], axis=0) for net in network_names}

subject_ids = np.array(subject_ids)
session_ids = np.array(session_ids)

# Relative time bins (closest to event = 0s): [0, -1, -2, ..., -19]
bin_times = -np.arange(n_bins) * bin_size

# For convenience (if you still like having explicit variables)
Y_DNa   = Y['DNa']
Y_DNb   = Y['DNb']
Y_dATNa = Y['dATNa']
Y_dATNb = Y['dATNb']
Y_FPCNa = Y['FPCNa']
Y_FPCNb = Y['FPCNb']
Y_SAL   = Y['SAL']

print(
    f'Final shapes:\n'
    f'X       = {X.shape}\n'
    f'Y_DNa   = {Y_DNa.shape}\n'
    f'Y_DNb   = {Y_DNb.shape}\n'
    f'Y_dATNa = {Y_dATNa.shape}\n'
    f'Y_dATNb = {Y_dATNb.shape}\n'
    f'Y_FPCNa = {Y_FPCNa.shape}\n'
    f'Y_FPCNb = {Y_FPCNb.shape}\n'
    f'Y_SAL   = {Y_SAL.shape}\n'
    f'subject_ids = {subject_ids.shape}\n'
    f'session_ids = {session_ids.shape}'
)

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------
processed_dir = os.path.join(data_path, 'processed_data')
os.makedirs(processed_dir, exist_ok=True)
save_path = os.path.join(processed_dir, 'eeg_fmri_data_7nets.npz')

np.savez_compressed(
    save_path,
    X=X,
    Y_DNa=Y_DNa,
    Y_DNb=Y_DNb,
    Y_dATNa=Y_dATNa,
    Y_dATNb=Y_dATNb,
    Y_FPCNa=Y_FPCNa,
    Y_FPCNb=Y_FPCNb,
    Y_SAL=Y_SAL,
    subject_ids=subject_ids,
    session_ids=session_ids,
    bin_times=bin_times
)

print(f'Data saved to {save_path}')
