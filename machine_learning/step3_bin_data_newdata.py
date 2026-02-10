"""
This script prepares EEG–fMRI data for predictive modeling by temporally and spectrally
aggregating EEG time–frequency representations. Specifically, it converts 1-s EEG
time–frequency power estimates into non-overlapping 2-s time bins over the 0–20 s
pre-BOLD window, averages power within canonical frequency bands (delta–gamma),
and outputs a compact feature tensor of shape (samples × channels × bands × time bins).
The resulting dataset is saved for subsequent XGBoost regression and model interpretation
analyses.
"""

import os
import numpy as np

# -------- Paths --------
data_root = os.path.join(os.getcwd(), 'eegfmri_data_11212025', 'processed_data')

# Use the 7-network file
load_path = os.path.join(data_root, 'eeg_fmri_data_7nets.npz')

# New output filename to avoid overwriting old 4-net file
save_path = os.path.join(
    data_root,
    'eeg_fmri_data_binned_2s_0to20s_canonicalbands_7nets.npz'
)

# -------- Load --------
data = np.load(load_path, allow_pickle=True)

X = data['X']                     # (n_samples, n_channels, n_freqs, n_bins)

# 7 subnetworks
Y_DNa   = data['Y_DNa']
Y_DNb   = data['Y_DNb']
Y_dATNa = data['Y_dATNa']
Y_dATNb = data['Y_dATNb']
Y_FPCNa = data['Y_FPCNa']
Y_FPCNb = data['Y_FPCNb']
Y_SAL   = data['Y_SAL']

subject_ids = data['subject_ids']
session_ids = data['session_ids']
bin_times   = data['bin_times']   # e.g., [0, -1, -2, ..., -19] (float or int)

n_samples, n_channels, n_freqs, n_bins = X.shape
epsilon = 1e-10

# If the frequency axis wasn't saved, assume the common 1–40 Hz grid used during extraction.
freqs = data['freqs'].astype(float) if 'freqs' in data.files else np.linspace(1, 40, n_freqs)

# -------- Log-transform (before aggregation across freqs) --------
# Log first; then average within bands (≈ geometric mean behavior).
X_log = np.log10(X + epsilon)

# -------- Define 2-s time bins within 0–20 s (pre-event) --------
# Pairs: [0,-1], [-2,-3], ..., [-18,-19]  -> 10 bins
def idx_of(t, times, atol=1e-6):
    hits = np.where(np.isclose(times, t, atol=atol))[0]
    if hits.size == 0:
        raise ValueError(f"Required time bin {t}s not found in bin_times.")
    return int(hits[0])

time_groups = []
time_bin_edges_sec = []
for start in range(0, 20, 2):          # 0,2,4,...,18
    a, b = start, start + 2            # numeric edges (seconds)
    # map to negative lag indices in bin_times
    t0 = 0 if a == 0 else - (a)        # left edge member index in 1-sec grid
    t1 = - (a + 1) if a != 0 else -1   # right member
    grp = [idx_of(t0, bin_times), idx_of(t1, bin_times)]
    time_groups.append(grp)
    time_bin_edges_sec.append((a, b))

time_bin_edges_sec = np.array(time_bin_edges_sec)                      # shape (10, 2)
time_bin_labels = [f"{a}\u2013{b}s" for (a, b) in time_bin_edges_sec]  # "0–2s", ..., "18–20s"

# Aggregate (mean) across the 1-s bins inside each 2-s bin
# Result shape after time rebin: (n_samples, n_channels, n_freqs, 10)
X_time_binned = np.stack(
    [X_log[..., grp].mean(axis=-1) for grp in time_groups],
    axis=-1
)

# -------- Canonical frequency bands --------
band_defs = {
    "delta": (1, 4),    # 1–4 Hz
    "theta": (4, 8),    # 4–8 Hz
    "alpha": (8, 13),   # 8–13 Hz
    "beta":  (13, 30),  # 13–30 Hz
    "gamma": (30, 40),  # 30–40 Hz
}
band_names = list(band_defs.keys())
band_bounds = np.array([band_defs[b] for b in band_names])

# Build index lists (left-inclusive, right-exclusive, last band right-inclusive)
band_indices = []
for i, (lo, hi) in enumerate(band_bounds):
    if i < len(band_bounds) - 1:
        idx = np.where((freqs >= lo) & (freqs < hi))[0]
    else:
        idx = np.where((freqs >= lo) & (freqs <= hi))[0]
    if idx.size == 0:
        raise ValueError(f"No frequency bins found for band {band_names[i]} [{lo},{hi}] Hz.")
    band_indices.append(idx)

# Aggregate (mean) within each band
# Output: (n_samples, n_channels, n_bands, 10)
X_binned = np.stack(
    [X_time_binned[:, :, idxs, :].mean(axis=2) for idxs in band_indices],
    axis=2
)

# -------- Clean printing --------
print(f"Original X shape: {X.shape}")
print(f"After time rebin (2s x 10): {X_time_binned.shape}")
print(f"Final X_binned shape (bands x time): {X_binned.shape}")
print("Bands: " + ", ".join(band_names))
print("Band bounds (Hz): " + ", ".join([f"[{lo},{hi}]" for (lo, hi) in band_bounds]))
print("Time bins: " + ", ".join(time_bin_labels))

# -------- Save --------
np.savez_compressed(
    save_path,
    X_binned=X_binned,                        # (n_samples, n_channels, n_bands, 10)
    bands=np.array(band_names),
    band_bounds=band_bounds,                  # [(lo,hi),...]
    time_bin_labels=np.array(time_bin_labels, dtype=object),
    time_bin_edges_sec=time_bin_edges_sec,    # numeric edges for programmatic use
    freqs=freqs,
    # 7 subnetworks
    Y_DNa=Y_DNa,
    Y_DNb=Y_DNb,
    Y_dATNa=Y_dATNa,
    Y_dATNb=Y_dATNb,
    Y_FPCNa=Y_FPCNa,
    Y_FPCNb=Y_FPCNb,
    Y_SAL=Y_SAL,
    subject_ids=subject_ids,
    session_ids=session_ids
)

print(f"Saved to {save_path}")
