import mne
import numpy as np
import matplotlib.pyplot as plt


def compute_pre_event_tfr_segments(eegfile, event_label='T  1', bin_size=2.0, n_bins=5,
                                   freqs=np.linspace(1, 40, 40), decim=10):
    """
    Compute multiple pre-event TFR segments of fixed size (bin_size) before each event.

    Parameters
    ----------
    eegfile : str
        Path to the EEGLAB .set file.
    event_label : str
        Event label to extract (default 'T  1').
    bin_size : float
        Duration in seconds of each time bin (default 2.0).
    n_bins : int
        Number of time bins to extract before each event (default 5).
    freqs : array-like
        Frequencies to compute TFR (default 1–40 Hz, 1 Hz step).
    decim : int
        Decimation factor for TFR computation (default 10).

    Returns
    -------
    tfr_segments : np.ndarray
        Array of shape (n_events, n_channels, n_freqs, n_bins), where each entry is
        the average power in a time bin before the event.
        Bin 0 = closest to the event, Bin n_bins-1 = farthest.
    """

    # Load EEG
    raw = mne.io.read_raw_eeglab(eegfile, preload=True)

    # Extract events
    events, event_id = mne.events_from_annotations(raw)
    if event_label not in event_id:
        raise ValueError(f"Event label '{event_label}' not found in the annotations.")
    
    idx_T = events[events[:, 2] == event_id[event_label], 0]
    idx_T_sec = raw.times[idx_T]

    # Compute TFR
    n_cycles = freqs / 2
    tfr = raw.compute_tfr(
        method='morlet',
        freqs=freqs,
        output='power',
        n_jobs=1,
        decim=decim,
        n_cycles=n_cycles,
        use_fft=True,
        zero_mean=True
    )

    sfreq = tfr.info['sfreq']
    times = tfr.times
    bin_samples = int(bin_size * sfreq)
    total_samples = bin_samples * n_bins

    # Store pre-event segments
    tfr_segments = []

    for ev_time in idx_T_sec:
        end_idx = np.argmin(np.abs(times - ev_time))
        start_idx = end_idx - total_samples

        # I want an error if out of bounds
        if start_idx < 0 or end_idx > len(times):
            raise ValueError(f"Event time {ev_time} is out of bounds for TFR data.")    

        bins = []
        for b in range(n_bins):
            b_start = start_idx + b * bin_samples
            b_end = b_start + bin_samples
            segment_bin = tfr.data[:, :, b_start:b_end].mean(axis=2)  # (n_channels, n_freqs)
            bins.append(segment_bin)

        # Stack bins along new axis → (n_channels, n_freqs, n_bins)
        bins = np.stack(bins, axis=-1)
        bins = bins[..., ::-1]  # <-- Flip bins: now 0 = closest to event
        tfr_segments.append(bins)

    # Final shape: (n_events, n_channels, n_freqs, n_bins)
    tfr_segments = np.array(tfr_segments)
    return tfr_segments




import numpy as np
import matplotlib.pyplot as plt
import mne
from functools import lru_cache

@lru_cache(maxsize=8)
def _load_info31(eegfile: str) -> mne.Info:
    """Load EEGLAB .set once and cache a 31-channel EEG-only Info."""
    raw = mne.io.read_raw_eeglab(eegfile, preload=False, verbose="ERROR")
    picks_all = mne.pick_types(raw.info, eeg=True, exclude="bads")
    if len(picks_all) < 31:
        raise ValueError(f"File has {len(picks_all)} EEG channels; need ≥31.")
    picks31 = picks_all[:31]  # deterministic first 31
    return mne.pick_info(raw.info, picks31)

@lru_cache(maxsize=8)
def _pos31(eegfile: str):
    """Compute and cache 2-D topomap coordinates for the 31-ch Info."""
    info31 = _load_info31(eegfile)
    picks = np.arange(31)
    # Try several options to be compatible across MNE versions
    try:
        from mne.viz.topomap import _get_pos
        pos, _ = _get_pos(info31, picks=picks, sphere="auto")
        return pos
    except Exception:
        pass
    try:
        from mne.channels.layout import _find_topomap_coords
        return _find_topomap_coords(info31, picks, ignore_overlap=True)
    except Exception:
        pass
    # Attach a standard montage then retry
    info_tmp = info31.copy()
    montage = mne.channels.make_standard_montage("standard_1020")
    info_tmp.set_montage(montage, match_case=False, on_missing="warn")
    try:
        from mne.viz.topomap import _get_pos as _gp2
        pos, _ = _gp2(info_tmp, picks=picks, sphere="auto")
        return pos
    except Exception:
        from mne.channels.layout import _find_topomap_coords as _ftc2
        return _ftc2(info_tmp, picks, ignore_overlap=True)

def plot_31ch_topomap(
    vals31,
    eegfile,
    title="Custom 31-channel values",
    *,
    ax=None,
    cmap=None,
    vmin=None,
    vmax=None,
    colorbar=False,
    sensors=True,
    contours=6,
    image_interp="linear",
):
    """
    Plot a 31-channel topomap compatible with grid/shared-colorbar scripts.
    Applies cmap/vmin/vmax AFTER drawing for compatibility with older MNE.
    """
    vals31 = np.asarray(vals31, dtype=float).ravel()
    if vals31.size != 31:
        raise ValueError(f"Expected 31 values, got {vals31.size}.")

    created = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(3.2, 3.2))
        created = True
    else:
        fig = ax.figure

    info31 = _load_info31(eegfile)

    # --- Try path 1: pass Info directly (preferred; supports sphere='auto')
    try:
        im, cn = mne.viz.plot_topomap(
            vals31, info31,
            axes=ax,
            sensors=sensors,
            contours=contours,
            image_interp=image_interp,
            outlines="head",
            sphere="auto",
            show=False,
        )
    except Exception:
        # --- Fallback path 2: pass 2-D positions and a numeric sphere
        pos = _pos31(eegfile)
        # Use a reasonable head radius in meters (avoid 'auto' with pos array)
        sphere_val = 0.095  # ~9.5 cm, close to MNE's EEG default
        im, cn = mne.viz.plot_topomap(
            vals31, pos,
            axes=ax,
            sensors=sensors,
            contours=contours,
            image_interp=image_interp,
            outlines="head",
            sphere=sphere_val,
            show=False,
        )

    # Now enforce color settings (older MNE may reject kwargs)
    try:
        if cmap is not None:
            im.set_cmap(cmap)
        if (vmin is not None) or (vmax is not None):
            lo = vmin if vmin is not None else np.nanmin(vals31)
            hi = vmax if vmax is not None else np.nanmax(vals31)
            if hi == lo:
                hi = lo + (1e-12 if hi == 0 else 1e-6 * abs(hi))
            im.set_clim(lo, hi)
    except Exception:
        pass

    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    if colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Never plt.show(); caller saves and closes.
    return fig
