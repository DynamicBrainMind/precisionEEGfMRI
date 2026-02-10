"""
This script visualizes SHAP-based feature importance results from XGBoost models
trained to predict network-specific BOLD activity from EEG features.

Specifically, for each target network, the script:
1) Plots band × time trajectories of mean absolute SHAP values (|SHAP|),
   highlighting time–frequency bins that are statistically significant based
   on permutation testing with FDR correction.
2) Generates channel × band topographic maps of time-averaged |SHAP| values.
3) Produces binary significance topomaps indicating channels that contribute
   significantly within each frequency band, based on permutation-derived
   p-values with per-band FDR correction.

Permutation testing is applied only to the band × time and channel × band
summaries of |SHAP| values (not the full channel × band × time tensor).
All figures are saved in publication-ready PNG and PDF formats.
"""


"""
shap_significance_time_series_and_topos.py

Plots:
  1) Band × time |SHAP| as time series (bands as lines, significant points highlighted)
  2) Channel × band |SHAP| as topomaps (time-averaged), plus binary significance topomaps
"""

import os
import numpy as np

# ---- Headless backend ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
plt.ioff()

import seaborn as sns

from utils import plot_31ch_topomap

# =========================
# Paths & config
# =========================
data_root   = os.path.join(os.getcwd(), "eegfmri_data_11212025")
results_dir = os.path.join(os.getcwd(), "results")
fig_root    = os.path.join(results_dir, "shap_significance_plots")
os.makedirs(fig_root, exist_ok=True)

# EEG file just for montage positions
example_sub = "sub-001"
example_ses = "001"
example_bld = "001"
eegfile = os.path.join(
    data_root, example_sub, f"ses-{example_ses}", "eeg",
    f"{example_sub}_ses-{example_ses}_bld{example_bld}_eeg_Bergen_CWreg_filt_ICA_rej.set"
)

# Targets (6 subnetworks, SAL removed)
TARGET_KEYS = [
    "Y_DNa",
    "Y_DNb",
    "Y_dATNa",
    "Y_dATNb",
    "Y_FPCNa",
    "Y_FPCNb",
]

# Human-readable names (keep dATNa / dATNb)
NAME_MAP = {
    "Y_DNa":   "DNa",
    "Y_DNb":   "DNb",
    "Y_dATNa": "dATNa",
    "Y_dATNb": "dATNb",
    "Y_FPCNa": "FPCNa",
    "Y_FPCNb": "FPCNb",
}

ALPHA   = 0.05   # FDR level
USE_FDR = True   # if False: use raw p <= ALPHA

sns.set(context="paper", style="whitegrid", font_scale=1.2)
TIME_SERIES_FIGSIZE = (7.4, 3.8)
TOPO_FIGHEIGHT      = 3.9
DPI = 300

ABS_CMAP      = "viridis"
SIG_TOPO_CMAP = "Greys"

# Subnetwork colors in R² script (for reference only, not used here):
# DNa:  #BB3738, DNb: #FE9386, FPCNa: #4F82B5, FPCNb: #A5DAF4, dATNa: #37773E, dATNb: #CEE0A4
# Use separate palette for frequency bands here.
BAND_COLOR_PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]
LINE_MARKERS = ["o", "s", "^", "D", "P", "X", "*", "v"]


# =========================
# Helpers
# =========================
def ensure_dir(p):
    os.makedirs(p, exist_ok=True)


def clean_band_names(band_arr):
    out = []
    for b in band_arr:
        if isinstance(b, (bytes, bytearray)):
            out.append(b.decode("utf-8"))
        else:
            out.append(str(b))
    return out


def bh_fdr(pvals, alpha=0.05):
    """
    Benjamini–Hochberg FDR for arbitrary-shaped p-value array.
    Returns:
        qvals : same shape as pvals
        sig   : boolean mask of same shape (qvals <= alpha)
    """
    p = np.asarray(pvals, dtype=float)
    shape = p.shape
    p_flat = p.ravel()

    # Treat NaNs as 1 (non-significant)
    nan_mask = ~np.isfinite(p_flat)
    p_flat = p_flat.copy()
    p_flat[nan_mask] = 1.0

    n = p_flat.size
    order = np.argsort(p_flat)
    ranked = np.arange(1, n + 1)
    p_sorted = p_flat[order]

    q_sorted = np.minimum.accumulate((p_sorted * n / ranked)[::-1])[::-1]
    q_flat = np.empty_like(p_flat)
    q_flat[order] = np.clip(q_sorted, 0, 1)

    qvals = q_flat.reshape(shape)
    sig = qvals <= alpha
    sig.reshape(-1)[nan_mask] = False
    return qvals, sig


def bh_fdr_1d(pvals, alpha=0.05):
    """
    Benjamini–Hochberg FDR for a 1D array of p-values.
    Used for per-band FDR across channels.
    Returns:
        qvals : same shape as pvals
        sig   : boolean mask (qvals <= alpha)
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size

    nan_mask = ~np.isfinite(p)
    p = p.copy()
    p[nan_mask] = 1.0

    order = np.argsort(p)
    ranked = np.arange(1, n + 1)
    p_sorted = p[order]

    q_sorted = np.minimum.accumulate((p_sorted * n / ranked)[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.clip(q_sorted, 0, 1)

    sig = q <= alpha
    sig[nan_mask] = False
    return q, sig


def time_from_bin_times(bin_times):
    """
    bin_times are centers in negative seconds (e.g., [-1,-3,...,-19]).
    Convert to positive seconds before event and sort ascending.
    """
    t_sec = -np.asarray(bin_times, dtype=float)
    order = np.argsort(t_sec)
    return t_sec[order], order


def pick_band_color(bname, i):
    # Dedicated band palette, distinct from subnetwork colors
    return BAND_COLOR_PALETTE[i % len(BAND_COLOR_PALETTE)]


def save_fig(fig, out_base, dpi=300):
    ensure_dir(os.path.dirname(out_base))
    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=dpi, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", dpi=dpi, bbox_inches="tight")
    print(f"Saved: {out_base}.png / .pdf")
    plt.close(fig)


def call_topomap_on_ax(ax, vals, title, cmap, vmin=None, vmax=None):
    kwargs = {"title": title, "ax": ax, "colorbar": False, "cmap": cmap}
    if vmin is not None:
        kwargs["vmin"] = vmin
    if vmax is not None:
        kwargs["vmax"] = vmax

    try:
        ret = plot_31ch_topomap(vals, eegfile, **kwargs)
    except TypeError:
        try:
            kwargs.pop("ax", None)
            ret = plot_31ch_topomap(vals, eegfile, **kwargs)
        except TypeError:
            ret = plot_31ch_topomap(vals, eegfile, title=title)
    return ret


def add_shared_cbar(fig, cax, cmap_name, vmin, vmax, label=""):
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap_name))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    if label:
        cb.set_label(label)
    return cb


# =========================
# First pass: global axes for time series
# =========================
global_ymax = 0.0
global_max_t = 0.0

for target_key in TARGET_KEYS:
    npz_path = os.path.join(results_dir, f"xgb_shap_alldata_{target_key}_0to20s.npz")
    if not os.path.exists(npz_path):
        continue

    D = np.load(npz_path, allow_pickle=True)

    if ("band_time_abs" not in D.files) or ("band_names" not in D.files):
        continue

    band_time_abs = D["band_time_abs"]  # (n_bands, n_bins)

    # Update global max |SHAP|
    if np.isfinite(band_time_abs).any():
        local_max = float(np.nanmax(band_time_abs))
        global_ymax = max(global_ymax, local_max)

    # Time axis for this target
    if "bin_times" in D.files:
        bin_times = D["bin_times"]
    else:
        n_bins_tmp = band_time_abs.shape[1]
        bin_times = -(np.arange(n_bins_tmp) * 2 + 1).astype(float)

    t_sec, _ = time_from_bin_times(bin_times)
    if t_sec.size:
        local_max_t = float(np.nanmax(t_sec))
        global_max_t = max(global_max_t, local_max_t)

# Fallbacks if something went wrong
if not np.isfinite(global_ymax) or global_ymax <= 0:
    global_ymax = 1.0
if not np.isfinite(global_max_t) or global_max_t <= 0:
    global_max_t = 1.0

# Add a bit of padding so points aren’t clipped
Y_PAD_FRAC = 0.05
GLOBAL_YMIN = 0.0
GLOBAL_YMAX = global_ymax + Y_PAD_FRAC * global_ymax

X_PAD = 0.5  # extend beyond last bin center (e.g., >19 s)
GLOBAL_TMAX = global_max_t + X_PAD


# =========================
# Main loop
# =========================
for target_key in TARGET_KEYS:
    npz_path = os.path.join(results_dir, f"xgb_shap_alldata_{target_key}_0to20s.npz")
    if not os.path.exists(npz_path):
        print(f"[WARN] Missing NPZ for {target_key}: {npz_path}")
        continue

    D = np.load(npz_path, allow_pickle=True)

    # Required base fields
    required = ["band_names", "band_time_abs", "imp_abs_cbb", "p_abs_perm_cb"]
    missing = [k for k in required if k not in D.files]
    if missing:
        print(f"[WARN] NPZ for {target_key} is missing: {missing}. Skipping.")
        continue

    band_names    = clean_band_names(D["band_names"])
    band_time_abs = D["band_time_abs"]        # (n_bands, n_bins)
    imp_abs_cbb   = D["imp_abs_cbb"]          # (n_channels, n_bands, n_bins)
    p_cb          = D["p_abs_perm_cb"]        # (n_channels, n_bands)

    # Band×time p-values: accept either new or old key name
    if "p_abs_perm_bt" in D.files:
        p_bt = D["p_abs_perm_bt"]
    elif "p_abs_perm" in D.files:
        p_bt = D["p_abs_perm"]                # legacy name from earlier script
    else:
        print(f"[WARN] NPZ for {target_key} has no band×time permutation p-values. Skipping.")
        continue

    # Time axis
    if "bin_times" in D.files:
        bin_times = D["bin_times"]
    else:
        n_bins_tmp = band_time_abs.shape[1]
        bin_times = -(np.arange(n_bins_tmp) * 2 + 1).astype(float)

    # Convert negative centers → positive seconds before event (ascending)
    t_sec, order_t = time_from_bin_times(bin_times)
    band_time_abs_ord = band_time_abs[:, order_t]
    p_bt_ord          = p_bt[:, order_t]

    n_bands, n_bins = band_time_abs_ord.shape
    n_channels      = imp_abs_cbb.shape[0]

    # -------------------------
    # FDR significance
    # -------------------------
    if USE_FDR:
        # Global FDR across all band×time cells
        q_bt, sig_bt = bh_fdr(p_bt_ord, alpha=ALPHA)   # (n_bands, n_bins)

        # Per-band FDR across channels for channel×band (31 tests per band)
        q_cb = np.zeros_like(p_cb)
        sig_cb = np.zeros_like(p_cb, dtype=bool)
        for b in range(p_cb.shape[1]):
            q_b, sig_b = bh_fdr_1d(p_cb[:, b], alpha=ALPHA)
            q_cb[:, b]   = q_b
            sig_cb[:, b] = sig_b
    else:
        sig_bt = p_bt_ord <= ALPHA
        sig_cb = p_cb     <= ALPHA

    # ====== CHECK SIGNIFICANCE FOR dATNb IN TIME SERIES ======
    if target_key == "Y_dATNb":
        any_sig = bool(sig_bt.any())
        print(f"[dATNb] Any significant band×time cells (FDR {ALPHA})? {any_sig}")

        if any_sig:
            for bi, bname in enumerate(band_names):
                sig_idx = np.where(sig_bt[bi])[0]   # indices of significant time bins for this band
                if sig_idx.size > 0:
                    t_sig = t_sec[sig_idx]
                    print(
                        f"[dATNb] Band {bname}: {sig_idx.size} significant bins "
                        f"at times (s): {np.round(t_sig, 2)}"
                    )


    # =========================
    # 1) Band × time time-series plot
    # =========================
    data_bt = band_time_abs_ord  # (bands, time)

    # Shared x ticks based on *un-padded* global_max_t
    xticks = np.arange(1, int(np.floor(global_max_t)) + 1, 2)

    fig_ts, ax_ts = plt.subplots(figsize=TIME_SERIES_FIGSIZE)
    marker_every = max(1, len(t_sec) // 10)

    for bi, bname in enumerate(band_names):
        color  = pick_band_color(bname, bi)
        marker = LINE_MARKERS[bi % len(LINE_MARKERS)]
        yvals  = data_bt[bi]

        ax_ts.plot(
            t_sec, yvals,
            lw=2.0, label=bname,
            color=color, marker=marker,
            markevery=marker_every, ms=5,
            zorder=2,
        )

        sig_mask_band = sig_bt[bi]  # (n_bins,) for this band
        if np.any(sig_mask_band):
            t_sig = t_sec[sig_mask_band]
            y_sig = yvals[sig_mask_band]
            ax_ts.scatter(
                t_sig, y_sig,
                s=80,
                facecolors="none",
                edgecolors="k",
                linewidths=1.5,
                zorder=5,
            )

    # Use global consistent axes with padding
    ax_ts.set_xlim(0.0, GLOBAL_TMAX)
    ax_ts.set_xticks(xticks)
    ax_ts.set_xticklabels([f"{int(x)}" for x in xticks])

    ax_ts.set_ylim(GLOBAL_YMIN, GLOBAL_YMAX)
    ax_ts.set_xlabel("Time before event (s)")
    ax_ts.set_ylabel("|SHAP| (channel-summed)")
    ax_ts.set_title(f"{NAME_MAP.get(target_key, target_key)} — |SHAP| over time by band")
    ax_ts.grid(True, alpha=0.25, linestyle="--", linewidth=0.8)

    ax_ts.legend(frameon=False, ncol=len(band_names),
                 loc="upper center", bbox_to_anchor=(0.5, -0.22))
    sns.despine(fig=fig_ts)

    out_ts = os.path.join(
        fig_root,
        f"{target_key}_absSHAP_time_series_bandMean_withSig"
    )
    save_fig(fig_ts, out_ts)

    # =========================
    # 2) Channel × band topomaps
    # =========================
    # Time-averaged |SHAP| per channel×band
    chan_band_abs = imp_abs_cbb.mean(axis=2)  # (channels, bands)

    vmax_abs = np.nanmax(chan_band_abs) if np.isfinite(chan_band_abs).any() else 1.0
    vmax_abs = 1e-12 if vmax_abs == 0 else vmax_abs
    vmin_abs = 0.0

    # ---- (a) unthresholded |SHAP| topomaps ----
    fig_mag = plt.figure(figsize=(3.4 * n_bands, TOPO_FIGHEIGHT))
    gs_mag = gridspec.GridSpec(2, n_bands, height_ratios=[20, 1], hspace=0.25, wspace=0.12)
    axes_mag = [fig_mag.add_subplot(gs_mag[0, i]) for i in range(n_bands)]
    cax_mag = fig_mag.add_subplot(gs_mag[1, :])

    for bi, (ax, bname) in enumerate(zip(axes_mag, band_names)):
        vals = chan_band_abs[:, bi]
        call_topomap_on_ax(
            ax, vals,
            title=f"{bname}",
            cmap=ABS_CMAP,
            vmin=vmin_abs,
            vmax=vmax_abs,
        )
        ax.set_xticks([])
        ax.set_yticks([])

    fig_mag.suptitle(
        f"{NAME_MAP.get(target_key, target_key)} — |SHAP| (time-avg) per band",
        y=0.98, fontsize=12
    )
    add_shared_cbar(
        fig_mag, cax_mag,
        cmap_name=ABS_CMAP,
        vmin=vmin_abs, vmax=vmax_abs,
        label="|SHAP| (time-avg)"
    )
    out_mag = os.path.join(
        fig_root,
        f"{target_key}_absSHAP_channelBand_topomap_grid"
    )
    save_fig(fig_mag, out_mag)

    # ---- (b) binary significance topomaps (0/1), FDR per band ----
    sig_cb_float = sig_cb.astype(float)  # (channels, bands), entries 0 or 1

    fig_sig = plt.figure(figsize=(3.4 * n_bands, TOPO_FIGHEIGHT))
    gs_sig = gridspec.GridSpec(2, n_bands, height_ratios=[20, 1], hspace=0.25, wspace=0.12)
    axes_sig = [fig_sig.add_subplot(gs_sig[0, i]) for i in range(n_bands)]
    cax_sig = fig_sig.add_subplot(gs_sig[1, :])

    for bi, (ax, bname) in enumerate(zip(axes_sig, band_names)):
        vals = sig_cb_float[:, bi]  # 0 or 1 per channel
        call_topomap_on_ax(
            ax, vals,
            title=f"{bname}",
            cmap=SIG_TOPO_CMAP,
            vmin=0.0, vmax=1.0,
        )
        ax.set_xticks([])
        ax.set_yticks([])

    fig_sig.suptitle(
        f"{NAME_MAP.get(target_key, target_key)} — significance (0–1) per band (FDR per band)",
        y=0.98, fontsize=12
    )
    add_shared_cbar(
        fig_sig, cax_sig,
        cmap_name=SIG_TOPO_CMAP,
        vmin=0.0, vmax=1.0,
        label="Significant (1) / Not (0)"
    )
    out_sig = os.path.join(
        fig_root,
        f"{target_key}_sig_channelBand_topomap_grid_binary"
    )
    save_fig(fig_sig, out_sig)

print("Done: time-series (band×time, global axes) and topomap (channel×band, FDR per band) significance plots.")
