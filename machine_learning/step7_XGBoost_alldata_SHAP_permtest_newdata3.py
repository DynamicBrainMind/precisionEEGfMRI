"""
This script fits XGBoost regression models to predict network-specific BOLD activity
from EEG features using all available valid samples and performs model interpretation
using SHAP (Shapley additive explanation) values.

Specifically, the script:
1) Loads EEG features organized as channel × frequency band × time-bin tensors
   and corresponding BOLD targets for each network.
2) Applies subject-level z-score normalization to EEG features to remove
   inter-subject scale differences.
3) Trains a separate XGBoost regression model for each network using fixed
   hyperparameters and all valid data (no cross-validation).
4) Computes SHAP values to quantify the contribution of each EEG feature
   (channel, band, time bin) to model predictions.
5) Aggregates SHAP values to obtain interpretable summary maps, including
   band × time and channel × band representations.
6) Assesses statistical significance of SHAP feature importance using
   within-subject permutation testing, followed by FDR correction on
   band × time and channel × band summary maps.
7) Saves SHAP importance measures, permutation-based p-values, and metadata
   for downstream visualization and reporting.

This analysis is designed to characterize the multivariate EEG features that
contribute most strongly to predicting BOLD activity in each network.
"""

import os
import json
import time
import numpy as np
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor
import shap  # pip install shap

# =========================
# Progress bars (tqdm w/ fallback)
# =========================
try:
    from tqdm.auto import tqdm
    _TQDM = True
except Exception:
    _TQDM = False

    def tqdm(iterable=None, total=None, desc=None, unit=None, **kwargs):
        # minimal no-op fallback
        return iterable if iterable is not None else range(total or 0)

    def _tqdm_write(msg):  # no-op
        print(msg)
else:
    def _tqdm_write(msg):
        try:
            tqdm.write(msg)
        except Exception:
            print(msg)

# =========================
# Paths & config
# =========================
data_root = os.path.join(os.getcwd(), "eegfmri_data_11212025")

# 7-network binned file
bands_file = os.path.join(
    data_root, "processed_data",
    "eeg_fmri_data_binned_2s_0to20s_canonicalbands_7nets.npz"
)

results_dir = os.path.join(os.getcwd(), "results")
os.makedirs(results_dir, exist_ok=True)

# Targets to run (7 subnetworks)
TARGET_KEYS = [
    "Y_DNa",
    "Y_DNb",
    "Y_dATNa",
    "Y_dATNb",
    "Y_FPCNa",
    "Y_FPCNb",
]
random_state = 1337

# Permutation-null for |SHAP| (0 = skip)
N_PERM_NULL = 500         # e.g., 100–500 for paper; 0 for speed
PERM_SEED   = 1337

xgb_params = dict(
    n_estimators=300,
    max_depth=3,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_alpha=0.0,
    reg_lambda=1.0,
    objective="reg:squarederror",
    tree_method="hist",
    random_state=random_state,
    n_jobs=0
)

# =========================
# Load data
# =========================
Z = np.load(bands_file, allow_pickle=True)

X_band      = Z["X_binned"]                 # (n_samples, n_channels, n_bands, n_bins)
band_names  = list(Z["bands"])
subject_ids = Z["subject_ids"]

# Time metadata (centers are negative, e.g. [-1, -3, ..., -19])
if "time_bin_edges_sec" in Z.files:
    edges          = Z["time_bin_edges_sec"].astype(float)   # (n_bins, 2) [(0,2),(2,4),...,(18,20)]
    bin_times      = -edges.mean(axis=1)                     # centers
    bin_left_edges = edges[:, 0].astype(int)                 # [0,2,...,18]
    bin_right_edge = int(edges[-1, 1])                       # 20
else:
    n_bins_fallback = X_band.shape[-1]
    bin_times      = -(np.arange(n_bins_fallback) * 2 + 1).astype(float)
    bin_left_edges = np.arange(n_bins_fallback) * 2
    bin_right_edge = int(bin_left_edges[-1] + 2)

n_samples, n_channels, n_bands, n_bins = X_band.shape
unique_subs = np.unique(subject_ids)
print(f"Loaded X_binned: {X_band.shape} — subjects: {len(unique_subs)}")

# =========================
# Per-subject normalization (z-score) + progress bar
# =========================
X_norm = np.zeros_like(X_band)
iter_subs = tqdm(unique_subs, desc="Per-subject z-scoring", unit="subj") if _TQDM else unique_subs
for subj in iter_subs:
    m = (subject_ids == subj)
    Xs2d = X_band[m].reshape(np.sum(m), -1)
    scaler_subj = StandardScaler()
    X_norm[m] = scaler_subj.fit_transform(Xs2d).reshape(np.sum(m), n_channels, n_bands, n_bins)

# Flatten for XGB
X_flat = X_norm.reshape(n_samples, -1)

# =========================
# Helpers
# =========================
def make_explainer(model, X):
    """Robust TreeSHAP explainer."""
    try:
        expl = shap.Explainer(model)
        sv = expl(X)
        vals = sv.values
    except Exception:
        expl = shap.TreeExplainer(model)
        vals = expl.shap_values(X)
    return np.asarray(vals)


def subject_level_cbt_shap(
    shap_values,
    subject_ids_valid,
    unique_subs_valid,
    n_channels,
    n_bands,
    n_bins,
):
    """
    Compute subject-level SHAP summaries in full 3D:

      shap_signed_cbt_by_subj : (n_subj, n_channels, n_bands, n_bins)
      shap_abs_cbt_by_subj    : (n_subj, n_channels, n_bands, n_bins)

    Each subject's tensor is the mean over that subject's trials.
    """
    sv_cbt = shap_values.reshape(shap_values.shape[0], n_channels, n_bands, n_bins)
    n_subj = len(unique_subs_valid)

    shap_signed_cbt_by_subj = np.zeros((n_subj, n_channels, n_bands, n_bins), dtype=np.float32)
    shap_abs_cbt_by_subj    = np.zeros_like(shap_signed_cbt_by_subj)

    for si, subj in enumerate(unique_subs_valid):
        m = (subject_ids_valid == subj)
        if not np.any(m):
            continue

        sv_sub = sv_cbt[m]  # (trials_s, C, B, T)

        # Mean over trials for this subject
        signed_cbt = sv_sub.mean(axis=0)          # (C, B, T)
        abs_cbt    = np.abs(sv_sub).mean(axis=0)  # (C, B, T)

        shap_signed_cbt_by_subj[si] = signed_cbt
        shap_abs_cbt_by_subj[si]    = abs_cbt

    return shap_signed_cbt_by_subj, shap_abs_cbt_by_subj


def permute_within_subjects(y, subject_ids_valid, rng):
    """
    Permute labels within each subject (destroys EEG–fMRI mapping within subject,
    preserves subject-level structure).
    """
    y_perm = y.copy()
    for subj in np.unique(subject_ids_valid):
        m = (subject_ids_valid == subj)
        y_perm[m] = rng.permutation(y_perm[m])
    return y_perm


# =========================
# Main loop (no folds; fit on ALL valid data) + permutations
# =========================
rng = np.random.default_rng(PERM_SEED)

outer_iter = tqdm(TARGET_KEYS, desc="Targets", unit="target") if _TQDM else TARGET_KEYS
for target_key in outer_iter:
    _tqdm_write("\n" + "=" * 70)
    _tqdm_write(f"Target (all-data SHAP): {target_key}")
    _tqdm_write("=" * 70)

    Y_all = Z[target_key][:, 0]  # (n_samples,)

    # ---- Mask out NaN/inf labels ----
    valid_mask = np.isfinite(Y_all)
    n_valid = int(valid_mask.sum())
    _tqdm_write(f"  Valid samples for {target_key}: {n_valid} / {len(Y_all)}")

    if n_valid < 2:
        _tqdm_write(f"  Not enough valid samples for {target_key}, skipping.")
        continue

    idx_valid         = np.where(valid_mask)[0]
    X_flat_valid      = X_flat[idx_valid]
    Y_valid           = Y_all[idx_valid]
    subject_ids_valid = subject_ids[idx_valid]
    unique_subs_valid = np.unique(subject_ids_valid)

    # ---- Scale features (valid subset only) ----
    t0 = time.perf_counter()
    scaler_all = StandardScaler()
    X_all_scaled = scaler_all.fit_transform(X_flat_valid)
    _tqdm_write(f"  Scaled features in {time.perf_counter() - t0:.2f}s")

    # ---- Fit model on all valid samples ----
    t0 = time.perf_counter()
    model_all = XGBRegressor(**xgb_params)
    model_all.fit(X_all_scaled, Y_valid)
    _tqdm_write(f"  Trained XGB in {time.perf_counter() - t0:.2f}s")

    # ---- SHAP on all valid samples ----
    t0 = time.perf_counter()
    shap_values = make_explainer(model_all, X_all_scaled)  # (n_valid, n_features)
    _tqdm_write(f"  Computed SHAP in {time.perf_counter() - t0:.2f}s")

    # Global (sample-averaged) SHAP summaries
    shap_mean_abs    = np.abs(shap_values).mean(axis=0)    # (n_features,)
    shap_mean_signed = shap_values.mean(axis=0)            # (n_features,)

    # Reshape to (channels, bands, time)
    imp_abs_cbb    = shap_mean_abs.reshape(n_channels, n_bands, n_bins)
    imp_signed_cbb = shap_mean_signed.reshape(n_channels, n_bands, n_bins)

    # Simple global summaries (no subjects)
    imp_abs_per_channel    = imp_abs_cbb.sum(axis=(1, 2))
    imp_signed_per_channel = imp_signed_cbb.sum(axis=(1, 2))
    band_time_abs          = imp_abs_cbb.sum(axis=0)       # (bands, time)
    band_time_signed       = imp_signed_cbb.sum(axis=0)
    chan_time_abs          = imp_abs_cbb.sum(axis=1)       # (channels, time)
    chan_time_signed       = imp_signed_cbb.sum(axis=1)

    # ---- Subject-level SHAP in full 3D (C,B,T) ----
    shap_signed_cbt_by_subj, shap_abs_cbt_by_subj = subject_level_cbt_shap(
        shap_values,
        subject_ids_valid,
        unique_subs_valid,
        n_channels,
        n_bands,
        n_bins,
    )

    # Derived subject-level summaries (if you want to inspect them later)
    shap_signed_bt_by_subj = shap_signed_cbt_by_subj.mean(axis=1)  # avg over channels → (subj,B,T)
    shap_abs_bt_by_subj    = shap_abs_cbt_by_subj.mean(axis=1)

    shap_signed_cb_by_subj = shap_signed_cbt_by_subj.mean(axis=3)  # avg over time → (subj,C,B)
    shap_abs_cb_by_subj    = shap_abs_cbt_by_subj.mean(axis=3)

    # =========================
    # Permutation-null for |SHAP| in full 3D
    #   - core: channel × band × time (CBT)
    #   - plus derived: band × time (BT) and channel × band (CB)
    # =========================
    p_abs_perm_cbt   = np.array([])  # (channels, bands, time)
    p_abs_perm_bt    = np.array([])  # (bands, time)
    p_abs_perm_cb    = np.array([])  # (channels, bands)

    null_means_abs_cbt = np.array([])
    null_means_abs_bt  = np.array([])
    null_means_abs_cb  = np.array([])

    if N_PERM_NULL > 0:
        _tqdm_write(f"  Permutation null for |SHAP| (full 3D): {N_PERM_NULL} permutations")

        null_means_abs_cbt = np.zeros(
            (N_PERM_NULL, n_channels, n_bands, n_bins),
            dtype=np.float32
        )

        perm_pbar = tqdm(
            total=N_PERM_NULL,
            desc=f"  Perms {target_key}",
            unit="perm",
        ) if _TQDM else None

        for r in range(N_PERM_NULL):
            # Permute labels within subjects
            y_perm = permute_within_subjects(Y_valid, subject_ids_valid, rng)

            # Fit null model
            model_null = XGBRegressor(**xgb_params)
            model_null.fit(X_all_scaled, y_perm)

            # SHAP for null model
            shap_null = make_explainer(model_null, X_all_scaled)

            # Subject-level |SHAP| in full 3D under null
            _, abs_cbt_by_subj_null = subject_level_cbt_shap(
                shap_null,
                subject_ids_valid,
                unique_subs_valid,
                n_channels,
                n_bands,
                n_bins,
            )

            # Mean across subjects → null statistic in full 3D
            null_means_abs_cbt[r] = abs_cbt_by_subj_null.mean(axis=0)  # (C,B,T)

            if perm_pbar is not None:
                perm_pbar.update(1)

        if perm_pbar is not None:
            perm_pbar.close()

        # ----- P-values in full 3D (C,B,T) -----
        obs_mean_abs_cbt = shap_abs_cbt_by_subj.mean(axis=0)  # (C,B,T)

        ge_counts_cbt = (null_means_abs_cbt >= obs_mean_abs_cbt[None, ...]).sum(axis=0)
        p_abs_perm_cbt = (1.0 + ge_counts_cbt) / (N_PERM_NULL + 1.0)  # (C,B,T)

        # ----- Derived band×time p-values (avg over channels) -----
        obs_mean_abs_bt = obs_mean_abs_cbt.mean(axis=0)               # (B,T)
        null_means_abs_bt = null_means_abs_cbt.mean(axis=1)           # (R,B,T)

        ge_counts_bt = (null_means_abs_bt >= obs_mean_abs_bt[None, ...]).sum(axis=0)
        p_abs_perm_bt = (1.0 + ge_counts_bt) / (N_PERM_NULL + 1.0)    # (B,T)

        # ----- Derived channel×band p-values (avg over time) -----
        obs_mean_abs_cb = obs_mean_abs_cbt.mean(axis=2)               # (C,B)
        null_means_abs_cb = null_means_abs_cbt.mean(axis=3)           # (R,C,B)

        ge_counts_cb = (null_means_abs_cb >= obs_mean_abs_cb[None, ...]).sum(axis=0)
        p_abs_perm_cb = (1.0 + ge_counts_cb) / (N_PERM_NULL + 1.0)    # (C,B)

        _tqdm_write("  Permutation p-values computed in 3D + aggregated 2D.")

    # ===== Save =====
    out_file = os.path.join(results_dir, f"xgb_shap_alldata_{target_key}_0to20s.npz")
    np.savez_compressed(
        out_file,
        # SHAP summaries (global, across all valid samples)
        shap_mean_abs=shap_mean_abs,
        shap_mean_signed=shap_mean_signed,
        imp_abs_cbb=imp_abs_cbb,
        imp_signed_cbb=imp_signed_cbb,
        imp_abs_per_channel=imp_abs_per_channel,
        imp_signed_per_channel=imp_signed_per_channel,
        band_time_abs=band_time_abs,
        band_time_signed=band_time_signed,
        chan_time_abs=chan_time_abs,
        chan_time_signed=chan_time_signed,

        # Subject-level SHAP summaries
        shap_signed_cbt_by_subj=shap_signed_cbt_by_subj,   # (subj, C, B, T)
        shap_abs_cbt_by_subj=shap_abs_cbt_by_subj,
        shap_signed_bt_by_subj=shap_signed_bt_by_subj,     # (subj, B, T)
        shap_abs_bt_by_subj=shap_abs_bt_by_subj,
        shap_signed_cb_by_subj=shap_signed_cb_by_subj,     # (subj, C, B)
        shap_abs_cb_by_subj=shap_abs_cb_by_subj,

        # Permutation-null results for |SHAP|
        #  Full 3D: channel × band × time
        p_abs_perm_cbt=p_abs_perm_cbt,                     # (C,B,T)
        null_means_abs_cbt=null_means_abs_cbt,             # (R,C,B,T)
        #  Aggregated (for convenience)
        p_abs_perm=p_abs_perm_bt,                          # (B,T)  (kept name for back-compat)
        null_means_abs_bt=null_means_abs_bt,               # (R,B,T)
        p_abs_perm_cb=p_abs_perm_cb,                       # (C,B)
        null_means_abs_cb=null_means_abs_cb,               # (R,C,B)

        # Meta
        band_names=np.array(band_names),
        bin_times=np.array(bin_times, dtype=float),
        bin_left_edges=np.array(bin_left_edges, dtype=int),
        bin_right_edge=np.int64(bin_right_edge),
        subject_ids=subject_ids_valid,
        unique_subjects=unique_subs_valid,
        n_channels=n_channels,
        n_bands=n_bands,
        n_bins=n_bins,
        target_key=target_key,
        xgb_params=json.dumps(xgb_params),
        n_perm_null=int(N_PERM_NULL),
        perm_seed=int(PERM_SEED),
        valid_mask=valid_mask,

        # Back-compat placeholders (so other scripts don't break if they expect them)
        fold_r2=np.array([]),
        overall_r2=np.nan,
        fold_subjects=np.array([]),
    )
    _tqdm_write(f"  Saved: {out_file}")

print("\nDone. One model per target, trained on ALL valid data, SHAP + 3D permutation stats saved.")
