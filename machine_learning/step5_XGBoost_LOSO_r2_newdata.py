"""
This script evaluates the predictive performance of EEG-derived features for
modeling fMRI BOLD activity using XGBoost regression.

EEG features consist of channel-level band-limited power values organized as
(channel × frequency band × time-bin) tensors derived from 2-second bins in the
20 seconds preceding each BOLD sample. Features are normalized within subjects
prior to modeling.

Prediction performance is assessed separately for each fMRI network using a
leave-one-subject-out (LOSO) cross-validation framework. In each fold, all data
from one subject are held out for testing while the model is trained on data
from the remaining subjects.

Model performance is quantified using the coefficient of determination (R²),
computed for each held-out subject. The script saves per-subject R² values,
out-of-fold predictions, and summary statistics for downstream statistical
analysis and visualization.
"""

import os
import json
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from xgboost import XGBRegressor

# =========================
# Paths & config
# =========================
data_root = os.path.join(os.getcwd(), "eegfmri_data_11212025")

# Use the 7-network binned file
bands_file = os.path.join(
    data_root, "processed_data",
    "eeg_fmri_data_binned_2s_0to20s_canonicalbands_7nets.npz"
)

results_dir = os.path.join(os.getcwd(), "results")
os.makedirs(results_dir, exist_ok=True)

# Run all seven targets
TARGET_KEYS = [
    "Y_DNa",
    "Y_DNb",
    "Y_dATNa",
    "Y_dATNb",
    "Y_FPCNa",
    "Y_FPCNb",
    "Y_SAL",
]
random_state = 1337

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
# Load X and metadata (shared across targets)
# =========================
Z = np.load(bands_file, allow_pickle=True)

X_band = Z["X_binned"]                               # (n_samples, n_channels, n_bands, n_bins)
band_names = list(Z["bands"])                        # ['delta','theta','alpha','beta','gamma']
subject_ids = Z["subject_ids"]

# Bin meta (edges + centers; centers are negative lag for naming)
if "time_bin_edges_sec" in Z.files:
    edges = Z["time_bin_edges_sec"].astype(float)    # (10,2): [(0,2),(2,4),...,(18,20)]
    bin_times = -edges.mean(axis=1)                  # centers: [-1,-3,...,-19]
    bin_left_edges = edges[:, 0].astype(int)         # [0,2,...,18]
    bin_right_edge = int(edges[-1, 1])               # 20
else:
    n_bins_fallback = X_band.shape[-1]
    bin_times = -(np.arange(n_bins_fallback) * 2 + 1).astype(float)
    bin_left_edges = np.arange(n_bins_fallback) * 2
    bin_right_edge = int(bin_left_edges[-1] + 2)

n_samples, n_channels, n_bands, n_bins = X_band.shape
unique_subs = np.unique(subject_ids)
print(f"Loaded X_binned: {X_band.shape} — subjects: {len(unique_subs)}")

# =========================
# Subject-wise normalization (shared, independent of target)
# =========================
X_norm = np.zeros_like(X_band)
for subj in unique_subs:
    m = (subject_ids == subj)
    Xs2d = X_band[m].reshape(np.sum(m), -1)
    scaler_subj = StandardScaler()
    X_norm[m] = scaler_subj.fit_transform(Xs2d).reshape(
        np.sum(m), n_channels, n_bands, n_bins
    )

X_flat = X_norm.reshape(n_samples, -1)

# =========================
# Loop over targets (R² only)
# =========================
summary_rows = []
for target_key in TARGET_KEYS:
    print("\n" + "="*70)
    print(f"Target: {target_key}")
    print("="*70)

    Y_all = Z[target_key][:, 0]  # (n_samples,)

    # ---- Mask out NaN/inf labels (XGBoost requirement) ----
    valid_mask = np.isfinite(Y_all)
    n_valid = valid_mask.sum()
    print(f"Valid samples for {target_key}: {n_valid} / {len(Y_all)}")

    if n_valid < 2:
        print(f"Not enough valid samples for {target_key}, skipping.")
        # still add a row to summary with NaNs
        summary_rows.append(dict(
            target_key=target_key,
            n_subjects=len(unique_subs),
            mean_fold_r2=np.nan,
            sd_fold_r2=np.nan,
            median_fold_r2=np.nan,
            overall_r2=np.nan,
            positives=0
        ))
        continue

    idx_valid = np.where(valid_mask)[0]
    X_flat_valid = X_flat[idx_valid]
    Y_valid = Y_all[idx_valid]
    subs_valid = subject_ids[idx_valid]
    unique_subs_valid = np.unique(subs_valid)

    # ========== LOSO (per-fold scaling) on valid subset ==========
    logo = LeaveOneGroupOut()
    # full-length arrays (aligned to original samples)
    oof_pred_full = np.full_like(Y_all, fill_value=np.nan, dtype=float)
    oof_true_full = Y_all.copy()
    fold_r2 = []
    fold_subjects = []

    for fold, (tr, te) in enumerate(
        logo.split(X_flat_valid, Y_valid, groups=subs_valid), start=1
    ):
        X_tr, X_te = X_flat_valid[tr], X_flat_valid[te]
        y_tr, y_te = Y_valid[tr], Y_valid[te]

        # In LOSO, te corresponds to one subject, but we should still be safe:
        if len(y_tr) < 2:
            print(f"  Fold {fold:02d}: not enough training samples after masking, skipping.")
            continue

        scaler_fold = StandardScaler()
        X_trs = scaler_fold.fit_transform(X_tr)
        X_tes = scaler_fold.transform(X_te)

        model = XGBRegressor(**xgb_params)
        model.fit(X_trs, y_tr)

        y_hat = model.predict(X_tes)
        # map back to full sample indices
        oof_pred_full[idx_valid[te]] = y_hat

        try:
            r2 = r2_score(y_te, y_hat)
        except ValueError:
            # This would happen if len(y_te) < 2, but we expect many trials per subject
            r2 = np.nan
        fold_r2.append(r2)
        fold_subjects.append(subs_valid[te][0])
        print(f"Fold {fold:02d}  subj={subs_valid[te][0]}  R^2: {r2:.3f}")

    # Overall R² only on valid samples
    overall_r2 = r2_score(Y_valid, oof_pred_full[idx_valid])
    print(f"Overall LOSO R^2 ({target_key}): {overall_r2:.3f}")

    # ========== Save per-target NPZ (R² only) ==========
    out_file = os.path.join(results_dir, f"xgb_loso_r2_{target_key}_0to20s.npz")
    np.savez_compressed(
        out_file,
        # OOF predictions and fold metrics
        oof_true_full=oof_true_full,          # includes NaNs
        oof_pred_full=oof_pred_full,          # NaN where label invalid or not tested
        valid_mask=valid_mask,                # to recover subset
        subject_ids=subject_ids,
        unique_subjects=unique_subs,
        unique_subjects_valid=unique_subs_valid,
        fold_r2=np.array(fold_r2, dtype=float),
        fold_subjects=np.array(fold_subjects),
        overall_r2=float(overall_r2),
        # Meta (useful for later alignment/plotting)
        band_names=np.array(band_names),
        bin_times=np.array(bin_times, dtype=float),           # centers (-1,-3,...,-19)
        bin_left_edges=np.array(bin_left_edges, dtype=int),   # [0,2,...,18]
        bin_right_edge=np.int64(bin_right_edge),
        n_channels=n_channels, n_bands=n_bands, n_bins=n_bins,
        target_key=target_key,
        xgb_params=json.dumps(xgb_params)
    )
    print(f"Saved: {out_file}")

    # Append to summary
    r2_arr = np.array(fold_r2, dtype=float)
    summary_rows.append(dict(
        target_key=target_key,
        n_subjects=len(unique_subs_valid),
        mean_fold_r2=float(np.nanmean(r2_arr)) if r2_arr.size > 0 else np.nan,
        sd_fold_r2=float(np.nanstd(r2_arr, ddof=1)) if r2_arr.size > 1 else np.nan,
        median_fold_r2=float(np.nanmedian(r2_arr)) if r2_arr.size > 0 else np.nan,
        overall_r2=float(overall_r2),
        positives=int(np.sum(r2_arr > 0))
    ))

# =========================
# Save a small TSV summary across targets
# =========================
summary_path = os.path.join(results_dir, "loso_r2_summary_0to20s_7nets.tsv")
with open(summary_path, "w") as f:
    cols = ["target_key","n_subjects","mean_fold_r2","sd_fold_r2",
            "median_fold_r2","overall_r2","positives"]
    f.write("\t".join(cols) + "\n")
    for row in summary_rows:
        f.write("\t".join(str(row[c]) for c in cols) + "\n")
print(f"Summary table: {summary_path}")
