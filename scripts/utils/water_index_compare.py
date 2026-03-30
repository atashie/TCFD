"""Statistical comparison engine for water index NetCDF files.

Pure functions for comparing two water index files (e.g., old RCP vs new SSP).
All functions accept 2D numpy arrays and return dicts with results.
"""

import numpy as np
from scipy import stats
from typing import Dict, Optional, Tuple


def spatial_correlation(old: np.ndarray, new: np.ndarray) -> Dict[str, float]:
    """Pearson and Spearman correlation between two 2D grids.

    Args:
        old: Reference 2D array (lat, lon).
        new: Comparison 2D array (lat, lon).

    Returns:
        Dict with pearson_r, pearson_p, spearman_r, spearman_p, n_valid.
    """
    mask = ~np.isnan(old) & ~np.isnan(new)
    n = mask.sum()
    if n < 3:
        return {"pearson_r": np.nan, "pearson_p": np.nan,
                "spearman_r": np.nan, "spearman_p": np.nan, "n_valid": int(n)}

    o, n_ = old[mask], new[mask]
    pr, pp = stats.pearsonr(o, n_)
    sr, sp = stats.spearmanr(o, n_)
    return {
        "pearson_r": float(pr), "pearson_p": float(pp),
        "spearman_r": float(sr), "spearman_p": float(sp),
        "n_valid": int(n),
    }


def summary_statistics(data: np.ndarray) -> Dict[str, float]:
    """Compute summary stats for a 2D array.

    Returns:
        Dict with mean, std, min, max, p5, p25, p50, p75, p95, n_valid, n_nan, coverage_pct.
    """
    valid = data[~np.isnan(data)]
    n_total = data.size
    n_valid = len(valid)
    if n_valid == 0:
        return {k: np.nan for k in ["mean", "std", "min", "max", "p5", "p25", "p50", "p75", "p95"]} | {
            "n_valid": 0, "n_nan": n_total, "coverage_pct": 0.0,
        }
    return {
        "mean": float(np.mean(valid)),
        "std": float(np.std(valid)),
        "min": float(np.min(valid)),
        "max": float(np.max(valid)),
        "p5": float(np.percentile(valid, 5)),
        "p25": float(np.percentile(valid, 25)),
        "p50": float(np.percentile(valid, 50)),
        "p75": float(np.percentile(valid, 75)),
        "p95": float(np.percentile(valid, 95)),
        "n_valid": int(n_valid),
        "n_nan": int(n_total - n_valid),
        "coverage_pct": float(n_valid / n_total * 100),
    }


def coverage_comparison(old: np.ndarray, new: np.ndarray) -> Dict[str, float]:
    """Compare NaN patterns between two grids.

    Returns:
        Dict with old_coverage, new_coverage, overlap_pct, old_only_pct, new_only_pct.
    """
    old_valid = ~np.isnan(old)
    new_valid = ~np.isnan(new)
    n_total = old.size

    old_count = old_valid.sum()
    new_count = new_valid.sum()
    overlap = (old_valid & new_valid).sum()
    old_only = (old_valid & ~new_valid).sum()
    new_only = (~old_valid & new_valid).sum()

    return {
        "old_coverage_pct": float(old_count / n_total * 100),
        "new_coverage_pct": float(new_count / n_total * 100),
        "overlap_pct": float(overlap / max(old_count, 1) * 100),
        "old_only_pct": float(old_only / n_total * 100),
        "new_only_pct": float(new_only / n_total * 100),
        "old_count": int(old_count),
        "new_count": int(new_count),
        "overlap_count": int(overlap),
    }


def detect_anomalous_differences(
    old: np.ndarray, new: np.ndarray, n_sigma: float = 3,
) -> Dict[str, object]:
    """Detect outlier differences between two grids.

    Args:
        old: Reference grid.
        new: Comparison grid.
        n_sigma: Threshold for outlier detection (default 3).

    Returns:
        Dict with diff_mean, diff_std, n_outliers, outlier_pct, outlier_locs.
    """
    mask = ~np.isnan(old) & ~np.isnan(new)
    if mask.sum() < 10:
        return {"diff_mean": np.nan, "diff_std": np.nan, "n_outliers": 0,
                "outlier_pct": 0.0, "outlier_locs": []}

    diff = new[mask] - old[mask]
    mu, sigma = np.mean(diff), np.std(diff)

    if sigma == 0:
        return {"diff_mean": float(mu), "diff_std": 0.0, "n_outliers": 0,
                "outlier_pct": 0.0, "outlier_locs": []}

    # Find outlier positions in original grid
    full_diff = np.full_like(old, np.nan)
    full_diff[mask] = new[mask] - old[mask]
    outlier_mask = np.abs(full_diff - mu) > n_sigma * sigma
    n_outliers = int(outlier_mask.sum())

    # Get up to 10 outlier locations (lat_idx, lon_idx)
    outlier_indices = np.argwhere(outlier_mask)[:10]
    outlier_locs = [(int(r[0]), int(r[1])) for r in outlier_indices]

    return {
        "diff_mean": float(mu),
        "diff_std": float(sigma),
        "n_outliers": n_outliers,
        "outlier_pct": float(n_outliers / mask.sum() * 100),
        "outlier_locs": outlier_locs,
    }


def ks_test(old: np.ndarray, new: np.ndarray) -> Dict[str, float]:
    """Kolmogorov-Smirnov test comparing distributions.

    Returns:
        Dict with ks_statistic, p_value.
    """
    o = old[~np.isnan(old)]
    n = new[~np.isnan(new)]
    if len(o) < 10 or len(n) < 10:
        return {"ks_statistic": np.nan, "p_value": np.nan}

    # Subsample for performance if very large
    max_sample = 100_000
    if len(o) > max_sample:
        rng = np.random.default_rng(42)
        o = rng.choice(o, max_sample, replace=False)
    if len(n) > max_sample:
        rng = np.random.default_rng(42)
        n = rng.choice(n, max_sample, replace=False)

    stat, pval = stats.ks_2samp(o, n)
    return {"ks_statistic": float(stat), "p_value": float(pval)}


def rmsd(old: np.ndarray, new: np.ndarray) -> float:
    """Root mean square difference between two grids."""
    mask = ~np.isnan(old) & ~np.isnan(new)
    if mask.sum() == 0:
        return np.nan
    diff = new[mask] - old[mask]
    return float(np.sqrt(np.mean(diff**2)))


def flag_suspicious_patterns(data: np.ndarray) -> Dict[str, bool]:
    """Check for suspicious patterns in a 2D grid.

    Returns:
        Dict of boolean flags: all_zero, all_nan, constant, has_inf, has_negative.
    """
    return {
        "all_nan": bool(np.all(np.isnan(data))),
        "all_zero": bool(np.nansum(np.abs(data)) == 0) if not np.all(np.isnan(data)) else False,
        "constant": bool(np.nanstd(data) == 0) if not np.all(np.isnan(data)) else False,
        "has_inf": bool(np.any(np.isinf(data))),
        "has_negative": bool(np.any(data[~np.isnan(data)] < 0)) if not np.all(np.isnan(data)) else False,
    }


def compare_seasonal_cycle(
    old_12: np.ndarray, new_12: np.ndarray,
) -> Dict[str, object]:
    """Compare 12-month seasonal cycles at a single location.

    Args:
        old_12: Array of 12 monthly values (old file).
        new_12: Array of 12 monthly values (new file).

    Returns:
        Dict with correlation, phase_shift, amplitude_ratio, rmsd.
    """
    mask = ~np.isnan(old_12) & ~np.isnan(new_12)
    n = mask.sum()

    if n < 6:
        return {"correlation": np.nan, "phase_shift_months": np.nan,
                "amplitude_ratio": np.nan, "rmsd": np.nan}

    o, n_ = old_12[mask], new_12[mask]

    # Correlation
    corr = float(np.corrcoef(o, n_)[0, 1]) if np.std(o) > 0 and np.std(n_) > 0 else np.nan

    # Amplitude ratio (range of new / range of old)
    old_range = np.ptp(o)
    new_range = np.ptp(n_)
    amp_ratio = float(new_range / old_range) if old_range > 0 else np.nan

    # Phase shift (cross-correlation peak)
    if n == 12 and np.std(o) > 0 and np.std(n_) > 0:
        cross_corr = np.correlate(
            (o - o.mean()) / o.std(),
            (n_ - n_.mean()) / n_.std(),
            mode="full",
        )
        # Peak position relative to center
        peak = np.argmax(cross_corr) - (len(cross_corr) // 2)
        phase_shift = int(peak)
    else:
        phase_shift = np.nan

    # RMSD
    rmsd_val = float(np.sqrt(np.mean((o - n_)**2)))

    return {
        "correlation": corr,
        "phase_shift_months": phase_shift,
        "amplitude_ratio": amp_ratio,
        "rmsd": rmsd_val,
    }
