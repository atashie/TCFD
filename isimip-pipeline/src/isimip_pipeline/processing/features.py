"""Feature extraction for ISIMIP climate data processing.

Implements statistical features matching the R workflow:
1. Kernel smoothing
2. Theil-Sen trend estimation
3. Spearman significance
4. Percentile ranking
5. Decadal aggregation
"""

from dataclasses import dataclass
from typing import Dict, List, Any

import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter1d


def kernel_smooth(
    x: np.ndarray,
    y: np.ndarray,
    bandwidth: float = 15,
) -> np.ndarray:
    """Apply kernel smoothing to time series.

    Uses Gaussian kernel smoothing similar to R's ksmooth.

    Args:
        x: Independent variable (e.g., years).
        y: Dependent variable (values to smooth).
        bandwidth: Smoothing bandwidth in same units as x.

    Returns:
        Smoothed y values.
    """
    # Calculate sigma from bandwidth (bandwidth ≈ 4*sigma for Gaussian)
    dx = np.mean(np.diff(x)) if len(x) > 1 else 1
    sigma = bandwidth / (4 * dx)

    smoothed = gaussian_filter1d(y.astype(float), sigma=sigma, mode="nearest")

    return smoothed


def theil_sen_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate Theil-Sen slope estimator.

    Robust alternative to OLS that is resistant to outliers.

    Args:
        x: Independent variable.
        y: Dependent variable.

    Returns:
        Theil-Sen slope estimate.
    """
    result = stats.theilslopes(y, x)
    return float(result.slope)


def spearman_significance(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate Spearman correlation p-value.

    Tests for monotonic relationship between x and y.

    Args:
        x: Independent variable.
        y: Dependent variable.

    Returns:
        P-value from Spearman correlation test.
    """
    _, p_value = stats.spearmanr(x, y)
    return float(p_value) if not np.isnan(p_value) else 1.0


def calculate_percentile_rank(
    value: float,
    historical: np.ndarray,
    bins: int = 100,
) -> int:
    """Calculate percentile rank of value against historical distribution.

    Args:
        value: Value to rank.
        historical: Historical distribution to compare against.
        bins: Number of percentile bins (default 100).

    Returns:
        Percentile rank from 1 to bins.
    """
    if len(historical) == 0:
        return 1

    percentile = stats.percentileofscore(historical, value, kind="weak")
    rank = max(1, min(bins, int(np.ceil(percentile * bins / 100))))

    return rank


def decadal_aggregation(
    years: np.ndarray,
    values: np.ndarray,
) -> Dict[int, np.ndarray]:
    """Aggregate values by decade.

    Args:
        years: Array of years (e.g., 2006-2099).
        values: Array of values corresponding to years.

    Returns:
        Dict mapping decade (10, 20, ..., 90) to values in that decade.
    """
    decades = {}

    base_year = years[0] - (years[0] % 10)  # Start of first decade

    for decade_num in range(10, 100, 10):
        # Calculate decade start and end years
        decade_start = base_year + decade_num - 6
        decade_end = decade_start + 10

        mask = (years >= decade_start) & (years < decade_end)
        if np.any(mask):
            decades[decade_num] = values[mask]

    return decades


@dataclass
class DecadalFeatures:
    """Features for a single decade."""

    smoothed_median: float
    percentile: int
    trend: float
    significance: float
    lower_bound: float
    upper_bound: float


class FeatureExtractor:
    """Extract statistical features from climate time series.

    Produces 6 feature types per decade:
    1. Smoothed median value
    2. Percentile rank
    3. Trend (Theil-Sen slope)
    4. Significance (Spearman p-value)
    5. Lower confidence bound
    6. Upper confidence bound
    """

    def __init__(
        self,
        bandwidth: float = 15,
        percentile_bins: int = 100,
    ):
        """Initialize feature extractor.

        Args:
            bandwidth: Kernel smoothing bandwidth.
            percentile_bins: Number of bins for percentile ranking.
        """
        self.bandwidth = bandwidth
        self.percentile_bins = percentile_bins
        self.historical_distribution: np.ndarray = np.array([])

    def process_timeseries(
        self,
        years: np.ndarray,
        values: np.ndarray,
    ) -> Dict[str, Any]:
        """Process a single time series to extract basic features.

        Args:
            years: Array of years.
            values: Array of values.

        Returns:
            Dict with smoothed, trend, and significance values.
        """
        smoothed = kernel_smooth(years, values, self.bandwidth)
        trend = theil_sen_slope(years, values)
        significance = spearman_significance(years, values)

        return {
            "smoothed": smoothed,
            "trend": trend,
            "significance": significance,
        }

    def calculate_decadal_features(
        self,
        years: np.ndarray,
        model_values: List[np.ndarray],
    ) -> Dict[int, Dict[str, Any]]:
        """Calculate all 6 feature types for each decade.

        Follows the R workflow pattern of computing median across models.

        Args:
            years: Array of years.
            model_values: List of value arrays, one per model.

        Returns:
            Dict mapping decade to dict of feature values.
        """
        # Smooth each model's values
        smoothed_models = [
            kernel_smooth(years, vals, self.bandwidth) for vals in model_values
        ]

        # Set up historical distribution from first decade if not set
        if len(self.historical_distribution) == 0:
            first_decade_values = []
            for vals in model_values:
                decades = decadal_aggregation(years, vals)
                if 10 in decades:
                    first_decade_values.extend(decades[10])
            self.historical_distribution = np.array(first_decade_values)

        decadal_features = {}
        decades_list = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        for decade in decades_list:
            # Get decade indices
            decade_start = years[0] - (years[0] % 10) + decade - 6
            decade_end = decade_start + 10
            decade_mask = (years >= decade_start) & (years < decade_end)

            if not np.any(decade_mask):
                continue

            # Calculate smoothed median across models
            decade_smoothed = []
            for sm in smoothed_models:
                decade_smoothed.extend(sm[decade_mask])
            smoothed_median = float(np.median(decade_smoothed))

            # Calculate percentile rank
            percentile = calculate_percentile_rank(
                smoothed_median,
                self.historical_distribution,
                self.percentile_bins,
            )

            # Calculate confidence bounds using IQR
            decade_raw = []
            for vals in model_values:
                decade_raw.extend(vals[decade_mask])
            q25, q50, q75 = np.percentile(decade_raw, [25, 50, 75])
            lower_bound = smoothed_median - abs(q50 - q25)
            upper_bound = smoothed_median + abs(q75 - q50)

            # Calculate trend up to this decade
            trend_mask = years < decade_end
            trend_years = years[trend_mask]
            trends = []
            significances = []
            for vals in model_values:
                trend_vals = vals[trend_mask]
                trends.append(theil_sen_slope(trend_years, trend_vals))
                significances.append(
                    spearman_significance(trend_years, trend_vals)
                )

            trend = float(np.median(trends))
            significance = float(np.median(significances))

            decadal_features[decade] = {
                "smoothed_median": smoothed_median,
                "percentile": percentile,
                "trend": trend,
                "significance": significance,
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
            }

        return decadal_features
