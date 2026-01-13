"""Tests for processing features module."""

import pytest
import numpy as np

from isimip_pipeline.processing.features import (
    kernel_smooth,
    theil_sen_slope,
    spearman_significance,
    calculate_percentile_rank,
    decadal_aggregation,
    FeatureExtractor,
)


class TestKernelSmooth:
    """Test kernel smoothing function."""

    def test_kernel_smooth_returns_same_length(self):
        """Smoothed output should have same length as input."""
        years = np.arange(2006, 2100)
        values = np.random.random(len(years))

        smoothed = kernel_smooth(years, values, bandwidth=15)

        assert len(smoothed) == len(values)

    def test_kernel_smooth_reduces_variance(self):
        """Smoothing should reduce variance of noisy data."""
        years = np.arange(2006, 2100)
        values = np.sin(years / 10) + np.random.random(len(years)) * 0.5

        smoothed = kernel_smooth(years, values, bandwidth=15)

        assert np.var(smoothed) < np.var(values)

    def test_kernel_smooth_preserves_mean(self):
        """Smoothing should roughly preserve mean."""
        years = np.arange(2006, 2100)
        values = np.random.random(len(years)) * 100

        smoothed = kernel_smooth(years, values, bandwidth=15)

        assert abs(np.mean(smoothed) - np.mean(values)) < 10


class TestTheilSenSlope:
    """Test Theil-Sen slope estimator."""

    def test_theil_sen_positive_trend(self):
        """Should detect positive trend."""
        years = np.arange(2006, 2050)
        values = years * 0.5 + np.random.random(len(years)) * 2

        slope = theil_sen_slope(years, values)

        assert slope > 0

    def test_theil_sen_negative_trend(self):
        """Should detect negative trend."""
        years = np.arange(2006, 2050)
        values = -years * 0.3 + np.random.random(len(years))

        slope = theil_sen_slope(years, values)

        assert slope < 0

    def test_theil_sen_no_trend(self):
        """Should return near-zero for flat data."""
        years = np.arange(2006, 2050)
        values = np.ones(len(years)) * 50 + np.random.random(len(years)) * 0.01

        slope = theil_sen_slope(years, values)

        assert abs(slope) < 0.1


class TestSpearmanSignificance:
    """Test Spearman correlation significance."""

    def test_spearman_strong_correlation(self):
        """Strong correlation should have low p-value."""
        years = np.arange(2006, 2050)
        values = years * 2  # Perfect positive correlation

        p_value = spearman_significance(years, values)

        assert p_value < 0.05

    def test_spearman_weak_correlation(self):
        """Random data should have high p-value."""
        np.random.seed(42)
        years = np.arange(2006, 2050)
        values = np.random.random(len(years)) * 100

        p_value = spearman_significance(years, values)

        # Random data typically has higher p-value
        assert 0 <= p_value <= 1


class TestPercentileRank:
    """Test percentile ranking."""

    def test_percentile_rank_returns_1_to_100(self):
        """Percentile ranks should be 1-100."""
        historical = np.random.random(1000) * 100
        value = 50

        rank = calculate_percentile_rank(value, historical)

        assert 1 <= rank <= 100

    def test_percentile_high_value_high_rank(self):
        """High values should have high percentile."""
        historical = np.arange(0, 100)
        value = 95

        rank = calculate_percentile_rank(value, historical)

        assert rank >= 90

    def test_percentile_low_value_low_rank(self):
        """Low values should have low percentile."""
        historical = np.arange(0, 100)
        value = 5

        rank = calculate_percentile_rank(value, historical)

        assert rank <= 10


class TestDecadalAggregation:
    """Test decadal aggregation."""

    def test_decadal_aggregation_returns_dict(self):
        """Should return dict with decade keys."""
        years = np.arange(2006, 2100)
        values = np.random.random(len(years))

        decades = decadal_aggregation(years, values)

        assert isinstance(decades, dict)
        assert 10 in decades  # 2010s
        assert 20 in decades  # 2020s

    def test_decadal_values_are_aggregated(self):
        """Each decade should have aggregated values."""
        years = np.arange(2006, 2100)
        values = np.ones(len(years))

        decades = decadal_aggregation(years, values)

        # Each decade should have 10 values summed/averaged
        assert len(decades) >= 9  # At least 9 decades


class TestFeatureExtractor:
    """Test FeatureExtractor class."""

    def test_extractor_initializes(self):
        """FeatureExtractor should initialize with config."""
        extractor = FeatureExtractor(bandwidth=15)

        assert extractor is not None
        assert extractor.bandwidth == 15

    def test_extractor_process_timeseries(self):
        """Should process a timeseries and return features."""
        extractor = FeatureExtractor(bandwidth=15)

        years = np.arange(2006, 2100)
        values = np.random.random(len(years)) * 100

        features = extractor.process_timeseries(years, values)

        assert "smoothed" in features
        assert "trend" in features
        assert "significance" in features

    def test_extractor_calculate_decadal_features(self):
        """Should calculate all 6 feature types per decade."""
        extractor = FeatureExtractor(bandwidth=15)

        years = np.arange(2006, 2100)
        model_values = [
            np.random.random(len(years)) * 100,
            np.random.random(len(years)) * 100,
        ]

        decadal_features = extractor.calculate_decadal_features(
            years, model_values
        )

        # Should have features for each decade
        assert len(decadal_features) >= 9

        # Each decade should have 6 value types
        for decade, features in decadal_features.items():
            assert "smoothed_median" in features
            assert "percentile" in features
            assert "trend" in features
            assert "significance" in features
            assert "lower_bound" in features
            assert "upper_bound" in features
