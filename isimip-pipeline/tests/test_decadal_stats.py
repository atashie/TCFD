"""Tests for the unified TCFD decadal-statistics contract.

The two properties that matter most and are cheapest to get wrong:
  * the boolean branch must be selected from VALUES, not names (GUARDRAILS §9);
  * the baseline decade's slopes must be NaN, never 0 -- a 0 there makes the whole
    ocean a finite zero and the QA report does not catch it.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.utils.decadal_stats import (  # noqa: E402
    decade_window_mask,
    expanding_slopes,
    is_boolean_field,
    pooled_decadal_stat,
)

YEARS = np.arange(2020, 2100)
BASELINE = 2020


# ---------------------------------------------------------------------------
# Field classification
# ---------------------------------------------------------------------------

def test_exact_binary_is_boolean():
    assert is_boolean_field(np.array([0.0, 1.0, 0.0, 1.0, np.nan]))


def test_single_valued_zero_field_is_boolean():
    assert is_boolean_field(np.zeros(100))


def test_continuous_fraction_is_not_boolean():
    """`let` is a continuous fraction in [0,1) -- must NOT take the mean/SD branch."""
    rng = np.random.RandomState(0)
    assert not is_boolean_field(rng.uniform(0, 1, 1000))


def test_zero_inflated_fraction_is_not_boolean():
    """Near-boolean is not boolean: 95% exact zeros plus a continuous tail."""
    vals = np.concatenate([np.zeros(950), np.random.RandomState(1).uniform(0.01, 0.9, 50)])
    assert not is_boolean_field(vals)


def test_all_nan_is_not_boolean():
    assert not is_boolean_field(np.full(10, np.nan))


def test_zero_and_one_hundred_is_not_boolean():
    """One-hot 0/100 presence flags (MC2 biome classes) are not {0,1}."""
    assert not is_boolean_field(np.array([0.0, 100.0, 0.0]))


# ---------------------------------------------------------------------------
# Central value + CI
# ---------------------------------------------------------------------------

def test_continuous_uses_median_and_iqr():
    # one cell, 1 member, 10 years, values 0..9 -> median 4.5, q25 2.25, q75 6.75
    annual = np.arange(10, dtype=np.float32).reshape(1, 10, 1)
    stat, lo, hi = pooled_decadal_stat(annual, np.arange(2020, 2030), 2020, boolean=False)
    assert stat[0] == pytest.approx(4.5)
    assert lo[0] == pytest.approx(2.25)
    assert hi[0] == pytest.approx(6.75)


def test_boolean_uses_mean_and_sd():
    # 6 ones, 4 zeros -> mean 0.6, population SD 0.4899
    vals = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.float32).reshape(1, 10, 1)
    stat, lo, hi = pooled_decadal_stat(vals, np.arange(2020, 2030), 2020, boolean=True)
    assert stat[0] == pytest.approx(0.6)
    sd = np.std(vals)
    assert lo[0] == pytest.approx(0.6 - sd, abs=1e-6)
    assert hi[0] == pytest.approx(0.6 + sd, abs=1e-6)


def test_pooling_is_over_year_and_member_jointly():
    """The sample is every (year, member) pair, not a two-stage collapse.

    Member A is flat 0, member B is flat 10. A mean-then-median or median-then-mean
    both give 5; the pooled median of {0 x10, 10 x10} is also 5, so use an odd,
    asymmetric layout that separates them: 3 members, values 0 / 0 / 12.
    Pooled median = 0. Two-stage (mean over members first) = 4.
    """
    annual = np.zeros((3, 10, 1), dtype=np.float32)
    annual[2, :, 0] = 12.0
    stat, _, _ = pooled_decadal_stat(annual, np.arange(2020, 2030), 2020, boolean=False)
    assert stat[0] == pytest.approx(0.0)


def test_all_nan_cell_stays_nan():
    annual = np.full((2, 10, 1), np.nan, dtype=np.float32)
    stat, lo, hi = pooled_decadal_stat(annual, np.arange(2020, 2030), 2020, boolean=False)
    assert np.isnan(stat[0]) and np.isnan(lo[0]) and np.isnan(hi[0])


def test_window_mask_is_ten_years_inclusive():
    m = decade_window_mask(YEARS, 2030)
    assert YEARS[m].min() == 2030 and YEARS[m].max() == 2039 and m.sum() == 10


# ---------------------------------------------------------------------------
# Slopes
# ---------------------------------------------------------------------------

def _linear(slope, n_members=3, noise=0.0, seed=0):
    """(member, year, cell) with a known per-YEAR slope."""
    rng = np.random.RandomState(seed)
    base = slope * (YEARS - YEARS[0])
    out = np.empty((n_members, YEARS.size, 1), np.float32)
    for m in range(n_members):
        out[m, :, 0] = base + m * 100.0 + (rng.normal(0, noise, YEARS.size) if noise else 0)
    return out


def test_baseline_decade_is_nan_not_zero():
    """The defect this contract exists to prevent."""
    res = expanding_slopes(_linear(0.5), YEARS, BASELINE, BASELINE)
    assert np.isnan(res.ols_slope[0])
    assert np.isnan(res.sen_slope[0])


def test_both_slopes_recover_a_clean_linear_trend():
    """Level offsets of +0/+100/+200 between members must not bias the slope."""
    res = expanding_slopes(_linear(0.5), YEARS, 2090, BASELINE)
    assert res.ols_slope[0] == pytest.approx(0.5, abs=1e-3)
    assert res.sen_slope[0] == pytest.approx(0.5, abs=1e-3)


def test_slope_units_are_per_year():
    res = expanding_slopes(_linear(2.0), YEARS, 2090, BASELINE)
    assert res.ols_slope[0] == pytest.approx(2.0, abs=1e-3)


def test_expanding_window_shortens_for_early_decades():
    """A series that only rises after 2060 has a smaller 2030s slope than 2090s."""
    annual = np.zeros((1, YEARS.size, 1), np.float32)
    late = YEARS >= 2060
    annual[0, late, 0] = (YEARS[late] - 2060) * 1.0
    early = expanding_slopes(annual, YEARS, 2030, BASELINE).ols_slope[0]
    full = expanding_slopes(annual, YEARS, 2090, BASELINE).ols_slope[0]
    assert abs(early) < abs(full)


def test_sen_matches_scipy_theilslopes():
    """Cross-check the vectorized implementation against scipy on stacked points."""
    annual = _linear(0.7, n_members=3, noise=2.0, seed=42)
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)

    mask = (YEARS >= BASELINE) & (YEARS <= 2099)
    y = annual[:, mask, 0].ravel()
    x = np.tile(YEARS[mask].astype(float), annual.shape[0])
    expected = stats.theilslopes(y, x)[0]

    assert res.sen_slope[0] == pytest.approx(expected, abs=1e-4)


def test_ols_matches_numpy_polyfit():
    annual = _linear(0.7, n_members=3, noise=2.0, seed=7)
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)

    mask = (YEARS >= BASELINE) & (YEARS <= 2099)
    y = annual[:, mask, 0].ravel()
    x = np.tile(YEARS[mask].astype(float), annual.shape[0])
    expected = np.polyfit(x, y, 1)[0]

    assert res.ols_slope[0] == pytest.approx(expected, abs=1e-4)


def test_same_year_ties_are_dropped_not_infinite():
    """12 members share every year; those pairs have dx=0 and must not produce inf."""
    annual = _linear(0.3, n_members=12, noise=1.0, seed=3)
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert np.isfinite(res.sen_slope[0])
    assert res.sen_slope[0] == pytest.approx(0.3, abs=0.05)


def test_sen_is_robust_to_outliers_where_ols_is_not():
    """The reason for carrying both fields."""
    annual = _linear(0.5, n_members=1, noise=0.0)
    annual[0, 5, 0] += 5000.0          # one wild point
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert res.sen_slope[0] == pytest.approx(0.5, abs=1e-2)
    assert abs(res.ols_slope[0] - 0.5) > 0.1


def test_sen_collapses_to_zero_on_zero_inflated_field():
    """Documented, expected degeneracy -- 91.3% of driedarea ssp126 cells.

    Pinned as a test so nobody 'fixes' the zero by changing the estimator: OLS is the
    field to read for these layers.
    """
    annual = np.zeros((1, YEARS.size, 1), np.float32)
    annual[0, -6:, 0] = 1.0            # activity only in the last few years
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert res.sen_slope[0] == 0.0
    assert res.ols_slope[0] > 0.0


def test_partially_masked_cell_still_gets_slopes():
    """Heterogeneous land masks must not strip trends.

    A cell present in 4 of 5 members has a well-defined median and CI; requiring a
    complete series would return NaN slopes there and silently blank ~13% of csoil
    land while `median` stayed finite.
    """
    annual = _linear(0.5, n_members=5, noise=0.0)
    annual[3, :, 0] = np.nan                      # one member masked out here
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert np.isfinite(res.ols_slope[0])
    assert np.isfinite(res.sen_slope[0])
    assert res.ols_slope[0] == pytest.approx(0.5, abs=1e-3)
    assert res.sen_slope[0] == pytest.approx(0.5, abs=1e-3)


def test_partial_masking_matches_the_unmasked_subset_exactly():
    """Dropping a member must give the same answer as never having had it."""
    full = _linear(0.37, n_members=5, noise=1.0, seed=21)
    masked = full.copy()
    masked[3, :, 0] = np.nan
    dropped = np.delete(full, 3, axis=0)

    a = expanding_slopes(masked, YEARS, 2090, BASELINE, max_pairs=None)
    b = expanding_slopes(dropped, YEARS, 2090, BASELINE, max_pairs=None)
    assert a.ols_slope[0] == pytest.approx(b.ols_slope[0], abs=1e-6)
    assert a.sen_slope[0] == pytest.approx(b.sen_slope[0], abs=1e-6)


def _flat_members(slope, n_members, offset):
    """Members with a common slope and a controlled level offset between them."""
    out = np.empty((n_members, YEARS.size, 1), np.float32)
    for m in range(n_members):
        out[m, :, 0] = slope * (YEARS - YEARS[0]) + m * offset
    return out


def test_scattered_missing_years_are_tolerated_when_members_share_a_level():
    """Gaps within a member's series, not just whole-member masking."""
    annual = _flat_members(0.5, 3, offset=0.0)
    annual[0, ::3, 0] = np.nan
    annual[2, 10:20, 0] = np.nan
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert res.ols_slope[0] == pytest.approx(0.5, abs=1e-3)
    assert res.sen_slope[0] == pytest.approx(0.5, abs=1e-3)


def test_unbalanced_members_bias_OLS_but_not_SEN():
    """The measured consequence of stacking members with level offsets.

    When member representation is UNBALANCED across years (scattered masking) and the
    members sit at different levels, the offset correlates with time-of-observation and
    OLS absorbs it as trend. Theil-Sen, being a median of pairwise slopes, does not.

    Measured here: true slope 0.5, offset 100 -> OLS reads +0.70 (a 40% overstatement),
    Sen reads +0.50 exactly. This is the documented 68.7x level-offset hazard on
    `csoil-total` showing up in the OLS field.

    Pinned so the pair is never collapsed to a single trend variable: the two estimators
    fail in OPPOSITE regimes -- Sen collapses to 0 on zero-inflated hazards (see
    test_sen_collapses_to_zero_on_zero_inflated_field) where OLS is fine, and OLS
    inflates on unbalanced ensembles where Sen is fine.
    """
    annual = _flat_members(0.5, 3, offset=100.0)
    annual[0, ::3, 0] = np.nan
    annual[2, 10:20, 0] = np.nan
    res = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=None)

    assert res.sen_slope[0] == pytest.approx(0.5, abs=1e-3)
    assert res.ols_slope[0] > 0.6          # biased high by the level offset
    assert res.ols_slope[0] == pytest.approx(0.7012, abs=5e-3)


def test_too_few_observations_gives_nan_not_a_spurious_trend():
    annual = np.full((1, YEARS.size, 1), np.nan, np.float32)
    annual[0, 0, 0] = 1.0
    annual[0, 5, 0] = 2.0                          # only 2 finite points < MIN_SLOPE_OBS
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert np.isnan(res.ols_slope[0]) and np.isnan(res.sen_slope[0])


def test_single_year_of_data_gives_nan_not_divide_by_zero():
    """All observations at one x -> zero variance in x; must be NaN, not inf."""
    annual = np.full((4, YEARS.size, 1), np.nan, np.float32)
    annual[:, 0, 0] = [1.0, 2.0, 3.0, 4.0]
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert np.isnan(res.ols_slope[0]) and np.isnan(res.sen_slope[0])


def test_large_magnitude_stock_with_small_trend_survives_precision():
    """A float32 cast before fitting would quantize this trend to zero.

    float32 spacing at 1e9 is ~64, so a stock at that level rising 1/yr collapses to
    repeated values. The module must not coerce the input cube.
    """
    base = 1.0e9
    annual = np.empty((2, YEARS.size, 1), np.float64)
    for m in range(2):
        annual[m, :, 0] = base + 1.0 * (YEARS - YEARS[0])
    res = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=None)
    assert res.ols_slope[0] == pytest.approx(1.0, rel=1e-3)
    assert res.sen_slope[0] == pytest.approx(1.0, rel=1e-3)


def test_float32_input_is_still_accepted():
    """Precision is preserved, not mandated -- float32 callers must keep working."""
    annual = _linear(0.5, n_members=3).astype(np.float32)
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert res.ols_slope[0] == pytest.approx(0.5, abs=1e-3)


def test_boolean_classifier_catches_a_single_rare_continuous_value():
    """A strided sample can step over a rare tail; the scan must not.

    Misclassifying a zero-inflated fraction as boolean silently swaps the published
    statistic from median/IQR to mean/SD.
    """
    vals = np.zeros(10_000_000, dtype=np.float32)
    vals[7_654_321] = 0.37                       # one non-boolean value, off-stride
    assert not is_boolean_field(vals)


def test_boolean_classifier_still_accepts_a_large_true_binary_field():
    vals = (np.arange(10_000_000) % 2).astype(np.float32)
    assert is_boolean_field(vals)


def test_subsample_error_is_not_spatially_uniform():
    """Per-chunk reseeding must decorrelate the approximation error across cells.

    With one shared pair sample every cell inherits the same error, which can read as
    geographic structure. Identical cells should NOT all land on the same wrong value.
    """
    one_cell = _linear(0.6, n_members=12, noise=1.5, seed=9)
    annual = np.repeat(one_cell, 400, axis=2)     # 400 identical cells
    res = expanding_slopes(annual, YEARS, 2090, BASELINE,
                           chunk_cells=50, max_pairs=20_000)
    assert np.isfinite(res.sen_slope).all()
    assert np.unique(res.sen_slope).size > 1      # errors differ across chunks


def test_results_are_reproducible_across_runs():
    """Decorrelated must not mean nondeterministic."""
    annual = np.repeat(_linear(0.6, n_members=12, noise=1.5, seed=9), 200, axis=2)
    a = expanding_slopes(annual, YEARS, 2090, BASELINE, chunk_cells=50, max_pairs=20_000)
    b = expanding_slopes(annual, YEARS, 2090, BASELINE, chunk_cells=50, max_pairs=20_000)
    assert np.array_equal(a.sen_slope, b.sen_slope)


def test_all_nan_cell_gives_nan_slopes():
    annual = np.full((2, YEARS.size, 1), np.nan, np.float32)
    res = expanding_slopes(annual, YEARS, 2090, BASELINE)
    assert np.isnan(res.ols_slope[0]) and np.isnan(res.sen_slope[0])


def test_chunking_does_not_change_results():
    annual = np.repeat(_linear(0.4, n_members=2, noise=1.0, seed=11), 300, axis=2)
    a = expanding_slopes(annual, YEARS, 2090, BASELINE, chunk_cells=17)
    b = expanding_slopes(annual, YEARS, 2090, BASELINE, chunk_cells=4096)
    assert np.allclose(a.sen_slope, b.sen_slope, equal_nan=True)
    assert np.allclose(a.ols_slope, b.ols_slope, equal_nan=True)


def test_default_is_exact_not_subsampled():
    """Production default must approximate nothing.

    Subsampling was measured to cost ~15% of the slope at csoil's 68x member-offset
    ratio, so DEFAULT_MAX_PAIRS is None. Pinned so a speedup is never quietly restored.
    """
    from scripts.utils.decadal_stats import DEFAULT_MAX_PAIRS
    assert DEFAULT_MAX_PAIRS is None

    annual = _linear(0.01, n_members=12, noise=0.5, seed=0)
    default = expanding_slopes(annual, YEARS, 2090, BASELINE)
    exact = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=None)
    assert default.sen_slope[0] == exact.sen_slope[0]


def test_subsampling_error_grows_with_member_level_offset():
    """The finding that disqualified subsampling as a default.

    Cross-member pairs dominate the pairwise-slope distribution and flatten it around
    the median, so the sample median degrades as the between-member level offset grows
    relative to interannual noise. Pinned so the accuracy table in the module docstring
    cannot silently rot.
    """
    def build(offset_sd, seed=0):
        rng = np.random.RandomState(seed)
        out = np.empty((12, YEARS.size, 1), np.float64)
        off = rng.normal(0, offset_sd, 12)
        for m in range(12):
            out[m, :, 0] = (0.01 * (YEARS - YEARS[0]) + off[m]
                            + rng.normal(0, 0.5, YEARS.size))
        return out

    def err(offset_sd):
        a = build(offset_sd)
        ex = expanding_slopes(a, YEARS, 2090, BASELINE, max_pairs=None)
        ap = expanding_slopes(a, YEARS, 2090, BASELINE, max_pairs=100_000)
        return abs(float(ap.sen_slope[0] - ex.sen_slope[0]))

    assert err(34.0) > 3 * err(0.0)


def test_ols_is_unaffected_by_the_pair_cap():
    """The cap is a Theil-Sen optimisation only; OLS uses every point regardless."""
    annual = _linear(0.01, n_members=12, noise=0.5, seed=0)
    a = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=20_000)
    b = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=None)
    assert a.ols_slope[0] == pytest.approx(b.ols_slope[0], abs=1e-9)


def test_pair_subsampling_is_deterministic_and_close():
    annual = _linear(0.6, n_members=4, noise=1.5, seed=5)
    exact = expanding_slopes(annual, YEARS, 2090, BASELINE)
    s1 = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=20_000)
    s2 = expanding_slopes(annual, YEARS, 2090, BASELINE, max_pairs=20_000)
    assert s1.sen_slope[0] == s2.sen_slope[0]          # deterministic
    assert s1.sen_slope[0] == pytest.approx(exact.sen_slope[0], abs=0.02)
