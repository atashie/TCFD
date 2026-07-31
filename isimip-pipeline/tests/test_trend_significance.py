"""Tests for scripts/utils/trend_significance.py.

The p-value gates a customer-facing significance claim, so the tests check it
against an INDEPENDENT textbook Mann-Kendall implementation written from the
formula rather than against the module under test. tau-b is checked against
scipy, where it is an exact algebraic identity with no approximation involved.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from utils.trend_significance import (  # noqa: E402
    MIN_OBS,
    mk_expanding,
    mk_pvalue,
    mk_stats,
    significance_definition,
)


# --------------------------------------------------------------------------
# Independent reference implementation (textbook MK, scalar, no vectorization)
# --------------------------------------------------------------------------
def ref_mk(y):
    """Two-sided tie-corrected MK p-value and tau-b for one series.

    Written directly from the standard formulation so it shares no code with the
    implementation being tested. x is the index, untied and increasing.
    """
    y = np.asarray([v for v in y if np.isfinite(v)], float)
    n = len(y)
    S = 0.0
    for i in range(n - 1):
        for j in range(i + 1, n):
            S += np.sign(y[j] - y[i])
    _, counts = np.unique(y, return_counts=True)
    ties = counts[counts > 1].astype(float)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(ties * (ties - 1) * (2 * ties + 5))) / 18.0
    n0 = n * (n - 1) / 2.0
    n_ties = np.sum(ties * (ties - 1) / 2.0)
    if var <= 0:
        return 1.0, 0.0, n
    z = (S - np.sign(S)) / np.sqrt(var)
    p = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    tau = S / np.sqrt(n0 * (n0 - n_ties))
    return p, tau, n


# --------------------------------------------------------------------------
# Agreement with the reference and with scipy
# --------------------------------------------------------------------------
@pytest.mark.parametrize("y", [
    np.arange(20.0),                                    # perfect monotone, untied
    np.arange(20.0)[::-1].copy(),                        # perfect decreasing
    np.array([0.0] * 10 + [1.0] * 10),                   # heavy ties, monotone
    np.array([0, 0, 1, 0, 1, 1, 0, 1, 1, 1.0] * 2),      # binary-like, noisy
    np.array([3.0, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8]),    # some ties, no clear trend
    np.linspace(0, 1, 30) + np.sin(np.arange(30)),       # noisy but rising
    np.array([5.0] * 12),                                # constant
])
def test_matches_independent_reference(y):
    p, tau, n = mk_pvalue(y)
    rp, rtau, rn = ref_mk(y)
    assert n == rn
    assert p == pytest.approx(rp, abs=1e-9)
    assert tau == pytest.approx(rtau, abs=1e-6)


@pytest.mark.parametrize("seed", range(12))
def test_tau_matches_scipy_taub(seed):
    """tau-b is an exact identity, so it must match scipy to float precision."""
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 5, size=25).astype(float)        # forces plenty of ties
    _, tau, _ = mk_pvalue(y)
    expected = stats.kendalltau(np.arange(len(y)), y, variant="b").statistic
    assert tau == pytest.approx(expected, abs=1e-6)


# --------------------------------------------------------------------------
# Documented behaviour at the edges
# --------------------------------------------------------------------------
def test_constant_series_is_p_one_not_nan():
    """A flat series has no evidence of a trend; p=1 is defined, NaN is not.

    This matters downstream: a NaN silently fails every `p < 0.05` comparison
    without ever being visible, whereas 1.0 renders and aggregates correctly.
    """
    p, tau, n = mk_pvalue([4.0] * 15)
    assert p == 1.0
    assert tau == 0.0
    assert n == 15


def test_perfect_monotone_is_highly_significant():
    p, tau, _ = mk_pvalue(np.arange(20.0))
    assert tau == pytest.approx(1.0)
    assert p < 1e-8


def test_sign_reversal_flips_tau_and_preserves_p():
    y = np.array([1.0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13])
    p_up, tau_up, _ = mk_pvalue(y)
    p_dn, tau_dn, _ = mk_pvalue(-y)
    assert p_up == pytest.approx(p_dn)
    assert tau_up == pytest.approx(-tau_dn)


def test_nan_years_are_excluded_and_counted():
    y = np.array([1.0, np.nan, 2.0, 3.0, np.nan, 4.0, 5.0])
    p, tau, n = mk_pvalue(y)
    assert n == 5
    assert tau == pytest.approx(1.0)
    p_dense, tau_dense, n_dense = mk_pvalue([1.0, 2.0, 3.0, 4.0, 5.0])
    assert (p, tau, n) == pytest.approx((p_dense, tau_dense, n_dense))


def test_all_nan_series_is_nan():
    p, tau, n = mk_pvalue([np.nan] * 10)
    assert np.isnan(p) and np.isnan(tau)
    assert n == 0


@pytest.mark.parametrize("n", range(0, MIN_OBS))
def test_below_min_obs_is_nan(n):
    p, tau, n_obs = mk_pvalue(list(np.arange(float(n))))
    assert np.isnan(p) and np.isnan(tau)
    assert n_obs == n


def test_p_is_bounded():
    rng = np.random.default_rng(7)
    for _ in range(200):
        y = rng.normal(size=rng.integers(4, 40))
        p, _, _ = mk_pvalue(y)
        assert 0.0 <= p <= 1.0


# --------------------------------------------------------------------------
# Vectorization must not change the answer
# --------------------------------------------------------------------------
def test_vectorized_matches_per_cell():
    rng = np.random.default_rng(3)
    y = rng.normal(size=(24, 5, 6)).cumsum(axis=0)
    y[:, 0, 0] = np.nan                                  # an all-NaN cell
    y[3:, 1, 1] = np.nan                                 # a ragged cell
    y[:, 2, 2] = 9.0                                     # a constant cell
    p, tau, n = mk_stats(y)
    assert p.shape == tau.shape == n.shape == (5, 6)
    for i in range(5):
        for j in range(6):
            ep, etau, en = mk_pvalue(y[:, i, j])
            assert n[i, j] == en
            if np.isnan(ep):
                assert np.isnan(p[i, j])
            else:
                assert p[i, j] == pytest.approx(ep, abs=1e-6)
                assert tau[i, j] == pytest.approx(etau, abs=1e-6)


def test_tie_terms_handle_a_single_group_spanning_the_series():
    """Regression guard: the run accumulator must flush the FINAL run.

    A series that is constant apart from its first value leaves one long tie run
    open at the end of the walk. Forgetting to close it understates the tie
    correction and inflates significance.
    """
    y = np.array([0.0] + [1.0] * 19)
    p, tau, _ = mk_pvalue(y)
    rp, rtau, _ = ref_mk(y)
    assert p == pytest.approx(rp, abs=1e-9)
    assert tau == pytest.approx(rtau, abs=1e-9)


# --------------------------------------------------------------------------
# Expanding windows
# --------------------------------------------------------------------------
DECADES = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]


def _rising(nyear=80, shape=(4, 3)):
    years = np.arange(2020, 2020 + nyear)
    ramp = np.linspace(0, 1, nyear).reshape(-1, 1, 1)
    return years, np.broadcast_to(ramp, (nyear,) + shape).copy()


def test_expanding_windows_have_expected_sample_sizes():
    years, annual = _rising()
    p, tau, n = mk_expanding(years, annual, DECADES)
    assert p.shape == (8, 4, 3)
    # The baseline decade has no elapsed period.
    assert np.all(np.isnan(p[0])) and np.all(np.isnan(tau[0]))
    assert np.all(n[0] == 0)
    for i, expected in enumerate([0, 20, 30, 40, 50, 60, 70, 80]):
        if i:
            assert np.all(n[i] == expected), f"decade {DECADES[i]} n={n[i,0,0]}"


def test_expanding_p_strengthens_as_the_record_lengthens():
    years, annual = _rising()
    p, tau, _ = mk_expanding(years, annual, DECADES)
    series = p[1:, 0, 0]
    assert np.all(np.diff(series) <= 0), series
    assert np.allclose(tau[1:, 0, 0], 1.0)


def test_expanding_baseline_decade_can_be_declared():
    years, annual = _rising()
    p, _, n = mk_expanding(years, annual, DECADES, baseline_decade=2020)
    assert np.all(np.isnan(p[0]))
    assert n[1, 0, 0] == 20


def test_expanding_rejects_a_non_increasing_year_axis():
    years, annual = _rising(nyear=10, shape=(2, 2))
    with pytest.raises(ValueError, match="strictly increasing"):
        mk_expanding(years[::-1].copy(), annual, [2020, 2030])
    with pytest.raises(ValueError, match="does not match|entries but"):
        mk_expanding(years[:5], annual, [2020, 2030])


def test_expanding_propagates_the_off_mask_pattern():
    years, annual = _rising()
    annual[:, 0, 0] = np.nan
    p, tau, n = mk_expanding(years, annual, DECADES)
    assert np.all(np.isnan(p[:, 0, 0]))
    assert np.all(n[1:, 0, 0] == 0)
    assert np.all(np.isfinite(p[1:, 1, 1]))


def test_definition_string_states_the_windows_and_the_caveat():
    text = significance_definition(DECADES)
    assert "2020..(decade+9)" in text
    assert "n=20" in text and "n=80" in text
    assert "NOT" in text and "inter-model agreement" in text
    assert "percentile_direction" in text


# --------------------------------------------------------------------------
# AnnualEnsembleMean — must reproduce the backfill's arithmetic exactly
# --------------------------------------------------------------------------
from utils.trend_significance import AnnualEnsembleMean  # noqa: E402


def _backfill_reference(years, members, shape):
    """The accumulation exactly as backfill_trend_significance does it."""
    total = np.zeros((len(years),) + shape, np.float64)
    count = np.zeros((len(years),) + shape, np.float32)
    for myears, mvals in members:
        aligned = np.full((len(years),) + shape, np.nan, np.float32)
        for k, y in enumerate(myears):
            idx = np.where(years == y)[0]
            if idx.size:
                aligned[idx[0]] = mvals[k]
        good = np.isfinite(aligned)
        total[good] += aligned[good]
        count += good
    return np.where(count > 0, total / np.maximum(count, 1), np.nan).astype(np.float32)


def test_accumulator_matches_the_backfill_reference_bitwise():
    """Same members in the same order must give bit-identical means.

    This is what licenses checking a processor's output against an
    already-published layer: float addition is not associative, so only identical
    arithmetic in identical order reproduces the published bits.
    """
    rng = np.random.default_rng(11)
    years = np.arange(2020, 2040)
    shape = (4, 5)
    members = []
    for m in range(6):
        my = years[m % 3:]                        # ragged start years
        vals = rng.normal(size=(len(my),) + shape).astype(np.float32)
        vals[vals < -1.5] = np.nan                # ragged coverage
        members.append((my, vals))

    acc = AnnualEnsembleMean(2020, 2039, shape)
    for my, vals in members:
        acc.add(my, vals)
    _, got = acc.result()
    expected = _backfill_reference(years, members, shape)
    assert np.array_equal(got, expected, equal_nan=True)
    assert acc.n_added == 6


def test_accumulator_ignores_years_outside_its_window():
    acc = AnnualEnsembleMean(2020, 2029, (2, 2))
    yrs = np.arange(2010, 2041)
    acc.add(yrs, np.ones((len(yrs), 2, 2), np.float32))
    years, mean = acc.result()
    assert years[0] == 2020 and years[-1] == 2029
    assert np.all(mean == 1.0)
    assert np.all(acc.members_per_year == 1)


def test_accumulator_counts_only_finite_cells():
    acc = AnnualEnsembleMean(2020, 2021, (2,))
    acc.add(np.array([2020, 2021]), np.array([[1.0, np.nan], [3.0, 5.0]], np.float32))
    acc.add(np.array([2020, 2021]), np.array([[3.0, 7.0], [np.nan, 9.0]], np.float32))
    _, mean = acc.result()
    assert mean[0, 0] == pytest.approx(2.0)     # (1+3)/2
    assert mean[0, 1] == pytest.approx(7.0)     # only the second member
    assert mean[1, 0] == pytest.approx(3.0)     # only the first member
    assert mean[1, 1] == pytest.approx(7.0)     # (5+9)/2


def test_accumulator_rejects_a_wrong_shaped_member():
    acc = AnnualEnsembleMean(2020, 2021, (3, 3))
    with pytest.raises(ValueError, match="does not match accumulator"):
        acc.add(np.array([2020, 2021]), np.zeros((2, 4, 4), np.float32))


def test_accumulator_all_nan_member_leaves_nan():
    acc = AnnualEnsembleMean(2020, 2021, (2,))
    acc.add(np.array([2020, 2021]), np.full((2, 2), np.nan, np.float32))
    _, mean = acc.result()
    assert np.all(np.isnan(mean))


# --------------------------------------------------------------------------
# theilsen_decadal — the shipped `trend`
# --------------------------------------------------------------------------
from utils.trend_significance import (  # noqa: E402
    TREND_MIN_OBS, theilsen, theilsen_decadal, trend_definition_decadal,
)


def test_decadal_sen_is_per_decade_not_per_year():
    """The x axis is the decade index, so a +1/decade ramp must return 1.0.

    Rescaling by window_years here would inflate every published trend tenfold.
    """
    dec = [2020, 2030, 2040, 2050]
    med = np.arange(4.0).reshape(4, 1, 1)          # +1 per decade
    out = theilsen_decadal(med, dec)
    assert np.isnan(out[0, 0, 0])                   # baseline panel
    assert out[1, 0, 0] == pytest.approx(1.0)
    assert out[3, 0, 0] == pytest.approx(1.0)


def test_decadal_sen_first_panel_equals_the_two_point_rate():
    """With 2 decades there is one pairwise slope, so it IS the anchored rate."""
    dec = [2020, 2030, 2040]
    med = np.array([2.0, 5.0, 11.0]).reshape(3, 1, 1)
    out = theilsen_decadal(med, dec)
    assert out[1, 0, 0] == pytest.approx((5.0 - 2.0) / 1)
    assert TREND_MIN_OBS == 2


DEC8 = [2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090]


def test_annual_sen_returns_zero_on_a_rising_rare_hazard():
    """Pins the pathology that moved `trend` off the annual series.

    A cell whose event rate doubles still reports a slope of EXACTLY zero, because
    the median pairwise slope over 80 mostly-zero years is 0. This is why
    `theilsen_decadal` exists; if this test ever stops failing to detect a trend,
    the annual series became viable and the choice can be revisited.
    """
    annual = np.zeros((80, 1, 1))
    annual[[5, 25], 0, 0] = 1.0                      # 2 events, first half
    annual[[45, 52, 61, 70], 0, 0] = 1.0             # 4 events, second half
    assert theilsen(annual)[0, 0] == 0.0


def test_decadal_sen_sees_a_trend_the_annual_version_misses():
    """On a realistically-resolved decadal series the slope is recovered.

    An ensemble mean over 10 years x many members takes many distinct values, so
    decadal ties are rare — that is the whole reason the decadal series works.
    """
    med = np.array([0.02, 0.03, 0.05, 0.06, 0.09, 0.11, 0.14, 0.17]).reshape(8, 1, 1)
    out = theilsen_decadal(med, DEC8)
    assert out[-1, 0, 0] > 0.0
    assert np.isnan(out[0, 0, 0])
    # and the annual series for the same cell would have been mostly zeros
    annual = np.zeros((80, 1, 1))
    annual[::9, 0, 0] = 1.0
    assert theilsen(annual)[0, 0] == 0.0


def test_decadal_sen_is_a_reduction_not_a_cure_when_decades_are_tied():
    """Honest limit: heavy quantization can still zero the decadal slope.

    With only 8 panels quantized to multiples of 0.1, enough pairs are exactly
    equal that the median slope is 0 even though the series rises. Measured on real
    data this residual is 13.7% of driedarea ssp126 cells (down from 91.3% on the
    annual series) — reduced, not eliminated. Do not document it as a cure.
    """
    med = np.array([0.1, 0.0, 0.1, 0.0, 0.1, 0.1, 0.1, 0.1]).reshape(8, 1, 1)
    out = theilsen_decadal(med, DEC8)
    assert out[-1, 0, 0] == 0.0


def test_decadal_sen_all_nan_and_constant_cells():
    dec = [2020, 2030, 2040]
    med = np.full((3, 1, 3), np.nan)
    med[:, 0, 1] = 7.0                              # constant
    med[:, 0, 2] = [1.0, 2.0, 3.0]                  # rising
    out = theilsen_decadal(med, dec)
    assert np.isnan(out[-1, 0, 0])
    assert out[-1, 0, 1] == 0.0
    assert out[-1, 0, 2] == pytest.approx(1.0)


def test_decadal_sen_rejects_a_mismatched_decade_axis():
    with pytest.raises(ValueError, match="decades but medians"):
        theilsen_decadal(np.zeros((3, 2, 2)), [2020, 2030])


def test_decadal_trend_definition_records_the_zero_inflation_reason():
    text = trend_definition_decadal([2020, 2030, 2040])
    assert "DECADAL MEDIAN series" in text
    assert "91%" in text and "exactly 0" in text
    assert "per decade" in text
