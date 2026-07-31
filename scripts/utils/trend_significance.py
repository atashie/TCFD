"""Mann-Kendall trend significance for TCFD/CDP layers.

The published ``trend`` value class is a BASELINE-ANCHORED RATE
``(median[decade] - median[2020s]) / elapsed decades`` (GUARDRAILS S10). That is
deliberately not a fitted slope, so **no p-value is derivable from it** — a rate
built from two numbers has no residual and no degrees of freedom. This module
supplies the missing significance test.

What is tested
--------------
The **ensemble-mean ANNUAL series**: for each cell, the members' annualized
values are averaged within each year, and Mann-Kendall is run over the years from
the baseline decade's first year up to the end of the target decade. So the 2030s
panel tests 2020-2039 (n=20), the 2090s panel tests 2020-2099 (n=80), and the
baseline decade itself has no elapsed period and is left NaN.

Why the ensemble mean rather than stacking members as observations
-----------------------------------------------------------------
Stacking every member-year into one sample looks more powerful (n = years x
members) but is not, because the sample is then dominated by BETWEEN-MODEL LEVEL
OFFSETS rather than by time. Measured on ``csoil-total`` ssp585, the
between-member level SD is **68.7x** the within-member interannual SD, and the
consequence is severe — share of land reaching p<0.05 over 2020-2039:

===========================  ======  ==============
pooling                      n       land p<0.05
===========================  ======  ==============
stack raw member values         240           14.4%
de-mean each member, stack      240           78.4%
ensemble-mean annual series       20           83.0%
===========================  ======  ==============

On a variable whose members share one scale (``driedarea``, binary, level ratio
0.56x) all three agree within a few points, so the choice only bites on the
continuous layers — where it bites hard. The ensemble mean also matches the
published ``median`` and ``trend``, which are themselves ensemble means, so the
sign of ``trend_tau`` cannot contradict the sign of ``trend``. Stacking members
would break that guarantee.

Why asymptotic rather than exact
--------------------------------
An exact tie-aware null distribution is enumerable only for very small n (the
null depends solely on the tie pattern, i.e. a partition of n, so n<=8 needs just
65 tables). At n=20-80 enumeration is impossible and unnecessary: the
tie-corrected normal approximation is well calibrated from about n>=10. Ties are
still handled exactly, through the standard variance correction — they are
common here (65.5% of ``driedarea`` cells and 43.8% of ``burntarea`` cells have a
tied decadal series), so ignoring them would misstate p on exactly the layers
that need it most.

What the p-value does and does not mean
---------------------------------------
It measures whether the ENSEMBLE-MEAN trajectory rises or falls monotonically.
It is **not** a statement of inter-model agreement: ``csoil-total`` has a sign
that is contested across models yet reaches p<0.05 on ~80% of land at ssp585,
because the ensemble mean can trend cleanly while its members disagree. Model
disagreement lives in ``lower_ci``/``upper_ci`` and ``n_models``. A report that
conflates the two will overclaim consensus.

Risk direction is NOT encoded here. ``trend_tau`` carries the direction of the
VALUE; whether a rising value is good or bad is the layer's
``percentile_direction`` attribute, which stays the single place that decision
lives.
"""

from __future__ import annotations

import warnings
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy import special

#: Fewest finite years for which a p-value is emitted. Below this the statistic
#: is degenerate (n=2 always yields p=1 regardless of the data), so it is NaN
#: rather than a number that looks like a result.
MIN_OBS = 3

#: Recorded verbatim into every backfilled layer as `significance_method`.
METHOD = ("mann_kendall_tie_corrected_asymptotic_two_sided_on_ensemble_mean_"
          "annual_series")


def _tie_terms(y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tie corrections for the MK variance and for tau-b, per cell.

    Args:
        y: ``(n, C)`` values with NaN for missing years.

    Returns:
        ``(n_finite, sum_u u(u-1)(2u+5), sum_u u(u-1))`` over the tie groups of
        each column's finite values, as float arrays of shape ``(C,)``.
    """
    n, C = y.shape
    finite = np.isfinite(y)
    n_fin = finite.sum(axis=0).astype(np.float64)

    # Sort with NaN pushed to the end so a tie run is a contiguous block; then
    # walk the axis maintaining the current run length. n is at most ~80, so a
    # Python loop over time with vector ops over cells is cheap and readable.
    s = np.sort(y, axis=0)
    s_fin = np.isfinite(s)

    t1 = np.zeros(C, dtype=np.float64)   # sum u(u-1)(2u+5)
    t3 = np.zeros(C, dtype=np.float64)   # sum u(u-1)
    run = s_fin[0].astype(np.float64)    # 0 where the column is entirely NaN

    def _flush(where: np.ndarray) -> None:
        """Close the current run on the selected columns and bank its terms."""
        u = np.where(where & (run > 1), run, 0.0)
        t1_add = u * (u - 1.0) * (2.0 * u + 5.0)
        t3_add = u * (u - 1.0)
        np.add(t1, t1_add, out=t1)
        np.add(t3, t3_add, out=t3)

    for j in range(1, n):
        same = s_fin[j] & (s[j] == s[j - 1])
        _flush(~same)
        run = np.where(same, run + 1.0, s_fin[j].astype(np.float64))
    _flush(np.ones(C, dtype=bool))

    return n_fin, t1, t3


def _kendall_S(y: np.ndarray) -> np.ndarray:
    """Kendall's S per cell for a strictly increasing, untied x axis.

    ``S = sum_{i<j} sign(y_j - y_i)``; NaN pairs contribute nothing. Because x is
    the year and is untied and increasing, ``sign(x_j - x_i)`` is always +1 and
    drops out of the sum.

    Args:
        y: ``(n, C)`` values with NaN for missing years.

    Returns:
        ``(C,)`` float array.
    """
    n, C = y.shape
    S = np.zeros(C, dtype=np.float64)
    for i in range(n - 1):
        d = y[i + 1:] - y[i]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            S += np.where(np.isfinite(d), np.sign(d), 0.0).sum(axis=0)
    return S


def mk_stats(y: np.ndarray, min_obs: int = MIN_OBS
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tie-corrected two-sided Mann-Kendall over axis 0, vectorized over cells.

    The x axis is implicit: consecutive, strictly increasing, untied (calendar
    years). Only y ties therefore enter the variance, which collapses the general
    correction to its first term.

    Args:
        y: ``(n, ...)`` values with NaN for missing years. Any trailing shape is
            flattened internally and restored on the way out.
        min_obs: Fewest finite values for which p and tau are emitted.

    Returns:
        ``(p, tau, n_obs)``, each float64 with y's trailing shape. ``p`` is the
        two-sided p-value with Kendall's continuity correction; ``tau`` is tau-b;
        ``n_obs`` is the count of finite years used. A perfectly constant series
        yields ``p=1.0, tau=0.0`` (defined, and honest: no evidence of a trend),
        not NaN. Cells with fewer than ``min_obs`` finite values yield NaN for p
        and tau, and their true count in ``n_obs``. Kept in float64 so the scalar
        path used for regional aggregation does not lose precision; callers
        writing a layer downcast to float32 to match the other value classes.
    """
    y = np.asarray(y, dtype=np.float64)
    n = y.shape[0]
    tail = y.shape[1:]
    if n == 0:
        # No time axis at all: nothing to rank. Return the right shapes rather
        # than letting reshape(0, -1) raise on an ambiguous dimension.
        return (np.full(tail, np.nan), np.full(tail, np.nan), np.zeros(tail))
    flat = y.reshape(n, -1)

    S = _kendall_S(flat)
    n_fin, t1, t3 = _tie_terms(flat)

    with np.errstate(invalid="ignore", divide="ignore"):
        var = (n_fin * (n_fin - 1.0) * (2.0 * n_fin + 5.0) - t1) / 18.0
        # tau-b: x is untied, so only the y-tie term reduces the denominator.
        n0 = n_fin * (n_fin - 1.0) / 2.0
        tau = S / np.sqrt(n0 * (n0 - t3 / 2.0))
        z = (S - np.sign(S)) / np.sqrt(var)
        p = special.erfc(np.abs(z) / np.sqrt(2.0))

    # var<=0 means every finite value is tied: S is exactly 0, there is no
    # ranking information, and the honest two-sided p is 1.
    degenerate = ~(var > 0)
    p = np.where(degenerate, 1.0, p)
    tau = np.where(degenerate, 0.0, tau)

    short = n_fin < min_obs
    p = np.where(short, np.nan, np.clip(p, 0.0, 1.0))
    tau = np.where(short, np.nan, tau)

    return p.reshape(tail), tau.reshape(tail), n_fin.reshape(tail)


def mk_pvalue(values: Sequence[float], min_obs: int = MIN_OBS
              ) -> Tuple[float, float, int]:
    """Scalar entry point: MK over one 1-D series of consecutive annual values.

    Used for the per-cell case in tests and, importantly, for REGIONAL
    aggregation: a polygon's p-value must be recomputed from its area-weighted
    mean annual series. Averaging per-cell p-values is meaningless, so downstream
    extraction must call this rather than aggregate the gridded ``trend_pvalue``.

    Args:
        values: Annual values in year order; NaN for missing years.
        min_obs: Fewest finite values for which p and tau are emitted.

    Returns:
        ``(p, tau, n_obs)`` as Python scalars.
    """
    arr = np.asarray(values, dtype=np.float64).reshape(-1, 1)
    p, tau, n = mk_stats(arr, min_obs=min_obs)
    return float(p[0]), float(tau[0]), int(n[0])


def mk_expanding(
    years: np.ndarray,
    annual: np.ndarray,
    decades: Sequence[int],
    window_years: int = 10,
    baseline_decade: Optional[int] = None,
    min_obs: int = MIN_OBS,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """MK over expanding windows: baseline start -> each decade's last year.

    Args:
        years: ``(n_year,)`` calendar years matching ``annual``'s first axis.
        annual: ``(n_year, lat, lon)`` ensemble-mean annual values, NaN off-mask.
        decades: Decade start years in ascending order, e.g. ``[2020, ..., 2090]``.
        window_years: Decade length, so decade ``d`` ends at ``d + window - 1``.
        baseline_decade: First decade of the tested record; defaults to
            ``decades[0]``. Its own panel is left NaN — no elapsed period.
        min_obs: Fewest finite years for which p and tau are emitted.

    Returns:
        ``(p, tau, n_obs)``, each ``(len(decades), lat, lon)`` float32. The
        baseline decade's panel is NaN in p and tau and 0 in n_obs.

    Raises:
        ValueError: If ``years`` is not strictly increasing (the untied,
            increasing x axis is what licenses the simplified variance), or if
            its length does not match ``annual``.
    """
    years = np.asarray(years)
    if years.shape[0] != annual.shape[0]:
        raise ValueError(
            f"years has {years.shape[0]} entries but annual has "
            f"{annual.shape[0]} time steps")
    if years.size > 1 and not np.all(np.diff(years) > 0):
        raise ValueError("years must be strictly increasing and untied")

    base = decades[0] if baseline_decade is None else baseline_decade
    shape = (len(decades),) + annual.shape[1:]
    p = np.full(shape, np.nan, np.float32)
    tau = np.full(shape, np.nan, np.float32)
    nobs = np.zeros(shape, np.float32)

    for i, d in enumerate(decades):
        if d == base:
            continue                      # no elapsed period to test
        sel = (years >= base) & (years <= d + window_years - 1)
        if not sel.any():
            continue
        p[i], tau[i], nobs[i] = mk_stats(annual[sel], min_obs=min_obs)
    return p, tau, nobs


def significance_definition(decades: Sequence[int], window_years: int = 10,
                            baseline_decade: Optional[int] = None) -> str:
    """Prose definition for the layer's ``significance_definition`` attribute."""
    base = decades[0] if baseline_decade is None else baseline_decade
    last = decades[-1] + window_years - 1
    return (
        f"trend_pvalue[decade] = two-sided Mann-Kendall p-value on the "
        f"ENSEMBLE-MEAN ANNUAL series over years {base}..(decade+{window_years - 1}), "
        f"i.e. an EXPANDING window anchored at the {base}s baseline: the "
        f"{base + 10}s panel tests {base}-{base + 19} (n=20) and the "
        f"{decades[-1]}s panel tests {base}-{last} (n={last - base + 1}). "
        f"Members are averaged WITHIN each year before the test, matching how "
        f"median and trend are built, so trend_tau cannot contradict trend; "
        f"stacking member-years instead would let between-model level offsets "
        f"dominate the sample (68.7x the interannual SD on csoil-total). "
        f"Variance is tie-corrected (asymptotic normal with Kendall's continuity "
        f"correction; ties are common and exact enumeration is infeasible at "
        f"n>=20). trend_tau is Kendall tau-b on the same series; trend_n_obs is "
        f"the number of finite years used. The {base}s panel is NaN (no elapsed "
        f"period). A perfectly constant series gives p=1.0, tau=0.0. "
        f"NOTE: this tests monotonicity of the ENSEMBLE-MEAN trajectory, NOT "
        f"inter-model agreement -- read it with lower_ci/upper_ci and n_models. "
        f"Risk direction requires combining sign(trend_tau) with the layer's "
        f"percentile_direction."
    )


# --- Theil-Sen slope on the same series -----------------------------------------
#
# Added 2026-07-30 (user decision) so that `trend` can be a FITTED SLOPE over the
# ensemble-mean annual series rather than the baseline-anchored two-point rate of
# GUARDRAILS S10. The point is coherence: trend, trend_tau and trend_pvalue then
# describe the *same* series over the *same* expanding window, so sign(trend) can
# only disagree with sign(trend_tau) if the trajectory is genuinely non-monotonic,
# and the p-value actually refers to the published slope.
#
# NOTE the divergence this creates: the nine layers published before this date use
# the S10 baseline-anchored rate. Until they are backfilled, `trend` does not mean
# the same thing across all layers -- see WORKFLOW-ISSUES.md 2026-07-30.


def theilsen(y: np.ndarray, min_obs: int = MIN_OBS) -> np.ndarray:
    """Theil-Sen slope over axis 0, vectorized over cells, NaN-aware.

    The median of all pairwise slopes ``(y_j - y_i) / (j - i)`` for ``i < j``. The x
    axis is implicit and consecutive (calendar years), so the denominators are just
    index gaps. Robust to outliers and to the heavy zero-inflation of count layers,
    where a least-squares slope is dragged by a few extreme years.

    Args:
        y: ``(n, ...)`` values, NaN for missing years.
        min_obs: Fewest finite values for which a slope is emitted.

    Returns:
        Slope in y-units **per single x step (per year)**, float64 with y's trailing
        shape. NaN where fewer than ``min_obs`` finite values exist. A constant
        series gives exactly 0.0.

    Memory: pairwise slopes are materialized in cell blocks, so peak use is bounded
    by ``block * n(n-1)/2`` floats regardless of how many cells are passed.
    """
    y = np.asarray(y, dtype=np.float64)
    n = y.shape[0]
    tail = y.shape[1:]
    flat = y.reshape(n, -1)
    ncell = flat.shape[1]
    out = np.full(ncell, np.nan)

    i_idx, j_idx = np.triu_indices(n, k=1)
    gaps = (j_idx - i_idx).astype(np.float64)
    npair = i_idx.size
    if npair == 0:
        return out.reshape(tail)

    # keep each block's pairwise buffer near ~200 MB
    block = max(1, int(2.5e7 // max(npair, 1)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN slices
        for s in range(0, ncell, block):
            e = min(s + block, ncell)
            chunk = flat[:, s:e]
            finite = np.isfinite(chunk).sum(axis=0)
            slopes = (chunk[j_idx, :] - chunk[i_idx, :]) / gaps[:, None]
            med = np.nanmedian(slopes, axis=0)
            med[finite < min_obs] = np.nan
            out[s:e] = med
            del slopes
    return out.reshape(tail)


def theilsen_expanding(
    years: np.ndarray,
    annual: np.ndarray,
    decades: Sequence[int],
    window_years: int = 10,
    baseline_decade: Optional[int] = None,
    min_obs: int = MIN_OBS,
    per_decade: bool = True,
) -> np.ndarray:
    """REJECTED for the published `trend` — see GUARDRAILS S10. Do not use for a layer.

    Fitting the ANNUAL series was attempted and rejected on 2026-07-30: Theil-Sen is a
    MEDIAN of pairwise slopes, so wherever more than half of year-pairs are tied the
    slope is **exactly 0** regardless of how much the quantity moved. On an annual
    hazard series most year-pairs are 0->0, so a low-frequency event whose frequency
    genuinely doubles (a 1-in-10-year event becoming 1-in-5) still reports zero trend.
    The windy-days layers are 49-97% zeros depending on threshold, which is precisely
    the failure case.

    Use :func:`theilsen_decadal` instead: it fits the DECADAL median series, where each
    point is already averaged over 10 years x every member, so ties are rare and a
    frequency change survives.

    Kept only because it is a faithful annual-series implementation and may be useful for
    diagnostics or for a genuinely continuous, non-zero-inflated field. Never wire it into
    a published `trend`. Note ``layer_publish.py``'s warning text still names this function;
    that message is stale relative to the rewritten S10.

    Theil-Sen slope over the SAME expanding windows as :func:`mk_expanding`.

    Window semantics are identical by construction — baseline decade's first year
    through the target decade's last year — so ``trend``, ``trend_tau`` and
    ``trend_pvalue`` are three statistics of one series over one window, and the
    baseline decade's own panel is NaN (no elapsed period).

    Args:
        years: ``(n_year,)`` calendar years matching ``annual``'s first axis.
        annual: ``(n_year, lat, lon)`` ensemble-mean annual values, NaN off-mask.
        decades: Decade start years, ascending, e.g. ``[2020, ..., 2090]``.
        window_years: Decade length; decade ``d`` ends at ``d + window_years - 1``.
        baseline_decade: First decade of the tested record; defaults to ``decades[0]``.
        min_obs: Fewest finite years for which a slope is emitted.
        per_decade: Return the slope per DECADE (per-year slope x window_years),
            matching the ``trend_units`` convention of the other layers. Set False
            for per-year units.

    Returns:
        ``(len(decades), lat, lon)`` float64 slopes, NaN in the baseline panel.
    """
    years = np.asarray(years)
    annual = np.asarray(annual, dtype=np.float64)
    base = decades[0] if baseline_decade is None else baseline_decade
    start = int(base)
    scale = float(window_years) if per_decade else 1.0

    out = np.full((len(decades), *annual.shape[1:]), np.nan)
    for k, d in enumerate(decades):
        if int(d) == int(base):
            continue                      # no elapsed period, same rule as mk_expanding
        end = int(d) + window_years - 1
        sel = (years >= start) & (years <= end)
        if sel.sum() < min_obs:
            continue
        out[k] = theilsen(annual[sel], min_obs=min_obs) * scale
    return out


def trend_definition_theilsen(decades: Sequence[int], window_years: int = 10,
                              baseline_decade: Optional[int] = None) -> str:
    """Provenance string for a Theil-Sen `trend`, for the layer's global attrs."""
    base = decades[0] if baseline_decade is None else baseline_decade
    last = int(decades[-1]) + window_years - 1
    return (
        f"Theil-Sen slope (median of pairwise slopes) of the ENSEMBLE-MEAN ANNUAL "
        f"series over an expanding window anchored at {int(base)}: each decade's "
        f"panel fits years {int(base)}-(decade end), so the final panel fits "
        f"{int(base)}-{last}. Units are value per decade. The {int(base)}s panel is "
        f"NaN (no elapsed period). Computed on the SAME series and windows as "
        f"trend_pvalue / trend_tau, so the three are mutually consistent. NOTE: this "
        f"is a FITTED SLOPE, not the GUARDRAILS S10 baseline-anchored two-point rate "
        f"used by layers published before 2026-07-30."
    )


# --- Building the ensemble-mean annual series inside a processor -----------------


class AnnualEnsembleMean:
    """Accumulate the ensemble-mean annual series one member at a time.

    A processor already streams each member once to build its decadal maps, but it
    discards the annual detail on the way — and `trend` (Theil-Sen, §10) and
    `trend_pvalue` (§15) both need the annual series. Holding every member's annual
    grid would cost ~1.4 GB for a 17-member layer, so this keeps a running sum and
    count instead: peak memory is two grids regardless of ensemble size.

    Add members in **sorted filename order**, per scenario. Floating-point addition
    is not associative, so the order fixes the last bits of the result; sorted order
    is what `scripts/backfill_trend_significance.py` uses, and matching it is what
    lets a processor's output be checked against an already-published layer.

    Example, inside a processor's member loop::

        acc = {s: AnnualEnsembleMean(MIN_YEAR, MAX_YEAR, (LAT, LON)) for s in scenarios}
        for f in sorted(files):
            da = load_member(f)
            acc[scenario_of(f)].add(da.year.values, da.values)
        years, mean_annual = acc[s].result()
        trend = theilsen_expanding(years, mean_annual, DECADES, window_years=WINDOW)
        p, tau, n = mk_expanding(years, mean_annual, DECADES, window_years=WINDOW)
    """

    def __init__(self, year0: int, year1: int, shape: Tuple[int, ...]):
        self.years = np.arange(int(year0), int(year1) + 1)
        self._shape = (len(self.years),) + tuple(shape)
        self._total = np.zeros(self._shape, np.float64)
        self._count = np.zeros(self._shape, np.float32)
        self.n_added = 0

    def add(self, years: np.ndarray, values: np.ndarray) -> None:
        """Add one member, aligning its years onto the accumulator's grid.

        Args:
            years: ``(n_year,)`` calendar years for ``values``.
            values: ``(n_year, *shape)`` annual values; NaN is skipped per cell, so a
                member that covers only part of the grid or only some years
                contributes exactly where it has data and nowhere else.

        Raises:
            ValueError: If ``values``' trailing shape does not match the grid, which
                would otherwise broadcast a wrong-sized member in silently.
        """
        years = np.asarray(years).astype(np.int64)
        values = np.asarray(values)
        if values.shape[1:] != self._shape[1:]:
            raise ValueError(
                f"member grid {values.shape[1:]} does not match accumulator "
                f"{self._shape[1:]}")

        aligned = np.full(self._shape, np.nan, np.float32)
        pos = np.searchsorted(self.years, years)
        inside = (pos >= 0) & (pos < len(self.years))
        pos_c = np.clip(pos, 0, len(self.years) - 1)
        keep = inside & (self.years[pos_c] == years)
        if keep.any():
            aligned[pos_c[keep]] = values[keep].astype(np.float32)

        good = np.isfinite(aligned)
        self._total[good] += aligned[good]
        self._count += good
        self.n_added += 1

    def result(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return ``(years, mean_annual)`` with NaN where no member contributed."""
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(self._count > 0,
                            self._total / np.maximum(self._count, 1), np.nan)
        return self.years, mean.astype(np.float32)

    @property
    def members_per_year(self) -> np.ndarray:
        """``(n_year, *shape)`` count of contributing members — read with the trend."""
        return self._count


# --- Theil-Sen on the DECADAL series (the shipped definition) --------------------
#
# WHY NOT THE ANNUAL SERIES (measured 2026-07-30, before shipping):
# Theil-Sen is the MEDIAN of pairwise slopes, so on a zero-inflated hazard, where
# more than half of all year-pairs are 0 -> 0, the median is EXACTLY ZERO however
# much the frequency actually moved. Share of cells with an exactly-zero slope at
# the 2090s, annual series vs decadal series:
#
#     layer                annual      decadal    (old anchored rate)
#     driedarea ssp126      91.3%        13.7%        4.5%
#     driedarea ssp585      56.2%        10.3%        3.6%
#     burntarea ssp585      14.4%        10.8%        9.8%
#     csoil     ssp585       3.7%         3.7%        0.4%
#
# On the annual series, 25.1% of driedarea ssp585 cells reported p<0.05 -- a
# significant trend -- with a slope of exactly 0. A decadal mean averages 10 years
# x every member, so exact zeros become rare and the pathology largely disappears;
# the residual sits close to the anchored rate's own zero share, i.e. cells that
# genuinely do not move. Continuous layers (csoil, timber) are unaffected either
# way, so the choice is made by the zero-inflated hazards.
#
# COST, accepted: `trend` is fitted on the decadal series while `trend_pvalue` is
# tested on the annual one, so they are no longer the same series. Measured sign
# agreement with tau is 87-97% rather than 100%. A coherent number that is wrong on
# 91% of cells is worse than a slightly less coherent one that is right.

#: Recorded as each layer's `trend_method`.
TREND_METHOD = "theil_sen_on_decadal_median_series"

#: A 2-decade window (the first post-baseline panel) has exactly one pairwise slope.
#: That is a well-defined rate — and equals the superseded anchored rate there — so
#: requiring 3 points would blank the first panel for no gain.
TREND_MIN_OBS = 2


def theilsen_decadal(
    medians: np.ndarray,
    decades: Sequence[int],
    window_years: int = 10,
    baseline_decade: Optional[int] = None,
    min_obs: int = TREND_MIN_OBS,
) -> np.ndarray:
    """Theil-Sen slope of the DECADAL median series over expanding windows.

    The x axis is the decade index, so the returned slope is already **per decade**
    — do not rescale it by ``window_years``. Getting that wrong inflates the trend
    tenfold, which is why this wrapper exists rather than callers reaching for
    :func:`theilsen_expanding` with a decade axis.

    Args:
        medians: ``(n_decade, lat, lon)`` published decadal medians, NaN off-mask.
        decades: Decade start years, ascending, matching ``medians``' first axis.
        window_years: Decade length, used only to bound each expanding window.
        baseline_decade: First decade of the window; defaults to ``decades[0]``.
        min_obs: Fewest finite decades for a slope. Default 2 — see TREND_MIN_OBS.

    Returns:
        ``(n_decade, lat, lon)`` float64 slopes in value per decade, NaN in the
        baseline panel (a fitted slope needs an elapsed period).

    Raises:
        ValueError: If ``decades`` does not match ``medians``' first axis.
    """
    medians = np.asarray(medians, dtype=np.float64)
    if len(decades) != medians.shape[0]:
        raise ValueError(
            f"{len(decades)} decades but medians has {medians.shape[0]} panels")
    base = decades[0] if baseline_decade is None else baseline_decade
    b_idx = list(decades).index(base)

    out = np.full(medians.shape, np.nan)
    for k in range(len(decades)):
        if k == b_idx:
            continue
        lo, hi = min(b_idx, k), max(b_idx, k)
        window = medians[lo:hi + 1]
        if window.shape[0] < min_obs:
            continue
        out[k] = theilsen(window, min_obs=min_obs)
    return out


def trend_definition_decadal(decades: Sequence[int], window_years: int = 10,
                             baseline_decade: Optional[int] = None) -> str:
    """Provenance string for the shipped Theil-Sen `trend`."""
    base = decades[0] if baseline_decade is None else baseline_decade
    return (
        f"Theil-Sen slope (median of pairwise slopes) of the DECADAL MEDIAN series "
        f"over an expanding window anchored at {int(base)}s: each panel fits the "
        f"decadal medians from {int(base)}s to that decade, so the final panel fits "
        f"all {len(decades)} decades. Units are value per decade (the x axis is the "
        f"decade index, so no rescaling is applied). The {int(base)}s panel is NaN "
        f"(no elapsed period) and the first post-baseline panel is the single "
        f"two-decade slope. Fitted on the DECADAL series, not the annual one, "
        f"because Theil-Sen on annual values returns exactly 0 wherever more than "
        f"half of all year-pairs are 0-to-0: on driedarea that was 91% of cells at "
        f"ssp126, and 25% of ssp585 cells reported a significant p-value alongside "
        f"a zero slope. trend_pvalue / trend_tau are tested on the ANNUAL series "
        f"(better powered), so trend and its p-value describe different series; "
        f"measured sign agreement is 87-97%. Supersedes the baseline-anchored "
        f"two-point rate used before 2026-07-30."
    )
