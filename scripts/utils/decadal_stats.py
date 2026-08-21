"""Shared decadal statistics for the TCFD/CDP product.

One contract for every hazard layer, so layers are mutually comparable. Replaces the
two divergent families that grew up in this codebase:

  * family A (the retired ``process_qg.py`` lineage, removed 2026-08-21) -- mean over models then
    MEDIAN over years; CI = interquartile range pooled over member x year; trend = an
    OLS slope fitted WITHIN each decade window.
  * family B (``process_csoil_soilcarbon.py``, ``process_burntarea_fire.py``,
    ``process_let_cyclone.py``, ``process_led_drought.py``) -- MEAN over years then MEAN
    over members (despite the output variable being named ``median``); CI = mean +/- 1
    inter-member SD; trend = a baseline-anchored two-point rate.

The contract implemented here (user decision, see WORKFLOW-ISSUES.md):

  ============  =====================================  ==============================
  variable      continuous field                       boolean {0,1} field
  ============  =====================================  ==============================
  median        MEDIAN over pooled (year x member)     MEAN over pooled (year x member)
  lower/upper   IQR: 25th / 75th pct of that pool      MEAN -/+ 1 SD of that pool
  ols_slope     OLS slope, stacked year x member       same
  sen_slope     Theil-Sen slope, stacked year x member same
  ============  =====================================  ==============================

``percentile`` is unchanged and stays in the individual processors -- it ranks against
the shared 2020s baseline distribution and needs the layer's ``percentile_direction``
and zero-inflation tier, neither of which belongs here.

Slope window: EXPANDING from the baseline decade through the target decade. The 2030s
panel fits years 2020-2039, the 2090s panel fits 2020-2099. The baseline panel itself
has no elapsed period and is **NaN** -- never 0. Writing 0 there makes the entire ocean
a finite zero, which the QA report does not catch because it only asserts that *finite*
baseline trends equal zero (this defect is live in the family-B processors).

Member pooling: every (year, member) observation is an independent point in the fit
(user decision). Note the measured hazard this accepts -- on ``csoil-total`` the
between-model level offset is 68.7x the interannual SD, so a pair drawn from two
different models at two different years mixes that offset into the slope. Fitting per
member and taking the median of slopes would be immune to it; stacking is what was
asked for and is what this module does.
"""

from __future__ import annotations

import warnings

import numpy as np

__all__ = [
    "is_boolean_field",
    "decade_window_mask",
    "pooled_decadal_stat",
    "expanding_slopes",
    "SlopeResult",
]

# Pairs of observations sharing a year have dx == 0 and no defined slope; they are
# dropped rather than treated as infinite.
_TIE_EPS = 1e-12

#: Minimum finite (year, member) observations before a slope is fitted at all. Below
#: this a slope is noise dressed as a trend, so the cell gets NaN.
MIN_SLOPE_OBS = 3


# ---------------------------------------------------------------------------
# Field classification
# ---------------------------------------------------------------------------

def is_boolean_field(values: np.ndarray, chunk: int = 4_000_000) -> bool:
    """True when the finite values are drawn only from {0, 1}.

    GUARDRAILS §9: a field's nature is established from its VALUES, never from its
    name, its ``long_name``, its ``units``, or a sibling variable. ``led`` is a binary
    per-cell flag while ``let`` is a continuous fraction in [0,1) -- same family, same
    naming convention, different nature -- so this must be measured per layer.

    A field that is *near*-boolean (a zero-inflated continuous fraction) is NOT boolean
    and must take the median/IQR branch; only an exactly-{0,1} field takes mean/SD.

    Parameters
    ----------
    values
        Any array of the raw member data. NaNs are ignored.
    chunk
        Block size for the scan. Bounds peak memory only -- every finite value is
        still inspected.

    Notes
    -----
    This is EXHAUSTIVE by design. An earlier version inspected a strided sample, which
    can step over a rare continuous tail and misclassify a zero-inflated fraction as
    boolean -- silently switching the layer from median/IQR to mean/SD. Since the
    branch this selects changes the published statistic, a sampled answer is not good
    enough. The scan short-circuits on the first offending value, so a genuinely
    continuous field costs one block.
    """
    arr = np.asarray(values).ravel()
    if arr.size == 0:
        return False

    seen_any = False
    for start in range(0, arr.size, chunk):
        blk = arr[start:start + chunk]
        finite = blk[np.isfinite(blk)]
        if finite.size == 0:
            continue
        seen_any = True
        if not np.all((finite == 0.0) | (finite == 1.0)):
            return False
    return seen_any


# ---------------------------------------------------------------------------
# Decade windows
# ---------------------------------------------------------------------------

def decade_window_mask(years: np.ndarray, decade: int, window_years: int = 10) -> np.ndarray:
    """Boolean mask over ``years`` selecting one decade window."""
    years = np.asarray(years)
    return (years >= decade) & (years <= decade + window_years - 1)


def _expanding_mask(years: np.ndarray, baseline_decade: int, decade: int,
                    window_years: int = 10) -> np.ndarray:
    """Mask from the START of the baseline decade through the END of ``decade``."""
    years = np.asarray(years)
    return (years >= baseline_decade) & (years <= decade + window_years - 1)


# ---------------------------------------------------------------------------
# Decadal central value + CI
# ---------------------------------------------------------------------------

def pooled_decadal_stat(annual: np.ndarray, years: np.ndarray, decade: int,
                        boolean: bool = False, window_years: int = 10,
                        central: str | None = None):
    """Central value and CI for one decade, pooling year x member into one sample.

    Parameters
    ----------
    annual
        ``(member, year, cell)`` array of annual values. ``cell`` may be a flattened
        land-cell axis or any trailing shape.
    years
        ``(year,)`` integer years matching axis 1.
    decade
        Decade start, e.g. ``2030``.
    boolean
        From :func:`is_boolean_field`. Selects mean/SD instead of median/IQR. Ignored
        when ``central`` is given explicitly.
    window_years
        Window length; 10 for a standard decade.
    central
        ``"median"`` (median + IQR) or ``"mean"`` (mean +/- 1 SD), overriding the
        ``boolean`` inference. Exists for the third branch described below; a layer
        that passes it MUST declare ``decadal_statistic`` accordingly.

    Returns
    -------
    (stat, lower, upper)
        Each shaped like ``annual.shape[2:]``. All-NaN cells give NaN, not 0.

    Notes
    -----
    **The third branch: extreme zero-inflation.** The boolean/continuous split is a
    proxy for the real question -- is the decade pool's distribution degenerate at
    zero? For most fields the proxy holds. It fails when a *continuous* field is so
    zero-inflated that its median is 0 almost everywhere, and then the median/IQR
    branch does not merely lose precision, it erases the layer.

    Measured on ``let`` (tropical-cyclone exposure) 2026-08-11, 2020s panel after
    spatial smoothing:

        ==========================  ==============  ==================
        branch                      exposed cells   exact-zero land
        ==========================  ==============  ==================
        pooled median + IQR              2,684            96.07%
        pooled mean +/- 1 SD            15,122            77.84%
        ==========================  ==============  ==================

    93% of exposed land reads exactly 0 under the median, and the two-tier percentile
    then assigns tier 1 to 96% of land. ``let`` is 97.84% exact-zero at annual
    resolution -- against 29.2% for ``burntarea``, which the median branch handles
    fine -- so this is a regime, not a slippery slope.

    ``let`` is the continuous analogue of ``led``: both estimate an EXPECTED ANNUAL
    EXPOSED FRACTION, and the spec already grants ``led`` the mean because a median of
    Bernoulli draws is meaningless. Passing ``central="mean"`` extends that same
    treatment to the continuous case (user decision 2026-08-11).

    This is a DECLARED deviation, not a silent one: such layers set
    ``decadal_statistic: pooled_mean_zero_inflated`` so a reader can tell which
    estimator produced the panel. Do NOT reach for it to "improve contrast" on an
    ordinary field -- measure the median's exact-zero share first, and record it.
    """
    annual = np.asarray(annual)
    if central is None:
        central = "mean" if boolean else "median"
    if central not in ("mean", "median"):
        raise ValueError(f"central must be 'mean' or 'median'; got {central!r}")
    mask = decade_window_mask(years, decade, window_years)
    if not mask.any():
        nan = np.full(annual.shape[2:], np.nan, np.float32)
        return nan, nan.copy(), nan.copy()

    # Pool member and year into a single sample axis.
    win = annual[:, mask, ...]
    pooled = win.reshape(-1, *win.shape[2:])

    all_nan = np.all(~np.isfinite(pooled), axis=0)

    # All-NaN cells (ocean, or a member's masked-out periphery) are expected and are
    # restored to NaN below; numpy's "All-NaN slice" notice is not actionable here.
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        if central == "mean":
            stat = np.nanmean(pooled, axis=0)
            sd = np.nanstd(pooled, axis=0)
            lower, upper = stat - sd, stat + sd
        else:
            stat = np.nanmedian(pooled, axis=0)
            lower = np.nanpercentile(pooled, 25, axis=0)
            upper = np.nanpercentile(pooled, 75, axis=0)

    for arr in (stat, lower, upper):
        arr[all_nan] = np.nan

    return (stat.astype(np.float32), lower.astype(np.float32), upper.astype(np.float32))


# ---------------------------------------------------------------------------
# Slopes
# ---------------------------------------------------------------------------

class SlopeResult(dict):
    """``{'ols_slope': array, 'sen_slope': array}`` with attribute access."""

    __getattr__ = dict.__getitem__


def _pair_indices(x: np.ndarray):
    """Upper-triangle index pairs with dx != 0, plus their dx.

    Computed once per panel and reused across every cell chunk -- the x layout is
    identical for all cells, so this is the expensive part done exactly once.
    """
    n = x.size
    i, j = np.triu_indices(n, k=1)
    dx = x[j] - x[i]
    keep = np.abs(dx) > _TIE_EPS
    return i[keep], j[keep], dx[keep]


#: Theil-Sen cost is quadratic in stacked points, and it is the ONLY expensive part of
#: this module (median/IQR and OLS are linear and effectively free). A 12-member,
#: 80-year expanding window is 960 points = 455,040 pairs PER CELL; over ~67k land cells
#: x 24 panels that measured 237 min/layer single-core (~24 min on 10 cores, 8.3 h for
#: 21 layers).
#:
#: Pair subsampling was evaluated as a speedup and is NOT enabled by default, because
#: its accuracy collapses on exactly the ensembles this pipeline has.
#:
#: The sample median is an APPROXIMATION of the exact Theil-Sen median, not an unbiased
#: estimator: the median is non-linear, so on a discrete, tied, or strongly asymmetric
#: slope distribution its expectation differs from the full-population median, and it
#: can land the far side of a large mass point (zero-inflated hazards have a big atom at
#: exactly 0). It converges as ``max_pairs`` rises, NOT as the ensemble grows with the
#: cap fixed.
#:
#: Crucially the error scales with the BETWEEN-MEMBER LEVEL OFFSET, because cross-member
#: pairs dominate the slope distribution and flatten it around the median. Measured max
#: abs error, true slope 0.010, interannual noise SD 0.5:
#:
#:     offset/noise   400,000     200,000     100,000      50,000
#:              0x   5.5e-05     1.7e-04     3.0e-04     4.4e-04
#:              1x   6.6e-05     2.4e-04     3.2e-04     5.4e-04
#:             10x   2.3e-04     1.2e-03     1.1e-03     2.0e-03
#:             68x   3.1e-04     1.8e-03     1.5e-03     3.0e-03
#:
#: ``csoil-total``'s measured ratio is 68.7x, where a 100k cap costs ~15% of the slope
#: -- unacceptable -- and the only tolerable cap (400k of 455k) buys just 1.14x. An
#: early benchmark suggesting ~2.8% error was run WITHOUT member offsets and was
#: therefore unrepresentative.
#:
#: The exact path meanwhile measures 21.9 min/layer on 10 cores (219 min single-core;
#: 7.7 h for a 21-layer backfill), which is inside the operational budget. So: exact by
#: default. ``max_pairs`` remains available to speed up iteration, and is safe on
#: layers whose members share a level, but production runs should leave it ``None``.
DEFAULT_MAX_PAIRS = None


def expanding_slopes(annual: np.ndarray, years: np.ndarray, decade: int,
                     baseline_decade: int, window_years: int = 10,
                     chunk_cells: int = 512,
                     max_pairs: int | None = DEFAULT_MAX_PAIRS,
                     rng_seed: int = 0) -> SlopeResult:
    """OLS and Theil-Sen slopes over an expanding window, stacked year x member.

    Every (year, member) observation is one point. Units are value per YEAR; multiply
    by 10 for value per decade if the layer declares decade units.

    The baseline decade returns all-NaN (no elapsed period).

    Parameters
    ----------
    annual
        ``(member, year, n_cells)``. Must be 3-D -- flatten spatial dims first so the
        chunking is over a single axis.
    max_pairs
        Cap on Theil-Sen pairs; see :data:`DEFAULT_MAX_PAIRS`. Sampling is deterministic
        (fixed seed), so reruns reproduce bit-identically. ``None`` forces the exact
        computation. Windows with fewer pairs than the cap are always exact.
    chunk_cells
        Cells per block. Peak memory is roughly ``4 x n_pairs x chunk_cells x itemsize``
        -- the two advanced-indexed operands ``sub[jj]`` and ``sub[ii]``, their
        difference, and nanmedian's internal partition copy all coexist, plus the
        float64 OLS temporaries. At the default cap and chunk size that is ~800 MB per
        worker on float32 input, so 10 parallel workers want ~8 GB. Lower
        ``chunk_cells`` if memory-bound.
    """
    # Do NOT coerce to float32 here. float32 spacing at 1e9 is ~64, so a stock field of
    # that magnitude carrying a real trend of ~1/yr would quantize to repeated values
    # and both slopes would read 0. Preserve the caller's precision through the fit and
    # cast only the returned slopes.
    annual = np.asarray(annual)
    if annual.dtype.kind != "f":
        annual = annual.astype(np.float64)
    if annual.ndim != 3:
        raise ValueError(f"expected (member, year, n_cells); got shape {annual.shape}")
    n_cells = annual.shape[2]

    if decade == baseline_decade:
        nan = np.full(n_cells, np.nan, np.float32)
        return SlopeResult(ols_slope=nan, sen_slope=nan.copy())

    mask = _expanding_mask(years, baseline_decade, decade, window_years)
    n_years = int(mask.sum())
    if n_years < 3:
        nan = np.full(n_cells, np.nan, np.float32)
        return SlopeResult(ols_slope=nan, sen_slope=nan.copy())

    win = annual[:, mask, :]                       # (member, year, cell)
    yrs = np.asarray(years)[mask].astype(np.float64)

    # Stack member x year into one observation axis, carrying x with it.
    n_mem = win.shape[0]
    y_all = win.reshape(n_mem * n_years, n_cells)  # (obs, cell)
    x_all = np.tile(yrs, n_mem)                    # (obs,)

    ols = np.full(n_cells, np.nan, np.float32)
    sen = np.full(n_cells, np.nan, np.float32)

    ii_all, jj_all, dx_all = _pair_indices(x_all)
    subsample = max_pairs is not None and ii_all.size > max_pairs

    # Keep the OLS regressor in float64: a float32 matmul emits spurious divide/overflow
    # FPE warnings from its vectorized remainder lanes, and the centred cross-product is
    # cheap enough that double precision costs nothing here.
    xc = (x_all - x_all.mean()).astype(np.float64)
    var_x = float(np.dot(xc, xc))
    if not np.isfinite(var_x) or var_x <= 0.0:
        nan = np.full(n_cells, np.nan, np.float32)
        return SlopeResult(ols_slope=nan, sen_slope=nan.copy())

    x64 = x_all.astype(np.float64)

    for c_idx, s in enumerate(range(0, n_cells, chunk_cells)):
        e = min(s + chunk_cells, n_cells)
        sub = y_all[:, s:e]                        # (obs, chunk)

        # Redraw the pair sample per chunk. Reusing ONE sample for the whole grid makes
        # the approximation error spatially coherent -- every cell inherits the same
        # over/under-representation of year-lags and member pairings, so the error can
        # look like a geographic signal instead of noise. Reseeding per chunk (and per
        # decade, so panels decorrelate too) breaks that up while staying exactly
        # reproducible. Measured cost: 0.06 s across 131 chunks, against ~1.8 s of real
        # work each.
        if subsample:
            rs = np.random.default_rng((rng_seed, decade, c_idx))
            pick = rs.choice(ii_all.size, size=max_pairs, replace=False)
            ii, jj, dx = ii_all[pick], jj_all[pick], dx_all[pick]
        else:
            ii, jj, dx = ii_all, jj_all, dx_all

        # A cell need NOT have a complete series. Ensembles have heterogeneous land
        # masks (csoil: 58,714-67,647 cells across 5 models), so requiring every
        # member-year to be finite would return NaN slopes for cells where `median`
        # and the CI are perfectly well defined -- ~13% of csoil land. Both estimators
        # below therefore use each cell's own finite subset.
        m = np.isfinite(sub)
        n_obs = m.sum(axis=0)
        if not (n_obs >= MIN_SLOPE_OBS).any():
            continue

        # --- OLS over each cell's finite subset.
        # Uses the sum form  slope = (n*Sxy - Sx*Sy) / (n*Sxx - Sx^2)  rather than a
        # shared centred regressor, because the valid subset -- and therefore the mean
        # of x -- differs per cell.
        y0 = np.where(m, sub, 0.0).astype(np.float64)
        mf = m.astype(np.float64)
        # Accelerate BLAS raises spurious FPE flags on all-finite matmuls (see the
        # note in the module docstring); the inputs here are finite by construction.
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            Sy = y0.sum(axis=0)
            Sx = x64 @ mf
            Sxy = x64 @ y0
            Sxx = (x64 * x64) @ mf
            denom = n_obs * Sxx - Sx * Sx
            ols_sub = (n_obs * Sxy - Sx * Sy) / denom
        ols_sub[~np.isfinite(denom) | (denom <= 0) | (n_obs < MIN_SLOPE_OBS)] = np.nan
        del y0, mf, Sy, Sx, Sxy, Sxx, denom

        # --- Theil-Sen: median of pairwise slopes, ignoring pairs touching a NaN.
        # A pair with a missing endpoint yields NaN and is skipped by nanmedian, which
        # is exactly the per-cell finite subset again. Differences are taken in the
        # input dtype (float64 for large-magnitude stocks) to avoid cancelling a small
        # trend against a large level.
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            slopes = (sub[jj, :] - sub[ii, :]) / dx[:, None]
        all_nan = np.all(~np.isfinite(slopes), axis=0)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-NaN slice encountered")
            sen_sub = np.nanmedian(slopes, axis=0)
        sen_sub[all_nan | (n_obs < MIN_SLOPE_OBS)] = np.nan
        del slopes

        ols[s:e] = ols_sub.astype(np.float32)
        sen[s:e] = sen_sub.astype(np.float32)

    return SlopeResult(ols_slope=ols, sen_slope=sen)
