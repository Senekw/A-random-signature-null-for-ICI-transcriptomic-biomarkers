"""Cox proportional-hazards association of a score with overall survival.

Reported as a log-hazard-ratio per standard deviation of the score, so
that magnitudes are comparable across signatures whose raw scores are on
different scales.

Uses ``lifelines`` when available and falls back to a Newton-Raphson fit
of the Breslow partial likelihood otherwise -- the univariate case is
small enough that the fallback is exact to numerical tolerance, and it
keeps this step from being a hard dependency.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["cox_loghr_per_sd"]


def cox_loghr_per_sd(score, time, event) -> tuple:
    """Fit a univariate Cox model of ``score`` on OS.

    Parameters
    ----------
    score, time, event
        Pandas Series sharing a sample index. ``event`` is 1 for a death,
        0 for censored.

    Returns
    -------
    (log_hr_per_sd, se, p_value)
        NaNs when the model is not estimable (too few events, no variance
        in the score, no overlapping samples).
    """
    df = pd.DataFrame({
        "x": pd.to_numeric(pd.Series(score), errors="coerce"),
        "t": pd.to_numeric(pd.Series(time), errors="coerce"),
        "e": pd.to_numeric(pd.Series(event), errors="coerce"),
    }).dropna()

    df = df[df.t > 0]
    if df.shape[0] < 10 or df.e.sum() < 3 or df.x.std(ddof=1) == 0:
        return float("nan"), float("nan"), float("nan")

    df["x"] = (df.x - df.x.mean()) / df.x.std(ddof=1)   # per-SD units

    try:
        from lifelines import CoxPHFitter

        cph = CoxPHFitter()
        cph.fit(df[["t", "e", "x"]], duration_col="t", event_col="e")
        return (
            float(cph.params_["x"]),
            float(cph.standard_errors_["x"]),
            float(cph.summary.loc["x", "p"]),
        )
    except Exception:
        return _cox_newton(df.x.to_numpy(), df.t.to_numpy(), df.e.to_numpy())


def _cox_newton(x, t, e, max_iter: int = 50, tol: float = 1e-9) -> tuple:
    """Breslow partial-likelihood Newton-Raphson for a single covariate."""
    order = np.argsort(t)
    x, t, e = x[order], t[order], e[order]

    beta = 0.0
    for _ in range(max_iter):
        eta = beta * x
        w = np.exp(eta)

        # risk-set sums via reverse cumulative sums (ties handled by
        # Breslow's approximation)
        s0 = np.cumsum(w[::-1])[::-1]
        s1 = np.cumsum((w * x)[::-1])[::-1]
        s2 = np.cumsum((w * x * x)[::-1])[::-1]

        d = e == 1
        if not d.any() or np.any(s0[d] <= 0):
            return float("nan"), float("nan"), float("nan")

        grad = float((x[d] - s1[d] / s0[d]).sum())
        info = float(((s2[d] / s0[d]) - (s1[d] / s0[d]) ** 2).sum())
        if info <= 0:
            return float("nan"), float("nan"), float("nan")

        step = grad / info
        beta += step
        if abs(step) < tol:
            break

    se = float(np.sqrt(1.0 / info))
    from scipy.stats import norm
    p = float(2 * norm.sf(abs(beta / se)))
    return float(beta), se, p
