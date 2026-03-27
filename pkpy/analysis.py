"""
pkpy.analysis
-------------
Non-compartmental analysis (NCA) and compartmental model fitting.

NCA parameters computed:
  Cmax, Tmax, AUC_last, AUC_inf, t_half, clearance (CL), Vd (IV only)

Compartmental fitting uses scipy.optimize.curve_fit with sensible bounds
and returns best-fit parameters plus a summary DataFrame.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

from .models import (
    one_compartment_iv,
    one_compartment_oral,
    two_compartment_iv,
)


# ---------------------------------------------------------------------------
# Non-compartmental analysis
# ---------------------------------------------------------------------------

def nca(time, conc, dose=None, route="iv", extrapolate=True):
    """
    Perform non-compartmental analysis on concentration-time data.

    Parameters
    ----------
    time : array-like
        Time points (same units throughout).
    conc : array-like
        Observed concentrations (same units throughout).
    dose : float, optional
        Administered dose. Required for CL and Vd calculations.
    route : str
        'iv' or 'oral'. Affects Vd calculation.
    extrapolate : bool
        If True, extrapolate AUC to infinity using terminal slope.

    Returns
    -------
    dict
        NCA parameters with descriptive keys.
    """
    time = np.array(time, dtype=float)
    conc = np.array(conc, dtype=float)

    # Remove any NaN pairs
    mask = ~(np.isnan(time) | np.isnan(conc))
    time, conc = time[mask], conc[mask]

    # Sort by time
    order = np.argsort(time)
    time, conc = time[order], conc[order]

    results = {}

    # Cmax and Tmax
    idx_max = np.argmax(conc)
    results["Cmax"] = float(conc[idx_max])
    results["Tmax"] = float(time[idx_max])

    # AUC_last via linear-log trapezoidal
    auc_last = _auc_lin_log_trap(time, conc)
    results["AUC_last"] = auc_last

    # Terminal elimination rate constant (lambda_z) via log-linear regression
    # Use last 3+ points above LOQ (conc > 0)
    positive_mask = conc > 0
    if positive_mask.sum() >= 3:
        # Use the last half of positive points for terminal slope
        pos_idx = np.where(positive_mask)[0]
        terminal_idx = pos_idx[len(pos_idx) // 2:]
        t_term = time[terminal_idx]
        c_term = np.log(conc[terminal_idx])
        if len(t_term) >= 2:
            slope, intercept = np.polyfit(t_term, c_term, 1)
            lambda_z = -slope  # positive value
        else:
            lambda_z = np.nan
    else:
        lambda_z = np.nan

    results["lambda_z"] = lambda_z

    # t_half
    if lambda_z and lambda_z > 0:
        results["t_half"] = np.log(2) / lambda_z
    else:
        results["t_half"] = np.nan

    # AUC_inf
    if extrapolate and not np.isnan(lambda_z) and lambda_z > 0 and conc[-1] > 0:
        auc_extrap = conc[-1] / lambda_z
        results["AUC_inf"] = auc_last + auc_extrap
        results["AUC_extrap_pct"] = 100 * auc_extrap / results["AUC_inf"]
    else:
        results["AUC_inf"] = np.nan
        results["AUC_extrap_pct"] = np.nan

    # Clearance and Vd (require dose)
    if dose is not None and not np.isnan(results["AUC_inf"]) and results["AUC_inf"] > 0:
        if route == "iv":
            results["CL"] = dose / results["AUC_inf"]
            if not np.isnan(lambda_z) and lambda_z > 0:
                results["Vd"] = results["CL"] / lambda_z
            else:
                results["Vd"] = np.nan
        elif route == "oral":
            # F unknown — reports CL/F and Vd/F
            results["CL_F"] = dose / results["AUC_inf"]
            if not np.isnan(lambda_z) and lambda_z > 0:
                results["Vd_F"] = results["CL_F"] / lambda_z
            else:
                results["Vd_F"] = np.nan
    else:
        results["CL"] = np.nan
        results["Vd"] = np.nan

    return results


def _auc_lin_log_trap(time, conc):
    """
    Linear-log trapezoidal rule for AUC.
    Uses linear trapezoid when conc is increasing, log trapezoid when decreasing.
    """
    auc = 0.0
    for i in range(1, len(time)):
        dt = time[i] - time[i - 1]
        c1, c2 = conc[i - 1], conc[i]
        if c1 <= 0 or c2 <= 0:
            # Fall back to linear
            auc += dt * (c1 + c2) / 2
        elif c2 >= c1:
            # Increasing — linear
            auc += dt * (c1 + c2) / 2
        else:
            # Decreasing — log trapezoidal
            auc += dt * (c1 - c2) / np.log(c1 / c2)
    return auc


# ---------------------------------------------------------------------------
# Compartmental fitting
# ---------------------------------------------------------------------------

MODELS = {
    "1comp_iv": {
        "func": one_compartment_iv,
        "params": ["C0", "kel"],
        "p0": [1.0, 0.1],
        "bounds": ([0, 0], [np.inf, np.inf]),
    },
    "1comp_oral": {
        "func": one_compartment_oral,
        "params": ["F_dose_over_Vd", "ka", "kel"],
        "p0": [1.0, 1.0, 0.1],
        "bounds": ([0, 0, 0], [np.inf, np.inf, np.inf]),
    },
    "2comp_iv": {
        "func": two_compartment_iv,
        "params": ["A", "alpha", "B", "beta"],
        "p0": [1.0, 1.0, 0.5, 0.1],
        "bounds": ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]),
    },
}


def fit_model(time, conc, model="1comp_iv"):
    """
    Fit a compartmental PK model to concentration-time data.

    Parameters
    ----------
    time : array-like
        Time points.
    conc : array-like
        Observed concentrations.
    model : str
        One of '1comp_iv', '1comp_oral', '2comp_iv'.

    Returns
    -------
    dict with keys:
        'model'      : model name
        'params'     : dict of fitted parameter name → value
        'params_se'  : dict of parameter name → standard error
        'r_squared'  : goodness of fit
        'fitted'     : array of fitted concentrations at observed time points
        'time'       : original time array
        'observed'   : original concentration array
    """
    if model not in MODELS:
        raise ValueError(f"Unknown model '{model}'. Choose from: {list(MODELS.keys())}")

    spec = MODELS[model]
    time = np.array(time, dtype=float)
    conc = np.array(conc, dtype=float)

    # Remove NaNs
    mask = ~(np.isnan(time) | np.isnan(conc))
    time_clean, conc_clean = time[mask], conc[mask]

    try:
        popt, pcov = curve_fit(
            spec["func"],
            time_clean,
            conc_clean,
            p0=spec["p0"],
            bounds=spec["bounds"],
            maxfev=10000,
        )
    except RuntimeError as e:
        raise RuntimeError(f"Model fitting failed: {e}")

    perr = np.sqrt(np.diag(pcov))
    fitted = spec["func"](time_clean, *popt)

    # R-squared
    ss_res = np.sum((conc_clean - fitted) ** 2)
    ss_tot = np.sum((conc_clean - np.mean(conc_clean)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        "model": model,
        "params": dict(zip(spec["params"], popt)),
        "params_se": dict(zip(spec["params"], perr)),
        "r_squared": r2,
        "fitted": fitted,
        "time": time_clean,
        "observed": conc_clean,
    }


def fit_all_models(time, conc):
    """
    Fit all available models and return results ranked by R².

    Parameters
    ----------
    time : array-like
    conc : array-like

    Returns
    -------
    list of dicts, sorted by r_squared descending.
    """
    results = []
    for model_name in MODELS:
        try:
            result = fit_model(time, conc, model=model_name)
            results.append(result)
        except Exception as e:
            results.append({"model": model_name, "r_squared": np.nan, "error": str(e)})
    results.sort(key=lambda x: x.get("r_squared", -1) or -1, reverse=True)
    return results


def summary_table(nca_results, fit_result=None):
    """
    Return a formatted pandas DataFrame summarizing NCA and optional fit results.

    Parameters
    ----------
    nca_results : dict
        Output from nca().
    fit_result : dict, optional
        Output from fit_model().

    Returns
    -------
    pd.DataFrame
    """
    rows = []
    units_map = {
        "Cmax": "conc units",
        "Tmax": "time units",
        "AUC_last": "conc·time",
        "AUC_inf": "conc·time",
        "AUC_extrap_pct": "%",
        "t_half": "time units",
        "lambda_z": "1/time",
        "CL": "volume/time",
        "CL_F": "volume/time",
        "Vd": "volume",
        "Vd_F": "volume",
    }
    for k, v in nca_results.items():
        rows.append({
            "Parameter": k,
            "Value": round(v, 4) if not np.isnan(v) else "NaN",
            "Units": units_map.get(k, ""),
            "Source": "NCA",
        })

    if fit_result and "params" in fit_result:
        for k, v in fit_result["params"].items():
            rows.append({
                "Parameter": k,
                "Value": round(v, 4),
                "Units": "",
                "Source": f"Fit ({fit_result['model']})",
            })
        rows.append({
            "Parameter": "R²",
            "Value": round(fit_result["r_squared"], 4),
            "Units": "",
            "Source": f"Fit ({fit_result['model']})",
        })

    return pd.DataFrame(rows)
