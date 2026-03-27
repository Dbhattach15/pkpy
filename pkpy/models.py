"""
pkpy.models
-----------
Pharmacokinetic compartmental models for curve fitting.

Supported models:
  - One-compartment IV bolus
  - One-compartment oral (first-order absorption)
  - Two-compartment IV bolus
"""

import numpy as np


# ---------------------------------------------------------------------------
# One-compartment models
# ---------------------------------------------------------------------------

def one_compartment_iv(t, C0, kel):
    """
    One-compartment IV bolus model.

    C(t) = C0 * exp(-kel * t)

    Parameters
    ----------
    t : array-like
        Time points.
    C0 : float
        Initial concentration (at t=0).
    kel : float
        Elimination rate constant (1/time).

    Returns
    -------
    np.ndarray
        Predicted concentrations.
    """
    return C0 * np.exp(-kel * np.array(t))


def one_compartment_oral(t, F_dose_over_Vd, ka, kel):
    """
    One-compartment first-order oral absorption model.

    C(t) = (F*Dose/Vd) * (ka / (ka - kel)) * (exp(-kel*t) - exp(-ka*t))

    Parameters
    ----------
    t : array-like
        Time points.
    F_dose_over_Vd : float
        F*Dose/Vd combined scaling parameter.
    ka : float
        Absorption rate constant (1/time).
    kel : float
        Elimination rate constant (1/time).

    Returns
    -------
    np.ndarray
        Predicted concentrations.
    """
    t = np.array(t)
    if abs(ka - kel) < 1e-9:
        ka += 1e-9  # avoid division by zero
    return F_dose_over_Vd * (ka / (ka - kel)) * (np.exp(-kel * t) - np.exp(-ka * t))


# ---------------------------------------------------------------------------
# Two-compartment model
# ---------------------------------------------------------------------------

def two_compartment_iv(t, A, alpha, B, beta):
    """
    Two-compartment IV bolus model (bi-exponential).

    C(t) = A * exp(-alpha * t) + B * exp(-beta * t)

    Parameters
    ----------
    t : array-like
        Time points.
    A : float
        Coefficient for distribution phase.
    alpha : float
        Distribution rate constant (1/time).
    B : float
        Coefficient for elimination phase.
    beta : float
        Terminal elimination rate constant (1/time).

    Returns
    -------
    np.ndarray
        Predicted concentrations.
    """
    t = np.array(t)
    return A * np.exp(-alpha * t) + B * np.exp(-beta * t)
