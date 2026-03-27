"""
pkpy — Open-Source Pharmacokinetic Modeling & NCA Toolkit
==========================================================

A lightweight Python package for non-compartmental analysis (NCA) and
compartmental PK modeling — an accessible alternative to WinNonlin for
small biotech labs and academic researchers.

Quick start
-----------
>>> from pkpy import nca, fit_model, plot_concentration_time
>>> results = nca(time, conc, dose=10, route='iv')
>>> fit = fit_model(time, conc, model='1comp_iv')
>>> plot_concentration_time(time, conc, fit_result=fit, nca_results=results)
"""

from .analysis import nca, fit_model, fit_all_models, summary_table
from .plotting import plot_concentration_time, plot_semilog, plot_residuals
from .models import one_compartment_iv, one_compartment_oral, two_compartment_iv

__version__ = "0.1.0"
__author__ = "Deb Bhattacharyya"
__all__ = [
    "nca",
    "fit_model",
    "fit_all_models",
    "summary_table",
    "plot_concentration_time",
    "plot_semilog",
    "plot_residuals",
    "one_compartment_iv",
    "one_compartment_oral",
    "two_compartment_iv",
]
