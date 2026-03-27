"""
pkpy.plotting
-------------
Publication-quality concentration-time plots for PK analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


PKPY_STYLE = {
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.frameon": False,
}


def plot_concentration_time(
    time,
    conc,
    fit_result=None,
    title="Concentration-Time Profile",
    xlabel="Time",
    ylabel="Concentration",
    log_y=False,
    nca_results=None,
    save_path=None,
    figsize=(8, 5),
):
    """
    Plot observed concentration-time data with optional fitted curve.

    Parameters
    ----------
    time : array-like
    conc : array-like
    fit_result : dict, optional
        Output from analysis.fit_model(). Overlays fitted curve.
    title : str
    xlabel : str
    ylabel : str
    log_y : bool
        If True, use log scale on y-axis (semi-log plot).
    nca_results : dict, optional
        If provided, annotates Cmax and Tmax on the plot.
    save_path : str, optional
        File path to save figure (e.g. 'pk_plot.png').
    figsize : tuple

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    with plt.rc_context(PKPY_STYLE):
        fig, ax = plt.subplots(figsize=figsize)

        time = np.array(time)
        conc = np.array(conc)

        # Observed data
        ax.scatter(
            time, conc,
            color="#2E86AB", s=60, zorder=5,
            label="Observed", edgecolors="white", linewidths=0.5
        )
        ax.plot(time, conc, color="#2E86AB", alpha=0.4, linewidth=1, linestyle="--")

        # Fitted curve
        if fit_result and "params" in fit_result:
            from .models import one_compartment_iv, one_compartment_oral, two_compartment_iv
            model_funcs = {
                "1comp_iv": one_compartment_iv,
                "1comp_oral": one_compartment_oral,
                "2comp_iv": two_compartment_iv,
            }
            func = model_funcs.get(fit_result["model"])
            if func:
                t_smooth = np.linspace(time.min(), time.max(), 300)
                c_smooth = func(t_smooth, *fit_result["params"].values())
                ax.plot(
                    t_smooth, c_smooth,
                    color="#E84855", linewidth=2,
                    label=f"Fitted ({fit_result['model']}, R²={fit_result['r_squared']:.3f})"
                )

        # Annotate Cmax / Tmax
        if nca_results:
            cmax = nca_results.get("Cmax")
            tmax = nca_results.get("Tmax")
            if cmax and tmax:
                ax.axhline(cmax, color="#6B4226", linestyle=":", linewidth=1, alpha=0.7)
                ax.axvline(tmax, color="#6B4226", linestyle=":", linewidth=1, alpha=0.7)
                ax.annotate(
                    f"Cmax = {cmax:.2f}",
                    xy=(tmax, cmax),
                    xytext=(tmax + (time.max() * 0.05), cmax * 1.05),
                    fontsize=9, color="#6B4226",
                    arrowprops=dict(arrowstyle="->", color="#6B4226", lw=0.8),
                )

        if log_y:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(ticker.ScalarFormatter())

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, pad=12)
        ax.legend()
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig, ax


def plot_semilog(time, conc, fit_result=None, **kwargs):
    """Convenience wrapper for semi-log concentration-time plot."""
    return plot_concentration_time(time, conc, fit_result=fit_result, log_y=True, **kwargs)


def plot_residuals(fit_result, save_path=None, figsize=(7, 4)):
    """
    Plot weighted residuals from a compartmental fit.

    Parameters
    ----------
    fit_result : dict
        Output from analysis.fit_model().
    save_path : str, optional
    figsize : tuple

    Returns
    -------
    fig, ax
    """
    if "fitted" not in fit_result:
        raise ValueError("fit_result must contain 'fitted' key.")

    obs = fit_result["observed"]
    pred = fit_result["fitted"]
    time = fit_result["time"]
    residuals = obs - pred

    with plt.rc_context(PKPY_STYLE):
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(time, residuals, color="#F4A261", s=55, edgecolors="white", linewidths=0.5)
        ax.axhline(0, color="#333", linewidth=1)
        ax.set_xlabel("Time")
        ax.set_ylabel("Residual (Observed − Predicted)")
        ax.set_title(f"Residuals — {fit_result['model']}", pad=10)
        fig.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig, ax
