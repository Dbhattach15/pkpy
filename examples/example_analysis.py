"""
example_analysis.py
-------------------
Full worked example using pkpy to analyze simulated in vitro ADME data.

Simulates a realistic oral PK profile, runs NCA, fits a 1-compartment
oral model, and generates publication-quality figures.

Run from the project root:
    python examples/example_analysis.py
"""

import numpy as np
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pkpy

# ---------------------------------------------------------------------------
# Simulate realistic oral PK data (e.g. CNS compound, rat PO dosing)
# ---------------------------------------------------------------------------
np.random.seed(42)

time_points = np.array([0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0])
dose = 10.0  # mg/kg

# True parameters: ka=1.2/h, kel=0.15/h, F*D/Vd=8.5
true_conc = pkpy.one_compartment_oral(time_points, 8.5, 1.2, 0.15)

# Add realistic log-normal noise (~15% CV)
noise = np.random.lognormal(mean=0, sigma=0.15, size=len(time_points))
observed_conc = true_conc * noise
observed_conc[0] = 0.0  # pre-dose = 0

print("=" * 60)
print("  pkpy — Pharmacokinetic Analysis Example")
print("  Compound: Simulated CNS oral compound")
print(f"  Dose: {dose} mg/kg PO")
print("=" * 60)

# ---------------------------------------------------------------------------
# Load into a DataFrame (mimics real ADME data import)
# ---------------------------------------------------------------------------
df = pd.DataFrame({"time_h": time_points, "conc_ng_mL": observed_conc})
print("\nRaw data:")
print(df.to_string(index=False))

# ---------------------------------------------------------------------------
# Non-compartmental analysis
# ---------------------------------------------------------------------------
print("\n--- Non-Compartmental Analysis (NCA) ---")
nca_results = pkpy.nca(
    df["time_h"].values,
    df["conc_ng_mL"].values,
    dose=dose,
    route="oral",
)
for k, v in nca_results.items():
    print(f"  {k:<20} {v:.4f}" if not np.isnan(float(v)) else f"  {k:<20} NaN")

# ---------------------------------------------------------------------------
# Compartmental model fitting
# ---------------------------------------------------------------------------
print("\n--- Compartmental Model Fitting ---")
fit = pkpy.fit_model(df["time_h"].values, df["conc_ng_mL"].values, model="1comp_oral")

print(f"  Model: {fit['model']}")
print(f"  R²:    {fit['r_squared']:.4f}")
print("  Parameters:")
for p, v in fit["params"].items():
    se = fit["params_se"][p]
    print(f"    {p:<20} {v:.4f}  ± {se:.4f}")

# ---------------------------------------------------------------------------
# Compare all models
# ---------------------------------------------------------------------------
print("\n--- Model Comparison ---")
all_fits = pkpy.fit_all_models(df["time_h"].values, df["conc_ng_mL"].values)
for r in all_fits:
    r2 = r.get("r_squared")
    r2_str = f"{r2:.4f}" if r2 is not None and not np.isnan(r2) else "failed"
    print(f"  {r['model']:<20} R² = {r2_str}")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
print("\n--- Summary Table ---")
table = pkpy.summary_table(nca_results, fit)
print(table.to_string(index=False))

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
print("\n--- Generating plots ---")

fig1, ax1 = pkpy.plot_concentration_time(
    df["time_h"].values,
    df["conc_ng_mL"].values,
    fit_result=fit,
    nca_results=nca_results,
    title="Simulated CNS Compound — Oral PK Profile",
    xlabel="Time (h)",
    ylabel="Concentration (ng/mL)",
    save_path="pk_profile.png",
)
print("  Saved: pk_profile.png")

fig2, ax2 = pkpy.plot_semilog(
    df["time_h"].values,
    df["conc_ng_mL"].values,
    fit_result=fit,
    title="Simulated CNS Compound — Semi-Log Plot",
    xlabel="Time (h)",
    ylabel="Concentration (ng/mL)",
    save_path="pk_semilog.png",
)
print("  Saved: pk_semilog.png")

fig3, ax3 = pkpy.plot_residuals(fit, save_path="pk_residuals.png")
print("  Saved: pk_residuals.png")

print("\nDone. pkpy analysis complete.")
