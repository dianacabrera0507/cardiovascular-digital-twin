"""Generate publication-ready figures and tables from validated local outputs.

This editorial revision preserves all frozen numerical results from analysis-v1.0.
It only improves labels, accessibility, and cross-feature comparability.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FEATURE_LABELS = {
    "median_mmhg": "Median ABP",
    "minimum_mmhg": "Minimum ABP",
    "maximum_mmhg": "Maximum ABP",
    "iqr_mmhg": "ABP IQR",
    "standard_deviation_mmhg": "ABP SD",
    "robust_amplitude_mmhg": "Robust amplitude",
    "map_slope_mmhg_per_minute": "ABP slope",
    "heart_rate_median_bpm": "Heart rate",
    "pulse_pressure_median_mmhg": "Pulse pressure",
    "tau_acceptance_fraction": "Accepted tau fraction",
    "tau_median_s": "Effective tau",
}
CONTROL_COLOR = "#0072B2"
CASE_COLOR = "#D55E00"
EFFECT_COLOR = "#1F4E79"


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required outputs: " + "; ".join(missing))


def save_figure(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def cohort_flow_figure(output_dir: Path) -> pd.DataFrame:
    rows = [
        ("Patients with physical ABP segments", 64),
        ("Patients with >=1 ABP segment lasting >=25 min", 54),
        ("Patients reviewed in event audit", 52),
        ("Patients with eligible hypotension event", 26),
        ("Patients with non-overlapping case-control pair", 23),
        ("Pairs with effective tau in both windows", 15),
    ]
    table = pd.DataFrame(rows, columns=["stage", "patients"])
    fig, ax = plt.subplots(figsize=(9.2, 6.7))
    ax.axis("off")
    colors = plt.cm.Blues(np.linspace(0.35, 0.85, len(table)))
    y_positions = np.arange(len(table))[::-1]
    widths = 0.42 + 0.48 * table["patients"].to_numpy(float) / table["patients"].max()
    for index, (_, row) in enumerate(table.iterrows()):
        y = y_positions[index]
        width = widths[index]
        rectangle = plt.Rectangle(
            (0.5 - width / 2, y - 0.32), width, 0.64,
            transform=ax.transData, facecolor=colors[index], edgecolor="black",
            linewidth=0.8,
        )
        ax.add_patch(rectangle)
        ax.text(0.5, y + 0.07, row["stage"], ha="center", va="center", fontsize=9)
        ax.text(0.5, y - 0.13, f"n = {int(row['patients'])}", ha="center", va="center", fontsize=10, fontweight="bold")
        if index < len(table) - 1:
            ax.annotate("", xy=(0.5, y - 0.66), xytext=(0.5, y - 0.34), arrowprops={"arrowstyle": "->", "lw": 1.0})
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, len(table) - 0.2)
    ax.set_title("Study cohort and paired analysis flow", fontsize=13, pad=12)
    save_figure(fig, output_dir / "figure_1_cohort_flow")
    return table


def paired_tau_figure(features: pd.DataFrame, output_dir: Path) -> None:
    paired = features.pivot(index="patient_code", columns="window_type", values="tau_median_s").dropna()
    fig, ax = plt.subplots(figsize=(6.4, 6.1))
    for _, row in paired.iterrows():
        ax.plot([0, 1], [row["control"], row["case"]], color="0.70", linewidth=1, alpha=0.85)
        ax.scatter([0, 1], [row["control"], row["case"]], color=[CONTROL_COLOR, CASE_COLOR], s=28, zorder=3)
    medians = [paired["control"].median(), paired["case"].median()]
    ax.plot([0, 1], medians, color="black", linewidth=3, marker="o", markersize=7, label="Group median")
    ax.set_xticks([0, 1], ["Control", "Pre-hypotension"])
    ax.set_ylabel("Effective diastolic tau (s)")
    ax.set_title(f"Paired effective diastolic time constant before hypotension (n = {len(paired)})")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.10))
    fig.subplots_adjust(bottom=0.18)
    save_figure(fig, output_dir / "figure_2_paired_tau")


def standardized_effect_figure(statistics: pd.DataFrame, output_dir: Path) -> None:
    plot = statistics.copy()
    plot["label"] = plot["feature"].map(FEATURE_LABELS).fillna(plot["feature"])
    plot = plot.sort_values("rank_biserial_correlation").reset_index(drop=True)
    values = plot["rank_biserial_correlation"].to_numpy(float)
    colors = np.where(values < 0, CASE_COLOR, CONTROL_COLOR)
    fig, ax = plt.subplots(figsize=(8.7, 6.8))
    y = np.arange(len(plot))
    ax.barh(y, values, color=colors, alpha=0.90)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlim(-1, 1)
    ax.set_yticks(y, plot["label"])
    ax.set_xlabel("Paired rank-biserial correlation")
    ax.set_title("Standardized paired effects: pre-hypotension versus control")
    ax.grid(axis="x", alpha=0.25)
    ax.text(-0.98, len(plot) - 0.15, "Lower in pre-hypotension", ha="left", va="bottom", fontsize=9, color=CASE_COLOR)
    ax.text(0.98, len(plot) - 0.15, "Higher in pre-hypotension", ha="right", va="bottom", fontsize=9, color=CONTROL_COLOR)
    for yi, value in zip(y, values):
        offset = 0.025 if value >= 0 else -0.025
        ax.text(value + offset, yi, f"{value:.2f}", ha="left" if value >= 0 else "right", va="center", fontsize=8)
    save_figure(fig, output_dir / "figure_3_standardized_effects")


def sensitivity_figure(sensitivity: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
    for color, (pressure, group) in zip(palette, sensitivity.groupby("venous_pressure_mmhg")):
        ordered = group.sort_values("fit_start_fraction")
        ax.plot(
            ordered["fit_start_fraction"] * 100,
            ordered["patient_tau_median_s"],
            marker="o", linewidth=1.8, color=color,
            label=f"Venous pressure = {pressure:g} mmHg",
        )
    ax.set_xlabel("Diastolic fit start (%)")
    ax.set_ylabel("Median patient-level effective tau (s)")
    ax.set_title("Sensitivity of effective tau to venous pressure and fit-start assumptions")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    save_figure(fig, output_dir / "figure_s1_tau_sensitivity")


def create_tables(statistics: pd.DataFrame, event_segments: pd.DataFrame, event_table: pd.DataFrame, features: pd.DataFrame, quality: pd.DataFrame, output_dir: Path) -> None:
    main = statistics.copy()
    main.insert(1, "feature_label", main["feature"].map(FEATURE_LABELS).fillna(main["feature"]))
    main.to_csv(output_dir / "table_1_paired_feature_results.csv", index=False)
    safe = event_table[event_table["leakage_safe_history"].fillna(False).astype(bool)]
    summary = pd.DataFrame([
        {"metric": "Patients reviewed in event audit", "value": int(event_segments["patient_code"].nunique())},
        {"metric": "Segments reviewed", "value": int(len(event_segments))},
        {"metric": "Valid hypotension events", "value": int(len(event_table))},
        {"metric": "Leakage-safe events", "value": int(len(safe))},
        {"metric": "Patients with leakage-safe event", "value": int(safe["patient_code"].nunique())},
        {"metric": "Paired patients", "value": int(features["patient_code"].nunique())},
        {"metric": "Pairs with tau in both windows", "value": int(features.pivot(index="patient_code", columns="window_type", values="tau_median_s").dropna().shape[0])},
        {"metric": "Median valid MAP fraction", "value": float(event_segments["finite_map_fraction"].median())},
    ])
    summary.to_csv(output_dir / "table_2_cohort_event_summary.csv", index=False)
    quality.to_csv(output_dir / "table_3_quality_exclusions.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate editorially revised final figures and tables.")
    parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    paths = {
        "statistics": root / "outputs" / "paired_statistical_analysis" / "paired_feature_statistics.csv",
        "features": root / "outputs" / "hypotension_case_control_features" / "paired_abp_features.csv",
        "sensitivity": root / "outputs" / "abp_sensitivity" / "sensitivity_configuration_summary.csv",
        "segments": root / "outputs" / "hypotension_event_audit" / "hypotension_segment_audit.csv",
        "events": root / "outputs" / "hypotension_event_audit" / "hypotension_events.csv",
        "quality": root / "outputs" / "abp_patient_statistics" / "pv05_start35" / "quality_exclusion_categories.csv",
    }
    require_files(list(paths.values()))
    statistics = pd.read_csv(paths["statistics"])
    features = pd.read_csv(paths["features"])
    sensitivity = pd.read_csv(paths["sensitivity"])
    event_segments = pd.read_csv(paths["segments"])
    event_table = pd.read_csv(paths["events"])
    quality = pd.read_csv(paths["quality"])
    output_dir = root / "outputs" / "final_results_artifacts_editorial"
    output_dir.mkdir(parents=True, exist_ok=True)
    flow = cohort_flow_figure(output_dir)
    paired_tau_figure(features, output_dir)
    standardized_effect_figure(statistics, output_dir)
    sensitivity_figure(sensitivity, output_dir)
    create_tables(statistics, event_segments, event_table, features, quality, output_dir)
    flow.to_csv(output_dir / "figure_1_cohort_flow_data.csv", index=False)
    print("Artefactos editoriales generados.")
    print("Figuras PNG:", len(list(output_dir.glob("*.png"))))
    print("Figuras PDF:", len(list(output_dir.glob("*.pdf"))))
    print("Tablas CSV:", len(list(output_dir.glob("*.csv"))))
    print("Resultados locales:", output_dir)


if __name__ == "__main__":
    main()
