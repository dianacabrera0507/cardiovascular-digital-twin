"""Generate reproducible patient-level ABP statistics from local sensitivity outputs.

This script reads the prespecified baseline configuration (Pv=5 mmHg,
fit start=35%), summarizes effective diastolic tau at patient level, computes
a percentile bootstrap confidence interval for the median, calculates a Wilson
confidence interval for estimability, and describes quality-based exclusions.
All generated outputs remain under outputs/, which is excluded from Git.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


QUALITY_COLUMNS = [
    "finite_fraction",
    "median_mmhg",
    "p05_mmhg",
    "p95_mmhg",
    "robust_amplitude_mmhg",
]


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Return the Wilson score interval for a binomial proportion."""
    if total <= 0:
        raise ValueError("total must be positive")
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    half_width = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / total
            + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return center - half_width, center + half_width


def bootstrap_median(
    values: np.ndarray,
    repetitions: int,
    random_seed: int,
) -> tuple[float, float, float]:
    """Return percentile CI limits and bootstrap standard error for the median."""
    if values.size == 0:
        raise ValueError("at least one patient value is required")
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    rng = np.random.default_rng(random_seed)
    sampled = rng.choice(values, size=(repetitions, values.size), replace=True)
    bootstrap_values = np.median(sampled, axis=1)
    lower, upper = np.quantile(bootstrap_values, [0.025, 0.975])
    standard_error = np.std(bootstrap_values, ddof=1)
    return float(lower), float(upper), float(standard_error)


def exclusion_categories(rejected: pd.DataFrame) -> pd.DataFrame:
    """Summarize overlapping quality-failure categories by window and patient."""
    definitions = {
        "missingness": rejected["finite_fraction"] < 0.99,
        "flat_or_near_flat": rejected["robust_amplitude_mmhg"] < 10.0,
        "low_or_zero_pressure": (
            (rejected["median_mmhg"] < 35.0) | (rejected["p05_mmhg"] < 15.0)
        ),
        "extreme_high_pressure": (
            (rejected["median_mmhg"] > 180.0) | (rejected["p95_mmhg"] > 250.0)
        ),
        "excessive_amplitude": rejected["robust_amplitude_mmhg"] > 120.0,
    }
    rows = []
    for name, mask in definitions.items():
        rows.append(
            {
                "category": name,
                "windows_affected": int(mask.sum()),
                "window_percent": float(100.0 * mask.mean()) if len(mask) else 0.0,
                "patients_affected": int(rejected.loc[mask, "patient_code"].nunique()),
            }
        )
    return pd.DataFrame(rows)


def summarize(
    input_dir: Path,
    output_dir: Path,
    bootstrap_repetitions: int,
    random_seed: int,
) -> None:
    """Create patient-level summaries and uncertainty estimates."""
    segment_file = input_dir / "abp_segments_summary.csv"
    window_file = input_dir / "abp_windows_summary.csv"
    beat_file = input_dir / "abp_beats_tau.csv"
    for path in (segment_file, window_file, beat_file):
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path.name}")

    segments = pd.read_csv(segment_file)
    windows = pd.read_csv(window_file)
    beats = pd.read_csv(beat_file)
    required_segment = {"patient_code", "segment_code"}
    required_window = {
        "patient_code",
        "segment_code",
        "analysis_status",
        "tau_median_s",
        *QUALITY_COLUMNS,
    }
    if not required_segment.issubset(segments.columns):
        raise ValueError("segment summary does not contain required columns")
    if not required_window.issubset(windows.columns):
        raise ValueError("window summary does not contain required columns")

    tau_windows = windows[windows["analysis_status"].eq("tau_available")].copy()
    patient_tau = (
        tau_windows.groupby("patient_code", as_index=False)
        .agg(
            windows_with_tau=("tau_median_s", "size"),
            tau_median_s=("tau_median_s", "median"),
            tau_q25_s=("tau_median_s", lambda x: x.quantile(0.25)),
            tau_q75_s=("tau_median_s", lambda x: x.quantile(0.75)),
            tau_min_s=("tau_median_s", "min"),
            tau_max_s=("tau_median_s", "max"),
        )
    )
    patient_tau["tau_range_s"] = patient_tau["tau_max_s"] - patient_tau["tau_min_s"]

    all_patients = pd.DataFrame(
        {"patient_code": sorted(segments["patient_code"].dropna().unique())}
    )
    patient_summary = all_patients.merge(patient_tau, on="patient_code", how="left")
    patient_summary["has_acceptable_tau"] = patient_summary["tau_median_s"].notna()

    total_patients = int(patient_summary["patient_code"].nunique())
    patients_with_tau = int(patient_summary["has_acceptable_tau"].sum())
    tau_values = patient_summary.loc[
        patient_summary["has_acceptable_tau"], "tau_median_s"
    ].to_numpy(dtype=float)
    ci_low, ci_high, bootstrap_se = bootstrap_median(
        tau_values,
        repetitions=bootstrap_repetitions,
        random_seed=random_seed,
    )
    wilson_low, wilson_high = wilson_interval(patients_with_tau, total_patients)

    patients_without_tau = set(
        patient_summary.loc[~patient_summary["has_acceptable_tau"], "patient_code"]
    )
    missing_windows = windows[windows["patient_code"].isin(patients_without_tau)].copy()
    rejected = missing_windows[missing_windows["analysis_status"].eq("quality_rejected")].copy()
    quality_exclusions = exclusion_categories(rejected)
    status_counts = (
        missing_windows["analysis_status"]
        .value_counts(dropna=False)
        .rename_axis("analysis_status")
        .reset_index(name="windows")
    )

    accepted_beats = beats[
        beats.get("acceptable_fit", pd.Series(False, index=beats.index))
        .fillna(False)
        .astype(bool)
    ]
    metrics = {
        "baseline_configuration": "pv05_start35",
        "patients_reviewed": total_patients,
        "patients_with_tau": patients_with_tau,
        "patients_without_tau": total_patients - patients_with_tau,
        "estimability_percent": 100.0 * patients_with_tau / total_patients,
        "estimability_wilson_95_low_percent": 100.0 * wilson_low,
        "estimability_wilson_95_high_percent": 100.0 * wilson_high,
        "patient_tau_median_s": float(np.median(tau_values)),
        "patient_tau_q25_s": float(np.quantile(tau_values, 0.25)),
        "patient_tau_q75_s": float(np.quantile(tau_values, 0.75)),
        "patient_tau_bootstrap_95_low_s": ci_low,
        "patient_tau_bootstrap_95_high_s": ci_high,
        "patient_tau_bootstrap_standard_error_s": bootstrap_se,
        "bootstrap_repetitions": bootstrap_repetitions,
        "random_seed": random_seed,
        "windows_total": int(len(windows)),
        "windows_with_tau": int(len(tau_windows)),
        "accepted_beats": int(len(accepted_beats)),
        "patients_with_internal_range_over_1s": int(
            (patient_tau["tau_range_s"] > 1.0).sum()
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    patient_summary.to_csv(output_dir / "patient_tau_summary.csv", index=False)
    status_counts.to_csv(output_dir / "patients_without_tau_status_counts.csv", index=False)
    quality_exclusions.to_csv(output_dir / "quality_exclusion_categories.csv", index=False)
    pd.DataFrame([metrics]).to_csv(output_dir / "patient_level_metrics.csv", index=False)
    (output_dir / "patient_level_metrics.json").write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("Resumen estadistico por paciente terminado.")
    print("Pacientes revisados:", total_patients)
    print("Pacientes con tau:", patients_with_tau)
    print("Pacientes sin tau:", total_patients - patients_with_tau)
    print("Mediana de tau:", round(metrics["patient_tau_median_s"], 3), "s")
    print(
        "IC 95% bootstrap:",
        round(ci_low, 3),
        "-",
        round(ci_high, 3),
        "s",
    )
    print(
        "Estimabilidad:",
        round(metrics["estimability_percent"], 1),
        "%",
    )
    print(
        "IC 95% Wilson:",
        round(metrics["estimability_wilson_95_low_percent"], 1),
        "-",
        round(metrics["estimability_wilson_95_high_percent"], 1),
        "%",
    )
    print("Resultados locales:", output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize patient-level ABP tau statistics."
    )
    parser.add_argument(
        "--input-label",
        default="pv05_start35",
        help="Sensitivity output label to summarize.",
    )
    parser.add_argument("--bootstrap-repetitions", type=int, default=10000)
    parser.add_argument("--random-seed", type=int, default=2026)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.bootstrap_repetitions <= 0:
        raise ValueError("--bootstrap-repetitions must be positive")
    if not args.input_label.replace("_", "").replace("-", "").isalnum():
        raise ValueError("--input-label contains unsupported characters")
    project_root = Path(__file__).resolve().parents[1]
    input_dir = project_root / "outputs" / "abp_sensitivity" / args.input_label
    output_dir = project_root / "outputs" / "abp_patient_statistics" / args.input_label
    summarize(
        input_dir=input_dir,
        output_dir=output_dir,
        bootstrap_repetitions=args.bootstrap_repetitions,
        random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()
