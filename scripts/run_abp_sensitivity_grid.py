"""Run the remaining ABP sensitivity grid and consolidate patient-level results."""
from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path

import pandas as pd


VENOUS_PRESSURES = (0, 5, 10, 15)
FIT_START_FRACTIONS = (0.25, 0.35, 0.45, 0.55)
BASELINE = (5, 0.35)


def label_for(venous_pressure: int, fit_fraction: float) -> str:
    return f"pv{venous_pressure:02d}_start{int(round(100 * fit_fraction)):02d}"


def run_configuration(project_root: Path, venous_pressure: int, fit_fraction: float) -> None:
    label = label_for(venous_pressure, fit_fraction)
    output_dir = project_root / "outputs" / "abp_sensitivity" / label
    required = (
        output_dir / "abp_segments_summary.csv",
        output_dir / "abp_windows_summary.csv",
        output_dir / "abp_beats_tau.csv",
    )
    if all(path.exists() and path.stat().st_size > 0 for path in required):
        print(f"OMITIDA {label}: resultados existentes")
        return
    command = [
        sys.executable,
        str(project_root / "scripts" / "run_abp_multwindow_analysis.py"),
        "--maximum-segments", "100",
        "--maximum-segments-per-patient", "2",
        "--venous-pressure-mmhg", str(venous_pressure),
        "--fit-start-fraction", str(fit_fraction),
        "--output-label", label,
    ]
    print(f"EJECUTANDO {label}")
    subprocess.run(command, cwd=project_root, check=True)


def consolidate(project_root: Path) -> None:
    root = project_root / "outputs" / "abp_sensitivity"
    configuration_rows: list[dict[str, object]] = []
    patient_rows: list[pd.DataFrame] = []
    for venous_pressure, fit_fraction in itertools.product(
        VENOUS_PRESSURES, FIT_START_FRACTIONS
    ):
        label = label_for(venous_pressure, fit_fraction)
        directory = root / label
        segment_file = directory / "abp_segments_summary.csv"
        window_file = directory / "abp_windows_summary.csv"
        beat_file = directory / "abp_beats_tau.csv"
        if not all(path.exists() for path in (segment_file, window_file, beat_file)):
            raise FileNotFoundError(f"Resultados incompletos para {label}")
        segments = pd.read_csv(segment_file)
        windows = pd.read_csv(window_file)
        beats = pd.read_csv(beat_file)
        tau_windows = windows[windows["analysis_status"].eq("tau_available")].copy()
        patients = (
            tau_windows.groupby("patient_code", as_index=False)["tau_median_s"]
            .median()
            .rename(columns={"tau_median_s": "patient_tau_median_s"})
        )
        patients["configuration"] = label
        patients["venous_pressure_mmhg"] = venous_pressure
        patients["fit_start_fraction"] = fit_fraction
        patient_rows.append(patients)
        accepted = beats[beats["acceptable_fit"].fillna(False).astype(bool)]
        patient_tau = patients["patient_tau_median_s"]
        configuration_rows.append({
            "configuration": label,
            "venous_pressure_mmhg": venous_pressure,
            "fit_start_fraction": fit_fraction,
            "patients_reviewed": int(segments["patient_code"].nunique()),
            "patients_with_tau": int(len(patients)),
            "windows_total": int(len(windows)),
            "windows_with_tau": int(len(tau_windows)),
            "accepted_beats": int(len(accepted)),
            "patient_tau_median_s": float(patient_tau.median()) if len(patient_tau) else float("nan"),
            "patient_tau_q25_s": float(patient_tau.quantile(0.25)) if len(patient_tau) else float("nan"),
            "patient_tau_q75_s": float(patient_tau.quantile(0.75)) if len(patient_tau) else float("nan"),
        })
    summary = pd.DataFrame(configuration_rows)
    patient_table = pd.concat(patient_rows, ignore_index=True)
    baseline_label = label_for(*BASELINE)
    baseline = (
        patient_table[patient_table["configuration"].eq(baseline_label)]
        [["patient_code", "patient_tau_median_s"]]
        .rename(columns={"patient_tau_median_s": "baseline_tau_median_s"})
    )
    paired = patient_table.merge(baseline, on="patient_code", how="inner")
    paired["absolute_change_s"] = (
        paired["patient_tau_median_s"] - paired["baseline_tau_median_s"]
    )
    paired["absolute_percent_change"] = (
        paired["absolute_change_s"].abs() / paired["baseline_tau_median_s"] * 100.0
    )
    stability = (
        paired.groupby("configuration", as_index=False)
        .agg(
            paired_patients=("patient_code", "nunique"),
            median_absolute_change_s=("absolute_change_s", lambda value: value.abs().median()),
            median_absolute_percent_change=("absolute_percent_change", "median"),
        )
    )
    summary = summary.merge(stability, on="configuration", how="left")
    summary.to_csv(root / "sensitivity_configuration_summary.csv", index=False)
    patient_table.to_csv(root / "sensitivity_patient_summary.csv", index=False)
    paired.to_csv(root / "sensitivity_paired_to_baseline.csv", index=False)
    print("\nRESUMEN CONSOLIDADO")
    print(summary.round(3).to_string(index=False))
    print("\nArchivos consolidados guardados en:", root)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    for venous_pressure, fit_fraction in itertools.product(
        VENOUS_PRESSURES, FIT_START_FRACTIONS
    ):
        run_configuration(project_root, venous_pressure, fit_fraction)
    consolidate(project_root)


if __name__ == "__main__":
    main()
