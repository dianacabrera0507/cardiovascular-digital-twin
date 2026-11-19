"""Extract ABP and effective-tau features for paired predictor windows."""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
import wfdb

from cardio_twin.abp import (
    characterize_beats,
    detect_systolic_peaks,
    estimate_effective_diastolic_tau,
    summarize_abp_window,
)

ABP_NAMES = {"ABP", "ART", "ARTERIAL", "AO"}
FEATURE_COLUMNS = [
    "patient_code", "segment_code", "window_type", "start_seconds", "end_seconds",
    "sampling_hz", "finite_fraction", "median_mmhg", "minimum_mmhg",
    "maximum_mmhg", "p05_mmhg", "p25_mmhg", "p75_mmhg", "p95_mmhg",
    "iqr_mmhg", "standard_deviation_mmhg", "robust_amplitude_mmhg",
    "map_slope_mmhg_per_minute", "n_peaks", "n_plausible_beats",
    "heart_rate_median_bpm", "pulse_pressure_median_mmhg", "n_tau_accepted",
    "tau_acceptance_fraction", "tau_median_s", "analysis_status", "error_type",
]


def code(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]


def patient_code(header: Path, root: Path) -> str:
    parts = header.relative_to(root).parts
    candidates = [p for p in parts if p.lower().startswith("p") and p[1:].isdigit()]
    return code(candidates[-1] if candidates else parts[0])


def segment_code(header: Path, root: Path) -> str:
    return code(header.relative_to(root).with_suffix("").as_posix())


def find_abp_index(names: list[str]) -> int | None:
    for index, name in enumerate(names):
        if name.strip().upper() in ABP_NAMES:
            return index
    return None


def locate_records(root: Path, wanted: set[str]) -> dict[str, tuple[Path, int, float]]:
    found: dict[str, tuple[Path, int, float]] = {}
    for header_path in sorted(root.rglob("*.hea")):
        current_code = segment_code(header_path, root)
        if current_code not in wanted:
            continue
        header = wfdb.rdheader(str(header_path.with_suffix("")))
        names = [str(name).strip() for name in header.sig_name]
        index = find_abp_index(names)
        if index is None:
            continue
        found[current_code] = (header_path.with_suffix(""), index, float(header.fs))
        if len(found) == len(wanted):
            break
    return found


def quality_ok(summary: dict[str, float]) -> bool:
    return bool(
        summary["finite_fraction"] >= 0.99
        and 35.0 <= summary["median_mmhg"] <= 180.0
        and 10.0 <= summary["robust_amplitude_mmhg"] <= 120.0
        and summary["p05_mmhg"] >= 15.0
        and summary["p95_mmhg"] <= 250.0
    )


def slope_per_minute(values: np.ndarray, sampling_hz: float) -> float:
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return float("nan")
    time_minutes = np.arange(values.size, dtype=float) / sampling_hz / 60.0
    return float(np.polyfit(time_minutes[finite], values[finite], 1)[0])


def analyze(
    signal: np.ndarray,
    sampling_hz: float,
    patient: str,
    segment: str,
    window_type: str,
    start_seconds: float,
    end_seconds: float,
    venous_pressure: float,
    fit_fraction: float,
) -> dict[str, object]:
    summary = summarize_abp_window(signal)
    finite_values = signal[np.isfinite(signal)]
    row: dict[str, object] = {
        "patient_code": patient,
        "segment_code": segment,
        "window_type": window_type,
        "start_seconds": start_seconds,
        "end_seconds": end_seconds,
        "sampling_hz": sampling_hz,
        "finite_fraction": summary["finite_fraction"],
        "median_mmhg": summary["median_mmhg"],
        "minimum_mmhg": float(np.min(finite_values)) if finite_values.size else np.nan,
        "maximum_mmhg": float(np.max(finite_values)) if finite_values.size else np.nan,
        "p05_mmhg": summary["p05_mmhg"],
        "p25_mmhg": float(np.quantile(finite_values, 0.25)) if finite_values.size else np.nan,
        "p75_mmhg": float(np.quantile(finite_values, 0.75)) if finite_values.size else np.nan,
        "p95_mmhg": summary["p95_mmhg"],
        "iqr_mmhg": float(np.quantile(finite_values, 0.75) - np.quantile(finite_values, 0.25)) if finite_values.size else np.nan,
        "standard_deviation_mmhg": float(np.std(finite_values, ddof=1)) if finite_values.size > 1 else np.nan,
        "robust_amplitude_mmhg": summary["robust_amplitude_mmhg"],
        "map_slope_mmhg_per_minute": slope_per_minute(signal, sampling_hz),
        "n_peaks": 0,
        "n_plausible_beats": 0,
        "heart_rate_median_bpm": np.nan,
        "pulse_pressure_median_mmhg": np.nan,
        "n_tau_accepted": 0,
        "tau_acceptance_fraction": np.nan,
        "tau_median_s": np.nan,
        "analysis_status": "quality_rejected",
        "error_type": "",
    }
    if not quality_ok(summary):
        return row
    if not np.isfinite(signal).all():
        row["analysis_status"] = "missing_values"
        return row
    try:
        peaks = detect_systolic_peaks(signal, sampling_hz)
        row["n_peaks"] = len(peaks)
        if len(peaks) < 3:
            row["analysis_status"] = "insufficient_peaks"
            return row
        beats = characterize_beats(signal, peaks, sampling_hz)
        plausible = [beat for beat in beats if bool(beat["plausible"])]
        row["n_plausible_beats"] = len(plausible)
        if not plausible:
            row["analysis_status"] = "no_plausible_beats"
            return row
        row["heart_rate_median_bpm"] = float(np.median([b["heart_rate_bpm"] for b in plausible]))
        row["pulse_pressure_median_mmhg"] = float(np.median([b["pulse_pressure_mmhg"] for b in plausible]))
        accepted: list[float] = []
        for beat in plausible:
            try:
                fit = estimate_effective_diastolic_tau(
                    values=signal,
                    sampling_hz=sampling_hz,
                    start_peak=int(beat["start_peak"]),
                    diastolic_index=int(beat["diastolic_index"]),
                    venous_pressure_mmhg=venous_pressure,
                    fit_start_fraction=fit_fraction,
                )
                if bool(fit["acceptable_fit"]):
                    accepted.append(float(fit["tau_seconds"]))
            except ValueError:
                continue
        row["n_tau_accepted"] = len(accepted)
        row["tau_acceptance_fraction"] = len(accepted) / len(plausible)
        if accepted:
            row["tau_median_s"] = float(np.median(accepted))
            row["analysis_status"] = "tau_available"
        else:
            row["analysis_status"] = "no_acceptable_tau"
    except (ValueError, RuntimeError) as error:
        row["analysis_status"] = "analysis_failed"
        row["error_type"] = type(error).__name__
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract paired ABP predictor features.")
    parser.add_argument("--venous-pressure-mmhg", type=float, default=5.0)
    parser.add_argument("--fit-start-fraction", type=float, default=0.35)
    args = parser.parse_args()
    root_value = os.environ.get("MIMIC_WAVEFORM_ROOT")
    if not root_value:
        raise EnvironmentError("Debe configurar MIMIC_WAVEFORM_ROOT.")
    project_root = Path(__file__).resolve().parents[1]
    waves_root = Path(root_value).expanduser() / "waves"
    pairs = pd.read_csv(
        project_root / "outputs" / "hypotension_case_control" / "case_control_pairs.csv"
    )
    records = locate_records(waves_root, set(pairs["segment_code"]))
    rows: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    for pair in pairs.itertuples(index=False):
        if pair.segment_code not in records:
            errors.append({"patient_code": pair.patient_code, "segment_code": pair.segment_code, "error_type": "record_not_found"})
            continue
        record_base, channel, sampling_hz = records[pair.segment_code]
        for window_type, start_seconds, end_seconds in (
            ("case", pair.case_predictor_start_seconds, pair.case_predictor_end_seconds),
            ("control", pair.control_predictor_start_seconds, pair.control_predictor_end_seconds),
        ):
            start = int(round(float(start_seconds) * sampling_hz))
            stop = int(round(float(end_seconds) * sampling_hz))
            try:
                record = wfdb.rdrecord(str(record_base), sampfrom=start, sampto=stop, channels=[channel], physical=True)
                signal = np.asarray(record.p_signal[:, 0], dtype=float)
                rows.append(analyze(signal, sampling_hz, pair.patient_code, pair.segment_code, window_type, float(start_seconds), float(end_seconds), args.venous_pressure_mmhg, args.fit_start_fraction))
            except Exception as error:
                errors.append({"patient_code": pair.patient_code, "segment_code": pair.segment_code, "error_type": type(error).__name__})
    features = pd.DataFrame(rows).reindex(columns=FEATURE_COLUMNS)
    error_table = pd.DataFrame(errors, columns=["patient_code", "segment_code", "error_type"])
    output_dir = project_root / "outputs" / "hypotension_case_control_features"
    output_dir.mkdir(parents=True, exist_ok=True)
    features.to_csv(output_dir / "paired_abp_features.csv", index=False)
    error_table.to_csv(output_dir / "feature_extraction_errors.csv", index=False)
    complete = features.pivot(index="patient_code", columns="window_type", values="analysis_status") if not features.empty else pd.DataFrame()
    complete_pairs = int(complete.dropna().shape[0]) if not complete.empty else 0
    tau_pairs = int(((complete.get("case") == "tau_available") & (complete.get("control") == "tau_available")).sum()) if {"case", "control"}.issubset(complete.columns) else 0
    print("Extraccion de caracteristicas terminada.")
    print("Ventanas extraidas:", len(features))
    print("Pares completos extraidos:", complete_pairs)
    print("Pares con tau en caso y control:", tau_pairs)
    print("Errores registrados:", len(error_table))
    print("Resultados locales:", output_dir)


if __name__ == "__main__":
    main()
