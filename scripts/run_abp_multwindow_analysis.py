"""Exploratory multi-window analysis of WFDB arterial pressure signals.

The script stores anonymized summaries under outputs/abp_multwindow. It does
not estimate vascular resistance and arterial compliance separately.
"""
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
POSITIONS = (0.05, 0.25, 0.50, 0.75)

SEGMENT_COLUMNS = [
    "patient_code", "segment_code", "sampling_hz", "total_samples",
    "duration_minutes", "abp_channel_name", "windows_requested",
]
WINDOW_COLUMNS = [
    "patient_code", "segment_code", "window_number", "start_seconds", "sampling_hz",
    "finite_fraction", "median_mmhg", "p05_mmhg", "p95_mmhg",
    "robust_amplitude_mmhg", "quality_candidate", "n_peaks",
    "n_plausible_beats", "n_tau_accepted", "tau_median_s",
    "analysis_status", "error_type",
]
BEAT_COLUMNS = [
    "patient_code", "segment_code", "window_number", "beat", "heart_rate_bpm",
    "systolic_mmhg", "diastolic_mmhg", "pulse_pressure_mmhg",
    "tau_seconds", "r_squared", "rmse_mmhg", "samples_used", "success",
    "acceptable_fit", "fit_error",
]


def anonymous_code(header_path: Path, waves_root: Path) -> str:
    relative_name = header_path.relative_to(waves_root).as_posix()
    return hashlib.sha256(relative_name.encode("utf-8")).hexdigest()[:12]



def patient_code(header_path: Path, waves_root: Path) -> str:
    """Create a stable anonymous patient-level grouping code."""
    relative_parts = header_path.relative_to(waves_root).parts
    patient_parts = [
        part
        for part in relative_parts
        if part.lower().startswith("p") and part[1:].isdigit()
    ]
    patient_part = patient_parts[-1] if patient_parts else relative_parts[0]
    return hashlib.sha256(patient_part.encode("utf-8")).hexdigest()[:12]


def find_abp_index(signal_names: list[str]) -> int | None:
    for index, name in enumerate(signal_names):
        if name.strip().upper() in ABP_NAMES:
            return index
    return None


def calculate_starts(total_samples: int, window_samples: int) -> list[int]:
    if total_samples <= window_samples:
        return [0]
    maximum_start = total_samples - window_samples
    return sorted({int(round(position * maximum_start)) for position in POSITIONS})


def quality_is_acceptable(summary: dict[str, float]) -> bool:
    return bool(
        summary["finite_fraction"] >= 0.99
        and 35.0 <= summary["median_mmhg"] <= 180.0
        and 10.0 <= summary["robust_amplitude_mmhg"] <= 120.0
        and summary["p05_mmhg"] >= 15.0
        and summary["p95_mmhg"] <= 250.0
    )


def analyze_window(
    signal: np.ndarray,
    sampling_hz: float,
    patient_code_value: str,
    segment_code: str,
    window_number: int,
    start_sample: int,
    venous_pressure: float,
    fit_fraction: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    summary = summarize_abp_window(signal)
    window_row: dict[str, object] = {
        "patient_code": patient_code_value,
        "segment_code": segment_code,
        "window_number": window_number,
        "start_seconds": start_sample / sampling_hz,
        "sampling_hz": sampling_hz,
        "finite_fraction": summary["finite_fraction"],
        "median_mmhg": summary["median_mmhg"],
        "p05_mmhg": summary["p05_mmhg"],
        "p95_mmhg": summary["p95_mmhg"],
        "robust_amplitude_mmhg": summary["robust_amplitude_mmhg"],
        "quality_candidate": False,
        "n_peaks": 0,
        "n_plausible_beats": 0,
        "n_tau_accepted": 0,
        "tau_median_s": np.nan,
        "analysis_status": "quality_rejected",
        "error_type": "",
    }
    if not quality_is_acceptable(summary):
        return window_row, []
    window_row["quality_candidate"] = True
    if not np.isfinite(signal).all():
        window_row["analysis_status"] = "missing_values"
        return window_row, []
    try:
        peaks = detect_systolic_peaks(signal, sampling_hz)
    except ValueError as error:
        window_row["analysis_status"] = "peak_detection_failed"
        window_row["error_type"] = type(error).__name__
        return window_row, []
    window_row["n_peaks"] = len(peaks)
    if len(peaks) < 3:
        window_row["analysis_status"] = "insufficient_peaks"
        return window_row, []
    beats = characterize_beats(signal, peaks, sampling_hz)
    plausible_beats = [beat for beat in beats if bool(beat["plausible"])]
    window_row["n_plausible_beats"] = len(plausible_beats)
    beat_rows: list[dict[str, object]] = []
    for beat in plausible_beats:
        beat_row: dict[str, object] = {
            "patient_code": patient_code_value,
            "segment_code": segment_code,
            "window_number": window_number,
            "beat": beat["beat"],
            "heart_rate_bpm": beat["heart_rate_bpm"],
            "systolic_mmhg": beat["systolic_mmhg"],
            "diastolic_mmhg": beat["diastolic_mmhg"],
            "pulse_pressure_mmhg": beat["pulse_pressure_mmhg"],
            "tau_seconds": np.nan,
            "r_squared": np.nan,
            "rmse_mmhg": np.nan,
            "samples_used": 0,
            "success": False,
            "acceptable_fit": False,
            "fit_error": "",
        }
        try:
            fit = estimate_effective_diastolic_tau(
                values=signal,
                sampling_hz=sampling_hz,
                start_peak=int(beat["start_peak"]),
                diastolic_index=int(beat["diastolic_index"]),
                venous_pressure_mmhg=venous_pressure,
                fit_start_fraction=fit_fraction,
            )
            beat_row.update(fit)
        except ValueError as error:
            beat_row["fit_error"] = str(error)
        beat_rows.append(beat_row)
    accepted = [row for row in beat_rows if bool(row["acceptable_fit"])]
    window_row["n_tau_accepted"] = len(accepted)
    if accepted:
        tau_values = np.asarray([row["tau_seconds"] for row in accepted], dtype=float)
        window_row["tau_median_s"] = float(np.median(tau_values))
        window_row["analysis_status"] = "tau_available"
    elif plausible_beats:
        window_row["analysis_status"] = "no_acceptable_tau"
    else:
        window_row["analysis_status"] = "no_plausible_beats"
    return window_row, beat_rows


def run_analysis(
    root: Path,
    output_dir: Path,
    maximum_segments: int,
    maximum_segments_per_patient: int,
    window_seconds: float,
    venous_pressure: float,
    fit_fraction: float,
) -> None:
    waves_root = root / "waves"
    if not waves_root.exists():
        raise FileNotFoundError("No se encontro la carpeta waves.")
    output_dir.mkdir(parents=True, exist_ok=True)
    segments: list[dict[str, object]] = []
    windows: list[dict[str, object]] = []
    beats: list[dict[str, object]] = []
    reviewed = loaded = skipped_headers = errors = 0
    patient_segment_counts: dict[str, int] = {}
    for header_path in sorted(waves_root.rglob("*.hea")):
        if reviewed >= maximum_segments:
            break
        try:
            record_base = header_path.with_suffix("")
            header = wfdb.rdheader(str(record_base))
            if not header.sig_name:
                skipped_headers += 1
                continue
            signal_names = [str(name).strip() for name in header.sig_name]
            abp_index = find_abp_index(signal_names)
            if abp_index is None:
                continue
            if header.sig_len is None or header.fs is None:
                skipped_headers += 1
                continue
            total_samples = int(header.sig_len)
            sampling_hz = float(header.fs)
            if total_samples <= 0 or sampling_hz <= 0:
                skipped_headers += 1
                continue
            requested_samples = int(round(window_seconds * sampling_hz))
            window_samples = min(total_samples, requested_samples)
            starts = calculate_starts(total_samples, window_samples)
            segment_code = anonymous_code(header_path, waves_root)
            patient_code_value = patient_code(header_path, waves_root)
            current_patient_count = patient_segment_counts.get(patient_code_value, 0)
            if current_patient_count >= maximum_segments_per_patient:
                continue
            patient_segment_counts[patient_code_value] = current_patient_count + 1
            reviewed += 1
            segments.append({
                "patient_code": patient_code_value,
                "segment_code": segment_code,
                "sampling_hz": sampling_hz,
                "total_samples": total_samples,
                "duration_minutes": total_samples / sampling_hz / 60.0,
                "abp_channel_name": signal_names[abp_index],
                "windows_requested": len(starts),
            })
            for number, start in enumerate(starts, start=1):
                stop = min(start + window_samples, total_samples)
                try:
                    record = wfdb.rdrecord(
                        str(record_base), sampfrom=start, sampto=stop,
                        channels=[abp_index], physical=True,
                    )
                    signal = np.asarray(record.p_signal[:, 0], dtype=float)
                    loaded += 1
                    window_row, beat_rows = analyze_window(
                        signal=signal,
                        sampling_hz=sampling_hz,
                        patient_code_value=patient_code_value,
                        segment_code=segment_code,
                        window_number=number,
                        start_sample=start,
                        venous_pressure=venous_pressure,
                        fit_fraction=fit_fraction,
                    )
                    windows.append(window_row)
                    beats.extend(beat_rows)
                except Exception as error:
                    errors += 1
                    windows.append({
                        "patient_code": patient_code_value,
                        "segment_code": segment_code,
                        "window_number": number,
                        "start_seconds": start / sampling_hz,
                        "sampling_hz": sampling_hz,
                        "analysis_status": "read_error",
                        "error_type": type(error).__name__,
                    })
        except Exception:
            errors += 1
    segment_table = pd.DataFrame(segments, columns=SEGMENT_COLUMNS)
    window_table = pd.DataFrame(windows).reindex(columns=WINDOW_COLUMNS)
    beat_table = pd.DataFrame(beats).reindex(columns=BEAT_COLUMNS)
    segment_table.to_csv(output_dir / "abp_segments_summary.csv", index=False)
    window_table.to_csv(output_dir / "abp_windows_summary.csv", index=False)
    beat_table.to_csv(output_dir / "abp_beats_tau.csv", index=False)
    quality_count = int(window_table["quality_candidate"].fillna(False).astype(bool).sum())
    tau_window_count = int((window_table["analysis_status"] == "tau_available").sum())
    accepted = beat_table[beat_table["acceptable_fit"].fillna(False).astype(bool)]
    print("Analisis terminado.")
    print("Pacientes anonimos revisados:", segment_table["patient_code"].nunique())
    print("Segmentos revisados:", reviewed)
    print("Maximo de segmentos por paciente:", maximum_segments_per_patient)
    print("Ventanas cargadas:", loaded)
    print("Encabezados omitidos:", skipped_headers)
    print("Errores registrados:", errors)
    print("Ventanas candidatas por calidad:", quality_count)
    print("Ventanas con tau disponible:", tau_window_count)
    print("Ajustes de tau aceptados:", len(accepted))
    if not accepted.empty:
        tau = accepted["tau_seconds"].astype(float)
        print("Tau mediana global:", round(float(tau.median()), 3), "s")
        print(
            "Rango intercuartil global:",
            round(float(tau.quantile(0.25)), 3), "-",
            round(float(tau.quantile(0.75)), 3), "s",
        )
    print("Resultados:", output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analisis multiventana ABP.")
    parser.add_argument("--maximum-segments", type=int, default=25)
    parser.add_argument("--maximum-segments-per-patient", type=int, default=2)
    parser.add_argument("--window-seconds", type=float, default=30.0)
    parser.add_argument("--venous-pressure-mmhg", type=float, default=5.0)
    parser.add_argument("--fit-start-fraction", type=float, default=0.35)
    parser.add_argument(
        "--output-label",
        type=str,
        default="baseline",
        help="Etiqueta segura para separar resultados locales.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.maximum_segments <= 0:
        raise ValueError("--maximum-segments debe ser positivo.")
    if args.maximum_segments_per_patient <= 0:
        raise ValueError("--maximum-segments-per-patient debe ser positivo.")
    if args.window_seconds <= 0:
        raise ValueError("--window-seconds debe ser positivo.")
    if not 0.0 <= args.fit_start_fraction < 1.0:
        raise ValueError("--fit-start-fraction debe estar entre 0 y 1.")
    if not args.output_label.replace("_", "").replace("-", "").isalnum():
        raise ValueError("--output-label solo admite letras, numeros, guion y guion bajo.")
    root_value = os.environ.get("MIMIC_WAVEFORM_ROOT")
    if not root_value:
        raise EnvironmentError("Debe configurar MIMIC_WAVEFORM_ROOT.")
    root = Path(root_value).expanduser()
    if not root.exists():
        raise FileNotFoundError("MIMIC_WAVEFORM_ROOT no existe.")
    project_root = Path(__file__).resolve().parents[1]
    output_dir = project_root / "outputs" / "abp_sensitivity" / args.output_label
    run_analysis(
        root=root,
        output_dir=output_dir,
        maximum_segments=args.maximum_segments,
        maximum_segments_per_patient=args.maximum_segments_per_patient,
        window_seconds=args.window_seconds,
        venous_pressure=args.venous_pressure_mmhg,
        fit_fraction=args.fit_start_fraction,
    )


if __name__ == "__main__":
    main()
