"""Audit hypotension events in long invasive-ABP WFDB segments.

The audit uses anonymized patient and segment codes, balances selection to a
maximum number of segments per patient, derives one-second mean arterial
pressure from the invasive ABP waveform, and applies the prespecified event
and leakage-safe timing rules. Generated CSV files remain under outputs/.
"""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
import wfdb

from cardio_twin.events import prediction_intervals, threshold_episodes

ABP_NAMES = {"ABP", "ART", "ARTERIAL", "AO"}

SEGMENT_COLUMNS = [
    "patient_code", "segment_code", "duration_minutes", "sampling_hz",
    "one_second_points", "finite_map_fraction", "events_total",
    "events_with_leakage_safe_history",
]
EVENT_COLUMNS = [
    "patient_code", "segment_code", "event_number",
    "event_start_seconds", "event_end_seconds", "event_duration_seconds",
    "hypotensive_seconds", "minimum_map_mmhg", "median_map_mmhg",
    "leakage_safe_history", "predictor_start_seconds",
    "predictor_end_seconds", "outcome_start_seconds",
]
ERROR_COLUMNS = [
    "patient_code", "segment_code", "error_type", "error_message",
]


def anonymous_code(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]


def patient_identifier(header_path: Path, waves_root: Path) -> str:
    parts = header_path.relative_to(waves_root).parts
    candidates = [
        part
        for part in parts
        if part.lower().startswith("p") and part[1:].isdigit()
    ]
    return candidates[-1] if candidates else parts[0]


def find_abp_index(names: list[str]) -> int | None:
    for index, name in enumerate(names):
        if name.strip().upper() in ABP_NAMES:
            return index
    return None


def one_second_map(
    signal: np.ndarray,
    sampling_hz: float,
    minimum_finite_fraction: float = 0.80,
    minimum_map_mmhg: float = 20.0,
    maximum_map_mmhg: float = 250.0,
    minimum_pulse_amplitude_mmhg: float = 10.0,
    maximum_pulse_amplitude_mmhg: float = 120.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Aggregate ABP into quality-controlled one-second MAP estimates.

    A second is invalidated when too many waveform samples are missing, its
    mean pressure is outside the prespecified physiological range, or its
    robust within-second pulse amplitude is implausibly small or large.
    """
    samples_per_second = int(round(sampling_hz))
    if samples_per_second <= 0:
        raise ValueError("invalid sampling frequency")
    complete_seconds = signal.size // samples_per_second
    if complete_seconds <= 0:
        return np.asarray([], dtype=float), np.asarray([], dtype=float)
    matrix = signal[: complete_seconds * samples_per_second].reshape(
        complete_seconds, samples_per_second
    )
    finite = np.isfinite(matrix)
    finite_count = finite.sum(axis=1)
    finite_fraction = finite_count / samples_per_second
    safe_values = np.where(finite, matrix, 0.0)
    map_values = np.divide(
        safe_values.sum(axis=1),
        finite_count,
        out=np.full(complete_seconds, np.nan, dtype=float),
        where=finite_count > 0,
    )
    p05 = np.full(complete_seconds, np.nan, dtype=float)
    p95 = np.full(complete_seconds, np.nan, dtype=float)
    rows_with_data = finite_count > 0
    if rows_with_data.any():
        p05[rows_with_data] = np.nanpercentile(
            matrix[rows_with_data], 5.0, axis=1
        )
        p95[rows_with_data] = np.nanpercentile(
            matrix[rows_with_data], 95.0, axis=1
        )
    pulse_amplitude = p95 - p05
    valid = (
        (finite_fraction >= minimum_finite_fraction)
        & (map_values >= minimum_map_mmhg)
        & (map_values <= maximum_map_mmhg)
        & (pulse_amplitude >= minimum_pulse_amplitude_mmhg)
        & (pulse_amplitude <= maximum_pulse_amplitude_mmhg)
    )
    map_values[~valid] = np.nan
    time_seconds = np.arange(complete_seconds, dtype=float)
    return time_seconds, map_values


def select_segments(
    waves_root: Path,
    maximum_segments: int,
    maximum_segments_per_patient: int,
    minimum_duration_minutes: float,
) -> tuple[list[dict[str, object]], int]:
    selected: list[dict[str, object]] = []
    patient_counts: dict[str, int] = {}
    skipped_headers = 0
    for header_path in sorted(waves_root.rglob("*.hea")):
        if len(selected) >= maximum_segments:
            break
        try:
            record_base = header_path.with_suffix("")
            header = wfdb.rdheader(str(record_base))
            if not header.sig_name or header.sig_len is None or header.fs is None:
                skipped_headers += 1
                continue
            names = [str(name).strip() for name in header.sig_name]
            abp_index = find_abp_index(names)
            if abp_index is None:
                continue
            total_samples = int(header.sig_len)
            sampling_hz = float(header.fs)
            if total_samples <= 0 or sampling_hz <= 0:
                skipped_headers += 1
                continue
            duration_minutes = total_samples / sampling_hz / 60.0
            if duration_minutes < minimum_duration_minutes:
                continue
            raw_patient = patient_identifier(header_path, waves_root)
            patient_code = anonymous_code(raw_patient)
            current_count = patient_counts.get(patient_code, 0)
            if current_count >= maximum_segments_per_patient:
                continue
            patient_counts[patient_code] = current_count + 1
            relative_record = header_path.relative_to(waves_root).with_suffix("").as_posix()
            selected.append(
                {
                    "record_base": record_base,
                    "patient_code": patient_code,
                    "segment_code": anonymous_code(relative_record),
                    "abp_index": abp_index,
                    "sampling_hz": sampling_hz,
                    "total_samples": total_samples,
                    "duration_minutes": duration_minutes,
                }
            )
        except Exception:
            skipped_headers += 1
    return selected, skipped_headers


def run_audit(args: argparse.Namespace) -> None:
    root_value = os.environ.get("MIMIC_WAVEFORM_ROOT")
    if not root_value:
        raise EnvironmentError("Debe configurar MIMIC_WAVEFORM_ROOT.")
    project_root = Path(__file__).resolve().parents[1]
    waves_root = Path(root_value).expanduser() / "waves"
    if not waves_root.exists():
        raise FileNotFoundError("No se encontro la carpeta waves.")

    selected, skipped_headers = select_segments(
        waves_root=waves_root,
        maximum_segments=args.maximum_segments,
        maximum_segments_per_patient=args.maximum_segments_per_patient,
        minimum_duration_minutes=args.minimum_duration_minutes,
    )
    segment_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    error_rows: list[dict[str, object]] = []
    read_errors = 0
    required_history_seconds = (
        args.predictor_window_minutes
        + args.exclusion_gap_minutes
        + args.prediction_horizon_minutes
    ) * 60.0

    for item in selected:
        try:
            record = wfdb.rdrecord(
                str(item["record_base"]),
                channels=[int(item["abp_index"])],
                physical=True,
            )
            signal = np.asarray(record.p_signal[:, 0], dtype=float)
            time_seconds, map_values = one_second_map(signal, float(item["sampling_hz"]))
            episodes = threshold_episodes(
                time_seconds,
                map_values,
                threshold=args.threshold_mmhg,
                minimum_duration=args.minimum_duration_seconds,
                merge_gap=args.merge_gap_seconds,
            )
            eligible_events = 0
            for event_number, (start, end) in enumerate(episodes, start=1):
                event_start = float(time_seconds[start])
                event_end = float(time_seconds[end - 1] + 1.0)
                leakage_safe = event_start >= required_history_seconds
                intervals: dict[str, tuple[float, float]] | None = None
                if leakage_safe:
                    intervals = prediction_intervals(
                        event_start_seconds=event_start,
                        predictor_window_seconds=args.predictor_window_minutes * 60.0,
                        exclusion_gap_seconds=args.exclusion_gap_minutes * 60.0,
                        prediction_horizon_seconds=args.prediction_horizon_minutes * 60.0,
                    )
                    eligible_events += 1
                event_map = map_values[start:end]
                hypotensive_map = event_map[
                    np.isfinite(event_map) & (event_map < args.threshold_mmhg)
                ]
                event_rows.append(
                    {
                        "patient_code": item["patient_code"],
                        "segment_code": item["segment_code"],
                        "event_number": event_number,
                        "event_start_seconds": event_start,
                        "event_end_seconds": event_end,
                        "event_duration_seconds": event_end - event_start,
                        "hypotensive_seconds": int(hypotensive_map.size),
                        "minimum_map_mmhg": float(np.min(hypotensive_map)),
                        "median_map_mmhg": float(np.median(hypotensive_map)),
                        "leakage_safe_history": leakage_safe,
                        "predictor_start_seconds": intervals["predictor"][0] if intervals else np.nan,
                        "predictor_end_seconds": intervals["predictor"][1] if intervals else np.nan,
                        "outcome_start_seconds": intervals["outcome_horizon"][0] if intervals else np.nan,
                    }
                )
            segment_rows.append(
                {
                    "patient_code": item["patient_code"],
                    "segment_code": item["segment_code"],
                    "duration_minutes": item["duration_minutes"],
                    "sampling_hz": item["sampling_hz"],
                    "one_second_points": len(map_values),
                    "finite_map_fraction": float(np.isfinite(map_values).mean()) if len(map_values) else 0.0,
                    "events_total": len(episodes),
                    "events_with_leakage_safe_history": eligible_events,
                }
            )
        except Exception as error:
            read_errors += 1
            error_rows.append(
                {
                    "patient_code": item["patient_code"],
                    "segment_code": item["segment_code"],
                    "error_type": type(error).__name__,
                    "error_message": str(error)[:500],
                }
            )
            segment_rows.append(
                {
                    "patient_code": item["patient_code"],
                    "segment_code": item["segment_code"],
                    "duration_minutes": item["duration_minutes"],
                    "sampling_hz": item["sampling_hz"],
                    "one_second_points": 0,
                    "finite_map_fraction": np.nan,
                    "events_total": 0,
                    "events_with_leakage_safe_history": 0,
                }
            )

    segments = pd.DataFrame(segment_rows, columns=SEGMENT_COLUMNS)
    events = pd.DataFrame(event_rows, columns=EVENT_COLUMNS)
    errors = pd.DataFrame(error_rows, columns=ERROR_COLUMNS)
    output_dir = project_root / "outputs" / "hypotension_event_audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    segments.to_csv(output_dir / "hypotension_segment_audit.csv", index=False)
    events.to_csv(output_dir / "hypotension_events.csv", index=False)
    errors.to_csv(output_dir / "hypotension_read_errors.csv", index=False)

    patients_reviewed = int(segments["patient_code"].nunique()) if not segments.empty else 0
    patients_with_events = int(events["patient_code"].nunique()) if not events.empty else 0
    eligible = events[events["leakage_safe_history"].fillna(False).astype(bool)] if not events.empty else events
    patients_with_eligible = int(eligible["patient_code"].nunique()) if not eligible.empty else 0
    print("Auditoria de eventos de hipotension terminada.")
    print("Pacientes revisados:", patients_reviewed)
    print("Segmentos revisados:", len(segments))
    print("Encabezados omitidos:", skipped_headers)
    print("Errores de lectura:", read_errors)
    print("Eventos totales:", len(events))
    print("Pacientes con eventos:", patients_with_events)
    print("Eventos con historia sin fuga:", len(eligible))
    print("Pacientes con eventos elegibles:", patients_with_eligible)
    print("Resultados locales:", output_dir)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit prespecified hypotension events.")
    parser.add_argument("--maximum-segments", type=int, default=100)
    parser.add_argument("--maximum-segments-per-patient", type=int, default=2)
    parser.add_argument("--minimum-duration-minutes", type=float, default=25.0)
    parser.add_argument("--threshold-mmhg", type=float, default=65.0)
    parser.add_argument("--minimum-duration-seconds", type=float, default=60.0)
    parser.add_argument("--merge-gap-seconds", type=float, default=60.0)
    parser.add_argument("--predictor-window-minutes", type=float, default=5.0)
    parser.add_argument("--exclusion-gap-minutes", type=float, default=5.0)
    parser.add_argument("--prediction-horizon-minutes", type=float, default=15.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    positive = (
        args.maximum_segments,
        args.maximum_segments_per_patient,
        args.minimum_duration_minutes,
        args.minimum_duration_seconds,
        args.predictor_window_minutes,
        args.prediction_horizon_minutes,
    )
    if any(value <= 0 for value in positive):
        raise ValueError("positive arguments must be greater than zero")
    if args.merge_gap_seconds < 0 or args.exclusion_gap_minutes < 0:
        raise ValueError("gap arguments must be non-negative")
    run_audit(args)


if __name__ == "__main__":
    main()
