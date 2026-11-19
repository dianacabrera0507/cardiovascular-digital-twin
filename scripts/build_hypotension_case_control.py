"""Build one leakage-safe hypotension case-control pair per patient."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PAIR_COLUMNS = [
    "patient_code", "segment_code", "patient_segment_rank",
    "case_event_number", "case_event_start_seconds",
    "case_predictor_start_seconds", "case_predictor_end_seconds",
    "case_outcome_start_seconds", "control_anchor_seconds",
    "control_predictor_start_seconds", "control_predictor_end_seconds",
    "control_outcome_start_seconds", "control_outcome_end_seconds",
    "control_distance_from_case_seconds",
]
EXCLUSION_COLUMNS = ["patient_code", "segment_code", "patient_segment_rank", "reason"]


def intervals_overlap(a_start: float, a_end: float, b_start: float, b_end: float) -> bool:
    return a_start < b_end and b_start < a_end


def add_segment_ranks(segments: pd.DataFrame) -> pd.DataFrame:
    ranked = segments.reset_index(drop=False).rename(columns={"index": "audit_order"})
    ranked["patient_segment_rank"] = (
        ranked.groupby("patient_code")["audit_order"].rank(method="first").astype(int)
    )
    return ranked


def find_control_anchors(
    segment_duration_seconds: float,
    segment_events: pd.DataFrame,
    case_start_seconds: float,
    predictor_seconds: float,
    gap_seconds: float,
    horizon_seconds: float,
    step_seconds: float,
) -> list[float]:
    block_seconds = predictor_seconds + gap_seconds + horizon_seconds
    case_block_start = case_start_seconds - block_seconds
    case_block_end = case_start_seconds
    anchors = np.arange(
        block_seconds,
        segment_duration_seconds + np.finfo(float).eps,
        step_seconds,
        dtype=float,
    )
    valid: list[float] = []
    for anchor in anchors:
        control_start = anchor - block_seconds
        control_end = anchor
        if intervals_overlap(
            control_start, control_end, case_block_start, case_block_end
        ):
            continue
        event_overlap = any(
            intervals_overlap(
                control_start,
                control_end,
                float(row.event_start_seconds),
                float(row.event_end_seconds),
            )
            for row in segment_events.itertuples(index=False)
        )
        if not event_overlap:
            valid.append(float(anchor))
    return valid


def build_pairs(
    segments: pd.DataFrame,
    events: pd.DataFrame,
    predictor_minutes: float,
    gap_minutes: float,
    horizon_minutes: float,
    step_seconds: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ranked = add_segment_ranks(segments)
    events_ranked = events.merge(
        ranked[["patient_code", "segment_code", "patient_segment_rank", "duration_minutes"]],
        on=["patient_code", "segment_code"],
        how="left",
        validate="many_to_one",
    )
    eligible = events_ranked[
        events_ranked["leakage_safe_history"].fillna(False).astype(bool)
    ].sort_values(
        ["patient_code", "patient_segment_rank", "event_start_seconds", "event_number"]
    )
    primary = eligible.groupby("patient_code", as_index=False).first()

    predictor_seconds = predictor_minutes * 60.0
    gap_seconds = gap_minutes * 60.0
    horizon_seconds = horizon_minutes * 60.0
    pair_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []

    for case in primary.itertuples(index=False):
        segment_events = events_ranked[events_ranked["segment_code"].eq(case.segment_code)]
        anchors = find_control_anchors(
            segment_duration_seconds=float(case.duration_minutes) * 60.0,
            segment_events=segment_events,
            case_start_seconds=float(case.event_start_seconds),
            predictor_seconds=predictor_seconds,
            gap_seconds=gap_seconds,
            horizon_seconds=horizon_seconds,
            step_seconds=step_seconds,
        )
        if not anchors:
            excluded_rows.append({
                "patient_code": case.patient_code,
                "segment_code": case.segment_code,
                "patient_segment_rank": int(case.patient_segment_rank),
                "reason": "no_nonoverlapping_control_in_case_segment",
            })
            continue
        anchor = min(
            anchors,
            key=lambda value: (abs(value - float(case.event_start_seconds)), value),
        )
        control_outcome_start = anchor - horizon_seconds
        control_predictor_end = control_outcome_start - gap_seconds
        control_predictor_start = control_predictor_end - predictor_seconds
        pair_rows.append({
            "patient_code": case.patient_code,
            "segment_code": case.segment_code,
            "patient_segment_rank": int(case.patient_segment_rank),
            "case_event_number": int(case.event_number),
            "case_event_start_seconds": float(case.event_start_seconds),
            "case_predictor_start_seconds": float(case.predictor_start_seconds),
            "case_predictor_end_seconds": float(case.predictor_end_seconds),
            "case_outcome_start_seconds": float(case.outcome_start_seconds),
            "control_anchor_seconds": anchor,
            "control_predictor_start_seconds": control_predictor_start,
            "control_predictor_end_seconds": control_predictor_end,
            "control_outcome_start_seconds": control_outcome_start,
            "control_outcome_end_seconds": anchor,
            "control_distance_from_case_seconds": abs(
                anchor - float(case.event_start_seconds)
            ),
        })

    pairs = pd.DataFrame(pair_rows, columns=PAIR_COLUMNS)
    exclusions = pd.DataFrame(excluded_rows, columns=EXCLUSION_COLUMNS)
    primary_columns = [
        "patient_code", "segment_code", "patient_segment_rank", "event_number",
        "event_start_seconds", "event_end_seconds", "event_duration_seconds",
        "hypotensive_seconds", "minimum_map_mmhg", "median_map_mmhg",
        "predictor_start_seconds", "predictor_end_seconds", "outcome_start_seconds",
    ]
    return pairs, exclusions, primary[primary_columns]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build one non-overlapping case-control pair per patient."
    )
    parser.add_argument("--predictor-window-minutes", type=float, default=5.0)
    parser.add_argument("--exclusion-gap-minutes", type=float, default=5.0)
    parser.add_argument("--prediction-horizon-minutes", type=float, default=15.0)
    parser.add_argument("--control-step-seconds", type=float, default=60.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if (
        args.predictor_window_minutes <= 0
        or args.prediction_horizon_minutes <= 0
        or args.control_step_seconds <= 0
        or args.exclusion_gap_minutes < 0
    ):
        raise ValueError("invalid timing argument")
    root = Path(__file__).resolve().parents[1]
    audit_dir = root / "outputs" / "hypotension_event_audit"
    segments = pd.read_csv(audit_dir / "hypotension_segment_audit.csv")
    events = pd.read_csv(audit_dir / "hypotension_events.csv")
    pairs, exclusions, primary = build_pairs(
        segments,
        events,
        args.predictor_window_minutes,
        args.exclusion_gap_minutes,
        args.prediction_horizon_minutes,
        args.control_step_seconds,
    )
    output_dir = root / "outputs" / "hypotension_case_control"
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output_dir / "case_control_pairs.csv", index=False)
    exclusions.to_csv(output_dir / "case_control_exclusions.csv", index=False)
    primary.to_csv(output_dir / "primary_cases.csv", index=False)
    print("Construccion caso-control terminada.")
    print("Casos principales:", len(primary))
    print("Pares completos:", len(pairs))
    print("Pacientes excluidos sin control:", len(exclusions))
    print("Resultados locales:", output_dir)


if __name__ == "__main__":
    main()
