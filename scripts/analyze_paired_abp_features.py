"""Paired statistical analysis of case-control ABP predictor features."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon

FEATURES = [
    "median_mmhg",
    "minimum_mmhg",
    "maximum_mmhg",
    "iqr_mmhg",
    "standard_deviation_mmhg",
    "robust_amplitude_mmhg",
    "map_slope_mmhg_per_minute",
    "heart_rate_median_bpm",
    "pulse_pressure_median_mmhg",
    "tau_acceptance_fraction",
    "tau_median_s",
]


def bootstrap_median_difference(
    differences: np.ndarray,
    repetitions: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    sampled = rng.choice(
        differences,
        size=(repetitions, differences.size),
        replace=True,
    )
    estimates = np.median(sampled, axis=1)
    lower, upper = np.quantile(estimates, [0.025, 0.975])
    return float(lower), float(upper), float(np.std(estimates, ddof=1))


def rank_biserial(differences: np.ndarray) -> float:
    nonzero = differences[differences != 0]
    if nonzero.size == 0:
        return 0.0
    ranks = rankdata(np.abs(nonzero), method="average")
    positive = float(ranks[nonzero > 0].sum())
    negative = float(ranks[nonzero < 0].sum())
    denominator = positive + negative
    return (positive - negative) / denominator if denominator else 0.0


def benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    adjusted = pd.Series(np.nan, index=p_values.index, dtype=float)
    finite = p_values.dropna().sort_values()
    count = len(finite)
    if count == 0:
        return adjusted
    raw = finite.to_numpy(dtype=float)
    ranks = np.arange(1, count + 1, dtype=float)
    corrected = raw * count / ranks
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    adjusted.loc[finite.index] = np.minimum(corrected, 1.0)
    return adjusted


def analyze_feature(
    features: pd.DataFrame,
    feature: str,
    repetitions: int,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    paired = features.pivot(
        index="patient_code",
        columns="window_type",
        values=feature,
    ).dropna(subset=["case", "control"])
    paired["difference_case_minus_control"] = paired["case"] - paired["control"]
    differences = paired["difference_case_minus_control"].to_numpy(dtype=float)
    row: dict[str, object] = {
        "feature": feature,
        "n_pairs": len(paired),
        "case_median": float(paired["case"].median()) if len(paired) else np.nan,
        "case_q25": float(paired["case"].quantile(0.25)) if len(paired) else np.nan,
        "case_q75": float(paired["case"].quantile(0.75)) if len(paired) else np.nan,
        "control_median": float(paired["control"].median()) if len(paired) else np.nan,
        "control_q25": float(paired["control"].quantile(0.25)) if len(paired) else np.nan,
        "control_q75": float(paired["control"].quantile(0.75)) if len(paired) else np.nan,
        "median_difference": float(np.median(differences)) if len(paired) else np.nan,
        "bootstrap_95_low": np.nan,
        "bootstrap_95_high": np.nan,
        "bootstrap_standard_error": np.nan,
        "wilcoxon_statistic": np.nan,
        "wilcoxon_p_value": np.nan,
        "rank_biserial_correlation": np.nan,
    }
    if len(paired):
        rng = np.random.default_rng(seed)
        low, high, standard_error = bootstrap_median_difference(
            differences,
            repetitions,
            rng,
        )
        row["bootstrap_95_low"] = low
        row["bootstrap_95_high"] = high
        row["bootstrap_standard_error"] = standard_error
        row["rank_biserial_correlation"] = rank_biserial(differences)
        if np.any(differences != 0):
            result = wilcoxon(
                paired["case"],
                paired["control"],
                zero_method="wilcox",
                correction=False,
                alternative="two-sided",
                method="auto",
            )
            row["wilcoxon_statistic"] = float(result.statistic)
            row["wilcoxon_p_value"] = float(result.pvalue)
    detail = paired.reset_index()
    detail.insert(1, "feature", feature)
    return row, detail


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze paired case-control ABP predictor features."
    )
    parser.add_argument("--bootstrap-repetitions", type=int, default=10000)
    parser.add_argument("--random-seed", type=int, default=2026)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.bootstrap_repetitions <= 0:
        raise ValueError("--bootstrap-repetitions must be positive")
    project_root = Path(__file__).resolve().parents[1]
    input_file = (
        project_root
        / "outputs"
        / "hypotension_case_control_features"
        / "paired_abp_features.csv"
    )
    if not input_file.exists():
        raise FileNotFoundError("paired ABP feature table is missing")
    features = pd.read_csv(input_file)
    missing = [feature for feature in FEATURES if feature not in features.columns]
    if missing:
        raise ValueError(f"missing feature columns: {missing}")

    rows: list[dict[str, object]] = []
    details: list[pd.DataFrame] = []
    for offset, feature in enumerate(FEATURES):
        row, detail = analyze_feature(
            features,
            feature,
            args.bootstrap_repetitions,
            args.random_seed + offset,
        )
        rows.append(row)
        details.append(detail)
    summary = pd.DataFrame(rows)
    summary["benjamini_hochberg_p_value"] = benjamini_hochberg(
        summary["wilcoxon_p_value"]
    )
    summary["significant_fdr_0_05"] = (
        summary["benjamini_hochberg_p_value"] < 0.05
    )
    paired_details = pd.concat(details, ignore_index=True)

    output_dir = project_root / "outputs" / "paired_statistical_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "paired_feature_statistics.csv", index=False)
    paired_details.to_csv(output_dir / "paired_feature_differences.csv", index=False)
    tau_summary = summary[summary["feature"].eq("tau_median_s")]
    tau_summary.to_csv(output_dir / "tau_primary_statistics.csv", index=False)

    print("Analisis estadistico pareado terminado.")
    print("Caracteristicas analizadas:", len(summary))
    print("Pares para tau:", int(tau_summary["n_pairs"].iloc[0]))
    print(
        "Diferencia mediana tau:",
        round(float(tau_summary["median_difference"].iloc[0]), 3),
        "s",
    )
    print(
        "IC 95% bootstrap tau:",
        round(float(tau_summary["bootstrap_95_low"].iloc[0]), 3),
        "-",
        round(float(tau_summary["bootstrap_95_high"].iloc[0]), 3),
        "s",
    )
    print(
        "p Wilcoxon tau:",
        round(float(tau_summary["wilcoxon_p_value"].iloc[0]), 4),
    )
    print(
        "p FDR tau:",
        round(float(tau_summary["benjamini_hochberg_p_value"].iloc[0]), 4),
    )
    print("Resultados locales:", output_dir)


if __name__ == "__main__":
    main()
