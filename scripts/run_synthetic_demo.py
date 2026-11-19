#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import yaml
from cardio_twin.simulation import make_synthetic_record
from cardio_twin.estimation import estimate_resistance_compliance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/base.yaml")
    parser.add_argument("--synthetic-config", default="configs/synthetic.yaml")
    args = parser.parse_args()

    with open(args.config, encoding="utf-8") as handle:
        base = yaml.safe_load(handle)
    with open(args.synthetic_config, encoding="utf-8") as handle:
        synthetic = yaml.safe_load(handle)["simulation"]

    record = make_synthetic_record(synthetic, seed=int(base["project"]["seed"]))
    estimate = estimate_resistance_compliance(
        record["time"], record["inflow"], record["pressure_observed"],
        venous_pressure=float(synthetic["venous_pressure"]),
    )
    output_dir = Path(base["project"]["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "truth": {"resistance": synthetic["resistance"], "compliance": synthetic["compliance"]},
        "estimate": estimate,
        "warning": "Synthetic demonstration only; not clinical evidence."
    }
    path = output_dir / "synthetic_demo_metrics.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
