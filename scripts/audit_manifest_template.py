#!/usr/bin/env python3
"""Create a metadata manifest without reading or exporting signal content."""
from pathlib import Path
import csv, argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, default=Path("outputs/file_manifest.csv"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in args.root.rglob("*"):
        if path.is_file():
            rows.append({"relative_path": str(path.relative_to(args.root)), "bytes": path.stat().st_size, "suffix": path.suffix})
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "suffix"])
        writer.writeheader(); writer.writerows(rows)
    print(f"Wrote {len(rows)} entries to {args.output}")


if __name__ == "__main__":
    main()
