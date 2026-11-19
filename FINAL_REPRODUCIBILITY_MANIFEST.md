# Final Reproducibility Manifest

## Analysis freeze

- Final analysis commit before manifest: `a9c8e85`
- Python: `3.11.18`
- NumPy: `2.4.6`
- pandas: `3.1.1`
- SciPy: `1.17.1`
- WFDB: `4.3.1`
- Matplotlib: `3.11.2`
- Dependency validation: `python -m pip check` returned no broken requirements.
- Automated validation: 17 tests passed.

## Privacy and repository audit

- No private Windows paths, local waveform directory names, tested patient identifiers, or tested record identifiers were found in tracked files.
- `local_private/` is ignored by Git.
- The only tracked files reported under `data/` and `outputs/` were documentation placeholders:
  - `data/interim/README.md`
  - `data/processed/README.md`
  - `data/raw/README.md`
  - `data/synthetic/README.md`
  - `outputs/README.md`
- Raw data, processed data, private local files, and generated result artifacts are not tracked.

## Final artifact inventory

- Files: 12
- Total bytes: 663488
- PNG figures: 4
- PDF figures: 4
- CSV files: 4
- Empty files: 0

## SHA-256 checksums

| SHA-256 | File |
|---|---|
| `8eb820e95afcbf0eca5d9b3f20a91a7326ff964da75fe55a1816bd629d4c0c79` | `figure_1_cohort_flow.pdf` |
| `d6cc6cc09bc00c9a8d202437ad8943eaa24170c8eee9a87142f48267fc8950c1` | `figure_1_cohort_flow.png` |
| `0502eb42cac6daab223f4559c88566d2b73efb47697d113afb66b13c31c9cbe7` | `figure_1_cohort_flow_data.csv` |
| `415152ace665707e6dd387848c3409b1736a89799ed0d64d026dfa0367d25667` | `figure_2_paired_tau.pdf` |
| `6f99db03439928ab404966d45cfbe1cae8f06b783cd0dc992950c6d024a07d22` | `figure_2_paired_tau.png` |
| `40a2392803fa42c544d0fc9a5b9cc4c5eb025c4fbe3452e1cd1eebb067487646` | `figure_3_feature_forest.pdf` |
| `788099e9d93d4ad00c6a1b12cf4e40e0700ddcf14a6b977b0a96abd447261ecc` | `figure_3_feature_forest.png` |
| `905e1afdee40a4721fb694bbf93b08869dfbe28ae07ea7dc9ef6602763406f37` | `figure_4_tau_sensitivity.pdf` |
| `16df96994189916f5223512780851a8be781f47b6202f257414064722d351f93` | `figure_4_tau_sensitivity.png` |
| `df4a3cfaef045ce0353184d36f57948c7aba5691da63625c2f37f63f83675e28` | `table_1_paired_feature_results.csv` |
| `6297c1f892a8cb1b42119e21f91c702b21dc734203c8a72347b3fc7be08b7841` | `table_2_cohort_event_summary.csv` |
| `e6c308b149dd5bec2968ca0231f4bd1d0308fced86046d715bbc27a4b08dd004` | `table_3_quality_exclusions.csv` |

## Frozen primary results

- Patients reviewed in event audit: 52
- Segments reviewed: 100
- Valid hypotension events: 276
- Leakage-safe events: 243
- Patients with leakage-safe events: 26
- Non-overlapping paired patients: 23
- Pairs with effective tau available in both windows: 15
- Median valid MAP fraction: 0.988668
- Effective tau median paired difference, case minus control: -0.0262 s
- Bootstrap 95% CI: -0.1204 to 0.0198 s
- Wilcoxon p-value: 0.1876
- Benjamini-Hochberg adjusted p-value: 0.5469

## Interpretation freeze

The paired analysis did not provide statistically sufficient evidence that effective diastolic tau differed before hypotension. The estimate was slightly lower in case windows, but the bootstrap confidence interval included zero, the Wilcoxon test was not significant, and the result did not survive false-discovery-rate correction. Results are exploratory because only 15 paired patients had acceptable tau estimates in both windows.
