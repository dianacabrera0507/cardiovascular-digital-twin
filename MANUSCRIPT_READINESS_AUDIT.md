# Manuscript-Readiness Audit

## Scope freeze

This work is an exploratory, retrospective, within-patient matched case-control analysis of arterial blood pressure waveforms. It is not a validated clinical prediction model and does not constitute a complete cardiovascular digital twin. The modeled quantity is an effective diastolic time constant derived under Windkessel assumptions.

## Data source

- Dataset: MIMIC-IV Waveform Database
- Version: 0.1.0
- Provider: PhysioNet
- DOI: 10.13026/a2mw-f949
- License: Open Data Commons Open Database License (ODbL)
- Data used: waveform data only
- Clinical MIMIC-IV tables: not used
- Original waveform redistribution: not permitted through this repository
- Patient-level derived outputs: not published through this repository

## Ethics status

The source database is deidentified and was created under the governance described by the MIMIC investigators. The present analysis involved no participant contact, intervention, or reidentification attempt. A formal Universidad Técnica del Norte determination on exemption or non-human-subject secondary analysis remains pending and must be obtained before manuscript submission.

## STROBE case-control readiness matrix

| Item | Requirement | Current status | Required manuscript action |
|---|---|---|---|
| 1 | Design in title/abstract; balanced abstract | Partial | Use “retrospective matched case-control waveform study”; write structured abstract after manuscript completion. |
| 2 | Scientific background and rationale | Pending | Explain hypotension relevance, Windkessel physiology, and uncertainty around effective tau. |
| 3 | Objectives and prespecified hypotheses | Partial | State primary hypothesis for tau and exploratory hypotheses for ABP features. |
| 4 | Key design elements | Ready | Describe retrospective within-patient matched case-control design early in Methods. |
| 5 | Setting, location, dates | Partial | Name MIMIC-IV Waveform v0.1.0, BIDMC ICU source, database publication/version, and analysis dates. |
| 6 | Eligibility, case ascertainment, control selection, matching | Ready | Report segment duration, signal quality, hypotension rules, first eligible event, one control per case, same patient/segment, and >=25-min separation. |
| 7 | Outcomes, exposures, predictors, confounders | Partial | Define hypotension, effective tau, ABP features; explicitly state unavailable clinical confounders. |
| 8 | Data sources and measurement | Ready | Report WFDB extraction, sampling, quality thresholds, beat detection, and tau fitting equation. |
| 9 | Bias | Partial | Discuss signal-quality selection, tau estimability, event-selection, temporal confounding, and lack of clinical covariates. |
| 10 | Study size | Ready with caveat | State convenience sample determined by available eligible waveforms; do not report retrospective power. |
| 11 | Quantitative variables | Ready | Report continuous treatment, medians/IQR, slopes, thresholds, and no data-driven categorization. |
| 12 | Statistical methods | Ready | Wilcoxon, 10,000 paired bootstrap resamples, rank-biserial correlation, Benjamini-Hochberg FDR, paired complete-case availability. |
| 13 | Participant flow | Ready | Use Figure 1 and report 64, 54, 52, 26, 23, and 15 at relevant stages. |
| 14 | Descriptive characteristics and missingness | Major gap | No linked clinical table; report waveform characteristics and missingness, and state demographic/clinical characteristics unavailable. |
| 15 | Outcome/exposure summary | Ready | Report event counts, eligible events, paired windows, and tau availability. |
| 16 | Main estimates and precision | Ready | Report paired medians, differences, bootstrap 95% CIs, Wilcoxon p, FDR-adjusted p, and effect sizes. |
| 17 | Other analyses | Partial | Report tau model-assumption sensitivity; label remaining features exploratory. |
| 18 | Key results | Pending writing | Summarize no statistically sufficient evidence and distinguish exploratory MAP trend. |
| 19 | Limitations | Pending writing | Address sample size, 15 tau pairs, no clinical covariates, signal selection, no external validation, and model assumptions. |
| 20 | Interpretation | Pending writing | Use cautious non-causal language; do not claim prediction or clinical utility. |
| 21 | Generalisability | Pending writing | Limit inference to similar ICU invasive ABP waveform settings. |
| 22 | Funding and sponsor role | Pending institutional input | Authors must provide actual funding/support and sponsor role; do not invent. |

## Blocking items before submission

1. Obtain written UTN ethics/exemption determination.
2. Confirm author list, affiliations, corresponding author, and CRediT roles.
3. Confirm funding, conflicts of interest, and acknowledgments.
4. Decide whether to add linked clinical information in a new analysis version or explicitly retain waveform-only scope.
5. Complete visual review of all four final figures.
6. Correct or split the mixed-unit forest plot before manuscript submission.
7. Complete the STROBE checklist with final manuscript page and line numbers.
8. Create a GitHub release from analysis-v1.0 and archive it with a DOI if available.
9. Add an explicit code license without implying redistribution rights over MIMIC data.
10. Select a target journal and adapt formatting, word count, abstract, and supplementary files.

## Frozen primary result

Among 15 paired patients with acceptable effective tau estimates in both windows, the median paired difference (case minus control) was -0.0262 s (bootstrap 95% CI -0.1204 to 0.0198 s; Wilcoxon p=0.1876; Benjamini-Hochberg adjusted p=0.5469). This does not provide statistically sufficient evidence of a pre-hypotension difference in effective tau.
