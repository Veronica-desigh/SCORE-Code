# S03 UKHLS — integrated 01–10 code package (2026-11-13)

This directory contains **only S03/UKHLS** code. It does not call S04 or ELSA scripts. Each numbered `.R` file contains its own functions and main execution.

## Order

- `00_Run_all_S03.R`: single entrypoint; configure `THREADS <- 5L`, `RUN_STEPS <- 1:10`, `M <- 20L`, `MAXIT <- 10L`, paths and `SMOKING_MODE` at top.
- 01–03: original interview, nurse, and blood cleaning (original operational definitions retained).
- 04: eligibility and score-complete cohort, preserving RM missingness.
- 05: SCORE/SCORE2 and fixed S0 treatment states, using the pre-existing M01 scoring function **only**; checks its MD5.
- 06: S0 six transition groups and risk-level descriptions.
- 07: RM-only MICE with fixed observed clinical predictors and eligible other RM source variables. Cross-RM prediction excludes variables having unresolved structural missing cells, preventing unobserved values from being silently encoded as negative. BMI, waist and hsCRP are imputed on their continuous scale, then the RM binaries are derived. It generates predictor edges and multi-imputation spread QC.
- 08: strict / expanded RM sets; validates observed SCORE and input hash; unresolved expanded RM can be skipped without suppressing strictly identifiable strict RM results. No RM unknown state is forced to zero. The `08_RM_set_run_status.csv` and `08_RM_executed_sets.csv` indicate which set has been computed.
- 09: risk classes, H/VH, treatment eligibility, treatment transitions, risk and six-group matrices, fixed original S0 group contribution additivity, and inter-imputation variability checks.
- 10: risk reclassification and treatment group transition vector PDFs.

## Launch / resume

Store all eleven files together under `LJW_code/Sensitivity_code/S03_UKHLS/` and run the **entire** `00_Run_all_S03.R` by Source or Rscript. In the RStudio Console (executing selected lines), specify `SCRIPT_DIR` explicitly. `PROJECT_ROOT` is the parent of `LJW_code`, `DATA`, and `LJW_outcome`. For resumed work use `RUN_STEPS <- 7:10` or `8:10` with matching existing files. Re-run 07 before 08 if the MICE predictive structure is revised.

## Statistical qualifications

- The default is 20 MICE chains, 10 iterations, 5 PSOCK workers, not 5 MICE chains. Adjust only in the master header. The default `MIN_DONORS` and `MIN_BINARY_EVENTS` are minimum numeric safeguards, **not** evidence of a well-identified model.
- Wave 3 physical activity with no comparable adult questionnaire is deliberately **not imputed**; if this prevents identifying expanded RM decisions, expanded is marked as not runnable, while strict may be retained. This is not an assumption of inactivity or activity.
- MI-mean counts are descriptive expected counts, potentially noninteger. No outcomes, risk difference, or confidence intervals are estimated.
- `Function_07_Guideline_strategies_S03` and `Function_08_Summarize_risk_transitions_S03` are old, no longer required. The only external source dependency remaining is M01 SCORE engine in 05.
- Encodings: R sources UTF-8; CSV reports UTF-8 with BOM for Windows Excel.

## Mandatory verification before publication

The R interpreter and the UKHLS raw/intermediate data were not available to execute an end-to-end analysis here. Run `parse()`/smoke test on your R installation; compare 05–06 to the previous S0 output; review `07_RM_predictor_matrix_S03.csv`, `07_RM_predictor_edges_S03.csv`, chained trace / events, `08_RM_set_run_status.csv` and 09 additivity checks. A successful R job alone does not establish the validity of the MAR assumption, model convergence or transportability.
