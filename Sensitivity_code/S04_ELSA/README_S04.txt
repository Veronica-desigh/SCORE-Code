S04 / ELSA — integrated workflow (2026-11-13)
============================================================

PURPOSE
  One substantive analysis per R script. Each script contains its own
  configuration, helper functions, analysis, and export/QA logic.
  Only the existing, validated M01 SCORE/SCORE2 algorithm is an
  external SCIENTIFIC dependency for Step 02; Step 03 no longer needs
  Function_08_Summarize_risk_transitions_S03.R.

ORDER
  01_Data_quality_S04.R                  source descriptions only
  02_Guideline_strategies_S04.R          freeze complete-case score cohort;
                                            SCORE/SCORE2, fixed S0 treatment
  03_Risk_transitions_S04.R             fixed S0 risk/6-group summaries
  04_RM_imputation_distribution_S04.R   ELSA RM MICE, distribution/diagnostics
  05_RM_strategies_S04.R                 strict/expanded RM; S0/SN/ST/SNT
  06_RM_tables_matrices_contributions_S04.R descriptive counts/matrices/contributions
  07_RM_vector_matrices_S04.R           vector PDF matrices
  00_Run_RM_05_07_S04.R                optional runner for 05-07 only

DATA
  Input directory: <PROJECT_ROOT>/DATA/ELSA/
    ELSA_data.rds, ELSA_variable.csv
  ELSA ID is idauniq everywhere outside Step 02's internal adapter.

FIXED INPUT/OUTPUT CONTRACTS
  02 CSV: LJW_outcome/data_intervention/S04_ELSA/complete_case/
          02_score_strategy_participants_S04.csv
  04 RDS: LJW_outcome/data_intervention/S04_ELSA/rm_imputation/
          04_RM_imputed_cohorts_S04.rds
  04 QC : same rm_imputation/quality_control/04_RM_run_metadata_S04.csv
  05/06/07: LJW_outcome/data_intervention/S04_ELSA/external_RM_analyses/
          01_individual_states; 02_per_imputation; 03_pooled;
          04_vector_figures; quality_control; logs

DEFAULTS
  M=5, MAXIT=5, THREADS=5, 04 seed=20261113.
  Override: S04_RM_M, S04_RM_MAXIT, S04_RM_THREADS,
            S04_PROJECT_ROOT, S04_INPUT_DIR, S04_ANALYSIS_DIR.
  Step 02 workers: S04_WORKERS=5; batch size: S04_BATCH_SIZE=1000.
  Step 02 M01 engine override: S04_M01_SCORING_FUNCTION=<absolute .R path>.
  Optional 05 RM mapping override: S04_RM_MAP=<CSV with domain,column>.

IMPORTANT RM RULES
  Only Wave 6 hs-CRP >2 mg/L; no Wave 4 or persistence.
  Fixed observed clinical fields are not imputed.
  Other ELSA observed/evaluable RM sources support MICE prediction.
  Female menopause is a separate conditional MICE model.
  A male's menopause field is not changed: for the *RM union only*,
  this sex-inapplicable modifier is a non-trigger (neutral 0).
  RM domains with no observed support are omitted ONLY from the
  ELSA-ASSESSABLE scenario domain set, and explicitly logged;
  they are not represented as observed negative.
  If an assessable domain remains unresolved for any baseline-moderate
  participant, 05 writes coverage QA and STOPS; it does not silently
  impute unsupported data, change denominators, or mark NA as negative.
  Strict nonlaboratory RM = obesity only; expanded = all ELSA assessable
  nonlaboratory RM. These are sensitivity scenarios, not official
  guideline-prescribed deterministic rules.
  2019 ST does not exist; comparisons explicitly use 2019 S0.
  SNT in 2019 = SN.

TEST MODE
  For 04, S04_RM_TEST_N=<positive integer> selects a reproducible
  sample within the fixed Step 02 cohort. Its RDS and metadata are in
  rm_imputation/test_N/. For 05-07, set S04_RM_TEST_N to the same N.
  05 matches the test MI IDs back to the FULL Step 02 score CSV;
  05-07 test reports are stored under external_RM_analyses/test_N/.
  Never combine test and full cohort MI files.

EVIDENCE/CAUTIONS
  06 pooled counts and percentages are arithmetic means across five
  imputed datasets, NOT Rubin inferential estimates or confidence intervals.
  01 audits source, while 04 audits the fixed scoring cohort; denominators
  can differ. ELSA RM proxy definitions do not restore missing UKB domains.
  05 validates the frozen S0 LDL-C treatment rules before applying RM.
  02 requires the canonical M01 scoring script to be present.

VALIDATION STATUS
  Source-level interface reconciliation and delimiter checks performed.
  No R interpreter installed in this environment: scripts are NOT
  certified by R parse() or end-to-end numeric runs. Run locally and
  check the logs and quality_control outputs before using research results.
