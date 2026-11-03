# Dorsal and Ventral White Matter Correlates of Word-Finding Difficulty in Neurologically Intact Adults

`analysis.py` computes age-adjusted Spearman partial correlations between word-finding difficulty and FA, MD, RD and AD in ten language tracts (one Benjamini–Hochberg correction across all 40 tests) and grey matter volume in 28 language regions, together with the age correlations, all 40 tract tests with added covariates and with naming response time as a covariate, and robustness checks of the significant tract results (naming response time as the outcome, permutation, winsorized and trimmed outcome). It writes the aggregate statistics (n, ρ, p, FDR q) to one CSV.

## Data availability

The data come from the Aging Brain Cohort (ABC) repository at the University of South Carolina. They are not openly available because of institutional restrictions and can be requested at https://abc.sc.edu/abc-repository-data-requests/.

## Input columns

One row per participant:

- `Age`, `sex` (`Male` or `Female`), `education` (ordinal)
- `Visual - TOT`: word-finding difficulty score
- `Visual - Mean RT`: mean naming response time (rows with 0 are excluded from the response-time analyses)
- `proc_speed_osd`: NIH Toolbox Oral Symbol Digit score
- `<tract>_dti_fa`, `<tract>_md`, `<tract>_rd`, `<tract>_ad` for `AF_L`, `AF_R`, `IFOF_L`, `IFOF_R`, `ILF_L`, `ILF_R`, `SLF_L`, `SLF_R`, `UF_L`, `UF_R`
- `<region> left volume %` and `<region> right volume %` for the 14 regions in `GM_REGIONS` in `analysis.py`

## Run

    pip install -r requirements.txt
    python analysis.py input.csv output.csv
