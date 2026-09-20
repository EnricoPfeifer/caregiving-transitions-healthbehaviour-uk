# Caregiving and Health Behaviours over the Lifecourse

Replication code for:

> Pfeifer, E., Lacey, R., Xue, B., Pikhart, H., & McMunn, A. (in press). Health behaviour changes around transitions into and out of unpaid caregiving in the UK: a longitudinal study. *SSM – Population Health.*

[Add the DOI here once the article is online.]

---

## Overview

This repository contains the Stata code used to produce all analyses reported in the paper and its supplementary material. The study examines longitudinal associations between transitions into and out of unpaid caregiving and four health behaviours: physical inactivity, fruit and vegetable consumption, problematic alcohol consumption, and smoking.

Two complementary strategies are used. Within-person fixed effects models are the primary specification. Piecewise growth curve models estimated on a propensity score matched sample are reported as a complementary robustness check.

The analyses draw on the UK Household Longitudinal Study (UKHLS, Understanding Society). Smoking is analysed using waves 5 to 13; physical inactivity, diet and alcohol consumption are analysed using waves 7, 9, 11 and 13, as those measures were collected biennially from wave 7 onwards.

---

## Data access

The UKHLS data are not included in this repository and cannot be redistributed by the authors. They are available to registered researchers from the UK Data Service (https://ukdataservice.ac.uk) under standard access conditions:

> University of Essex, Institute for Social and Economic Research. (2023). *Understanding Society: Waves 1–13, 2009–2022 and Harmonised BHPS: Waves 1–18, 1991–2009* [data collection]. UK Data Service. SN: 6614.
> https://doi.org/10.5255/UKDA-SN-6614-19

To reproduce these analyses, download the Stata (`.dta`) files and set the `ukhls` global at the top of each do-file to the directory holding them.

---

## Software

Stata 17 or later. Earlier versions may work but have not been tested.

The following user-written packages are required:

```stata
ssc install kmatch,               replace
ssc install ebalance,             replace
ssc install grstyle,              replace
ssc install palettes,             replace
ssc install colrspace,            replace
ssc install winsor2,              replace
ssc install mplotoffset,          replace
ssc install combomarginsplot,     replace
ssc install mdesc,                replace
ssc install how_many_imputations, replace
```

The distributional checks in `08_sensitivity.do` additionally require `prcounts`, which is distributed with SPost (Long and Freese) rather than through SSC. [Add the installation command you used.]

---

## Repository structure

```
code/
  00_graph_settings.do      Graph style settings for all publication figures
  01_extract_wide.do        Extract variables from UKHLS and build wide dataset
  02_clean.do               Clean variables and construct analysis dataset
  03_fixed_effects.do       Fixed effects models (entry and exit)
  04_psm_entry.do           PSM growth curve models for caregiving entry
  05_psm_entry_strata.do    PSM entry models stratified by sex and age group
  06_psm_exit.do            PSM growth curve models for caregiving exit
  07_psm_exit_strata.do     PSM exit models stratified by sex and age group
  08_sensitivity.do         Continuous and standardised re-estimation, and
                            distributional checks for the AUDIT-C outcome
  09_multiple_imputation.do Multiple imputation on the matched sample

data/                       Data access instructions (no data files stored here)
output/                     Output directory (no results files stored here)
```

---

## Configuration

Each do-file has a short path block near the top. Before running a file, set the path to the folder holding your UKHLS data:

```stata
global ukhls "/path/to/UKHLS/stata/stata13_se/ukhls"
```

The same line appears in every do-file, so set it once per file or set it once at the start of a Stata session and it will apply to all of them. Figures and margins files are written to the working directory unless you change it.

---

## How to reproduce the analyses

Each do-file is run by setting a small number of local macros at the top. The workflow for a single outcome is as follows. Repeat for each of the four outcomes: `pa`, `alc`, `diet`, `smok`.

1. Graph style, once per session

```stata
do 00_graph_settings.do
```

2. Extract wide-format data

```stata
* Set: local outcome "pa"   // pa | alc | diet | smok
do 01_extract_wide.do
```

3. Clean and construct analysis variables

```stata
* Set: local outcome "pa"
do 02_clean.do
```

4. Fixed effects models

```stata
* Set: local outcome "pa"
do 03_fixed_effects.do
```

5. PSM growth curve models, entry

```stata
* Set: local outcome "pa"
do 04_psm_entry.do
```

6. PSM entry, sex and age strata

```stata
* Set: local outcome "pa"
*      local mode    "sex_strat"   // sex_int | sex_strat | age_int | age_strat
*      local strat   1             // sex: 0=men, 1=women; age: 0=16-29, 1=30-49, 2=50-64, 3=65+
do 05_psm_entry_strata.do
```

7. PSM growth curve models, exit

```stata
* Set: local outcome "pa"
*      local control "nocare"      // nocare | alwayscare
do 06_psm_exit.do
```

8. PSM exit, sex and age strata

```stata
* Set: local outcome "pa"
*      local control "nocare"
*      local mode    "sex_strat"
*      local strat   1
do 07_psm_exit_strata.do
```

9. Sensitivity analyses (supplementary material)

```stata
* Set: local block   "fe"      // fe | psm | dist
*      local outcome "pa"      // pa | alc
*      local sample  "entry"   // entry | exit_nocare | exit_alwayscare
do 08_sensitivity.do
```

10. Multiple imputation (supplementary material)

```stata
* Set: local outcome "pa"      // pa | alc | diet | smok
*      local sample  "entry"   // entry | exit_nocare | exit_alwayscare
*      local M       10
do 09_multiple_imputation.do
```

---

## Outcome variable summary

| Local macro | Outcome | Waves used | Model type |
|---|---|---|---|
| `pa` | Physical inactivity (binary) | g, i, k, m (7, 9, 11, 13) | `xtlogit` / `logit` |
| `alc` | Problematic drinking (binary, AUDIT-C ≥3 women / ≥4 men) | g, i, k, m | `xtlogit` / `logistic` |
| `diet` | Mean daily portions fruit and vegetables (continuous) | g, i, k, m | `xtreg` / `regress` |
| `smok` | Current smoker (binary) | e–m (5–13) | `xtlogit` / `logistic` |

---

## Key analytical decisions

Fixed effects models (`03_fixed_effects.do`). Within-person fixed effects logit for the binary outcomes and linear models for diet, with factor wave dummies. Entry into and exit from caregiving are modelled separately, each on its own analytic sample. Sex and age-group interactions are tested jointly with `testparm`.

PSM models (`04_psm_entry.do`, `06_psm_exit.do`). Kernel propensity score matching (`kmatch ps`), 1:3 nearest-neighbour matching with replacement, exact matching on sex, age at first observation and wave of first observation. Entropy balancing (`ebalance`) is applied after matching. Piecewise growth curves (`mkspline`) are centred on the transition wave, with spline knots at waves 2 and 3 for the four-wave outcomes and at waves 6 and 7 for smoking.

Stratified analyses (`05_psm_entry_strata.do`, `07_psm_exit_strata.do`). Sex or age is removed from matching and balancing when the corresponding interaction is tested. For stratified models the sample is restricted to the stratum of interest before matching.

Exit analysis (`06_psm_exit.do`, `07_psm_exit_strata.do`). Two comparison groups are used: never-caregivers (`nocare`) and continuing caregivers who never exited (`alwayscare`). Caregiving hours and place of care are taken from the wave immediately prior to exit using a forward-fill procedure.

Sensitivity analyses (`08_sensitivity.do`). Physical activity and AUDIT-C are dichotomised in the main analysis. Both are re-estimated on their underlying continuous scale and again after standardisation, in the fixed effects models and in the matched-sample growth curves. Physical activity minutes are winsorised at the 99th percentile before standardisation; AUDIT-C is not. Standardisation uses `egen std()` on the sample in memory, so coefficients are in standard deviations of the sample used in that model. Observed and predicted distributions for AUDIT-C are compared across Poisson, negative binomial and zero-inflated specifications.

Multiple imputation (`09_multiple_imputation.do`). Listwise deletion at the matching stage removes a non-trivial share of the sample, whereas the fixed effects models are structurally robust to missingness on time-invariant covariates. Imputation is therefore applied to the matched analysis. Following the "within" approach of Leyrat et al. (2019), matching, balancing and estimation are repeated inside each imputed dataset and the estimates pooled with Rubin's rules. Only baseline covariates and the baseline outcome enter the imputation model; repeated wave-specific outcomes do not. Confidence intervals use the normal approximation rather than a Barnard-Rubin degrees of freedom correction. The number of imputations used is reported in the supplementary material.

---

## Notes on reproducibility

- Seeds are set in the matching and imputation files. See the individual do-files for the values used.
- Results may differ slightly across Stata versions owing to numerical precision in iterative estimation.
- Combined publication figures spanning several strata are wrapped in comment blocks and should be run only once all individual strata runs are complete and their margins files exist on disk.
- The multiple imputation routine repeats matching and balancing for every imputation, so it is computationally heavy. Run it with a small number of imputations first to confirm it executes.

---

## Citation

If you use this code, please cite both the paper and the repository:

> Pfeifer, E., Lacey, R., Xue, B., Pikhart, H., & McMunn, A. (in press). Health behaviour changes around transitions into and out of unpaid caregiving in the UK: a longitudinal study. *SSM – Population Health.*

> Pfeifer, E. (2026). *Replication code: health behaviour changes around caregiving transitions in the UK* (Version 1.0.0) [Computer software]. Zenodo. [DOI]

---

## Funding

This project is funded by the UK Economic and Social Research Council (ESRC) through the UBEL-DTP (UCL, Bloomsbury and East London Doctoral Training Partnership), grant reference ES/P000592/1. The funder was not involved in project design, analysis or write up of findings.

AM and BX were supported by the Joint Programming Initiative More Years, Better Lives (JPI MYBL) through the ESRC (ES/W001454/1). AM, BX and EP are supported by ESRC funding for Equalise: ESRC Centre for Lifecourse Health Equity (ES/Z504270/1).

---

## Licence

This code is released under the [MIT Licence](LICENSE). The underlying UKHLS data are subject to the UK Data Service licence conditions and cannot be redistributed.

---

## Contact

Enrico Pfeifer, Institute of Epidemiology and Healthcare, University College London, 1-19 Torrington Place, London WC1E 7HB, United Kingdom.

For questions about the code, please open an issue on this repository.
