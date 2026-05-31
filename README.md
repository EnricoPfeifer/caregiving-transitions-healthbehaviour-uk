# Caregiving and Health Behaviours over the Lifecourse

## Acknowledgements

This code was developed by Enrico Pfeifer. Portions of the analysis 
code were adapted from code originally written by Rebecca Lacey and 
Baowen Xue (UCL Department of Epidemiology and Public Health), to 
whom the authors are grateful.

**Replication code for:**
Pfeifer, E. et al. (under review). Health behaviour changes around transitions into and out of unpaid caregiving in the UK: a longitudinal study. *SSM – Population Health.*

---

## Overview

This repository contains the Stata code used to produce all analyses in the above paper. The study examines the longitudinal associations between transitions into and out of informal caregiving and four health behaviours — physical inactivity, fruit and vegetable consumption, problematic alcohol consumption, and smoking — across the adult lifecourse in the UK.

The analyses draw on data from the **UK Household Longitudinal Study (UKHLS; Understanding Society)**, waves 1–13 (2009–2023). Methods include within-person fixed effects models and propensity score matching with piecewise growth curve models.

---

## Data access

The UKHLS data are not included in this repository. They are available under a standard end-user licence from the **UK Data Service**:

> University of Essex, Institute for Social and Economic Research. (2023). *Understanding Society: Waves 1–14, 2009–2023 and Harmonised BHPS: Waves 1–18, 1991–2009* [data collection]. UK Data Service. SN: 6614.
> https://doi.org/10.5255/UKDA-SN-6614-19

To reproduce these analyses, download the Stata (`.dta`) files and place them in the directory you assign to `$ukhls` in each do-file (see **Configuration** below).

---

## Software

- **Stata** 17 or later (earlier versions may work but have not been tested)
- The following user-written packages are required and can be installed from SSC:

```stata
ssc install kmatch,      replace
ssc install ebalance,    replace
ssc install grstyle,     replace
ssc install palettes,    replace
ssc install colrspace,   replace
ssc install winsor2,     replace
ssc install mplotoffset, replace
ssc install combomarginsplot, replace
ssc install fre,         replace
ssc install twopm,       replace
```

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

data/
  README.md                 Data access instructions (no data files stored here)

output/
  README.md                 Output directory (no results files stored here)
```
