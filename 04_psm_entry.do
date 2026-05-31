/*===========================================================================
  04_psm_entry.do

  PURPOSE: Propensity score matching (PSM) and piecewise growth curve models
           for caregiving ENTRY. Matches individuals transitioning into
           caregiving to non-caregivers, then estimates the trajectory of
           the health behaviour outcome before and after caregiving onset.

  PREREQUISITE: 02_clean.do must have been run for this outcome.

  OUTCOMES SUPPORTED:
    pa    -- physical inactivity   (logit;    waves g i k m; spline knots 2/3)
    alc   -- problematic drinking  (logistic; waves g i k m; spline knots 2/3)
    diet  -- fruit/veg portions    (regress;  waves g i k m; spline knots 2/3)
    smok  -- smoking               (logistic; waves e-m;     spline knots 6/7)

  MATCHING METHOD: kmatch ps (kernel propensity score matching), nn(3),
                   with replacement (wor), exact match on sex, age at
                   first observation, and wave of first observation.
                   Balance weights generated with ebalance.

  NOTE: Paths below use the author's local machine. Replace with your own
        paths before running. These will be removed prior to publication.

  AUTHOR: Enrico Pfeifer
===========================================================================*/

*---------------------------------------------------------------------------
* SET OUTCOME HERE
* Options: pa | alc | diet | smok
*---------------------------------------------------------------------------
local outcome "pa"


*---------------------------------------------------------------------------
* PATHS  --  replace with your own directory structure
*---------------------------------------------------------------------------
clear all
set more off
set maxvar 10000
set cformat %4.2f

cd "/Users/Enrico/Documents/UCL/UBEL-DTP/Data/UKHLS/stata/stata13_se/ukhls"
global ukhls "/Users/Enrico/Documents/UCL/UBEL-DTP/Data/UKHLS/stata/stata13_se/ukhls/"


*===========================================================================
* SECTION 1: LOAD DATA AND SAMPLE SELECTION
*===========================================================================

if "`outcome'" == "smok" {
    use "$ukhls/smok_wide_beforeps", clear
}
else {
    use "$ukhls/`outcome'_wide_cleaned", clear
}

keep if ecare != .
drop if care_fiob == 1          // exclude baseline caregivers

gen wave_fiob = fiob


*---------------------------------------------------------------------------
* Last wave of observation (used to check pre/post requirement)
*---------------------------------------------------------------------------
if "`outcome'" == "smok" {
    gen lastwave = 9 if m_carebi != .
    replace lastwave = 8 if l_carebi != . & m_carebi == . & lastwave == .
    replace lastwave = 7 if k_carebi != . & l_carebi == . & lastwave == .
    replace lastwave = 6 if j_carebi != . & k_carebi == . & lastwave == .
    replace lastwave = 5 if i_carebi != . & j_carebi == . & lastwave == .
    replace lastwave = 4 if h_carebi != . & i_carebi == . & lastwave == .
    replace lastwave = 3 if g_carebi != . & h_carebi == . & lastwave == .
    replace lastwave = 2 if f_carebi != . & g_carebi == . & lastwave == .
    replace lastwave = 1 if e_carebi != . & f_carebi == . & lastwave == .
}
else {
    gen lastwave = 4 if m_carebi != .
    replace lastwave = 3 if k_carebi != . & m_carebi == . & lastwave == .
    replace lastwave = 2 if i_carebi != . & k_carebi == . & lastwave == .
    replace lastwave = 1 if g_carebi != . & i_carebi == . & lastwave == .
}

*--- Keep only caregivers with at least one wave before and one after onset ---
gen samplepiece = 1 if (wave_fiob < ficare & ficare != .) & (lastwave >= ficare)
drop if samplepiece != 1 & ecare == 1

*--- Keep only controls with at least two waves of observation ---
drop if nwaves == 1 & ecare == 0


*===========================================================================
* SECTION 2: BASELINE VARIABLES
*===========================================================================

*--- Baseline household ID ---
if "`outcome'" == "smok" {
    gen hhid_fiob = e_hidp if fiob == 1
    replace hhid_fiob = f_hidp if fiob == 2
    replace hhid_fiob = g_hidp if fiob == 3
    replace hhid_fiob = h_hidp if fiob == 4
    replace hhid_fiob = i_hidp if fiob == 5
    replace hhid_fiob = j_hidp if fiob == 6
    replace hhid_fiob = k_hidp if fiob == 7
    replace hhid_fiob = l_hidp if fiob == 8
    replace hhid_fiob = m_hidp if fiob == 9
}
else {
    gen hhid_fiob = g_hidp if fiob == 1
    replace hhid_fiob = i_hidp if fiob == 2
    replace hhid_fiob = k_hidp if fiob == 3
    replace hhid_fiob = m_hidp if fiob == 4
}
sort hhid_fiob pidp
by hhid_fiob: gen nucarehh = _N

*--- Outcome at first observation ---
if "`outcome'" == "pa" {
    gen pabi_fiob  = g_pabi  if fiob == 1
    replace pabi_fiob  = i_pabi  if fiob == 2
    replace pabi_fiob  = k_pabi  if fiob == 3
    replace pabi_fiob  = m_pabi  if fiob == 4
    gen patot_fiob = g_patot if fiob == 1
    replace patot_fiob = i_patot if fiob == 2
    replace patot_fiob = k_patot if fiob == 3
    replace patot_fiob = m_patot if fiob == 4
    gen watot_fiob = g_watot if fiob == 1
    replace watot_fiob = i_watot if fiob == 2
    replace watot_fiob = k_watot if fiob == 3
    replace watot_fiob = m_watot if fiob == 4
}
else if "`outcome'" == "alc" {
    gen alcbi_fiob    = g_alcbi    if fiob == 1
    replace alcbi_fiob    = i_alcbi    if fiob == 2
    replace alcbi_fiob    = k_alcbi    if fiob == 3
    replace alcbi_fiob    = m_alcbi    if fiob == 4
    gen sumaudit_fiob = g_sumaudit if fiob == 1
    replace sumaudit_fiob = i_sumaudit if fiob == 2
    replace sumaudit_fiob = k_sumaudit if fiob == 3
    replace sumaudit_fiob = m_sumaudit if fiob == 4
}
else if "`outcome'" == "diet" {
    gen diet_fiob     = g_diet     if fiob == 1
    replace diet_fiob     = i_diet     if fiob == 2
    replace diet_fiob     = k_diet     if fiob == 3
    replace diet_fiob     = m_diet     if fiob == 4
    gen meandiet_fiob = g_meandiet if fiob == 1
    replace meandiet_fiob = i_meandiet if fiob == 2
    replace meandiet_fiob = k_meandiet if fiob == 3
    replace meandiet_fiob = m_meandiet if fiob == 4
}
else if "`outcome'" == "smok" {
    gen smok_fiob = e_smok if fiob == 1
    replace smok_fiob = f_smok if fiob == 2
    replace smok_fiob = g_smok if fiob == 3
    replace smok_fiob = h_smok if fiob == 4
    replace smok_fiob = i_smok if fiob == 5
    replace smok_fiob = j_smok if fiob == 6
    replace smok_fiob = k_smok if fiob == 7
    replace smok_fiob = l_smok if fiob == 8
    replace smok_fiob = m_smok if fiob == 9
}

*--- Caregiving hours and place at onset ---
if "`outcome'" == "smok" {
    capture drop carehrs_fiob
    gen carehrs_fiob = e_carehrs if ficare == 1
    replace carehrs_fiob = f_carehrs if ficare == 2
    replace carehrs_fiob = g_carehrs if ficare == 3
    replace carehrs_fiob = h_carehrs if ficare == 4
    replace carehrs_fiob = i_carehrs if ficare == 5
    replace carehrs_fiob = j_carehrs if ficare == 6
    replace carehrs_fiob = k_carehrs if ficare == 7
    replace carehrs_fiob = l_carehrs if ficare == 8
    replace carehrs_fiob = m_carehrs if ficare == 9
    gen carecat_fiob = f_carecat if ficare == 2
    replace carecat_fiob = g_carecat if ficare == 3
    replace carecat_fiob = h_carecat if ficare == 4
    replace carecat_fiob = i_carecat if ficare == 5
    replace carecat_fiob = j_carecat if ficare == 6
    replace carecat_fiob = k_carecat if ficare == 7
    replace carecat_fiob = l_carecat if ficare == 8
    replace carecat_fiob = m_carecat if ficare == 9
}
else {
    gen carehrs_fiob = i_carehrs if ficare == 2
    replace carehrs_fiob = k_carehrs if ficare == 3
    replace carehrs_fiob = m_carehrs if ficare == 4
    gen carecat_fiob = i_carecat if ficare == 2
    replace carecat_fiob = k_carecat if ficare == 3
    replace carecat_fiob = m_carecat if ficare == 4
}
label val carehrs_fiob carehval
recode carecat_fiob (0=.)
label val carecat_fiob carecatval

*--- Intensity binary (low / high) ---
gen intense = carehrs_fiob
recode intense (1=0) (2=1)
label define intenseval 0 "Low (<20 hrs/week)" 1 "High (>=20 hrs/week)"
label val intense intenseval

*--- Ever smoked / ever not smoked (smoking only) ---
if "`outcome'" == "smok" {
    egen eversmok    = anycount(*_smok), values(1)
    recode eversmok 1/max=1
    egen smokmiss    = rowmiss(*_smok)
    replace eversmok = . if smokmiss == 9
    label var eversmok "Ever smoked"
    egen everNOsmoke = anycount(*_smok), values(0)
    recode everNOsmoke 1/max=1
    replace everNOsmoke = . if smokmiss == 9
    label var everNOsmoke "Ever did not smoke"
}

sort pidp


*===========================================================================
* SECTION 3: PROPENSITY SCORE MATCHING (kmatch)
*===========================================================================
* Matches non-caregivers (ecare=0) to caregivers (ecare=1)
* Exact match on: sex, age at first observation, wave of first observation
* Caliper matching on: sociodemographic covariates
* Outcome at baseline included as a conditioning variable (except smoking)
* nn(3): 1:3 matching; wor: with replacement

if "`outcome'" == "pa" {
    kmatch ps ecare i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob ///
        sf_fiob i.children_fiob (pabi_fiob), ///
        ematch(sex wave_fiob agefiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else if "`outcome'" == "alc" {
    kmatch ps ecare i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob ///
        (alcbi_fiob), ///
        ematch(sex wave_fiob agefiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else if "`outcome'" == "diet" {
    kmatch ps ecare i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob ///
        (meandiet_fiob), ///
        ematch(sex agemin agemax) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else if "`outcome'" == "smok" {
    kmatch ps ecare i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob, ///
        ematch(sex agefiob wave_fiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}

tab kpscore
gen matchwt = _KM_mw
gen pscore  = _KM_ps


*===========================================================================
* SECTION 4: ASSIGN ONSET WAVE TO MATCHED CONTROLS
* Each matched control receives the onset wave of their matched carer
* (ficare_all). This defines the time scale for growth curves.
*===========================================================================

*--- Pass caregiver onset wave to matched controls (up to 3 controls) ---
forval k = 1/3 {
    gen double id`k' = controlid`k'
    replace id`k' = pidp if kpscore == 0
    gen ficare_tmp`k' = ficare if kpscore == 1
    sort id`k' ficare_tmp`k'
    by id`k': replace ficare_tmp`k' = ficare_tmp`k'[1]
    gen ficare_g`k' = ficare_tmp`k' if kpscore == 0
    drop ficare_tmp`k'
}

*--- Resolve onset wave for controls using matched carers ---
egen occ_carer = anycount(*_carebi), values(1)
replace ficare = 0 if occ_carer == 0

gen ficare_all = ficare

*   Assign matched carer onset; check pre/post requirement; cascade to
*   next control if condition not met
replace ficare_all = ficare_g1 if ficare_g1 != . & ficare_all == 0
gen preafter0 = 1 if (wave_fiob < ficare_all) & (lastwave >= ficare_all) & ficare_all != .
replace ficare_all = ficare_g2 if preafter0 != 1
gen preafter  = 1 if (wave_fiob < ficare_all) & (lastwave >= ficare_all) & ficare_all != .
replace ficare_all = ficare_g3 if preafter != 1

*--- Final usesample flag ---
gen usesample = 1 if (wave_fiob < ficare_all) & (lastwave >= ficare_all) & ficare_all != .
keep if usesample == 1


*===========================================================================
* SECTION 5: COVARIATE BALANCE (ebalance)
* Generates balance weights for the main model (balancewt) and for the
* intensity sub-analysis (intensewt)
*===========================================================================

if "`outcome'" == "pa" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity sf_fiob i.children_fiob (pabi_fiob), gen(balancewt)
    ebalance intense agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity sf_fiob i.children_fiob (pabi_fiob), gen(intensewt)
}
else if "`outcome'" == "alc" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity (alcbi_fiob), gen(balancewt)
    ebalance intense agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity (alcbi_fiob), gen(intensewt)
}
else if "`outcome'" == "diet" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob workbi_fiob ///
        i.ethnicity (meandiet_fiob), gen(balancewt)
    ebalance intense agefiob i.sex i.cohab_fiob nwaves ///
        iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity (meandiet_fiob), gen(intensewt)
}
else if "`outcome'" == "smok" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity (smok_fiob), gen(balancewt)
    ebalance intense agefiob i.sex i.cohab_fiob nwaves ///
        iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity (smok_fiob), gen(intensewt)
}

sort pidp

*--- Save wide post-PSM dataset ---
save "$ukhls/`outcome'_wide_afterpsm", replace


*===========================================================================
* SECTION 6: RESHAPE TO LONG FORMAT
*===========================================================================

use "$ukhls/`outcome'_wide_afterpsm", clear
sort pidp

if "`outcome'" == "smok" {
    rename e_* *5
    rename f_* *6
    rename g_* *7
    rename h_* *8
    rename i_* *9
    rename j_* *10
    rename k_* *11
    rename l_* *12
    rename m_* *13
    reshape long smok ciggy age_dv carebi carehrs carecat, i(pidp) j(waves)
    gen wave = waves - 4         // wave runs 1-9
}
else {
    rename g_* *1
    rename i_* *2
    rename k_* *3
    rename m_* *4
    if "`outcome'" == "pa" {
        reshape long age_dv patot pabi watot motot vigtot carebi carehrs carecat, ///
            i(pidp) j(wave)
    }
    else if "`outcome'" == "alc" {
        reshape long age_dv sumaudit alcbi carebi carehrs, i(pidp) j(wave)
    }
    else if "`outcome'" == "diet" {
        reshape long age_dv diet carebi meandiet, i(pidp) j(wave)
        winsor2 meandiet, suffix(_trim) cut(0 99) trim label
    }
}

by pidp (wave), sort: gen countid = _n
by pidp (wave): gen count = _N

*--- Time scale: waves relative to caregiving onset ---
gen wavetoset = wave - ficare_all

if "`outcome'" == "smok" {
    gen poyrtoset = wavetoset + 7
    recode poyrtoset (-1=.)
    label var poyrtoset "Waves before/after caregiving onset (add 7)"
    mkspline part1 6 part2 7 part3 = poyrtoset, marginal
}
else {
    gen poyrtoset = wavetoset + 3
    label var poyrtoset "Waves before/after caregiving onset (add 3)"
    mkspline part1 2 part2 3 part3 = poyrtoset, marginal
}

xtset pidp wave

*--- Save long post-PSM dataset ---
save "$ukhls/`outcome'_long_afterpsm", replace


*===========================================================================
* SECTION 7: GROWTH CURVE MODELS
*
* For each analysis:
*   [MAIN]    Main model: caregiver vs non-caregiver trajectory
*   [HOURS]   By caregiving intensity (hours) at onset
*   [INTENSE] By caregiving intensity: matched sub-sample (intensewt)
*   [PLACE]   By place of caregiving (inside/outside/dual)
*   [INT]     Interaction test with kpscore or modifier
*   [SPLINE]  Piecewise spline interaction test
*===========================================================================

use "$ukhls/`outcome'_long_afterpsm", clear
xtset pidp wave

*--- Select model command ---
if "`outcome'" == "diet" {
    local cmd    "regress"
    local depvar "meandiet_trim"
    local wave_term ""             // diet models do not include i.wave
}
else if "`outcome'" == "smok" {
    local cmd    "logistic"
    local depvar "smok"
    local wave_term "i.wave"      // smoking includes wave adjustment
}
else if "`outcome'" == "pa" {
    local cmd    "logit"
    local depvar "pabi"
    local wave_term ""
}
else if "`outcome'" == "alc" {
    local cmd    "logistic"
    local depvar "alcbi"
    local wave_term ""
}


*--- [MAIN] Caregiver vs non-caregiver ---
`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_psm1, replace)

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_psm0, replace)

combomarginsplot `outcome'_psm0 `outcome'_psm1, ///
    labels("Non-caregiver" "Entered caregiving") ///
    savefile(`outcome'_intro, replace)

*--- [INT] Overall interaction test ---
`cmd' `depvar' i.poyrtoset##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm poyrtoset#kpscore

*--- [SPLINE] Piecewise interaction test ---
`cmd' `depvar' part1 c.part2##i.kpscore c.part3##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)


*--- [HOURS] By caregiving hours at onset ---
replace carehrs_fiob = 0 if kpscore == 0    // controls get 0 hours

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_h0, replace)

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 & carehrs_fiob==1 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_h1, replace)

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 & carehrs_fiob==2 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_h2, replace)

combomarginsplot `outcome'_h0 `outcome'_h1 `outcome'_h2, ///
    labels("Non-caregiver" "Low intensity (<20 hrs/week)" "High intensity (>=20 hrs/week)") ///
    savefile(`outcome'_intro_hours, replace)

`cmd' `depvar' i.poyrtoset##i.carehrs_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm poyrtoset#carehrs_fiob

`cmd' `depvar' part1 c.part2##i.carehrs_fiob c.part3##i.carehrs_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm carehrs_fiob#c.part2
testparm carehrs_fiob#c.part3


*--- [INTENSE] Intensity matched sub-sample ---
`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 & intense==0 [pw=intensewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_i0, replace)

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 & intense==1 [pw=intensewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_i1, replace)

combomarginsplot `outcome'_i0 `outcome'_i1, ///
    labels("Low (<20 hrs/week)" "High (>=20 hrs/week)") ///
    savefile(`outcome'_intro_intense, replace)

`cmd' `depvar' i.poyrtoset##i.intense `wave_term' [pw=intensewt], vce(cluster hhid_fiob)
testparm poyrtoset#intense

`cmd' `depvar' part1 c.part2##i.intense c.part3##i.intense `wave_term' [pw=intensewt], vce(cluster hhid_fiob)
testparm intense#c.part2
testparm intense#c.part3


*--- [PLACE] By place of caregiving ---
replace carecat_fiob = 0 if kpscore == 0    // controls get 0 (no caregiving)

`cmd' `depvar' i.poyrtoset `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoset, saving(`outcome'_cat0, replace)

forval cat = 1/3 {
    `cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 & carecat_fiob==`cat' [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoset, saving(`outcome'_cat`cat', replace)
}

combomarginsplot `outcome'_cat0 `outcome'_cat1 `outcome'_cat2 `outcome'_cat3, ///
    labels("Non-caregiver" ///
           "Outside household" ///
           "Inside household" ///
           "Dual (inside and outside)") ///
    savefile(`outcome'_intro_carecat, replace)

`cmd' `depvar' i.poyrtoset##i.carecat_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm poyrtoset#carecat_fiob

`cmd' `depvar' part1 c.part2##i.carecat_fiob c.part3##i.carecat_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm carecat_fiob#c.part2
testparm carecat_fiob#c.part3

di as result "===== 04_psm_entry.do completed for outcome: `outcome' ====="
