/*===========================================================================
  06_psm_exit.do

  PURPOSE: Propensity score matching and piecewise growth curve models for
           caregiving EXIT. Two control group options:
             nocare    -- non-caregivers (never caregivers)
             alwayscare -- continuing caregivers (caregivers at baseline
                           who never exited)

  Key differences from 04_psm_entry.do:
    - Outcome is exitcarebi (exited caregiving) not ecare
    - Time scale is wavetoexit / poyrtoexit centred on exit wave
    - careprior and catprior: care hours and place in the wave prior to exit
      (derived via a forward-fill reshape before the main reshape)
    - nocare: pool of controls is never-caregivers; alwayscare: pool is
      baseline caregivers who never exited

  PREREQUISITE: 02_clean.do must have been run for this outcome.

  OUTCOMES SUPPORTED:
    pa    -- physical inactivity   (logit;    waves g i k m)
    alc   -- problematic drinking  (logistic; waves g i k m)
    diet  -- fruit/veg portions    (regress;  waves g i k m)
    smok  -- smoking               (logistic; waves e-m; includes i.wave)

  NOTE: Paths below use the author's local machine. Replace with your own
        paths before running. These will be removed prior to publication.

  AUTHOR: Enrico Pfeifer
===========================================================================*/

*---------------------------------------------------------------------------
* SET OUTCOME AND CONTROL GROUP HERE
* outcome options : pa | alc | diet | smok
* control options : nocare | alwayscare
*---------------------------------------------------------------------------
local outcome "pa"
local control "nocare"


*---------------------------------------------------------------------------
* PATHS  --  replace with your own directory structure
*---------------------------------------------------------------------------
clear all
set more off
set maxvar 10000
set cformat %4.2f
set seed 12345

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

*--- Control group restriction ---
if "`control'" == "nocare" {
    * Exclude always-caregivers; controls are never-caregivers
    egen occ_noncarer = anycount(*_carebi), values(0)
    drop if occ_noncarer == 0
}
else {
    * Restrict to baseline caregivers; controls are those who never exited
    keep if ecare == 1
    keep if care_fiob == 1
}

gen wave_fiob = fiob


*---------------------------------------------------------------------------
* Last wave of observation
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


*---------------------------------------------------------------------------
* Exit caregiving wave and indicator
*---------------------------------------------------------------------------
gen exitcarewave = .

if "`outcome'" == "smok" {
    replace exitcarewave = 2 if f_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 3 if g_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 4 if h_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 5 if i_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 6 if j_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 7 if k_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 8 if l_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 9 if m_carebi == 0 & care_fiob == 1 & exitcarewave == .
}
else {
    replace exitcarewave = 1 if g_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 2 if i_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 3 if k_carebi == 0 & care_fiob == 1 & exitcarewave == .
    replace exitcarewave = 4 if m_carebi == 0 & care_fiob == 1 & exitcarewave == .
}

gen exitcarebi = exitcarewave
recode exitcarebi 1/max=1

if "`control'" == "nocare" {
    replace exitcarebi = 0 if ecare == 0 & exitcarewave == .
    label define exitval 0 "Never caregiver" 1 "Exited caregiving"
}
else {
    replace exitcarebi = 0 if ecare == 1 & exitcarewave == .
    label define exitval 0 "Remained caregiver" 1 "Exited caregiving"
}
label val exitcarebi exitval
label var exitcarewave "Wave in which caregiving was exited"
label var exitcarebi   "Ever exited caregiving"

*--- samplepiece: at least one wave before AND one wave after exit ---
gen samplepiece = 1 if (wave_fiob < exitcarewave) & (lastwave >= exitcarewave) & exitcarewave != .
drop if samplepiece != 1 & exitcarebi == 1
drop if nwaves == 1 & exitcarebi == 0

*--- Age at exit (smoking only, used in nocare matching) ---
if "`outcome'" == "smok" & "`control'" == "nocare" {
    gen agexit = .
    forval w = 1/9 {
        local wl : word `w' of e f g h i j k l m
        replace agexit = `wl'_age_dv if exitcarewave == `w' & agexit == .
    }
    replace agexit = 0 if exitcarebi == 0
}


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
    recode smok_fiob (min/-1=.)
}

sort pidp


*===========================================================================
* SECTION 3: PROPENSITY SCORE MATCHING
* Treatment: exitcarebi (1=exited, 0=control group)
*===========================================================================

if "`outcome'" == "pa"   local outvar "(pabi_fiob)"
if "`outcome'" == "alc"  local outvar "(alcbi_fiob)"
if "`outcome'" == "diet" local outvar "(meandiet_fiob)"
if "`outcome'" == "smok" local outvar "(smok_fiob)"

if "`outcome'" == "pa" {
    kmatch ps exitcarebi i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob ///
        sf_fiob i.children_fiob `outvar', ///
        ematch(sex wave_fiob agefiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else if "`outcome'" == "smok" {
    kmatch ps exitcarebi i.oclass_fiob i.cohab_fiob i.hhgroup_fiob nwaves ///
        i.edu i.workbi_fiob i.iwealth_fiob i.ethnicity ///
        ghq_fiob i.ghealth_fiob i.children_fiob `outvar', ///
        ematch(sex agefiob wave_fiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else {
    * alc and diet
    kmatch ps exitcarebi i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves ///
        i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob ///
        i.children_fiob `outvar', ///
        ematch(sex wave_fiob agefiob) nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}

tab kpscore
gen matchwt = _KM_mw
gen pscore  = _KM_ps


*===========================================================================
* SECTION 4: ASSIGN EXIT WAVE TO MATCHED CONTROLS
*===========================================================================

forval k = 1/3 {
    gen double id`k' = controlid`k'
    replace id`k' = pidp if kpscore == 0
    gen exitcare_tmp`k' = exitcarewave if kpscore == 1
    sort id`k' exitcare_tmp`k'
    by id`k': replace exitcare_tmp`k' = exitcare_tmp`k'[1]
    gen exitcare_g`k' = exitcare_tmp`k' if kpscore == 0
    drop exitcare_tmp`k'
}

gen exitcare_all = exitcarewave

if "`control'" == "nocare" {
    replace exitcare_all = exitcare_g1 if exitcare_g1 != . & exitcare_all == 0
}
else {
    replace exitcare_all = exitcare_g1 if exitcare_g1 != . & exitcare_all == .
}

gen preafter0 = 1 if (wave_fiob < exitcare_all) & (lastwave >= exitcare_all) & exitcare_all != .
replace exitcare_all = exitcare_g2 if preafter0 != 1
gen preafter  = 1 if (wave_fiob < exitcare_all) & (lastwave >= exitcare_all) & exitcare_all != .
replace exitcare_all = exitcare_g3 if preafter != 1

gen usesample = 1 if (wave_fiob < exitcare_all) & (lastwave >= exitcare_all) & exitcare_all != .
keep if usesample == 1


*===========================================================================
* SECTION 5: CARE HOURS AND PLACE PRIOR TO EXIT
* Forward-fill reshape to get the last observed care hours/place before exit
*===========================================================================

save "$ukhls/`outcome'_wide_afterps_prior_`control'", replace

sort pidp

*--- Rename to numeric suffixes for reshape ---
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
    local carewaves "5 6 7 8 9 10 11 12 13"
    local priorwaves "5 6 7 8 9 10 11 12"
    keep carehrs5-carehrs13 carecat5-carecat13 pidp
}
else {
    rename g_* *7
    rename i_* *9
    rename k_* *11
    rename m_* *13
    local carewaves "7 9 11 13"
    local priorwaves "7 9 11"
    keep carehrs7-carehrs13 carecat7-carecat13 pidp
}

reshape long carehrs carecat, i(pidp) j(wave)
sort pidp wave

*--- Forward-fill to carry last observed care hours/place forward ---
gen carehrs_dummy = carehrs
gen carecat_dummy = carecat
by pidp: replace carehrs_dummy = carehrs_dummy[_n-1] if missing(carehrs_dummy)
by pidp: replace carecat_dummy = carecat_dummy[_n-1] if missing(carecat_dummy)
drop carehrs carecat
reshape wide carehrs_dummy carecat_dummy, i(pidp) j(wave)

merge 1:1 pidp using "$ukhls/`outcome'_wide_afterps_prior_`control'"

*--- Extract care hours/place at the wave immediately prior to exit ---
gen prior = exitcarewave - 1
gen careprior = .
gen catprior  = .

foreach w of local priorwaves {
    replace careprior = carehrs_dummy`w' if prior == `w' & careprior == . & carehrs_dummy`w' != .
    replace catprior  = carecat_dummy`w' if prior == `w' & catprior  == . & carecat_dummy`w' != .
}
recode catprior (0=.)


*===========================================================================
* SECTION 6: COVARIATE BALANCE (ebalance)
*===========================================================================

if "`outcome'" == "pa" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity sf_fiob i.children_fiob (pabi_fiob), gen(balancewt)
}
else if "`outcome'" == "alc" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob workbi_fiob ///
        i.ethnicity i.children_fiob (alcbi_fiob), gen(balancewt)
}
else if "`outcome'" == "diet" {
    ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
        i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob i.workcom_fiob ///
        i.ethnicity i.children_fiob (meandiet_fiob), gen(balancewt)
}
else if "`outcome'" == "smok" {
    if "`control'" == "nocare" {
        ebalance exitcarebi agefiob i.sex i.cohab_fiob ///
            hhincome_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
            ghq_fiob i.hhgroup_fiob i.ethnicity ///
            nwaves i.workbi_fiob i.children_fiob (smok_fiob), gen(balancewt)
    }
    else {
        ebalance exitcarebi agefiob i.sex i.cohab_fiob ///
            i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
            ghq_fiob i.hhgroup_fiob i.ethnicity nwaves workbi_fiob ///
            i.children_fiob, gen(balancewt)
    }
}

sort pidp


*===========================================================================
* SECTION 7: RESHAPE TO LONG FORMAT
*===========================================================================

if "`outcome'" == "smok" {
    * already renamed above; re-rename from wide_afterps_prior merge result
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
    gen wave = waves - 4
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
        reshape long age_dv alcbi sumaudit carebi carehrs carecat, i(pidp) j(wave)
    }
    else if "`outcome'" == "diet" {
        reshape long age_dv diet carebi meandiet, i(pidp) j(wave)
        winsor2 meandiet, suffix(_trim) cut(0 99) trim label
    }
}

by pidp (wave), sort: gen countid = _n
by pidp (wave): gen count = _N

gen wavetoexit = wave - exitcare_all

if "`outcome'" == "smok" {
    gen poyrtoexit = wavetoexit + 7
    recode poyrtoexit (min/-1=.) (13=.)
    label var poyrtoexit "Waves before/after caregiving exit (add 7)"
    mkspline part1 6 part2 7 part3 = poyrtoexit, marginal
}
else {
    gen poyrtoexit = wavetoexit + 3
    label var poyrtoexit "Waves before/after caregiving exit (add 3)"
    mkspline part1 2 part2 3 part3 = poyrtoexit, marginal
}

xtset pidp wave

save "$ukhls/`outcome'_exit_`control'_long_afterPSM", replace


*===========================================================================
* SECTION 8: OUTCOME-SPECIFIC GRAPH SETTINGS
*===========================================================================

if "`outcome'" == "pa" {
    local ylab   "0.45(0.05)0.65"
    local yscale "range(0.45 0.65)"
    local ytitle "Probability of physical inactivity"
}
else if "`outcome'" == "alc" {
    local ylab   "0.2(0.1)0.7"
    local yscale "range(0.2 0.7)"
    local ytitle "Probability of problematic drinking"
}
else if "`outcome'" == "diet" {
    local ylab   "2(0.5)5"
    local yscale "range(2 5)"
    local ytitle "Portions fruit and vegetables per day"
}
else if "`outcome'" == "smok" {
    local ylab   "0.0(0.1)0.4"
    local yscale "range(0.0 0.4)"
    local ytitle "Probability of smoking"
}

if "`outcome'" == "smok" {
    local xline1 6
    local xline2 7
    local xlab `"0 "-7" 1 "-6" 2 "-5" 3 "-4" 4 "-3" 5 "-2" 6 "-1" 7 "0" 8 "1" 9 "2" 10 "3" 11 "4" 12 "5" 13 "6" 14 "7""'
}
else {
    local xline1 2
    local xline2 3
    local xlab `"0 "-6" 1 "-4" 2 "-2" 3 "0" 4 "2" 5 "4""'
}

*--- Control group label for graphs ---
if "`control'" == "nocare" {
    local ctrl_label "Non-caregiver"
}
else {
    local ctrl_label "Continued caregiving"
}


*===========================================================================
* SECTION 9: GROWTH CURVE MODELS
*===========================================================================

use "$ukhls/`outcome'_exit_`control'_long_afterPSM", clear
xtset pidp wave

*--- Model command ---
if "`outcome'" == "diet" {
    local cmd    "regress"
    local depvar "meandiet_trim"
    local wave_term ""
}
else if "`outcome'" == "smok" {
    local cmd    "logistic"
    local depvar "smok"
    local wave_term "i.wave"
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

local tag "`outcome'_exit_`control'"


*--- [MAIN] Exited vs control group ---
`cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==1 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoexit, saving(`tag'_1, replace)

`cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoexit, saving(`tag'_0, replace)

combomarginsplot `tag'_0 `tag'_1, ///
    labels("`ctrl_label'" "Exited caregiving") ///
    savefile(`tag', replace)

mplotoffset using `tag', ///
    xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
    xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
    title("`ytitle'" ///
          "– exited caregiving vs. `ctrl_label'", size(large)) ///
    xtitle("Years centred on caregiving exit year", size(large)) ///
    ytitle("Predictive AME", size(medium)) ///
    xlabel(`xlab', labsize(large)) ///
    ylabel(`ylab', labsize(large)) ///
    yscale(`yscale') ///
    legend(size(medium))

*--- [SPLINE] Overall interaction test ---
`cmd' `depvar' i.poyrtoexit##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm poyrtoexit#kpscore

`cmd' `depvar' part1 c.part2##i.kpscore c.part3##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm kpscore#c.part2
testparm kpscore#c.part3


*--- [HOURS] By caregiving intensity prior to exit ---
replace careprior = 0 if kpscore == 0

`cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==1 & careprior==1 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoexit, saving(`tag'_h1, replace)

`cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==1 & careprior==2 [pw=balancewt], vce(cluster hhid_fiob)
margins i.poyrtoexit, saving(`tag'_h2, replace)

combomarginsplot `tag'_0 `tag'_h1 `tag'_h2, ///
    labels("`ctrl_label'" ///
           "Exited caregiving – low (<20 hrs/week)" ///
           "Exited caregiving – high (>=20 hrs/week)") ///
    savefile(`tag'_hours, replace)

mplotoffset using `tag'_hours, ///
    xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
    xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
    title("`ytitle'" ///
          "– by caregiving intensity prior to exit", size(large)) ///
    xtitle("Years centred on caregiving exit year", size(large)) ///
    ytitle("Predictive AME", size(medium)) ///
    xlabel(`xlab', labsize(large)) ///
    ylabel(`ylab', labsize(large)) ///
    yscale(`yscale') ///
    legend(size(medium))

`cmd' `depvar' part1 c.part2##i.careprior c.part3##i.careprior `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm careprior#c.part2
testparm careprior#c.part3


*--- [PLACE] By place of caregiving prior to exit ---
replace catprior = 0 if kpscore == 0

forval cat = 1/3 {
    `cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==1 & catprior==`cat' [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoexit, saving(`tag'_c`cat', replace)
}

combomarginsplot `tag'_0 `tag'_c1 `tag'_c2 `tag'_c3, ///
    labels("`ctrl_label'" ///
           "Exited caregiving – outside household" ///
           "Exited caregiving – inside household" ///
           "Exited caregiving – dual location") ///
    savefile(`tag'_place, replace)

mplotoffset using `tag'_place, ///
    xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
    xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
    title("`ytitle'" ///
          "– by place of caregiving prior to exit", size(large)) ///
    xtitle("Years centred on caregiving exit year", size(large)) ///
    ytitle("Predictive AME", size(medium)) ///
    xlabel(`xlab', labsize(large)) ///
    ylabel(`ylab', labsize(large)) ///
    yscale(`yscale') ///
    legend(size(small))

`cmd' `depvar' part1 c.part2##i.catprior c.part3##i.catprior `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
testparm catprior#c.part2
testparm catprior#c.part3


*===========================================================================
* COMBINED THREE-GROUP GRAPH
* Requires BOTH control groups to have been run first.
* Comment out until both nocare and alwayscare runs are complete.
*===========================================================================

/*

combomarginsplot ///
    `outcome'_exit_alwayscare_0 ///
    `outcome'_exit_nocare_0     ///
    `outcome'_exit_alwayscare_1, ///
    labels("Continued caregiving" ///
           "Non-caregiver" ///
           "Exited caregiving (vs. continued)") ///
    savefile(`outcome'_exit_3group, replace)

mplotoffset using `outcome'_exit_3group, ///
    xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
    xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
    title("`ytitle' – caregiving status comparison", size(large)) ///
    xtitle("Years centred on caregiving exit year", size(large)) ///
    ytitle("Predictive AME", size(medium)) ///
    xlabel(`xlab', labsize(large)) ///
    ylabel(`ylab', labsize(large)) ///
    yscale(`yscale') ///
    legend(size(medium))

combomarginsplot ///
    `outcome'_exit_alwayscare_0 ///
    `outcome'_exit_alwayscare_1 ///
    `outcome'_exit_nocare_0     ///
    `outcome'_exit_nocare_1, ///
    labels("Continued caregiving" ///
           "Exited caregiving (vs. continued)" ///
           "Non-caregiver" ///
           "Exited caregiving (vs. non-caregiver)") ///
    savefile(`outcome'_exit_4group, replace)

mplotoffset using `outcome'_exit_4group, ///
    xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
    xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
    title("`ytitle' – all comparisons around exit", size(large)) ///
    xtitle("Years centred on caregiving exit year", size(large)) ///
    ytitle("Predictive AME", size(medium)) ///
    xlabel(`xlab', labsize(large)) ///
    ylabel(`ylab', labsize(large)) ///
    yscale(`yscale') ///
    legend(size(medium))

*/

di as result "===== 06_psm_exit.do completed: outcome=`outcome' control=`control' ====="
