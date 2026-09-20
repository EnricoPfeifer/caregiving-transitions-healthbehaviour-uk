/*===========================================================================
  03_fixed_effects.do

  PURPOSE: Fixed effects models for caregiving transitions (entry and exit)
           and health behaviours. Covers:
             - FE model (caregiving status)
             - Caregiving hours
             - Sex and age-group interactions
             - Stratified models by age group

  PREREQUISITE: 02_clean.do must have been run for this outcome.

  OUTCOMES SUPPORTED:
    pa    -- physical inactivity (binary; xtlogit; waves g i k m)
    alc   -- problematic alcohol consumption (binary; xtlogit; waves g i k m)
    diet  -- daily fruit and vegetable portions (continuous; xtreg; waves g i k m)
    smok  -- smoking status (binary; xtlogit; waves e-m)

  WAVE NUMBERING (after reshape):
    pa / alc / diet : wave = 1 (w7) 2 (w9) 3 (w11) 4 (w13)
    smok            : wave = 1 (w5) 2 (w6) ... 9 (w13)

  OUTPUT: results stored in memory; export to log or use estimates store
          as required

  NOTE: Update file paths to your own file paths

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

cd "[insert_own_filepath]"
global ukhls "[insert_own_filepath]"


*===========================================================================
* SECTION 1: LOAD DATA AND PREPARE WIDE-FORMAT VARIABLES
*===========================================================================

*---------------------------------------------------------------------------
* Load cleaned wide dataset
*---------------------------------------------------------------------------
if "`outcome'" == "smok" {
    use "$ukhls/smok_wide_beforeps", clear
}
else {
    use "$ukhls/`outcome'_wide_cleaned", clear
}

keep if ecare != .

*---------------------------------------------------------------------------
* Drop if outcome is missing in ALL waves
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    egen out_miss = rmiss(*pabi)
    drop if out_miss == 4
}
else if "`outcome'" == "alc" {
    egen out_miss = rmiss(*alcbi)
    drop if out_miss == 4
}
else if "`outcome'" == "diet" {
    egen out_miss = rmiss(*meandiet)
    drop if out_miss == 4
}
else if "`outcome'" == "smok" {
    egen out_miss = rmiss(*smok)
    drop if out_miss == 9
}
drop out_miss

*---------------------------------------------------------------------------
* Baseline household ID
*---------------------------------------------------------------------------
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

*---------------------------------------------------------------------------
* Outcome at first observation
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    gen pabi_fiob   = g_pabi  if fiob == 1
    replace pabi_fiob   = i_pabi  if fiob == 2
    replace pabi_fiob   = k_pabi  if fiob == 3
    replace pabi_fiob   = m_pabi  if fiob == 4
    gen patot_fiob  = g_patot if fiob == 1
    replace patot_fiob  = i_patot if fiob == 2
    replace patot_fiob  = k_patot if fiob == 3
    replace patot_fiob  = m_patot if fiob == 4
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

*---------------------------------------------------------------------------
* Caregiving hours at onset (smoking file uses ficare, not fiob)
*---------------------------------------------------------------------------
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
    label val carehrs_fiob carehval
}

*---------------------------------------------------------------------------
* Ever smoked / ever not smoked (smoking only)
*---------------------------------------------------------------------------
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

*---------------------------------------------------------------------------
* Outcome and exposure change flags (4-wave outcomes only)
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    egen min_out = rowmin(g_pabi  i_pabi  k_pabi  m_pabi)
    egen max_out = rowmax(g_pabi  i_pabi  k_pabi  m_pabi)
}
else if "`outcome'" == "alc" {
    egen min_out = rowmin(g_alcbi    i_alcbi    k_alcbi    m_alcbi)
    egen max_out = rowmax(g_alcbi    i_alcbi    k_alcbi    m_alcbi)
}
else if "`outcome'" == "diet" {
    egen min_out = rowmin(g_meandiet i_meandiet k_meandiet m_meandiet)
    egen max_out = rowmax(g_meandiet i_meandiet k_meandiet m_meandiet)
}

if "`outcome'" != "smok" {
    gen outcomechange  = (min_out  != max_out)  if min_out  != .
    egen min_carebi    = rowmin(g_carebi i_carebi k_carebi m_carebi)
    egen max_carebi    = rowmax(g_carebi i_carebi k_carebi m_carebi)
    gen exposurechange = (min_carebi != max_carebi) if min_carebi != .
    drop min_out max_out min_carebi max_carebi
}

*---------------------------------------------------------------------------
* Intensity variable (5-category caregiving hours)
*---------------------------------------------------------------------------
label define intensityval ///
    1 "0-4 hours"    ///
    2 "5-9 hours"    ///
    3 "10-19 hours"  ///
    4 "20-34 hours"  ///
    5 "35+ hours"

if "`outcome'" == "smok" {
    foreach w in e f g h i j k l m {
        gen `w'_intensity = `w'_aidhrs
        recode `w'_intensity (-10/-1=.) (1=1) (2=2) (3=3) (4=4) (5/7=5) (8=3) (9=5) (97=.)
        label val `w'_intensity intensityval
        label var `w'_intensity "Caregiving hours (5 categories)"
    }
}
else {
    foreach w in g i k m {
        gen `w'_intensity = `w'_aidhrs
        recode `w'_intensity (-10/-1=.) (1=1) (2=2) (3=3) (4=4) (5/7=5) (8=3) (9=5) (97=.)
        label val `w'_intensity intensityval
        label var `w'_intensity "Caregiving hours (5 categories)"
    }
}

sort pidp


*===========================================================================
* SECTION 2: RESHAPE TO LONG FORMAT
* All outcomes harmonised to wave = 1/2/3/4 (4-wave) or 1-9 (smoking)
*===========================================================================

if "`outcome'" == "smok" {
    rename e_* *1
    rename f_* *2
    rename g_* *3
    rename h_* *4
    rename i_* *5
    rename j_* *6
    rename k_* *7
    rename l_* *8
    rename m_* *9

    reshape long smok ciggy age_dv carebi carehrs carecat intensity ///
        oclass cohab hhgroup workbi iwealth ghq ghealth children, ///
        i(pidp) j(wave)
    * wave already runs 1-9 after rename
}
else if "`outcome'" == "pa" {
    rename g_* *1
    rename i_* *2
    rename k_* *3
    rename m_* *4

    reshape long pabi patot watot age_dv carebi carehrs carecat intensity ///
        oclass cohab hhgroup workbi iwealth ghq ghealth sf12p children, ///
        i(pidp) j(wave)
    * wave runs 1-4
}
else if "`outcome'" == "alc" {
    rename g_* *1
    rename i_* *2
    rename k_* *3
    rename m_* *4

    reshape long alcbi sumaudit age_dv carebi carehrs carecat intensity ///
        oclass cohab hhgroup workbi iwealth ghq ghealth children, ///
        i(pidp) j(wave)
    * wave runs 1-4
}
else if "`outcome'" == "diet" {
    rename g_* *1
    rename i_* *2
    rename k_* *3
    rename m_* *4

    reshape long age_dv diet meandiet carebi carehrs carecat ///
        oclass cohab hhgroup workbi iwealth ghq ghealth children, ///
        i(pidp) j(wave)
    * wave runs 1-4

    * Winsorise continuous diet outcome at 99th percentile for modelling
    winsor2 meandiet, suffix(_trim) cut(0 99) trim label
}

by pidp (wave), sort: gen countid = _n
by pidp (wave): gen count = _N

gen carebirev = carebi
recode carebirev (0=1) (1=0)

xtset pidp wave


*===========================================================================
* SECTION 3: DEFINE COVARIATES AND TRANSITION VARIABLES
*===========================================================================

if "`outcome'" == "pa" {
    global covariates age_dv i.oclass i.cohab i.hhgroup i.workbi i.iwealth ///
        ghq i.ghealth sf12p children
}
else {
    global covariates age_dv i.oclass i.cohab i.hhgroup i.workbi i.iwealth ///
        ghq i.ghealth children
}

*--- Entry and exit transitions (derived from lagged carebi) ---
bysort pidp (wave): gen caregiving_intro = (carebi[_n-1] == 0 & carebi == 1)
bysort pidp (wave): gen caregiving_exit  = (carebi[_n-1] == 1 & carebi == 0)

*--- Lag variable ---
gen carebi_lag = L.carebi
bysort pidp (wave): replace carebi_lag = . if _n == 1
gen start_care = (carebi == 1 & carebi_lag == 0)
gen stop_care  = (carebi == 0 & carebi_lag == 1)


*===========================================================================
* SECTION 4: FIXED EFFECTS MODELS
*
* Notation:
*   [ENTRY]  = models for caregiving entry / intro
*   [EXIT]   = models for caregiving exit
*   [INT]    = interaction test
*   [STRAT]  = stratified by age group
*===========================================================================

*---------------------------------------------------------------------------
* Select model command based on outcome type
* diet uses xtreg (continuous outcome); all others use xtlogit
*---------------------------------------------------------------------------
if "`outcome'" == "diet" {
    local cmd   "xtreg"
    local depvar "meandiet_trim"
    local opts  "fe i(pidp)"
}
else {
    local cmd   "xtlogit"
    local depvar "`outcome'bi" // pabi / alcbi / smok (note: smok uses smok not smokbi)
    if "`outcome'" == "smok" local depvar "smok"
    local opts  "fe i(pidp) or"
}

*--- LR test: linear vs factor wave ---
`cmd' `depvar' carebi wave,       `opts'
est store wave_linear
`cmd' `depvar' carebi i.wave,     `opts'
est store wave_factor
lrtest wave_linear wave_factor
/* Expected result: factor wave preferred (p<0.001) */


*===========================================================================
* [ENTRY] CAREGIVING HOURS MODELS
*===========================================================================

*--- caregiving status ---
`cmd' `depvar' carebi i.wave, `opts'

*--- Caregiving hours ---
`cmd' `depvar' i.carehrs i.wave, `opts'
testparm i.carehrs

*--- [INT] Sex interaction: entry ---
`cmd' `depvar' i.carebi##i.sex i.wave, `opts'
testparm carebi#sex

*--- [INT] Age interaction: entry ---
`cmd' `depvar' i.carebi##i.agegr_fiob i.wave, `opts'
testparm carebi#agegr_fiob

*--- [INT] Sex interaction: hours ---
`cmd' `depvar' i.carehrs##i.sex i.wave, `opts'
testparm carehrs#sex

*--- [INT] Age interaction: hours ---
`cmd' `depvar' i.carehrs##i.agegr_fiob i.wave, `opts'
testparm carehrs#agegr_fiob

*--- [STRAT] Entry stratified by age group ---
forval ag = 0/3 {
    `cmd' `depvar' i.carebi i.wave if agegr_fiob == `ag', `opts'
}

*--- [STRAT] Hours stratified by age group ---
forval ag = 0/3 {
    `cmd' `depvar' i.carehrs i.wave if agegr_fiob == `ag', `opts'
    testparm i.carehrs
}


*===========================================================================
* [EXIT] CAREGIVING EXIT MODELS
*===========================================================================

*--- [EXIT] Main exit model ---
`cmd' `depvar' caregiving_exit i.wave, `opts'
`cmd' `depvar' caregiving_exit i.wave $covariates, `opts'

*--- [EXIT] Restricted to baseline caregivers ---
`cmd' `depvar' caregiving_exit i.wave if care_fiob == 1, `opts'
`cmd' `depvar' caregiving_exit i.wave $covariates if care_fiob == 1, `opts'

*--- [INT] Sex interaction: exit ---
`cmd' `depvar' i.caregiving_exit##i.sex i.wave, `opts'
testparm caregiving_exit#sex

*--- [INT] Age interaction: exit ---
`cmd' `depvar' i.caregiving_exit##i.agegr_fiob i.wave, `opts'
testparm caregiving_exit#agegr_fiob

*--- [STRAT] Exit stratified by age group ---
forval ag = 0/3 {
    `cmd' `depvar' i.caregiving_exit i.wave if agegr_fiob == `ag', `opts'
}
