/*===========================================================================
  07_psm_exit_strata.do

  PURPOSE: Sex and age-group stratified PSM growth curve models for
           caregiving EXIT. Mirrors 05_psm_entry_strata.do but for exit,
           with the added `control` dimension (nocare or alwayscare).

  MODES:
    sex_int   -- sex interaction test (full sample; sex removed from
                 matching and balance)
    sex_strat -- stratified by sex (restrict before matching)
    age_int   -- age interaction test (full sample; age removed from
                 matching and balance)
    age_strat -- stratified by age group (restrict before matching)

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
* SET OUTCOME, CONTROL GROUP, MODE AND STRATUM HERE
*
* outcome options : pa | alc | diet | smok
* control options : nocare | alwayscare
* mode options    : sex_int | sex_strat | age_int | age_strat
* strat options   :
*   sex_strat -> 0 (men) | 1 (women)
*   age_strat -> 0 (16-29) | 1 (30-49) | 2 (50-64) | 3 (65+)
*   (ignored for sex_int and age_int)
*---------------------------------------------------------------------------
local outcome "diet"
local control "nocare"
local mode    "sex_strat"
local strat   1           // 1 = women


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
* SECTION 1: LOAD DATA AND RESTRICT SAMPLE
*===========================================================================

if "`outcome'" == "smok" {
    use "$ukhls/smok_wide_beforeps", clear
}
else {
    use "$ukhls/`outcome'_wide_cleaned", clear
}

keep if ecare != .

*--- Control group restriction (applied first) ---
if "`control'" == "nocare" {
    egen occ_noncarer = anycount(*_carebi), values(0)
    drop if occ_noncarer == 0
}
else {
    keep if ecare == 1
    keep if care_fiob == 1
}

*--- Stratum restriction (applied after control group restriction) ---
if "`mode'" == "sex_strat" {
    keep if sex == `strat'
}
else if "`mode'" == "age_strat" {
    keep if agegr_fiob == `strat'
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

gen samplepiece = 1 if (wave_fiob < exitcarewave) & (lastwave >= exitcarewave) & exitcarewave != .
drop if samplepiece != 1 & exitcarebi == 1
drop if nwaves == 1 & exitcarebi == 0


*===========================================================================
* SECTION 2: BASELINE VARIABLES
*===========================================================================

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

if "`outcome'" == "pa" {
    gen pabi_fiob = g_pabi if fiob == 1
    replace pabi_fiob = i_pabi if fiob == 2
    replace pabi_fiob = k_pabi if fiob == 3
    replace pabi_fiob = m_pabi if fiob == 4
}
else if "`outcome'" == "alc" {
    gen alcbi_fiob = g_alcbi if fiob == 1
    replace alcbi_fiob = i_alcbi if fiob == 2
    replace alcbi_fiob = k_alcbi if fiob == 3
    replace alcbi_fiob = m_alcbi if fiob == 4
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
*
* Sex modes:  sex removed from ematch; age kept
* Age modes:  agefiob removed from ematch; sex kept (via ematch(wave_fiob))
*===========================================================================

if "`outcome'" == "pa"   local outvar "(pabi_fiob)"
if "`outcome'" == "alc"  local outvar "(alcbi_fiob)"
if "`outcome'" == "diet" local outvar "(meandiet_fiob)"
if "`outcome'" == "smok" local outvar "(smok_fiob)"

local core_cov "i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob i.children_fiob"

if inlist("`mode'", "sex_int", "sex_strat") {
    local ematch_vars "wave_fiob agefiob"     // sex excluded
}
else {
    local ematch_vars "sex wave_fiob"         // agefiob excluded
}

if "`outcome'" == "smok" {
    local ematch_vars = subinstr("`ematch_vars'", "agefiob", "agefiob", .)
    // smoking uses agefiob in ematch for sex modes, same variable name
}

kmatch ps exitcarebi `core_cov' `outvar', ///
    ematch(`ematch_vars') nn(3) ///
    idgenerate(controlid) idvar(pidp) generate(kpscore) wor

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
*===========================================================================

save "$ukhls/`outcome'_wide_afterps_prior_strata", replace
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
    local priorwaves "5 6 7 8 9 10 11 12"
    keep carehrs5-carehrs13 carecat5-carecat13 pidp
}
else {
    rename g_* *7
    rename i_* *9
    rename k_* *11
    rename m_* *13
    local priorwaves "7 9 11"
    keep carehrs7-carehrs13 carecat7-carecat13 pidp
}

reshape long carehrs carecat, i(pidp) j(wave)
sort pidp wave
gen carehrs_dummy = carehrs
gen carecat_dummy = carecat
by pidp: replace carehrs_dummy = carehrs_dummy[_n-1] if missing(carehrs_dummy)
by pidp: replace carecat_dummy = carecat_dummy[_n-1] if missing(carecat_dummy)
drop carehrs carecat
reshape wide carehrs_dummy carecat_dummy, i(pidp) j(wave)
merge 1:1 pidp using "$ukhls/`outcome'_wide_afterps_prior_strata"

gen prior    = exitcarewave - 1
gen careprior = .
gen catprior  = .
foreach w of local priorwaves {
    replace careprior = carehrs_dummy`w' if prior == `w' & careprior == . & carehrs_dummy`w' != .
    replace catprior  = carecat_dummy`w' if prior == `w' & catprior  == . & carecat_dummy`w' != .
}
recode catprior (0=.)


*===========================================================================
* SECTION 6: COVARIATE BALANCE
*
* Sex modes:  age kept; sex excluded from ebalance
* Age modes:  sex kept; agefiob excluded from ebalance
*===========================================================================

if "`outcome'" == "pa"   local eb_out "(pabi_fiob)"
if "`outcome'" == "alc"  local eb_out "(alcbi_fiob)"
if "`outcome'" == "diet" local eb_out "(meandiet_fiob)"
if "`outcome'" == "smok" local eb_out "(smok_fiob)"

local eb_core "i.cohab_fiob nwaves i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ghq_fiob i.hhgroup_fiob workbi_fiob i.ethnicity i.children_fiob"

if inlist("`mode'", "sex_int", "sex_strat") {
    * Age kept; sex excluded
    ebalance kpscore agefiob `eb_core' `eb_out', gen(balancewt)
}
else {
    * Sex kept; age excluded
    ebalance kpscore i.sex `eb_core' `eb_out', gen(balancewt)
}

sort pidp


*===========================================================================
* SECTION 7: RESHAPE TO LONG FORMAT
*===========================================================================

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
recode wavetoexit (-4=.) if "`outcome'" != "smok"

if "`outcome'" == "smok" {
    gen poyrtoexit = wavetoexit + 7
    recode poyrtoexit (min/-1=.) (13=.)
    mkspline part1 6 part2 7 part3 = poyrtoexit, marginal
}
else {
    gen poyrtoexit = wavetoexit + 3
    mkspline part1 2 part2 3 part3 = poyrtoexit, marginal
}

*--- Save long dataset ---
if inlist("`mode'", "sex_int", "age_int") {
    save "$ukhls/`outcome'_exit_`control'_`mode'_long", replace
}
else {
    save "$ukhls/`outcome'_exit_`control'_`mode'`strat'_long", replace
}


*===========================================================================
* SECTION 8: GRAPH SETTINGS
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
    local ytitle "Daily portions of fruit and vegetables"
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

if "`control'" == "nocare" {
    local ctrl_label "Non-caregiver"
}
else {
    local ctrl_label "Continued caregiving"
}

*--- Stratum labels ---
if "`mode'" == "sex_strat" {
    if `strat' == 0 {
        local stratlabel "men"
        local strattitle "Male"
        local leg1 "Male: `ctrl_label'"
        local leg2 "Male: exited caregiving"
    }
    else {
        local stratlabel "women"
        local strattitle "Female"
        local leg1 "Female: `ctrl_label'"
        local leg2 "Female: exited caregiving"
    }
}
else if "`mode'" == "age_strat" {
    if `strat' == 0 {
        local stratlabel "age0"
        local strattitle "Early adulthood (16–29)"
        local leg1 "16–29: `ctrl_label'"
        local leg2 "16–29: exited caregiving"
    }
    else if `strat' == 1 {
        local stratlabel "age1"
        local strattitle "Early mid-adulthood (30–49)"
        local leg1 "30–49: `ctrl_label'"
        local leg2 "30–49: exited caregiving"
    }
    else if `strat' == 2 {
        local stratlabel "age2"
        local strattitle "Late mid-adulthood (50–64)"
        local leg1 "50–64: `ctrl_label'"
        local leg2 "50–64: exited caregiving"
    }
    else if `strat' == 3 {
        local stratlabel "age3"
        local strattitle "Late adulthood (65+)"
        local leg1 "65+: `ctrl_label'"
        local leg2 "65+: exited caregiving"
    }
}

local tag "`outcome'_exit_`control'_`mode'"
if inlist("`mode'", "sex_strat", "age_strat") local tag "`tag'`strat'"


*===========================================================================
* SECTION 9: MODELS AND GRAPHS
*===========================================================================

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


*---------------------------------------------------------------------------
* INTERACTION TEST (sex_int and age_int modes)
*---------------------------------------------------------------------------
if "`mode'" == "sex_int" {
    `cmd' `depvar' part1 c.part2##i.kpscore##i.sex c.part3##i.kpscore##i.sex `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#sex#c.part2
    testparm kpscore#sex#c.part3
}
else if "`mode'" == "age_int" {
    `cmd' `depvar' part1 c.part2##i.kpscore##i.agegr_fiob c.part3##i.kpscore##i.agegr_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#agegr_fiob#c.part2
    testparm kpscore#agegr_fiob#c.part3
}


*---------------------------------------------------------------------------
* STRATIFIED MODELS AND GRAPHS (sex_strat and age_strat modes)
*---------------------------------------------------------------------------
if inlist("`mode'", "sex_strat", "age_strat") {

    *--- Main trajectory ---
    `cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==1 [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoexit, saving(`tag'_1, replace)

    `cmd' `depvar' i.poyrtoexit `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoexit, saving(`tag'_0, replace)

    combomarginsplot `tag'_0 `tag'_1, ///
        labels("`leg1'" "`leg2'") ///
        savefile(`tag'_intro, replace)

    mplotoffset using `tag'_intro, ///
        xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
        xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
        title("`ytitle' – `strattitle'", size(large)) ///
        xtitle("Waves centred on caregiving exit wave", size(large)) ///
        ytitle("Predictive AME", size(medium)) ///
        xlabel(`xlab', labsize(large)) ///
        ylabel(`ylab', labsize(large)) ///
        yscale(`yscale') ///
        legend(size(medium))

    *--- Spline interaction test within stratum ---
    `cmd' `depvar' part1 c.part2##i.kpscore c.part3##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#c.part2
    testparm kpscore#c.part3

    *--- Publication panel (no labels) ---
    mplotoffset using `tag'_intro, ///
        xline(`xline1', lpattern(shortdash) lcolor(black)) ///
        xline(`xline2', lpattern(shortdash) lcolor(black)) ///
        xtitle("") ytitle("") ///
        xlabel(, nolabel) ylabel(, nolabel) ///
        yscale(`yscale') ///
        title("", size(vlarge)) ///
        legend(off) ///
        graphregion(color(white))
    graph save `outcome'_`control'_`stratlabel'_exit_panel.gph, replace

    *--- Publication panel with legend ---
    mplotoffset using `tag'_intro, ///
        xline(`xline1', lpattern(shortdash) lcolor(black)) ///
        xline(`xline2', lpattern(shortdash) lcolor(black)) ///
        title("`strattitle'", size(vhuge)) ///
        xtitle("") ytitle("") ///
        xlabel(, nolabel) ylabel(, nolabel) ///
        yscale(`yscale') ///
        legend(position(11) ring(1) rows(2) size(huge) ///
               label(1 "`leg1'" 2 "`leg2'") ///
               region(color(white))) ///
        graphregion(color(white))
    graph save `outcome'_`control'_`stratlabel'_exit_panel_legend.gph, replace
}

di as result "===== 07_psm_exit_strata.do completed: outcome=`outcome' control=`control' mode=`mode' strat=`strat' ====="


*===========================================================================
* COMBINED THREE-GROUP GRAPHS
* Requires BOTH control group strata runs to be complete.
* Run only after nocare AND alwayscare strata have been completed.
* The combined graph superimposes: always-care control, nocare control,
* and exiters (from the alwayscare run, which is the primary comparison).
*===========================================================================

/*

*--- Combined sex strata ---
if "`mode'" == "sex_strat" {
    combomarginsplot ///
        `outcome'_exit_alwayscare_`mode'`strat'_0 ///
        `outcome'_exit_nocare_`mode'`strat'_0     ///
        `outcome'_exit_alwayscare_`mode'`strat'_1, ///
        labels("`strattitle': `ctrl_label'" ///
               "`strattitle': non-caregiver" ///
               "`strattitle': exited caregiving") ///
        savefile(`outcome'_exit_3grp_`stratlabel', replace)

    mplotoffset using `outcome'_exit_3grp_`stratlabel', ///
        xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
        xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
        title("`ytitle' – `strattitle'", size(large)) ///
        xtitle("Waves centred on caregiving exit wave", size(large)) ///
        ytitle("Predictive AME", size(medium)) ///
        xlabel(`xlab', labsize(large)) ///
        ylabel(`ylab', labsize(large)) ///
        yscale(`yscale') ///
        legend(size(medium))
}

*--- Combined publication panels: all six strata without legend ---
graph combine ///
    `outcome'_nocare_men_exit_panel.gph    `outcome'_nocare_women_exit_panel.gph ///
    `outcome'_nocare_age0_exit_panel.gph   `outcome'_nocare_age1_exit_panel.gph  ///
    `outcome'_nocare_age2_exit_panel.gph   `outcome'_nocare_age3_exit_panel.gph, ///
    cols(6) xcommon ycommon ///
    imargin(0 0 0 0) ///
    graphregion(color(white)) ///
    title("", size(large)) ///
    xsize(15) ysize(3)

graph export "`outcome'_exit_nocare_strata.jpg", as(jpg) name("Graph") quality(90) replace

graph combine ///
    `outcome'_alwayscare_men_exit_panel.gph    `outcome'_alwayscare_women_exit_panel.gph ///
    `outcome'_alwayscare_age0_exit_panel.gph   `outcome'_alwayscare_age1_exit_panel.gph  ///
    `outcome'_alwayscare_age2_exit_panel.gph   `outcome'_alwayscare_age3_exit_panel.gph, ///
    cols(6) xcommon ycommon ///
    imargin(0 0 0 0) ///
    graphregion(color(white)) ///
    title("", size(large)) ///
    xsize(15) ysize(3)

graph export "`outcome'_exit_alwayscare_strata.jpg", as(jpg) name("Graph") quality(90) replace

*/
