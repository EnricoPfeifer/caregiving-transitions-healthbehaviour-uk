/*===========================================================================
  05_psm_entry_strata.do

  PURPOSE: Sex and age-group stratified PSM growth curve models for
           caregiving ENTRY. Runs in three modes:
             sex_int   -- sex interaction test (full sample; sex removed
                          from matching and balance)
             sex_strat -- separate models for men and women
             age_int   -- age-group interaction test (full sample; age
                          removed from matching and balance)
             age_strat -- separate model for one age group at a time

  The key analytical logic follows 04_psm_entry.do but with different
  sample restrictions and matching specifications.

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
* SET OUTCOME AND MODE HERE
*
* outcome options : pa | alc | diet | smok
* mode options    : sex_int | sex_strat | age_int | age_strat
* stratum options :
*   sex_strat  -> strat = 0 (men) | 1 (women)
*   age_strat  -> strat = 0 (16-29) | 1 (30-49) | 2 (50-64) | 3 (65+)
*   (ignored for sex_int and age_int)
*---------------------------------------------------------------------------
local outcome "pa"
local mode    "sex_strat"
local strat   1            // women; set to 0 for men


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
* SECTION 1: LOAD AND RESTRICT SAMPLE
*===========================================================================

if "`outcome'" == "smok" {
    use "$ukhls/smok_wide_beforeps", clear
}
else {
    use "$ukhls/`outcome'_wide_cleaned", clear
}

keep if ecare != .
drop if care_fiob == 1
gen wave_fiob = fiob

*--- Apply stratum restriction BEFORE matching ---
if "`mode'" == "sex_strat" {
    keep if sex == `strat'
}
else if "`mode'" == "age_strat" {
    keep if agegr_fiob == `strat'
}

*--- Sample eligibility ---
gen samplepiece = 1 if wave_fiob < ficare & ficare != .
drop if samplepiece != 1 & ecare == 1
drop if nwaves == 1 & ecare == 0


*===========================================================================
* SECTION 2: BASELINE VARIABLES
* (identical to 04_psm_entry.do; repeated here for self-contained running)
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
    gen pabi_fiob = g_pabi if fiob == 1
    replace pabi_fiob = i_pabi if fiob == 2
    replace pabi_fiob = k_pabi if fiob == 3
    replace pabi_fiob = m_pabi if fiob == 4
    gen patot_fiob = g_patot if fiob == 1
    replace patot_fiob = i_patot if fiob == 2
    replace patot_fiob = k_patot if fiob == 3
    replace patot_fiob = m_patot if fiob == 4
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
}

sort pidp


*===========================================================================
* SECTION 3: PROPENSITY SCORE MATCHING
*
* Sex modes:  sex removed from ematch; age kept
* Age modes:  age removed from ematch; sex kept
*===========================================================================

*--- Select outcome variable name for kmatch conditioning ---
if "`outcome'" == "pa"   local outvar "pabi_fiob"
if "`outcome'" == "alc"  local outvar "alcbi_fiob"
if "`outcome'" == "diet" local outvar "meandiet_fiob"
if "`outcome'" == "smok" local outvar ""   // smoking: no outcome conditioning

*--- Core covariates (same for all modes) ---
local core_cov "i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu nwaves i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob i.ghealth_fiob"

*--- ematch specification differs by mode ---
if inlist("`mode'", "sex_int", "sex_strat") {
    local ematch_vars "wave_fiob agefiob"    // sex excluded
}
else {
    local ematch_vars "sex wave_fiob"        // age excluded
}

if "`outvar'" != "" {
    kmatch ps ecare `core_cov' (`outvar'), ///
        ematch(`ematch_vars') nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}
else {
    kmatch ps ecare `core_cov', ///
        ematch(`ematch_vars') nn(3) ///
        idgenerate(controlid) idvar(pidp) generate(kpscore) wor
}

tab kpscore
gen matchwt = _KM_mw
gen pscore  = _KM_ps


*===========================================================================
* SECTION 4: ASSIGN ONSET WAVE TO MATCHED CONTROLS
*===========================================================================

forval k = 1/3 {
    gen double id`k' = controlid`k'
    replace id`k' = pidp if kpscore == 0
    gen ficare_tmp`k' = ficare if kpscore == 1
    sort id`k' ficare_tmp`k'
    by id`k': replace ficare_tmp`k' = ficare_tmp`k'[1]
    gen ficare_g`k' = ficare_tmp`k' if kpscore == 0
    drop ficare_tmp`k'
}

egen occ_carer = anycount(*_carebi), values(1)
replace ficare = 0 if occ_carer == 0

gen ficare_all = ficare
replace ficare_all = ficare_g1 if ficare_g1 != . & ficare_all == 0
gen preafter0 = 1 if (wave_fiob < ficare_all) & (ficare_all >= wave_fiob) & ficare_all != .
replace ficare_all = ficare_g2 if preafter0 != 1
gen preafter  = 1 if (wave_fiob < ficare_all) & (ficare_all >= wave_fiob) & ficare_all != .
replace ficare_all = ficare_g3 if preafter != 1

gen usesample = 1 if wave_fiob < ficare_all & ficare_all != .
keep if usesample == 1


*===========================================================================
* SECTION 5: COVARIATE BALANCE (ebalance)
*
* Sex modes:  sex and age included; no age group (sex already restricted or
*             being tested as interaction)
* Age modes:  sex included; age excluded (already restricted or being tested)
*===========================================================================

*--- Select outcome variable for ebalance conditioning ---
if "`outcome'" == "pa"   local eb_out "(pabi_fiob)"
if "`outcome'" == "alc"  local eb_out "(alcbi_fiob)"
if "`outcome'" == "diet" local eb_out "(meandiet_fiob)"
if "`outcome'" == "smok" local eb_out "(smok_fiob)"

if inlist("`mode'", "sex_int", "sex_strat") {
    * Age kept; sex excluded from ebalance
    ebalance kpscore agefiob i.cohab_fiob ///
        hhincome_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob `eb_out', gen(balancewt)
}
else {
    * Sex kept; age excluded from ebalance
    ebalance kpscore i.sex i.cohab_fiob ///
        hhincome_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
        ghq_fiob i.hhgroup_fiob `eb_out', gen(balancewt)
}

sort pidp


*===========================================================================
* SECTION 6: RESHAPE TO LONG FORMAT
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
        reshape long age_dv sumaudit alcbi carebi carehrs, i(pidp) j(wave)
    }
    else if "`outcome'" == "diet" {
        reshape long age_dv diet carebi meandiet, i(pidp) j(wave)
        winsor2 meandiet, suffix(_trim) cut(0 99) trim label
    }
}

by pidp (wave), sort: gen countid = _n
by pidp (wave): gen count = _N

gen wavetoset = wave - ficare_all

if "`outcome'" == "smok" {
    gen poyrtoset = wavetoset + 7
    recode poyrtoset (-1=.)
    mkspline part1 6 part2 7 part3 = poyrtoset, marginal
}
else {
    gen poyrtoset = wavetoset + 3
    mkspline part1 2 part2 3 part3 = poyrtoset, marginal
}

*--- Save long dataset (named by mode and stratum) ---
if inlist("`mode'", "sex_int", "age_int") {
    save "$ukhls/`outcome'_long_`mode'", replace
}
else {
    save "$ukhls/`outcome'_long_`mode'_`strat'", replace
}


*===========================================================================
* SECTION 7: GROWTH CURVE MODELS
*===========================================================================

*--- Select model command ---
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

*--- Shorthand label for saved margins files ---
if inlist("`mode'", "sex_int", "age_int") {
    local tag "`outcome'_`mode'"
}
else {
    local tag "`outcome'_`mode'`strat'"
}


*---------------------------------------------------------------------------
* INTERACTION TEST (sex_int and age_int modes)
*---------------------------------------------------------------------------
if "`mode'" == "sex_int" {
    `cmd' `depvar' i.poyrtoset##i.kpscore##i.sex `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm poyrtoset#kpscore#sex

    `cmd' `depvar' part1 c.part2##i.kpscore##i.sex c.part3##i.kpscore##i.sex `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#sex#c.part2
    testparm kpscore#sex#c.part3
}
else if "`mode'" == "age_int" {
    `cmd' `depvar' i.poyrtoset##i.kpscore##i.agegr_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm poyrtoset#kpscore#agegr_fiob

    `cmd' `depvar' part1 c.part2##i.kpscore##i.agegr_fiob c.part3##i.kpscore##i.agegr_fiob `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#agegr_fiob#c.part2
    testparm kpscore#agegr_fiob#c.part3
}


*---------------------------------------------------------------------------
* OUTCOME-SPECIFIC GRAPH SETTINGS
* Adjust ylabel/yscale per outcome to match the appropriate probability range
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    local ylab   "0.4(0.05)0.8"
    local yscale "range(0.4 0.8)"
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
    local ylab   "0(0.1)0.4"
    local yscale "range(0 0.4)"
    local ytitle "Probability of smoking"
}

*--- Spline xline positions differ by outcome ---
if "`outcome'" == "smok" {
    local xline1 6
    local xline2 7
    local xlab   `"0 "-7" 1 "-6" 2 "-5" 3 "-4" 4 "-3" 5 "-2" 6 "-1" 7 "0" 8 "1" 9 "2" 10 "3" 11 "4" 12 "5" 13 "6" 14 "7""'
}
else {
    local xline1 2
    local xline2 3
    local xlab   `"0 "-6" 1 "-4" 2 "-2" 3 "0" 4 "2" 5 "4""'
}

*--- Stratum labels for graph titles and legends ---
if "`mode'" == "sex_strat" {
    if `strat' == 0 {
        local stratlabel "men"
        local straттitle "Male"
        local leg1 "Male – Non-caregiver"
        local leg2 "Male – Entering caregiving"
    }
    else {
        local stratlabel "women"
        local strattitle "Female"
        local leg1 "Female – Non-caregiver"
        local leg2 "Female – Entering caregiving"
    }
}
else if "`mode'" == "age_strat" {
    if `strat' == 0 {
        local stratlabel "age0"
        local strattitle "Early adulthood (16–29)"
        local leg1 "16–29: Non-caregiver"
        local leg2 "16–29: Entering caregiving"
    }
    else if `strat' == 1 {
        local stratlabel "age1"
        local strattitle "Early mid-adulthood (30–49)"
        local leg1 "30–49: Non-caregiver"
        local leg2 "30–49: Entering caregiving"
    }
    else if `strat' == 2 {
        local stratlabel "age2"
        local strattitle "Late mid-adulthood (50–64)"
        local leg1 "50–64: Non-caregiver"
        local leg2 "50–64: Entering caregiving"
    }
    else if `strat' == 3 {
        local stratlabel "age3"
        local strattitle "Late adulthood (65+)"
        local leg1 "65+: Non-caregiver"
        local leg2 "65+: Entering caregiving"
    }
}


*---------------------------------------------------------------------------
* STRATIFIED MODELS AND GRAPHS (sex_strat and age_strat modes)
*---------------------------------------------------------------------------
if inlist("`mode'", "sex_strat", "age_strat") {

    *--- Main trajectory: caregiver vs non-caregiver ---
    `cmd' `depvar' i.poyrtoset `wave_term' if kpscore==1 [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoset, saving(`tag'_1, replace)

    `cmd' `depvar' i.poyrtoset `wave_term' if kpscore==0 [pw=balancewt], vce(cluster hhid_fiob)
    margins i.poyrtoset, saving(`tag'_0, replace)

    combomarginsplot `tag'_0 `tag'_1, ///
        labels("`leg1'" "`leg2'") ///
        savefile(`tag'_intro, replace)

    mplotoffset using `tag'_intro, ///
        xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
        xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
        title("`ytitle'" ///
              "`strattitle'", size(large)) ///
        xtitle("Years centred on caregiving onset year", size(large)) ///
        ytitle("Predictive AME", size(medium)) ///
        xlabel(`xlab', labsize(large)) ///
        ylabel(`ylab', labsize(large)) ///
        yscale(`yscale') ///
        legend(size(medium))

    *--- Spline interaction test within stratum ---
    `cmd' `depvar' part1 c.part2##i.kpscore c.part3##i.kpscore `wave_term' [pw=balancewt], vce(cluster hhid_fiob)
    testparm kpscore#c.part2
    testparm kpscore#c.part3

    *--- Save individual publication panel (no legend, no axis labels) ---
    mplotoffset using `tag'_intro, ///
        xline(`xline1', lpattern(shortdash) lcolor(black)) ///
        xline(`xline2', lpattern(shortdash) lcolor(black)) ///
        xtitle("") ytitle("") ///
        xlabel(, nolabel) ylabel(, nolabel) ///
        yscale(`yscale') ///
        title("", size(vlarge)) ///
        legend(off) ///
        graphregion(color(white))
    graph save `outcome'_`stratlabel'_panel.gph, replace

    *--- Save individual publication panel WITH legend ---
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
    graph save `outcome'_`stratlabel'_panel_legend.gph, replace
}

di as result "===== 05_psm_entry_strata.do completed: outcome=`outcome' mode=`mode' strat=`strat' ====="


*===========================================================================
* COMBINED PUBLICATION GRAPHS
* Run this section ONLY after all strata have been run and all individual
* panel .gph files and margins files exist on disk.
*
* Comment out or wrap in an if-block to prevent accidental re-execution.
*===========================================================================

/*

*--- Combined: all strata caregivers only (for within-stratum comparison) ---

* Sex strata
if "`mode'" == "sex_strat" & "`outcome'" != "smok" {
    combomarginsplot ///
        `outcome'_sex_strat0_0 `outcome'_sex_strat0_1 ///
        `outcome'_sex_strat1_0 `outcome'_sex_strat1_1, ///
        labels("Male – Non-caregiver"    "Male – Entering caregiving" ///
               "Female – Non-caregiver"  "Female – Entering caregiving") ///
        savefile(`outcome'_sex_combined, replace)

    mplotoffset using `outcome'_sex_combined, ///
        xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
        xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
        title("`ytitle' by sex", size(large)) ///
        xtitle("Years centred on caregiving onset year", size(large)) ///
        ytitle("Predictive AME", size(medium)) ///
        xlabel(`xlab', labsize(large)) ///
        ylabel(`ylab', labsize(large)) ///
        yscale(`yscale') ///
        legend(size(medium))
}

* Age strata: caregivers only across all four groups
if "`mode'" == "age_strat" {
    combomarginsplot ///
        `outcome'_age_strat0_1 `outcome'_age_strat1_1 ///
        `outcome'_age_strat2_1 `outcome'_age_strat3_1, ///
        labels("Early adulthood (16–29)"      ///
               "Early mid-adulthood (30–49)"  ///
               "Late mid-adulthood (50–64)"   ///
               "Late adulthood (65+)") ///
        savefile(`outcome'_ageall, replace)

    mplotoffset using `outcome'_ageall, ///
        xline(`xline1', lpattern(shortdash) lcolor(gray)) ///
        xline(`xline2', lpattern(shortdash) lcolor(gray)) ///
        title("`ytitle' among caregiving entrants, by age group", size(large)) ///
        xtitle("Years centred on caregiving onset year", size(large)) ///
        ytitle("Predictive AME", size(medium)) ///
        xlabel(`xlab', labsize(large)) ///
        ylabel(`ylab', labsize(large)) ///
        yscale(`yscale') ///
        legend(size(medium))
}

*--- Publication panel: all six strata (men, women, 4 age groups) without legend ---
graph combine ///
    `outcome'_men_panel.gph    `outcome'_women_panel.gph ///
    `outcome'_age0_panel.gph   `outcome'_age1_panel.gph ///
    `outcome'_age2_panel.gph   `outcome'_age3_panel.gph, ///
    cols(6) xcommon ycommon ///
    imargin(0 0 0 0) ///
    graphregion(color(white)) ///
    title("", size(large)) ///
    xsize(15) ysize(3)

graph export "`outcome'_intro_strata.jpg", as(jpg) name("Graph") quality(90) replace

*--- Publication panel: all six strata WITH legend ---
graph combine ///
    `outcome'_men_panel_legend.gph    `outcome'_women_panel_legend.gph ///
    `outcome'_age0_panel_legend.gph   `outcome'_age1_panel_legend.gph ///
    `outcome'_age2_panel_legend.gph   `outcome'_age3_panel_legend.gph, ///
    cols(6) xcommon ycommon ///
    imargin(0 0 0 0) ///
    graphregion(color(white)) ///
    title("", size(large)) ///
    xsize(15) ysize(3)

graph export "`outcome'_intro_strata_legend.jpg", as(jpg) name("Graph") quality(90) replace

*/
