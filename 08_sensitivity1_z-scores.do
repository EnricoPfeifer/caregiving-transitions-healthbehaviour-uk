/*===========================================================================
  08_sensitivity.do

  PURPOSE: Sensitivity analyses for the two dichotomised outcomes,
           reported in the Supplementary Material. Three blocks, selected
           with `block' below:

             fe    -- Fixed effects models re-estimated on the continuous
                      and standardised outcome scales
             psm   -- Piecewise growth curve models re-estimated on the
                      standardised outcome scale, in the matched samples
             dist  -- Distributional checks for the AUDIT-C count outcome
                      (Poisson, negative binomial, zero-inflated variants)

  RATIONALE: The main analysis dichotomises physical activity (inactive
           vs active, from weekly minutes of moderate and vigorous
           activity) and alcohol consumption (AUDIT-C at or above the
           sex-specific threshold). Both are re-estimated here on their
           underlying scale, and again standardised, to establish that
           the findings are not an artefact of the cut-points.

  OUTCOMES SUPPORTED: pa, alc

  SAMPLES SUPPORTED (psm block):
    entry             -- entry into caregiving vs non-caregivers
    exit_nocare       -- exit from caregiving vs non-caregivers
    exit_alwayscare   -- exit from caregiving vs continuing caregivers

  PREREQUISITES:
    fe   : 03_fixed_effects.do run for this outcome, saving the long
           dataset as "$ukhls/`outcome'_long_fe"
    psm  : 04_psm_entry.do (entry) or 06_psm_exit.do (exit) run for this
           outcome and comparison group
    dist : as for fe

  STANDARDISATION: egen std() is applied to the dataset in memory, so
           the fixed effects and matched-sample blocks standardise over
           different samples. Coefficients are therefore in standard
           deviations of the sample used in that model. This is stated in
           the Supplement.

  USER-WRITTEN COMMANDS: winsor2, combomarginsplot, mplotoffset, and
           prcounts. prcounts is distributed with SPost (Long and Freese)
           rather than SSC; see the README for installation.

  AUTHOR: Enrico Pfeifer
===========================================================================*/


*---------------------------------------------------------------------------
* SET BLOCK, OUTCOME AND SAMPLE HERE
*
*   block   : fe | psm | dist
*   outcome : pa | alc
*   sample  : entry | exit_nocare | exit_alwayscare   (psm block only)
*---------------------------------------------------------------------------
local block   "fe"
local outcome "pa"
local sample  "entry"


*---------------------------------------------------------------------------
* PATHS
* Set $ukhls to the directory holding the UKHLS Stata files and $out to a
* directory for figures. No absolute paths appear elsewhere in this file.
*---------------------------------------------------------------------------
clear all
set more off
set maxvar 10000
set cformat %4.2f

global ukhls ""      // e.g. "C:/data/ukhls/"
global out   ""      // e.g. "C:/output/"

if "$ukhls" == "" {
    display as error "Set the global macro ukhls before running this file."
    exit 198
}

do 00_graph_settings.do

if !inlist("`outcome'", "pa", "alc") {
    display as error "This file applies to the pa and alc outcomes only."
    exit 198
}


*---------------------------------------------------------------------------
* Outcome-specific settings
*
*   cvar : continuous counterpart of the dichotomised outcome
*   trim : whether the continuous outcome is winsorised before use
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    local cvar   "patot"
    local trim   1                  // winsorise at the 99th percentile
    local clab   "Weekly minutes of moderate and vigorous activity"
    local zlab   "Standardised physical activity minutes"
}
else if "`outcome'" == "alc" {
    local cvar   "sumaudit"
    local trim   0                  // AUDIT-C is bounded at 0-12
    local clab   "AUDIT-C score"
    local zlab   "Standardised AUDIT-C score"
}


*===========================================================================
*===========================================================================
* BLOCK 1: FIXED EFFECTS MODELS ON THE CONTINUOUS AND STANDARDISED SCALE
*===========================================================================
*===========================================================================

if "`block'" == "fe" {

    use "$ukhls/`outcome'_long_fe", clear
    xtset pidp wave

    *-----------------------------------------------------------------------
    * Transition indicators
    * 03_fixed_effects.do already creates caregiving_intro and
    * caregiving_exit. They are recreated here only if absent, so that this
    * file can also be run on a long dataset saved at an earlier stage.
    *-----------------------------------------------------------------------
    capture confirm variable caregiving_intro
    if _rc {
        bysort pidp (wave): gen caregiving_intro = (carebi[_n-1] == 0 & carebi == 1)
        label var caregiving_intro "Entry into caregiving"
    }
    capture confirm variable caregiving_exit
    if _rc {
        bysort pidp (wave): gen caregiving_exit = (carebi[_n-1] == 1 & carebi == 0)
        label var caregiving_exit "Exit from caregiving"
    }

    *-----------------------------------------------------------------------
    * Continuous and standardised outcome
    *-----------------------------------------------------------------------
    if `trim' == 1 {
        winsor2 `cvar', suffix(_trim) cut(0 99) trim label
        local dvc "`cvar'_trim"
    }
    else {
        local dvc "`cvar'"
    }

    summarize `dvc', detail
    egen `outcome'_z = std(`dvc')
    label var `outcome'_z "`zlab'"
    local dvz "`outcome'_z"

    *-----------------------------------------------------------------------
    * Models. Standard errors are clustered on the individual, matching the
    * main fixed effects specification.
    *-----------------------------------------------------------------------
    display _newline as text "{hline 75}"
    display as text "Fixed effects, continuous and standardised: `outcome'"
    display as text "{hline 75}"

    *--- [STATUS] Caregiving status, continuous scale ---
    xtreg `dvc' carebi i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_c_status

    *--- [STATUS] Caregiving status, standardised scale ---
    xtreg `dvz' carebi i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_z_status

    *--- [ENTRY] ---
    xtreg `dvc' caregiving_intro i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_c_entry
    xtreg `dvz' caregiving_intro i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_z_entry

    *--- [EXIT] ---
    xtreg `dvc' caregiving_exit i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_c_exit
    xtreg `dvz' caregiving_exit i.wave, fe i(pidp) vce(cluster pidp)
    est store fe_z_exit

    *--- Comparison table ---
    estimates table fe_c_status fe_z_status fe_c_entry fe_z_entry ///
                    fe_c_exit   fe_z_exit, ///
        b(%9.3f) se(%9.3f) stats(N) ///
        keep(carebi caregiving_intro caregiving_exit)

    *-----------------------------------------------------------------------
    * Count specification for the AUDIT-C score
    *
    * AUDIT-C is a bounded count, so a fixed effects Poisson model is
    * reported alongside the linear model as a specification check.
    *-----------------------------------------------------------------------
    if "`outcome'" == "alc" {
        xtpoisson `dvc' carebi i.wave, fe i(pidp) irr vce(robust)
        est store fe_pois_status

        xtpoisson `dvc' caregiving_exit i.wave, fe i(pidp) irr vce(robust)
        est store fe_pois_exit

        xtnbreg `dvc' carebi i.wave, fe irr
    }

    if "`outcome'" == "pa" {
        xtpoisson `dvc' carebi i.wave, fe i(pidp) irr vce(robust)
        est store fe_pois_status

        xtpoisson `dvc' caregiving_exit i.wave, fe i(pidp) irr vce(robust)
        est store fe_pois_exit
    }

}


*===========================================================================
*===========================================================================
* BLOCK 2: PIECEWISE GROWTH CURVES ON THE STANDARDISED SCALE
*
* Same matched sample, balance weights and spline structure as the main
* analysis; only the outcome scale changes.
*===========================================================================
*===========================================================================

if "`block'" == "psm" {

    *-----------------------------------------------------------------------
    * Sample-specific settings
    *   timevar : centred time variable created in 04 or 06
    *   lab0    : comparison group label
    *   lab1    : transition group label
    *-----------------------------------------------------------------------
    if "`sample'" == "entry" {
        local dta     "`outcome'_long_afterpsm"
        local timevar "poyrtoset"
        local suffix  ""
        local lab0    "Non-caregiver"
        local lab1    "Entered caregiving"
        local xtitle  "Years centred on caregiving entry"
    }
    else if "`sample'" == "exit_nocare" {
        local dta     "`outcome'_exit_nocare_long_afterPSM"
        local timevar "poyrtoexit"
        local suffix  "_exnocare"
        local lab0    "Non-caregiver"
        local lab1    "Exited caregiving"
        local xtitle  "Years centred on caregiving exit"
    }
    else if "`sample'" == "exit_alwayscare" {
        local dta     "`outcome'_exit_alwayscare_long_afterPSM"
        local timevar "poyrtoexit"
        local suffix  "_exalways"
        local lab0    "Continued caregiving"
        local lab1    "Exited caregiving"
        local xtitle  "Years centred on caregiving exit"
    }
    else {
        display as error "sample must be entry, exit_nocare or exit_alwayscare."
        exit 198
    }

    use "$ukhls/`dta'", clear

    *-----------------------------------------------------------------------
    * Standardised outcome
    *-----------------------------------------------------------------------
    if `trim' == 1 {
        winsor2 `cvar', suffix(_trim) cut(0 99) trim label
        egen `outcome'_z = std(`cvar'_trim)
    }
    else {
        egen `outcome'_z = std(`cvar')
    }
    label var `outcome'_z "`zlab'"
    local dvz "`outcome'_z"

    *-----------------------------------------------------------------------
    * Group-specific trajectories
    *-----------------------------------------------------------------------
    regress `dvz' i.`timevar' if kpscore == 0 [pw=balancewt]
    margins i.`timevar', saving("$out/`outcome'_standard0`suffix'", replace)

    regress `dvz' i.`timevar' if kpscore == 1 [pw=balancewt]
    margins i.`timevar', saving("$out/`outcome'_standard1`suffix'", replace)

    combomarginsplot "$out/`outcome'_standard0`suffix'" ///
                     "$out/`outcome'_standard1`suffix'", ///
        labels("`lab0'" "`lab1'") ///
        savefile("$out/`outcome'_stand`suffix'", replace)

    mplotoffset using "$out/`outcome'_stand`suffix'", ///
        xline(2, lpattern(shortdash) lcolor(black)) ///
        xline(3, lpattern(shortdash) lcolor(black)) ///
        title("`zlab'", size(large)) ///
        ytitle("z-scores", size(large)) ///
        xtitle("`xtitle'", size(large)) ///
        xlabel(0 "-6" 1 "-4" 2 "-2" 3 "0" 4 "2" 5 "4", labsize(large)) ///
        ylabel(-0.1(0.05)0.2, labsize(large)) ///
        yscale(range(-0.1 0.2)) ///
        legend(position(12) ring(1) rows(1) size(large) ///
               label(1 "`lab0'") label(2 "`lab1'") ///
               region(color(white))) ///
        graphregion(color(white))

    graph export "$out/`outcome'_stand`suffix'.jpg", as(jpg) quality(90) replace

    *-----------------------------------------------------------------------
    * Piecewise interaction test
    * The two terms reported in the Supplement are c.part2#1.kpscore
    * (during the transition) and c.part3#1.kpscore (after the transition).
    *-----------------------------------------------------------------------
    regress `dvz' part1 c.part2##i.kpscore c.part3##i.kpscore ///
        [pw=balancewt], vce(cluster hhid_fiob)
    est store psm_`outcome'`suffix'


    *-----------------------------------------------------------------------
    * Combined exit figure: both comparison groups on one plot
    *
    * Requires the exit_nocare and exit_alwayscare runs above to have been
    * completed, so that all three margins files exist. Run last.
    *-----------------------------------------------------------------------
    /*
    combomarginsplot "$out/`outcome'_standard0_exalways" ///
                     "$out/`outcome'_standard0_exnocare" ///
                     "$out/`outcome'_standard1_exalways", ///
        labels("Continued caregiving" ///
               "Non-caregiver" ///
               "Exited caregiving") ///
        savefile("$out/`outcome'_stand_exit", replace)

    mplotoffset using "$out/`outcome'_stand_exit", ///
        xline(2, lpattern(shortdash) lcolor(gray)) ///
        xline(3, lpattern(shortdash) lcolor(gray)) ///
        title("`zlab' by caregiving status around exit", size(large)) ///
        xtitle("Years centred on caregiving exit", size(large)) ///
        ytitle("z-scores", size(medium)) ///
        xlabel(0 "-6" 1 "-4" 2 "-2" 3 "0" 4 "2" 5 "4", labsize(large)) ///
        ylabel(-0.1(0.05)0.2, labsize(large)) ///
        yscale(range(-0.1 0.2)) ///
        legend(size(medium))

    graph export "$out/`outcome'_stand_exit.jpg", as(jpg) quality(90) replace
    */

}


*===========================================================================
*===========================================================================
* BLOCK 3: DISTRIBUTIONAL CHECKS FOR THE AUDIT-C COUNT OUTCOME
*
* AUDIT-C is a bounded count with a large number of zeros. Observed and
* predicted distributions are compared across four specifications to
* establish which distributional assumption fits best. These are pooled
* models used solely for distributional comparison; they are not the
* substantive estimates reported in the paper.
*
* Requires prcounts (SPost, Long and Freese).
*===========================================================================
*===========================================================================

if "`block'" == "dist" {

    if "`outcome'" != "alc" {
        display as error "The distributional checks apply to the alc outcome."
        exit 198
    }

    use "$ukhls/`outcome'_long_fe", clear
    xtset pidp wave

    local dvc "sumaudit"

    *--- Poisson ---
    poisson `dvc' carebi i.wave, irr
    capture drop numdiff*
    prcounts numdiff_, plot max(10)

    label var numdiff_val  "AUDIT-C score"
    label var numdiff_obeq "Observed proportion"
    label var numdiff_preq "Predicted proportion (Poisson)"

    twoway connected numdiff_preq numdiff_obeq numdiff_val, ///
        ytitle("Probability") ylabel(0(.1).8) xlabel(0(1)6) ///
        title("AUDIT-C: observed vs Poisson prediction")
    graph export "$out/alc_dist_poisson.jpg", as(jpg) quality(90) replace

    *--- Negative binomial ---
    nbreg `dvc' carebi i.wave
    capture drop numdiff*
    prcounts numdiff_, plot max(10)

    label var numdiff_val  "AUDIT-C score"
    label var numdiff_obeq "Observed proportion"
    label var numdiff_preq "Predicted proportion (negative binomial)"

    twoway connected numdiff_preq numdiff_obeq numdiff_val, ///
        ytitle("Probability") ylabel(0(.1).8) xlabel(0(1)6) ///
        title("AUDIT-C: observed vs negative binomial prediction")
    graph export "$out/alc_dist_nbreg.jpg", as(jpg) quality(90) replace

    *--- Zero-inflated Poisson ---
    zip `dvc' carebi i.wave, inflate(carebi i.wave)
    capture drop numdiff*
    prcounts numdiff_, plot max(10)

    label var numdiff_val  "AUDIT-C score"
    label var numdiff_obeq "Observed proportion"
    label var numdiff_preq "Predicted proportion (zero-inflated Poisson)"

    twoway connected numdiff_preq numdiff_obeq numdiff_val, ///
        ytitle("Probability") ylabel(0(.1).8) xlabel(0(1)6) ///
        title("AUDIT-C: observed vs zero-inflated Poisson prediction")
    graph export "$out/alc_dist_zip.jpg", as(jpg) quality(90) replace

    *--- Zero-inflated negative binomial ---
    zinb `dvc' carebi i.wave, inflate(carebi i.wave)
    capture drop numdiff*
    prcounts numdiff_, plot max(10)

    label var numdiff_val  "AUDIT-C score"
    label var numdiff_obeq "Observed proportion"
    label var numdiff_preq "Predicted proportion (zero-inflated negative binomial)"

    twoway connected numdiff_preq numdiff_obeq numdiff_val, ///
        ytitle("Probability") ylabel(0(.1).8) xlabel(0(1)6) ///
        title("AUDIT-C: observed vs zero-inflated negative binomial prediction")
    graph export "$out/alc_dist_zinb.jpg", as(jpg) quality(90) replace

}


*===========================================================================
* NOTES
*
* - The y-axis range in the standardised figures is fixed at -0.1 to 0.2 so
*   that entry and exit plots are directly comparable. Widen it if an
*   estimate falls outside that range for another outcome.
*===========================================================================
