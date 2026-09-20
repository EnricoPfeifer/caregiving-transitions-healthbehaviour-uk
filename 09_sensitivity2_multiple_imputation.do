/*===========================================================================
  09_multiple_imputation.do

  PURPOSE: Multiple imputation by chained equations applied to the matched
           analysis, reported in the Supplementary Material. This file
           consolidates the eleven outcome-by-sample scripts into one
           parameterised routine.

  RATIONALE: Listwise deletion at the matching stage removes a non-trivial
           share of the sample, whereas the fixed effects models are
           structurally robust to missingness on time-invariant
           covariates. Imputation is therefore most informative where the
           missingness is concentrated, which is the matched analysis.

  METHOD: The "within" or MIte approach of Leyrat et al. (2019). Matching,
           balancing and estimation are carried out separately inside each
           imputed dataset, and the estimates are pooled with Rubin's
           rules. Only baseline covariates and the baseline outcome enter
           the imputation model; repeated wave-specific outcomes do not.
           This trade-off is stated in the Supplement.

  PREREQUISITE: the pre-matching dataset for this outcome and sample, saved
           by 04_psm_entry.do or 06_psm_exit.do immediately before the
           matching step:

             04: save "$ukhls/`outcome'_prematch_entry", replace
             06: save "$ukhls/`outcome'_prematch_`sample'", replace

           Loading these rather than rebuilding the sample keeps a single
           definition of every variable across the repository, and removes
           the filename collisions present in the original scripts.

  OUTCOMES SUPPORTED: pa, alc, diet, smok
  SAMPLES SUPPORTED:  entry, exit_nocare, exit_alwayscare

  RUNTIME: Matching and balancing are repeated M times for each outcome and
           sample. Run with M = 5 first to confirm the pipeline executes,
           then set M to the number reported in the Supplement.

  AUTHOR: Enrico Pfeifer
===========================================================================*/


*---------------------------------------------------------------------------
* SET OUTCOME, SAMPLE AND NUMBER OF IMPUTATIONS HERE
*
*   outcome : pa | alc | diet | smok
*   sample  : entry | exit_nocare | exit_alwayscare
*   M       : number of imputations reported in the Supplement
*---------------------------------------------------------------------------
local outcome "pa"
local sample  "entry"
local M       10


*---------------------------------------------------------------------------
* PATHS
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


*---------------------------------------------------------------------------
* Sample-specific settings
*
*   treat    : treatment indicator
*   wavevar  : wave of the transition, for cases
*   allvar   : transition wave assigned to cases and matched controls
*   centred  : time centred on the transition
*---------------------------------------------------------------------------
if "`sample'" == "entry" {
    local treat   "ecare"
    local wavevar "ficare"
    local allvar  "ficare_all"
    local centred "poyrtoset"
    local lab0    "Non-caregiver"
    local lab1    "Entered caregiving"
    local xtitle  "Waves centred on caregiving entry"
}
else if inlist("`sample'", "exit_nocare", "exit_alwayscare") {
    local treat   "exitcarebi"
    local wavevar "exitcarewave"
    local allvar  "exitcare_all"
    local centred "poyrtoexit"
    local lab1    "Exited caregiving"
    local xtitle  "Waves centred on caregiving exit"
    if "`sample'" == "exit_nocare"     local lab0 "Non-caregiver"
    if "`sample'" == "exit_alwayscare" local lab0 "Continued caregiving"
}
else {
    display as error "sample must be entry, exit_nocare or exit_alwayscare."
    exit 198
}

local prematch "`outcome'_prematch_`sample'"
local tag      "`outcome'_`sample'"

*--- Spline knots and centring follow the wave structure of the outcome.
*    Smoking is observed annually across waves 5 to 13; the remaining
*    outcomes are observed at waves 7, 9, 11 and 13. ---
if "`outcome'" == "smok" {
    local knot1     6
    local knot2     7
    local shift     7
    local wave_term "i.wave"
}
else {
    local knot1     2
    local knot2     3
    local shift     3
    local wave_term ""
}


*---------------------------------------------------------------------------
* Outcome-specific settings
*
*   basevar : outcome at first observation, a matching covariate
*   depvar  : outcome in long format
*   cmd     : estimation command
*   imp_*   : incomplete variables, grouped by imputation method
*   sfterm  : SF-12 physical component, matched and balanced on for the
*             physical activity outcome only
*   rslong  : variables reshaped to long format
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {
    local basevar    "pabi_fiob"
    local depvar     "pabi"
    local cmd        "logit"
    local imp_pmm    "ghq_fiob sf_fiob"
    local imp_ologit "edu iwealth_fiob"
    local imp_logit  "pabi_fiob cohab_fiob ghealth_fiob"
    local imp_mlogit "oclass_fiob ethnicity workbi_fiob"
    local sfterm     "sf_fiob"
    local rslong     "age_dv patot pabi watot motot vigtot carebi carehrs carecat"
    local ytitle     "Probability of physical inactivity"
}
else if "`outcome'" == "alc" {
    local basevar    "alcbi_fiob"
    local depvar     "alcbi"
    local cmd        "logit"
    local imp_pmm    "ghq_fiob"
    local imp_ologit "edu iwealth_fiob"
    local imp_logit  "alcbi_fiob cohab_fiob ghealth_fiob"
    local imp_mlogit "oclass_fiob ethnicity workbi_fiob"
    local sfterm     ""
    local rslong     "age_dv sumaudit alcbi carebi carehrs carecat"
    local ytitle     "Probability of problematic drinking"
}
else if "`outcome'" == "diet" {
    local basevar    "meandiet_fiob"
    local depvar     "meandiet_trim"
    local cmd        "regress"
    local imp_pmm    "ghq_fiob meandiet_fiob"
    local imp_ologit "edu iwealth_fiob"
    local imp_logit  "cohab_fiob ghealth_fiob"
    local imp_mlogit "oclass_fiob ethnicity workbi_fiob"
    local sfterm     ""
    local rslong     "age_dv meandiet carebi carehrs carecat"
    local ytitle     "Portions of fruit and vegetables"
}
else if "`outcome'" == "smok" {
    local basevar    "smok_fiob"
    local depvar     "smok"
    local cmd        "logit"
    local imp_pmm    "ghq_fiob"
    local imp_ologit "edu iwealth_fiob"
    local imp_logit  "smok_fiob cohab_fiob ghealth_fiob"
    local imp_mlogit "oclass_fiob ethnicity workbi_fiob"
    local sfterm     ""
    local rslong     "smok ciggy age_dv carebi carehrs carecat"
    local ytitle     "Probability of current smoking"
}
else {
    display as error "outcome must be pa, alc, diet or smok."
    exit 198
}


*===========================================================================
* STAGE 1: IMPUTATION
*
* The pre-matching dataset already carries every sample restriction that
* matching impose
*
* The treatment indicator enters the imputation model as a complete
* predictor, so that the imputation model is compatible with the analysis
* model that follows.
*===========================================================================

use "$ukhls/`prematch'", clear

*--- Confirm which variables are incomplete before imputing ---
mdesc `imp_pmm' `imp_ologit' `imp_logit' `imp_mlogit' ///
      `treat' sex agefiob wave_fiob nwaves hhgroup_fiob children_fiob

xtset, clear
mi set wide

mi register imputed `imp_pmm' `imp_ologit' `imp_logit' `imp_mlogit'

*--- Complete predictors. Wave-specific outcomes are deliberately excluded;
*    only the baseline outcome enters, through `basevar'. ---
mi register regular `treat' sex agefiob wave_fiob nwaves ///
    hhid_fiob `wavevar' hhgroup_fiob children_fiob

mi impute chained ///
    (pmm, knn(5)) `imp_pmm' ///
    (ologit)      `imp_ologit' ///
    (logit)       `imp_logit' ///
    (mlogit)      `imp_mlogit' ///
    = i.`treat' i.sex agefiob i.wave_fiob nwaves ///
      i.hhgroup_fiob i.children_fiob ///
    , add(`M') rseed(20260602) augment dots

save "$ukhls/`tag'_imputed", replace

*--- Diagnostic: number of imputations implied by the target precision,
*    reported for transparency. The number actually used is set by `M'
*    above and is stated in the Supplement. ---
mi estimate, dots: `cmd' `basevar' i.`treat' i.sex agefiob i.wave_fiob ///
    nwaves i.oclass_fiob i.cohab_fiob i.edu i.iwealth_fiob ///
    i.workbi_fiob i.ethnicity i.hhgroup_fiob i.children_fiob ghq_fiob ///
    i.ghealth_fiob `sfterm'

capture how_many_imputations


*===========================================================================
* STAGE 2: REPEAT THE FULL PIPELINE INSIDE EACH IMPUTED DATASET
*
* The matched sample itself depends on the imputed covariates, so matching
* and balancing are repeated rather than carried over.
*===========================================================================

forvalues m = 1/`M' {
    quietly {
        use "$ukhls/`tag'_imputed", clear
        mi extract `m', clear

        *--- Propensity score matching ---
        kmatch ps `treat' i.oclass_fiob i.cohab_fiob i.hhgroup_fiob i.edu ///
            nwaves i.workbi_fiob i.iwealth_fiob i.ethnicity ghq_fiob ///
            i.ghealth_fiob `sfterm' i.children_fiob (`basevar'), ///
            ematch(sex wave_fiob agefiob) nn(3) ///
            idgenerate(controlid) idvar(pidp) generate(kpscore) wor

        gen matchwt = _KM_mw
        gen pscore  = _KM_ps

        *--- Assign the transition wave of the matched case to each control.
        *    Up to three matched cases are considered in turn, taking the
        *    first that leaves at least one observation on either side of
        *    the transition. ---
        forvalues j = 1/3 {
            gen double id`j' = controlid`j'
            replace id`j' = pidp if kpscore == 0
            gen tw`j' = `wavevar' if kpscore == 1
            sort id`j' tw`j'
            by id`j': replace tw`j' = tw`j'[1]
            gen `allvar'_g`j' = tw`j' if kpscore == 0
        }

        if "`sample'" == "entry" {

            *--- Non-caregivers have no onset wave of their own ---
            egen occ_carer = anycount(*_carebi), values(1)
            label var occ_carer "Number of occasions observed as a caregiver"
            replace `wavevar' = 0 if occ_carer == 0

            gen `allvar' = `wavevar'
            replace `allvar' = `allvar'_g1 if `allvar'_g1 != . & `allvar' == 0

            gen preafter0 = 1 if wave_fiob < `allvar' & `allvar' != .
            replace `allvar' = `allvar'_g2 if preafter0 != 1

            gen preafter = 1 if wave_fiob < `allvar' & `allvar' != .
            replace `allvar' = `allvar'_g3 if preafter != 1

            gen usesample = 1 if wave_fiob < `allvar' & `allvar' != .
        }
        else {

            gen `allvar' = `wavevar'
            replace `allvar' = `allvar'_g1 if `allvar'_g1 != . & `allvar' == .

            gen preafter0 = 1 if wave_fiob < `allvar' & lastwave >= `allvar' ///
                & `allvar' != .
            replace `allvar' = `allvar'_g2 if preafter0 != 1

            gen preafter = 1 if wave_fiob < `allvar' & lastwave >= `allvar' ///
                & `allvar' != .
            replace `allvar' = `allvar'_g3 if preafter != 1

            gen usesample = 1 if wave_fiob < `allvar' & lastwave >= `allvar' ///
                & `allvar' != .
        }

        keep if usesample == 1

        *--- Covariate balance ---
        ebalance kpscore agefiob i.sex i.cohab_fiob nwaves ///
            i.iwealth_fiob i.edu i.oclass_fiob i.ghealth_fiob ///
            ghq_fiob i.hhgroup_fiob i.workbi_fiob ///
            i.ethnicity `sfterm' i.children_fiob (`basevar'), gen(balancewt)

        *--- Long format ---
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
            reshape long `rslong', i(pidp) j(waves)
            gen wave = waves - 4
        }
        else {
            rename g_* *1
            rename i_* *2
            rename k_* *3
            rename m_* *4
            reshape long `rslong', i(pidp) j(wave)
        }

        if "`outcome'" == "diet" {
            winsor2 meandiet, suffix(_trim) cut(0 99) trim label
        }

        by pidp (wave), sort: gen countid = _n
        by pidp (wave): gen count = _N

        *--- Time centred on the transition ---
        gen wavetoset = wave - `allvar'
        gen `centred' = wavetoset + `shift'
        if "`outcome'" == "smok" recode `centred' (min/-1 = .)
        label var `centred' "Waves before and after the transition (plus `shift')"

        mkspline part1 `knot1' part2 `knot2' part3 = `centred', marginal

        *--- Piecewise model. Stored before margins, because margins, post
        *    overwrites e(). ---
        `cmd' `depvar' part1 c.part2##i.kpscore c.part3##i.kpscore ///
            `wave_term' [pw=balancewt], vce(cluster hhid_fiob)

        matrix b`m' = e(b)
        matrix V`m' = e(V)
        if `m' == 1 local pnames : colfullnames e(b)

        *--- Group-specific predicted trajectories on a fixed grid, so that
        *    every imputation returns the same cells in the same order ---
        `cmd' `depvar' i.`centred' `wave_term' if kpscore == 1 ///
            [pw=balancewt], vce(cluster hhid_fiob)
        margins i.`centred', post
        matrix mb`m'  = e(b)
        matrix mbV`m' = e(V)

        `cmd' `depvar' i.`centred' `wave_term' if kpscore == 0 ///
            [pw=balancewt], vce(cluster hhid_fiob)
        margins i.`centred', post
        matrix mc`m'  = e(b)
        matrix mcV`m' = e(V)
    }
    display as text "imputation `m' of `M' complete"
}


*===========================================================================
* STAGE 3: POOL WITH RUBIN'S RULES
*
*   Qbar = mean of the M estimates
*   Ubar = mean within-imputation variance
*   B    = between-imputation variance
*   T    = Ubar + (1 + 1/M) * B
*
* Confidence intervals use the normal approximation rather than a
* Barnard-Rubin degrees of freedom correction. The choice is recorded here
* so that it is visible to anyone reproducing the Supplement.
*===========================================================================

mata:
void rubin_pool(string scalar bstem, string scalar vstem, real scalar M,
                string scalar qout, string scalar seout)
{
    real matrix  b1, Qsum, Usum, Ball, bm, Vm, Qbar, Ubar, Bvar, d, T
    real scalar  k, m

    b1   = st_matrix(bstem + "1")
    k    = cols(b1)
    Qsum = J(1, k, 0)
    Usum = J(k, k, 0)
    Ball = J(M, k, .)

    for (m = 1; m <= M; m++) {
        bm   = st_matrix(bstem + strofreal(m))
        Vm   = st_matrix(vstem + strofreal(m))
        Qsum = Qsum + bm
        Usum = Usum + Vm
        Ball[m, .] = bm
    }

    Qbar = Qsum / M
    Ubar = Usum / M

    Bvar = J(k, k, 0)
    for (m = 1; m <= M; m++) {
        d    = Ball[m, .] :- Qbar
        Bvar = Bvar + d' * d
    }
    Bvar = Bvar / (M - 1)

    T = Ubar + (1 + 1/M) * Bvar

    st_matrix(qout, Qbar)
    st_matrix(seout, sqrt(diagonal(T))')
}
end

*--- Pooled piecewise coefficients ---
mata: rubin_pool("b", "V", `M', "Qbar", "SEp")

mata:
Q  = st_matrix("Qbar")
SE = st_matrix("SEp")
Z  = Q :/ SE
P  = 2 * normal(-abs(Z))
LB = Q :- 1.96 :* SE
UB = Q :+ 1.96 :* SE
st_matrix("results", (Q', LB', UB', P'))
end

matrix rownames results = `pnames'
matrix colnames results = Coef LL_95 UL_95 p

display _newline as text "Pooled piecewise estimates, M = `M'  [`tag']"
matrix list results, format(%9.4f)

*--- The two terms reported in the Supplement are the slope-by-group
*    interactions: c.part2#1.kpscore (during the transition) and
*    c.part3#1.kpscore (after the transition). ---


*===========================================================================
* STAGE 4: POOLED TRAJECTORY FIGURE
*===========================================================================

mata: rubin_pool("mb", "mbV", `M', "Qcare", "SEcare")
mata: rubin_pool("mc", "mcV", `M', "Qctrl", "SEctrl")

local k = colsof(Qcare)

clear
set obs `k'
gen tcell = _n
gen prob_care = .
gen se_care   = .
gen prob_ctrl = .
gen se_ctrl   = .

forvalues j = 1/`k' {
    quietly replace prob_care = Qcare[1, `j']  in `j'
    quietly replace se_care   = SEcare[1, `j'] in `j'
    quietly replace prob_ctrl = Qctrl[1, `j']  in `j'
    quietly replace se_ctrl   = SEctrl[1, `j'] in `j'
}

foreach g in care ctrl {
    gen lo_`g' = prob_`g' - 1.96 * se_`g'
    gen hi_`g' = prob_`g' + 1.96 * se_`g'
}

*--- Map the margin cells onto the centred wave axis used in the main
*    figures. Waves are biennial for pa, alc and diet, annual for smok. ---
if "`outcome'" == "smok" {
    gen year = tcell - `shift'
}
else {
    gen year = 2 * tcell - 8
}

*--- Small horizontal offset so the error bars do not overlap ---
gen year_ctrl = year - 0.15
gen year_care = year + 0.15

twoway ///
    (rcap lo_ctrl hi_ctrl year_ctrl, pstyle(p1)) ///
    (rcap lo_care hi_care year_care, pstyle(p2)) ///
    (connected prob_ctrl year_ctrl, pstyle(p1)) ///
    (connected prob_care year_care, pstyle(p2)) ///
    , ///
    ytitle("`ytitle'") ///
    xtitle("`xtitle'") ///
    xline(0, lpattern(dash) lcolor(gs8)) ///
    legend(order(3 "`lab0'" 4 "`lab1'") pos(3) col(1)) ///
    xsize(6.5) ysize(4.5)

graph export "$out/`tag'_mi_m`M'.jpg", as(jpg) quality(90) replace


