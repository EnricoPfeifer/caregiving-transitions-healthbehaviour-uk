/*===========================================================================
  02_clean.do
  
  PURPOSE: Clean the wide-format UKHLS dataset for a single outcome.
           Constructs caregiving variables, sociodemographic covariates,
           and the outcome-specific variables, then saves an analysis-ready
           wide dataset.
  
  PREREQUISITE: 01_extract_wide.do must have been run for this outcome.
  
  OUTCOMES SUPPORTED:
    pa    -- physical inactivity (binary; waves g i k m)
    alc   -- problematic alcohol consumption (binary; waves g i k m)
    diet  -- daily fruit and vegetable portions (continuous; waves g i k m)
    smok  -- smoking status (binary; waves e f g h i j k l m)
  
  OUTPUT: [outcome]_wide_cleaned.dta saved to $ukhls
  
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

cd "[insert_own_filepath]"
global ukhls "[insert_own_filepath]"


*---------------------------------------------------------------------------
* WAVE LOOP MACRO
* All loops that run across caregiving and covariate waves use this macro
*---------------------------------------------------------------------------
if "`outcome'" == "smok" {
    local waves "e f g h i j k l m"
    local owaves "e f g h i j k l m"   // waves with outcome variable
}
else {
    local waves "g i k m"
    local owaves "g i k m"
}


*===========================================================================
* LOAD DATA
*===========================================================================
use "$ukhls/`outcome'_wide", clear

label define bival 0 "No" 1 "Yes"

*--- Sex (constructed early as it is needed for outcome-specific cut-offs) ---
egen sex = rowmax(*_sex_dv)
recode sex (min/-1=.) (1=0) (2=1)
label define sexval 0 "Men" 1 "Women"
label val sex sexval


*===========================================================================
* SECTION 1: CAREGIVING VARIABLES
* Constructed first as they are needed for first-observation flags below
*===========================================================================

*--- Binary caregiving (inside and outside household) ---
foreach w in `waves' {
    recode `w'_aidhh  (-9/-1=.), gen(`w'_carehh)
    recode `w'_aidxhh (-10/-1=.), gen(`w'_careohh)

    gen `w'_carebi = 1 if `w'_carehh == 1 | `w'_careohh == 1
    replace `w'_carebi = 2 if `w'_carehh == 2 & `w'_careohh == 2
    replace `w'_carebi = 2 if `w'_carehh == 2 & `w'_careohh == .
    replace `w'_carebi = 2 if `w'_carehh == . & `w'_careohh == 2

    recode `w'_carebi 2=0
    recode `w'_careohh 2=0
    recode `w'_carehh 2=0
    replace `w'_carehh = 0 if `w'_hhsize == 1  // single-person HH cannot be HH carer

    label val `w'_carebi  bival
    label val `w'_carehh  bival
    label val `w'_careohh bival
}

*--- Caregiving location category ---
label define carecatval  ///
    0 "No caregiving"           ///
    1 "Non-residential"         ///
    2 "Household"               ///
    3 "HH + non-residential"

foreach w in `waves' {
    gen `w'_carecat = 0 if `w'_carebi != . | `w'_carehh != . | `w'_careohh != .
    replace `w'_carecat = 1 if `w'_careohh == 1
    replace `w'_carecat = 2 if `w'_carehh  == 1
    replace `w'_carecat = 3 if `w'_carehh  == 1 & `w'_careohh == 1
    label val `w'_carecat carecatval
    label var `w'_carecat "Caregiving location category"
}

*--- Caregiving hours (binary: <20 / 20+) ---
label define carehval ///
    0 "Not a caregiver"    ///
    1 "Less than 20 hours" ///
    2 "20 hours or more"

foreach w in `waves' {
    gen `w'_carehrs = `w'_aidhrs
    recode `w'_carehrs (-10/-1=.) (1/3=1) (4/7=2) (8=1) (9=2) (97=.)
    replace `w'_carehrs = 0 if `w'_carebi == 0
    label val `w'_carehrs carehval
    label var `w'_carehrs "Caregiving hours category"
}

*--- Restrict caregiving to adults aged 16+ ---
foreach w in `waves' {
    replace `w'_carehh  = . if `w'_age_dv < 16
    replace `w'_careohh = . if `w'_age_dv < 16
    replace `w'_carebi  = . if `w'_age_dv < 16
}


*===========================================================================
* SECTION 2: FIRST-OBSERVATION FLAG VARIABLES
* These derive baseline and onset information from the panel
*===========================================================================

*--- Ever a caregiver ---
egen ticarehh  = anycount(*_carehh),  values(1)
egen ticareohh = anycount(*_careohh), values(1)
egen ticarebi  = anycount(*_carebi),  values(1)

egen rowmisscare = rowmiss(*_carebi)

if "`outcome'" == "smok" {
    recode ticarebi (1/9=1), gen(ecare)
    replace ecare = . if rowmisscare == 9
}
else {
    recode ticarebi (1/4=1), gen(ecare)
    replace ecare = . if rowmisscare == 4
}
label var ecare "Caregiver in any observed wave"
label define ecareval 0 "No care at any wave" 1 "Care at at least one wave"
label val ecare ecareval

*--- Number of observed waves ---
egen nwaves = rownonmiss(*_carebi)

if "`outcome'" == "smok" {
    drop if nwaves == 0
}

*--- Wave when first observed caring ---
if "`outcome'" == "smok" {
    gen ficare = 1 if e_carebi == 1
    foreach w in f g h i j k l m {
        local n = strpos("efghijklm", "`w'")
        replace ficare = `n' if `w'_carebi == 1 & ficare == .
    }
    gen fiob = 1 if e_carebi != .
    foreach w in f g h i j k l m {
        local n = strpos("efghijklm", "`w'")
        replace fiob = `n' if `w'_carebi != . & fiob == .
    }
}
else {
    gen ficare = 1 if g_carebi == 1
    replace ficare = 2 if i_carebi == 1 & ficare == .
    replace ficare = 3 if k_carebi == 1 & ficare == .
    replace ficare = 4 if m_carebi == 1 & ficare == .
    gen fiob = 1 if g_carebi != .
    replace fiob = 2 if i_carebi != . & fiob == .
    replace fiob = 3 if k_carebi != . & fiob == .
    replace fiob = 4 if m_carebi != . & fiob == .
}

*--- Recode age to missing if below 16 ---
foreach w in `waves' {
    recode `w'_age_dv (min/15=.)
}

*--- Age at first observation ---
if "`outcome'" == "smok" {
    gen agefiob = e_age_dv if fiob == 1
    replace agefiob = f_age_dv if fiob == 2
    replace agefiob = g_age_dv if fiob == 3
    replace agefiob = h_age_dv if fiob == 4
    replace agefiob = i_age_dv if fiob == 5
    replace agefiob = j_age_dv if fiob == 6
    replace agefiob = k_age_dv if fiob == 7
    replace agefiob = l_age_dv if fiob == 8
    replace agefiob = m_age_dv if fiob == 9
    gen ageonset = .
    replace ageonset = e_age_dv if ficare == 1
    replace ageonset = f_age_dv if ficare == 2
    replace ageonset = g_age_dv if ficare == 3
    replace ageonset = h_age_dv if ficare == 4
    replace ageonset = i_age_dv if ficare == 5
    replace ageonset = j_age_dv if ficare == 6
    replace ageonset = k_age_dv if ficare == 7
    replace ageonset = l_age_dv if ficare == 8
    replace ageonset = m_age_dv if ficare == 9
}
else {
    gen agefiob = g_age_dv if fiob == 1
    replace agefiob = i_age_dv if fiob == 2
    replace agefiob = k_age_dv if fiob == 3
    replace agefiob = m_age_dv if fiob == 4
    gen ageonset = .
    replace ageonset = g_age_dv if ficare == 1
    replace ageonset = i_age_dv if ficare == 2
    replace ageonset = k_age_dv if ficare == 3
    replace ageonset = m_age_dv if ficare == 4
}
replace ageonset = 0 if ecare == 0

gen agegr_fiob = agefiob
recode agegr_fiob (min/15=.) (16/29=0) (30/49=1) (50/64=2) (65/max=3)
label define agegrval 0 "16-29" 1 "30-49" 2 "50-64" 3 "65+"
label val agegr_fiob agegrval

egen agemin = rowmin(*_age_dv)
egen agemax = rowmax(*_age_dv)

*--- Caregiving hours at onset ---
if "`outcome'" == "smok" {
    gen carehrs_ficare = e_carehrs if ficare == 1
    replace carehrs_ficare = f_carehrs if ficare == 2
    replace carehrs_ficare = g_carehrs if ficare == 3
    replace carehrs_ficare = h_carehrs if ficare == 4
    replace carehrs_ficare = i_carehrs if ficare == 5
    replace carehrs_ficare = j_carehrs if ficare == 6
    replace carehrs_ficare = k_carehrs if ficare == 7
    replace carehrs_ficare = l_carehrs if ficare == 8
    replace carehrs_ficare = m_carehrs if ficare == 9
}
else {
    gen carehrs_ficare = g_carehrs if ficare == 1
    replace carehrs_ficare = i_carehrs if ficare == 2
    replace carehrs_ficare = k_carehrs if ficare == 3
    replace carehrs_ficare = m_carehrs if ficare == 4
}

*--- Caregiving hours summary (midpoints) ---
foreach w in `waves' {
    recode `w'_aidhrs (-10/-1=.)
    gen `w'_hours = .
    replace `w'_hours = 2    if `w'_aidhrs == 1 & `w'_carebi != .
    replace `w'_hours = 7    if `w'_aidhrs == 2 & `w'_carebi != .
    replace `w'_hours = 14.5 if `w'_aidhrs == 3 & `w'_carebi != .
    replace `w'_hours = 27   if `w'_aidhrs == 4 & `w'_carebi != .
    replace `w'_hours = 42   if `w'_aidhrs == 5 & `w'_carebi != .
    replace `w'_hours = 74.5 if `w'_aidhrs == 6 & `w'_carebi != .
    replace `w'_hours = 100  if `w'_aidhrs == 7 & `w'_carebi != .
    replace `w'_hours = 10   if `w'_aidhrs == 8 & `w'_carebi != .
    replace `w'_hours = 34   if `w'_aidhrs == 9 & `w'_carebi != .
    replace `w'_hours = .    if `w'_aidhrs == 97
}
egen meanhour = rowmean(*_hours)
gen meanhour_g = 1 if meanhour < 5
replace meanhour_g = 2 if meanhour < 10  & meanhour_g == .
replace meanhour_g = 3 if meanhour < 20  & meanhour_g == .
replace meanhour_g = 4 if meanhour < 35  & meanhour_g == .
replace meanhour_g = 5 if meanhour < 50  & meanhour_g == .
replace meanhour_g = 6 if meanhour < 100 & meanhour_g == .
replace meanhour_g = 7 if meanhour == 100 & meanhour_g == .
replace meanhour_g = 0 if ecare == 0
recode meanhour_g (5/7=5)

*--- Caregiver at first observation (baseline caregiving status) ---
if "`outcome'" == "smok" {
    gen care_fiob = e_carebi if fiob == 1
    replace care_fiob = f_carebi if fiob == 2
    replace care_fiob = g_carebi if fiob == 3
    replace care_fiob = h_carebi if fiob == 4
    replace care_fiob = i_carebi if fiob == 5
    replace care_fiob = j_carebi if fiob == 6
    replace care_fiob = k_carebi if fiob == 7
    replace care_fiob = l_carebi if fiob == 8
    replace care_fiob = m_carebi if fiob == 9
}
else {
    gen care_fiob = g_carebi if fiob == 1
    replace care_fiob = i_carebi if fiob == 2
    replace care_fiob = k_carebi if fiob == 3
    replace care_fiob = m_carebi if fiob == 4
}
label var care_fiob "Caregiver at first observation"
label define carefival 0 "Non-caregiver at baseline" 1 "Caregiver at baseline"
label val care_fiob carefival

*--- Transition INTO caregiving ---
gen caretrans = .
replace caretrans = 1 if care_fiob == 0 & ecare == 1
replace caretrans = 0 if ecare == 0
label define transval 0 "Never transitioned" 1 "Transitioned INTO care"
label val caretrans transval

gen transhrs = .
replace transhrs = 0 if ecare == 0
replace transhrs = 1 if carehrs_ficare == 1
replace transhrs = 2 if carehrs_ficare == 2
label define transhrval 0 "No transition" 1 "INTO <20 hrs" 2 "INTO 20+ hrs"
label val transhrs transhrval

*--- Always/never caregiver summary ---
if "`outcome'" == "smok" {
    egen alwayscarer = anycount(*_carebi), values(0)
    recode alwayscarer (0=1) (1/9=0)
    egen nevercarer = anycount(*_carebi), values(1)
    recode nevercarer (0=1) (1/9=0)
}
else {
    egen alwayscarer = anycount(*_carebi), values(0)
    recode alwayscarer (0=1) (1/4=0)
    egen nevercarer = anycount(*_carebi), values(1)
    recode nevercarer (0=1) (1/4=0)
}
label define alwayscareval 0 "Non-carer in at least one wave" 1 "Always carer"
label define nevercareval  0 "Carer in at least one wave"    1 "Never carer"
label val alwayscarer alwayscareval
label val nevercarer  nevercareval

gen caresummary = .
replace caresummary = 0 if nevercarer == 1
replace caresummary = 1 if alwayscarer == 0 | nevercarer == 0
replace caresummary = 2 if alwayscarer == 1
label define caresumval 0 "Never caregiver" 1 "Transitioned" 2 "Always caregiver"
label val caresummary caresumval

gen caresum2 = caresummary
recode caresum2 2=3
replace caresum2 = 1 if care_fiob == 0 & nevercarer == 0
replace caresum2 = 2 if care_fiob == 1 & alwayscarer == 0
label define caresum2val 0 "Never" 1 "Trans INTO" 2 "EXIT care" 3 "Always"
label val caresum2 caresum2val


*===========================================================================
* SECTION 3: OUTCOME-SPECIFIC VARIABLES
*===========================================================================

*---------------------------------------------------------------------------
* PHYSICAL ACTIVITY (pa)
*---------------------------------------------------------------------------
if "`outcome'" == "pa" {

    * Vigorous PA
    foreach w in g i k m {
        recode `w'_vday  (-10/-1=.), gen(`w'_vigday)
        recode `w'_vdhrs (-10/-1=.), gen(`w'_vighrs)
        replace `w'_vighrs = 0 if `w'_vday == 0
        recode `w'_vdmin  (-10/-1=.), gen(`w'_vigmin)
        replace `w'_vigmin = 0 if `w'_vday == 0

        gen  `w'_vighm1a = `w'_vigday * `w'_vighrs * 60
        gen  `w'_vighm1b = `w'_vigday * `w'_vigmin
        egen `w'_vighm1  = rowtotal(`w'_vighm1a `w'_vighm1b)
        replace `w'_vighm1 = . if `w'_vighm1a == . & `w'_vighm1b == .

        gen  `w'_vigwh   = `w'_vwhrs
        recode `w'_vigwh -9/-1=.
        gen  `w'_vigwhm  = `w'_vigwh * 60
        gen  `w'_vigwm   = `w'_vwmin
        recode `w'_vigwm -9/-1=.
        egen `w'_vighm2  = rowtotal(`w'_vigwhm `w'_vigwm)
        replace `w'_vighm2 = . if `w'_vigwhm == . & `w'_vigwm == .

        egen `w'_vigtot = rowtotal(`w'_vighm1 `w'_vighm2)
        replace `w'_vigtot = . if `w'_vighm1 == . & `w'_vighm2 == .
        replace `w'_vigtot = `w'_vighm2 if `w'_vighm1 == .

        gen `w'_vigbi = `w'_vigtot
        recode `w'_vigbi (0/74=1) (75/max=0)
        label val `w'_vigbi bival
        label var `w'_vigbi "Met 75 min vigorous PA threshold"
    }

    * Moderate PA
    foreach w in g i k m {
        recode `w'_mday  (-10/-1=.), gen(`w'_moday)
        recode `w'_mdhrs (-10/-1=.), gen(`w'_mohrs)
        replace `w'_mohrs = 0 if `w'_moday == 0
        recode `w'_mdmin  (-10/-1=.), gen(`w'_momin)
        replace `w'_momin = 0 if `w'_moday == 0

        gen  `w'_mohm1a = `w'_moday * `w'_mohrs * 60
        gen  `w'_mohm1b = `w'_moday * `w'_momin
        egen `w'_mohm1  = rowtotal(`w'_mohm1a `w'_mohm1b)
        replace `w'_mohm1 = . if `w'_mohm1a == . & `w'_mohm1b == .

        gen  `w'_mowh  = `w'_mwhrs
        recode `w'_mowh -9/-1=.
        gen  `w'_mowhm = `w'_mowh * 60
        gen  `w'_mowm  = `w'_mdmin
        recode `w'_mowm -9/-1=.
        egen `w'_mohm2 = rowtotal(`w'_mowhm `w'_mowm)
        replace `w'_mohm2 = . if `w'_mowhm == . & `w'_mowm == .

        egen `w'_motot = rowtotal(`w'_mohm1 `w'_mohm2)
        replace `w'_motot = . if `w'_mohm1 == . & `w'_mohm2 == .

        gen `w'_mobi = `w'_motot
        recode `w'_mobi (0/149=1) (150/max=0)
        label val `w'_mobi bival
        label var `w'_mobi "Met 150 min moderate PA threshold"
    }

    * Combined PA outcome
    foreach w in g i k m {
        egen `w'_patot = rowtotal(`w'_vigtot `w'_motot)
        replace `w'_patot = . if `w'_vigtot == . & `w'_motot == .
        label var `w'_patot "Minutes moderate + vigorous PA per week"

        gen `w'_pabi = `w'_patot
        recode `w'_pabi (min/149=1) (150/max=0)
        replace `w'_pabi = 0 if `w'_vigbi == 0 | `w'_mobi == 0
        replace `w'_pabi = . if `w'_vigbi == . & `w'_mobi == .
        label val `w'_pabi bival
        label var `w'_pabi "Physically inactive"
    }
}


*---------------------------------------------------------------------------
* ALCOHOL (alc)
*---------------------------------------------------------------------------
if "`outcome'" == "alc" {

    label define auditcval1 ///
        0 "Never"             ///
        1 "Monthly or less"   ///
        2 "2-4 times a month" ///
        3 "2-3 times a week"  ///
        4 "4+ times per week"

    label define auditcval2 ///
        0 "1-2"   ///
        1 "3-4"   ///
        2 "5-6"   ///
        3 "7-9"   ///
        4 "10+"

    label define auditcval3 ///
        0 "Never"                          ///
        1 "Less than monthly"              ///
        2 "Monthly"                        ///
        3 "Weekly"                         ///
        4 "Daily or almost daily"

    foreach w in g i k m {
        recode `w'_auditc1 (-10/-1=.), gen(`w'_drinkyr)   // drank past 12 months
        recode `w'_auditc2 (-10/-1=.), gen(`w'_nodrink)   // always non-drinker
        recode `w'_auditc3 (-10/-1=.), gen(`w'_alcfre)    // drinking frequency
        recode `w'_auditc4 (-10/-1=.), gen(`w'_drinkd)    // drinks per occasion
        recode `w'_auditc5 (-10/-1=.), gen(`w'_binge)     // 6+ drinks frequency

        recode `w'_drinkyr (2=0)
        label var `w'_drinkyr "Drank in past 12 months"
        label val `w'_drinkyr bival

        recode `w'_alcfre (1=0) (2=1) (3=2) (4=3) (5=4)
        replace `w'_alcfre = 0 if `w'_drinkyr == 0
        label var `w'_alcfre "Alcohol frequency past 12 months"
        label val `w'_alcfre auditcval1

        recode `w'_drinkd (1=0) (2=1) (3=2) (4=3) (5=4)
        replace `w'_drinkd = 0 if `w'_drinkyr == 0
        label var `w'_drinkd "Standard drinks per drinking occasion"
        label val `w'_drinkd auditcval2

        recode `w'_binge (1=0) (2=1) (3=2) (4=3) (5=4)
        replace `w'_binge = 0 if `w'_drinkyr == 0
        label var `w'_binge "Frequency of 6+ drinks on one occasion"
        label val `w'_binge auditcval3

        * AUDIT-C total score (0-12)
        gen `w'_auditc = `w'_alcfre + `w'_drinkd + `w'_binge
        replace `w'_auditc = 0 if `w'_drinkyr == 0
        replace `w'_auditc = . if `w'_alcfre == . & `w'_drinkd == . & `w'_binge == .
        label var `w'_auditc "AUDIT-C score (0-12)"

        * Problematic drinking: sex-specific AUDIT-C thresholds
        * >= 3 for women, >= 4 for men
        * Note: sex must already exist at this point; it is constructed above
        gen `w'_probdrink = . 
        replace `w'_probdrink = (`w'_auditc >= 3) if `w'_auditc != . & sex == 1
        replace `w'_probdrink = (`w'_auditc >= 4) if `w'_auditc != . & sex == 0
        label val `w'_probdrink bival
        label var `w'_probdrink "Problematic drinking (AUDIT-C >=3 women, >=4 men)"
    }
}


*---------------------------------------------------------------------------
* DIET (diet)
*---------------------------------------------------------------------------
if "`outcome'" == "diet" {

    foreach w in g i k m {
        * Fruit
        recode `w'_wkfruit (-9/-1=.), gen(`w'_weekfruit)
        recode `w'_weekfruit (1=0) (2=2) (3=5) (4=7)
        label var `w'_weekfruit "Days per week eating fruit"

        recode `w'_fruitamt (-9/-1=.), gen(`w'_portionfruit)
        replace `w'_portionfruit = 0 if `w'_weekfruit == 0
        label var `w'_portionfruit "Portions of fruit per day"

        gen `w'_meanfruit = (`w'_weekfruit * `w'_portionfruit) / 7
        label var `w'_meanfruit "Mean daily portions of fruit"

        * Vegetables
        recode `w'_wkvege (-9/-1=.), gen(`w'_weekvege)
        recode `w'_weekvege (1=0) (2=2) (3=5) (4=7)
        label var `w'_weekvege "Days per week eating vegetables"

        recode `w'_vegeamt (-9/-1=.), gen(`w'_portionvege)
        replace `w'_portionvege = 0 if `w'_weekvege == 0
        label var `w'_portionvege "Portions of vegetables per day"

        gen `w'_meanvege = (`w'_weekvege * `w'_portionvege) / 7
        label var `w'_meanvege "Mean daily portions of vegetables"

        * Combined fruit and vegetable outcome
        gen `w'_fvtot = `w'_meanfruit + `w'_meanvege
        replace `w'_fvtot = . if `w'_meanfruit == . & `w'_meanvege == .
        label var `w'_fvtot "Mean daily portions of fruit and vegetables"
    }
}


*---------------------------------------------------------------------------
* SMOKING (smok)
*---------------------------------------------------------------------------
if "`outcome'" == "smok" {

    * Wave e uses smnow/smever (different question set from wave f onwards)
    gen e_smok = e_smnow
    recode e_smok (-8/-1=.) (2=0)
    replace e_smok = 0 if e_smever == 2
    label val e_smok bival
    label var e_smok "Current smoker (wave e)"

    foreach w in f g h i j k l m {
        recode `w'_smoker (-10/-1=.), gen(`w'_smok)
        recode `w'_smok 2=0
        label val `w'_smok bival
        label var `w'_smok "Current smoker"
    }

    * Number of cigarettes per day
    foreach w in e f g h i j k l m {
        gen `w'_ciggy = `w'_ncigs
        recode `w'_ciggy (-10/-1=.)
        replace `w'_ciggy = 0 if `w'_smok == 0
        label var `w'_ciggy "Cigarettes per day"
    }

    * Smoking summary across waves
    egen smokwaves = rowtotal(*_smok)

    gen alwaysmok = smokwaves - nwaves
    recode alwaysmok (min/-1=0) (0=1) (1/max=0)
    label define asmokval 0 "Not always smoker" 1 "Smoker every observed wave"
    label val alwaysmok asmokval
    label var alwaysmok "Smoker at every wave"

    gen neversmok = smokwaves
    recode neversmok (0=1) (1/max=0)
    label define nsmokval 0 "Smoked at least once" 1 "Never smoked"
    label val neversmok nsmokval
}


*===========================================================================
* SECTION 4: COVARIATES
*===========================================================================

*--- Number of children in household (merged from indresp) ---
save "$ukhls/`outcome'_wide_cleaned", replace

foreach w in `waves' {
    use pidp `w'_nchild_dv using "$ukhls/`w'_indresp", clear
    sort pidp
    merge 1:1 pidp using "$ukhls/`outcome'_wide_cleaned"
    drop _merge
    sort pidp
    save "$ukhls/`outcome'_wide_cleaned", replace
}

label define childval 0 "0" 1 "1" 2 "2" 3 "3+"
foreach w in `waves' {
    gen `w'_children = `w'_nchild_dv
    recode `w'_children (min/-1=.) (3/max=3)
    label var `w'_children "Number of children in household"
    label val `w'_children childval
}

if "`outcome'" == "smok" {
    gen children_fiob = e_children if fiob == 1
    replace children_fiob = f_children if fiob == 2
    replace children_fiob = g_children if fiob == 3
    replace children_fiob = h_children if fiob == 4
    replace children_fiob = i_children if fiob == 5
    replace children_fiob = j_children if fiob == 6
    replace children_fiob = k_children if fiob == 7
    replace children_fiob = l_children if fiob == 8
    replace children_fiob = m_children if fiob == 9
}
else {
    gen children_fiob = g_children if fiob == 1
    replace children_fiob = i_children if fiob == 2
    replace children_fiob = k_children if fiob == 3
    replace children_fiob = m_children if fiob == 4
}
label val children_fiob childval

*--- Equivalised household income (quintiles) ---
foreach w in `waves' {
    recode `w'_ieqmoecd_dv   (-9/-1=.)
    recode `w'_fihhmnnet1_dv (min/-0.0001=.)
    gen `w'_ehhincome = `w'_fihhmnnet1_dv / `w'_ieqmoecd_dv
    xtile `w'_iwealth = `w'_ehhincome, n(5)
}

if "`outcome'" == "smok" {
    gen hhincome_fiob = e_ehhincome if fiob == 1
    replace hhincome_fiob = f_ehhincome if fiob == 2
    replace hhincome_fiob = g_ehhincome if fiob == 3
    replace hhincome_fiob = h_ehhincome if fiob == 4
    replace hhincome_fiob = i_ehhincome if fiob == 5
    replace hhincome_fiob = j_ehhincome if fiob == 6
    replace hhincome_fiob = k_ehhincome if fiob == 7
    replace hhincome_fiob = l_ehhincome if fiob == 8
    replace hhincome_fiob = m_ehhincome if fiob == 9
}
else {
    gen hhincome_fiob = g_ehhincome if fiob == 1
    replace hhincome_fiob = i_ehhincome if fiob == 2
    replace hhincome_fiob = k_ehhincome if fiob == 3
    replace hhincome_fiob = m_ehhincome if fiob == 4
}
xtile iwealth_fiob = hhincome_fiob, n(5)

*--- Employment status ---
label define jbstatusval ///
    2  "Employed"                 ///
    3  "Unemployed"               ///
    4  "Retired"                  ///
    6  "Family care or home"      ///
    7  "Full-time student"        ///
    8  "LT sick or disabled"      ///
    9  "Government training"      ///
    97 "Something else"

foreach w in `waves' {
    gen `w'_jbstatus = `w'_jbstat
    recode `w'_jbstatus (-9/-1=.) (11=9) (10=97) (5=2) (12/13=3) (1=2) (9=7)
    label val `w'_jbstatus jbstatusval
}

if "`outcome'" == "smok" {
    gen jbstatus_fiob = e_jbstat if fiob == 1
    replace jbstatus_fiob = f_jbstat if fiob == 2
    replace jbstatus_fiob = g_jbstat if fiob == 3
    replace jbstatus_fiob = h_jbstat if fiob == 4
    replace jbstatus_fiob = i_jbstat if fiob == 5
    replace jbstatus_fiob = j_jbstat if fiob == 6
    replace jbstatus_fiob = k_jbstat if fiob == 7
    replace jbstatus_fiob = l_jbstat if fiob == 8
    replace jbstatus_fiob = m_jbstat if fiob == 9
}
else {
    gen jbstatus_fiob = g_jbstat if fiob == 1
    replace jbstatus_fiob = i_jbstat if fiob == 2
    replace jbstatus_fiob = k_jbstat if fiob == 3
    replace jbstatus_fiob = m_jbstat if fiob == 4
}
label val jbstatus_fiob jbstatusval
recode jbstatus_fiob (-9/-1=.) (11=9) (10=97) (5=2) (12/13=3) (1=2) (9=7)

*--- Occupational class (NS-SEC 3-category) ---
label define oclassval ///
    0 "Not employed"              ///
    1 "Managerial & professional" ///
    2 "Intermediate"              ///
    3 "Routine"

foreach w in `waves' {
    gen `w'_oclass = `w'_jbnssec3_dv
    recode `w'_oclass (-9/-1=.)
    replace `w'_oclass = 0 if `w'_jbstatus == 3 | `w'_jbstatus == 4
    replace `w'_oclass = 0 if `w'_jbstatus >= 6 & `w'_jbstatus < .
    label val `w'_oclass oclassval
    label var `w'_oclass "Occupational class"
}

if "`outcome'" == "smok" {
    gen oclass_fiob = e_jbnssec3_dv if fiob == 1
    replace oclass_fiob = f_jbnssec3_dv if fiob == 2
    replace oclass_fiob = g_jbnssec3_dv if fiob == 3
    replace oclass_fiob = h_jbnssec3_dv if fiob == 4
    replace oclass_fiob = i_jbnssec3_dv if fiob == 5
    replace oclass_fiob = j_jbnssec3_dv if fiob == 6
    replace oclass_fiob = k_jbnssec3_dv if fiob == 7
    replace oclass_fiob = l_jbnssec3_dv if fiob == 8
    replace oclass_fiob = m_jbnssec3_dv if fiob == 9
}
else {
    gen oclass_fiob = g_jbnssec3_dv if fiob == 1
    replace oclass_fiob = i_jbnssec3_dv if fiob == 2
    replace oclass_fiob = k_jbnssec3_dv if fiob == 3
    replace oclass_fiob = m_jbnssec3_dv if fiob == 4
}
replace oclass_fiob = 0 if jbstatus_fiob == 3 | jbstatus_fiob == 4
replace oclass_fiob = 0 if jbstatus_fiob >= 6 & jbstatus_fiob < .
recode oclass_fiob (min/-1=.)
label val oclass_fiob oclassval

*--- Employment (full-time / part-time) ---
label define workval 0 "Not in paid employment" 1 "Full-time" 2 "Part-time"
foreach w in `waves' {
    gen `w'_workbi = `w'_employ
    recode `w'_workbi (-9/-1=.) (2=0)
    replace `w'_workbi = 1 if `w'_jbft_dv == 1
    replace `w'_workbi = 2 if `w'_jbft_dv == 2
    label val `w'_workbi workval
    label var `w'_workbi "Employment type"
}

*--- Cohabitation status ---
label define cohabval 0 "Single, separated, or widowed" 1 "Married or cohabiting"
label define marriedval 1 "Single" 2 "Married" 3 "Separated" 4 "Living as couple" 6 "Widowed"

foreach w in `waves' {
    gen `w'_married = `w'_mastat_dv
    recode `w'_married (0=1) (3=2) (4/5=3) (7/8=3) (9/10=4)
    label val `w'_married marriedval
    label var `w'_married "Marital status"

    gen `w'_cohab = `w'_married
    recode `w'_cohab (-9/-1=.) (1=0) (2=1) (3=0) (4=1) (6=0)
    label val `w'_cohab bival
    label var `w'_cohab "Married or cohabiting"
}

if "`outcome'" == "smok" {
    gen cohab_fiob = e_mastat_dv if fiob == 1
    replace cohab_fiob = f_mastat_dv if fiob == 2
    replace cohab_fiob = g_mastat_dv if fiob == 3
    replace cohab_fiob = h_mastat_dv if fiob == 4
    replace cohab_fiob = i_mastat_dv if fiob == 5
    replace cohab_fiob = j_mastat_dv if fiob == 6
    replace cohab_fiob = k_mastat_dv if fiob == 7
    replace cohab_fiob = l_mastat_dv if fiob == 8
    replace cohab_fiob = m_mastat_dv if fiob == 9
}
else {
    gen cohab_fiob = g_mastat_dv if fiob == 1
    replace cohab_fiob = i_mastat_dv if fiob == 2
    replace cohab_fiob = k_mastat_dv if fiob == 3
    replace cohab_fiob = m_mastat_dv if fiob == 4
}
recode cohab_fiob (min/-1=.) (0=1) (3=2) (4/5=3) (7/8=3) (9/10=4)
gen tmp_cohab = cohab_fiob
recode tmp_cohab (1=0) (2=1) (3=0) (4=1) (6=0)
drop cohab_fiob
rename tmp_cohab cohab_fiob
label val cohab_fiob cohabval

*--- Ethnicity (time-invariant) ---
egen ethin = rowmax(*_ethn_dv)
recode ethin 56=15
gen ethnicity = 1 if ethin >= 1 & ethin <= 4
replace ethnicity = 2 if inlist(ethin, 5, 6, 14, 15, 16)
replace ethnicity = 3 if ethin == 9
replace ethnicity = 4 if ethin == 10
replace ethnicity = 5 if ethin == 11
replace ethnicity = 6 if inlist(ethin, 7, 8, 12, 13, 17, 97)
recode ethnicity (5=4) (6=5)
label define ethnicityval ///
    1 "White"                   ///
    2 "Black"                   ///
    3 "Indian"                  ///
    4 "Pakistani/Bangladeshi"   ///
    5 "Other Asian/Other"
label val ethnicity ethnicityval

*--- Education (highest qualification; time-invariant) ---
egen edu = rowmax(*_hiqual_dv)
recode edu (-9/-1=.) (9=0) (3/5=1) (1/2=2)
label define eduval ///
    0 "No qualification"                       ///
    1 "A-Level, GCSE, or other qualification"  ///
    2 "Degree or higher"
label var edu "Highest qualification"
label val edu eduval

*--- GHQ psychological distress ---
foreach w in `waves' {
    gen `w'_ghq = `w'_scghq1_dv
    recode `w'_ghq (min/-1=.)
    label var `w'_ghq "GHQ psychological distress score"
}

if "`outcome'" == "smok" {
    gen ghq_fiob = e_scghq1_dv if fiob == 1
    replace ghq_fiob = f_scghq1_dv if fiob == 2
    replace ghq_fiob = g_scghq1_dv if fiob == 3
    replace ghq_fiob = h_scghq1_dv if fiob == 4
    replace ghq_fiob = i_scghq1_dv if fiob == 5
    replace ghq_fiob = j_scghq1_dv if fiob == 6
    replace ghq_fiob = k_scghq1_dv if fiob == 7
    replace ghq_fiob = l_scghq1_dv if fiob == 8
    replace ghq_fiob = m_scghq1_dv if fiob == 9
}
else {
    gen ghq_fiob = g_scghq1_dv if fiob == 1
    replace ghq_fiob = i_scghq1_dv if fiob == 2
    replace ghq_fiob = k_scghq1_dv if fiob == 3
    replace ghq_fiob = m_scghq1_dv if fiob == 4
}
recode ghq_fiob (min/-1=.)

*--- Self-rated general health ---
label define ghval 0 "Excellent, very good, or good" 1 "Fair or poor"
foreach w in `waves' {
    gen `w'_ghealth = `w'_scsf1
    recode `w'_ghealth (-9/-1=.) (1/3=0) (4/5=1)
    label val `w'_ghealth ghval
    label var `w'_ghealth "Self-rated general health"
}

if "`outcome'" == "smok" {
    gen ghealth_fiob = e_ghealth if fiob == 1
    replace ghealth_fiob = f_ghealth if fiob == 2
    replace ghealth_fiob = g_ghealth if fiob == 3
    replace ghealth_fiob = h_ghealth if fiob == 4
    replace ghealth_fiob = i_ghealth if fiob == 5
    replace ghealth_fiob = j_ghealth if fiob == 6
    replace ghealth_fiob = k_ghealth if fiob == 7
    replace ghealth_fiob = l_ghealth if fiob == 8
    replace ghealth_fiob = m_ghealth if fiob == 9
}
else {
    gen ghealth_fiob = g_ghealth if fiob == 1
    replace ghealth_fiob = i_ghealth if fiob == 2
    replace ghealth_fiob = k_ghealth if fiob == 3
    replace ghealth_fiob = m_ghealth if fiob == 4
}

*--- Household size ---
if "`outcome'" == "smok" {
    gen hhsize_fiob = e_hhsize if fiob == 1
    replace hhsize_fiob = f_hhsize if fiob == 2
    replace hhsize_fiob = g_hhsize if fiob == 3
    replace hhsize_fiob = h_hhsize if fiob == 4
    replace hhsize_fiob = i_hhsize if fiob == 5
    replace hhsize_fiob = j_hhsize if fiob == 6
    replace hhsize_fiob = k_hhsize if fiob == 7
    replace hhsize_fiob = l_hhsize if fiob == 8
    replace hhsize_fiob = m_hhsize if fiob == 9
}
else {
    gen hhsize_fiob = g_hhsize if fiob == 1
    replace hhsize_fiob = i_hhsize if fiob == 2
    replace hhsize_fiob = k_hhsize if fiob == 3
    replace hhsize_fiob = m_hhsize if fiob == 4
}

label define hhgroupval 1 "1" 2 "2" 3 "3-4" 4 "5 or more"
gen hhgroup_fiob = hhsize_fiob
recode hhgroup_fiob (3/4=3) (5/max=4)
label val hhgroup_fiob hhgroupval

foreach w in `waves' {
    gen `w'_hhgroup = `w'_hhsize
    recode `w'_hhgroup (-9/-1=.) (1=0) (2=1) (3/4=2) (5/max=3)
    label var `w'_hhgroup "Household size group"
}

*--- Physical health (SF-12 physical component) ---
foreach w in `waves' {
    recode `w'_sf12pcs_dv (-10/-1=.), gen(`w'_sf12p)
    label var `w'_sf12p "SF-12 physical health score"
}

if "`outcome'" != "smok" {
    gen sf_fiob = g_sf12p if fiob == 1
    replace sf_fiob = i_sf12p if fiob == 2
    replace sf_fiob = k_sf12p if fiob == 3
    replace sf_fiob = m_sf12p if fiob == 4
}


*===========================================================================
* SAVE
*===========================================================================
save "$ukhls/`outcome'_wide_cleaned", replace

di as result "===== `outcome'_wide_cleaned.dta saved successfully ====="
