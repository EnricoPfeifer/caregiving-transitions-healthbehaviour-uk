/*===========================================================================
  01_extract_wide.do
  
  PURPOSE: Extract variables from UKHLS waves and build a wide-format dataset
           for a single outcome. Run this file once per outcome by setting
           the `outcome` local macro below.
  
  OUTCOMES SUPPORTED:
    pa    -- physical activity (waves g i k m)
    alc   -- alcohol consumption (waves g i k m)
    diet  -- fruit and vegetable intake (waves g i k m)
    smok  -- smoking (waves e f g h i j k l m)
  
  OUTPUT: [outcome]_wide.dta saved to $ukhls
  
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
cd "/Users/Enrico/Documents/UCL/UBEL-DTP/Data"

global ukhls   "/Users/Enrico/Documents/UCL/UBEL-DTP/Data/UKHLS/stata/stata13_se/ukhls/"
global outputpath "/Users/Enrico/Documents/UCL/UBEL-DTP/Data"

clear all
set more off


*---------------------------------------------------------------------------
* WAVE SELECTION
* smok uses all waves from e onwards; all other outcomes use 4 waves only
*---------------------------------------------------------------------------
if "`outcome'" == "smok" {
    local allWaves "a b c d e f g h i j k l m"
    local keepWaves "e f g h i j k l m"   // waves with smoking + caregiving
}
else {
    local allWaves "a b c d e f g h i j k l m"
    local keepWaves "g i k m"             // waves with outcome + caregiving
}


*---------------------------------------------------------------------------
* SHARED VARIABLE LISTS
* indvars differs by outcome (outcome-specific questions); all else shared
*---------------------------------------------------------------------------

local indallvars "age_dv country ethn_dv gor_dv hhsize hidp mastat_dv nchild_dv pidp pno psnen01_lw psnen01_xw psnen91_lw psnen91_xw psnen99_xw psnenbh_lw psnenbh_xw psnenub_lw psnenub_xw psnenui_lw psnenui_xw psnenus_lw psnenus_xw psu racel_dv sex_dv strata urban_dv"

local hhvars "country fihhmnnet1_dv gor_dv hhden01_xw hhden91_xw hhden99_xw hhdenbh_xw hhdenub_xw hhdenui_xw hhdenus_xw hhsize hhtype_dv hidp ieqmoecd_dv nkids_dv psu strata tenure_dv urban_dv"

local chvars "age_dv country gor_dv hidp pidp pno psnen01_lw psnen91_lw psnenbh_lw psnenbh_xw psnenub_lw psnenub_xw psnenui_lw psnenui_xw psnenus_lw psnenus_xw psu sex_dv strata urban_dv"

local youthvars "age_dv country ethn_dv gor_dv hidp pidp pno psu racel_dv sex_dv strata urban_dv ythscbh_xw ythscub_xw ythscui_xw ythscus_xw"

* Shared core indresp variables (present in all outcomes)
local indvars_core "age_dv aidhh aidhrs aidxhh country ethn_dv gor_dv hhsize hhtype_dv hidp ind5mus_lw ind5mus_xw indbd91_lw indbdub_lw indin01_lw indin01_xw indin91_lw indin91_xw indin99_lw indinbh_xw indinub_lw indinub_xw indinui_lw indinui_xw indinus_lw indinus_xw indns91_lw indnsub_lw indpxbh_xw indpxub_lw indpxub_xw indpxui_lw indpxui_xw indpxus_lw indpxus_xw indscbh_xw indscub_lw indscub_xw indscui_lw indscui_xw indscus_lw indscus_xw jbstat mastat_dv nchild_dv pidp pno psu racel_dv sex_dv smoker strata tenure_dv ukborn urban_dv jbnssec3_dv jbft_dv employ jbhrs scghq1_dv sf12pcs_dv scsf1 hiqual_dv"

* Outcome-specific indresp variables appended to the core list
if "`outcome'" == "pa" {
    local indvars_extra "vday vdhrs vdmin vwhrs vwmin mday mdhrs mdmin mwhrs mwmin wday wdhrs wdmin wwhrs wwmin"
}
else if "`outcome'" == "alc" {
    local indvars_extra "auditc1 auditc2 auditc3 auditc4 auditc5 scfalcdrnk sceverdrnk"
}
else if "`outcome'" == "diet" {
    local indvars_extra "wkfruit fruitamt wkvege vegeamt"
}
else if "`outcome'" == "smok" {
    local indvars_extra "smever smnow ncigs"
}

local indvars "`indvars_core' `indvars_extra'"
local outputfilename "UKHLS_wide_abcdefghijkl"


*===========================================================================
* DO NOT EDIT BELOW THIS LINE
* Wave extraction and merging logic
*===========================================================================

* Program: prefix variable names with wave letter
program define getVars, rclass
    version 14.0
    if ("`1'" != "") {
        local wavemyvars = " `1'"
        local wavemyvars = subinstr("`wavemyvars'"," "," `2'_",.)
        local wavemyvars = substr("`wavemyvars'",2,.)
    }
    else local wavemyvars = ""
    return local fixedVars "`wavemyvars'"
end

* Program: return only those variables that exist in the current dataset
program define getExistingVars, rclass
    version 14.0
    local all = ""
    foreach var in `1' {
        capture confirm variable `var'
        if !_rc {
            local all = "`all' `var'"
        }
    }
    return local existingVars "`all'"
end


*---------------------------------------------------------------------------
* LOOP THROUGH ALL WAVES: extract and save one temp file per wave
*---------------------------------------------------------------------------
foreach wave in `allWaves' {
    local waveno = strpos("abcdefghijklmnopqrstuvwxyz","`wave'")

    getVars "`hhvars'"      `wave' ; local wavehhvars    = "`r(fixedVars)'"
    getVars "`indvars'"     `wave' ; local waveindvars   = "`r(fixedVars)'"
    getVars "`indallvars'"  `wave' ; local waveindallvars = "`r(fixedVars)'"
    getVars "`chvars'"      `wave' ; local wavechvars    = "`r(fixedVars)'"
    getVars "`youthvars'"   `wave' ; local waveyouthvars = "`r(fixedVars)'"

    use "$ukhls/`wave'_hhresp", clear
    getExistingVars "`wave'_hidp `wavehhvars'"
    keep `r(existingVars)'

    if ("`indvars'" != "" | "`chvars'" != "" | "`youthvars'" != "") {
        merge 1:m `wave'_hidp using "$ukhls/`wave'_indall"
        drop _merge
        drop if (pidp == .)
        getExistingVars "pidp `wave'_hidp `wavehhvars' `waveindallvars'"
        keep `r(existingVars)'

        if ("`indvars'" != "") {
            merge 1:1 pidp using "$ukhls/`wave'_indresp"
            drop _merge
            getExistingVars "pidp `wave'_hidp `wavehhvars' `waveindvars' `waveyouthvars' `wavechvars' `waveindallvars'"
            keep `r(existingVars)'
        }
        if ("`waveyouthvars'" != "") {
            merge 1:1 pidp using "$ukhls/`wave'_youth"
            drop _merge
            getExistingVars "pidp `wave'_hidp `wavehhvars' `waveindvars' `waveyouthvars' `wavechvars' `waveindallvars'"
            keep `r(existingVars)'
        }
        if ("`wavechvars'" != "") {
            merge 1:1 pidp using "$ukhls/`wave'_child"
            drop _merge
            getExistingVars "pidp `wave'_hidp `wavehhvars' `waveindvars' `waveyouthvars' `wavechvars' `waveindallvars'"
            keep `r(existingVars)'
        }
    }

    save temp_`wave', replace
}


*---------------------------------------------------------------------------
* MERGE ALL WAVE FILES INTO ONE WIDE DATASET
*---------------------------------------------------------------------------
local firstWave = substr("`allWaves'", 1, 1)
local firstvar  = subinstr("`joinvars'","#","`firstWave'",.)
use pidp `firstvar' using "temp_`firstWave'", clear
sort pidp
save "$outputpath/`outputfilename'", replace

foreach wave in `allWaves' {
    use * using "temp_`wave'", clear
    sort pidp
    merge 1:1 pidp using "$outputpath/`outputfilename'"
    drop _merge
    sort pidp
    order pidp, first
    save "$outputpath/`outputfilename'", replace
}

* Remove temporary wave files
foreach w in `allWaves' {
    erase temp_`w'.dta
}


*---------------------------------------------------------------------------
* DROP WAVES NOT NEEDED FOR THIS OUTCOME
* Keep only the analysis waves; drop all earlier waves to reduce file size
*---------------------------------------------------------------------------
use "$outputpath/`outputfilename'", clear

if "`outcome'" == "smok" {
    * Keep waves e-m; drop a-d
    drop a_* b_* c_* d_*
}
else {
    * Keep waves g, i, k, m; drop all others
    drop a_* b_* c_* d_* e_* f_* h_* j_* l_*
}

save "$ukhls/`outcome'_wide", replace

di as result "===== `outcome'_wide.dta saved successfully ====="
