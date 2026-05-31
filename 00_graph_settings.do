/*===========================================================================
  00_graph_settings.do

  PURPOSE: Sets a consistent publication-ready graph style using grstyle.
           Run this file once at the start of any session that produces
           figures, before calling mplotoffset or any graph commands.

  REQUIRED PACKAGES (install once if not already present):
    ssc install grstyle, replace
    ssc install palettes, replace
    ssc install colrspace, replace

  SOURCE: Style approach adapted from:
    https://www.csae.ox.ac.uk/files/coderscornerttweek8fmpdf

  AUTHOR: Enrico Pfeifer
===========================================================================*/


*===========================================================================
* INITIALISE
*===========================================================================

set scheme s2color   // start from the default s2color base scheme
grstyle init         // initialise grstyle


*===========================================================================
* LAYOUT AND ORIENTATION
*===========================================================================

grstyle set horizontal        // horizontal tick labels
grstyle set compact           // compact overall layout
grstyle set size small:  subheading axis_title
grstyle set size vsmall: small_body


*===========================================================================
* LEGEND
*===========================================================================

grstyle set legend 6, nobox   // legend at 6 o'clock; no box


*===========================================================================
* LINE WIDTHS
*===========================================================================

grstyle set linewidth medthick: p         // plot lines
grstyle set linewidth medthin:  ci        // confidence interval spikes
grstyle set linewidth medthin:  axisline  // axis lines
grstyle set linewidth medthin:  tick      // tick marks
grstyle set linewidth vthin:    xyline    // reference lines (e.g. xline)
grstyle set linewidth thin:     major_grid


*===========================================================================
* SYMBOLS
*===========================================================================

grstyle set symbolsize medium
grstyle set symbol circle diamond triangle square plus X smcircle


*===========================================================================
* LINE PATTERNS
*===========================================================================

grstyle set lpattern solid dash shortdash dash_dot longdash dot


*===========================================================================
* BACKGROUND AND PLOT REGION (publication style: white throughout)
*===========================================================================

grstyle color background white
grstyle set color black*.7:  tick tick_label
grstyle set color White:     plotregion plotregion_line
grstyle set color black*.04: major_grid
grstyle set color black*.7:  small_body
grstyle set color white*.08, opacity(0): pbarline


*===========================================================================
* COLOUR PALETTE
* Uses the plottig palette (publication-friendly, colour-blind accessible).
* The Tableau and Wong palettes below are kept as alternatives.
*===========================================================================

*--- Active palette: plottig ---
grstyle set color plottig, order(1 2 3 4 5 6)

/*
*--- Alternative: Tableau palette ---
grstyle set color "78 121 167"  ///  blue
                  "225 87 89"   ///  red
                  "242 142 43"  ///  orange
                  "0 110 80"    ///  dark green
                  "102 45 145"  ///  purple
                  "237 201 73"  ///  yellow
                  "89 161 79"   ///  bright green
                  "175 122 161"     // mauve

*--- Alternative: Wong (maximally distinct, colour-blind safe) ---
grstyle set color "0 0 0" "0 114 178" "213 94 0" "0 158 115" "204 121 167" "86 180 233" "230 159 0"
*/


*===========================================================================
* CONFIRM SETTINGS ARE ACTIVE
*===========================================================================

di as result "===== 00_graph_settings.do applied successfully ====="
