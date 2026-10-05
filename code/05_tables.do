* Inputs: scalars and bin datasets. Outputs: LaTeX fragments and numeric macros.
* Captions and column headers belong to the report's LaTeX section files.
version 19.0 // Use the same Stata version as the estimation stages.
local llp_text : display %12.3fc ll_c // Format the class likelihood with three decimals and a thousands separator.
local llp_math = subinstr(strtrim("`llp_text'"), ",", "{,}", .) // Brace the comma so LaTeX treats it as a separator rather than mathematical punctuation.
file open fh using "${TABLES}/tab-estimates.tex", write replace // Open the summary-table fragment, replacing the previous run's version.
file write fh "Choice trials & " %9.0fc (n_i) " & " %9.0fc (pool_n) " \\" _n // Write the individual and class choice sample sizes as whole numbers.
file write fh "Choice sensitivity & " %8.6f (alpha_i) " & " %8.6f (alpha_c) " \\" _n // Write both choice-sensitivity estimates to six decimal places.
file write fh "Log likelihood & \(" %12.3f (ll_i) "\) & \(" "`llp_math'" "\) \\" _n // Write the two log likelihoods in math mode so their minus signs render correctly.
file write fh "Choice RMSE (percentage points) & " %6.2f (100*choice_rmse_i) " & " %6.2f (100*choice_rmse_c) " \\" _n // Convert choice RMSE from a frequency to percentage points and show two decimals.
file write fh "\midrule" _n // Separate choice results from response-time results with a booktabs rule.
file write fh "Response-time bins & 15 & 15 \\" _n // Record the 15 observations in each response-time regression.
file write fh "Time scale, $\widehat{A}$ & " %10.3fc (A_i) " & " %10.3fc (A_c) " \\" _n // Write the fitted time scales with three decimals and thousands separators.
file write fh "Nondecision time, $\widehat{t}_0$ (ms) & " %10.3fc (t0_i) " & " %10.3fc (t0_c) " \\" _n // Write the fitted nondecision times in milliseconds with three decimals.
file write fh "\(R^2\) of bin means & " %6.4f (r2_i) " & " %6.4f (r2_c) " \\" _n // Report R-squared for variation in the 15 bin means, to four decimal places.
file write fh "Response-time RMSE (ms) & " %8.2f (rmse_i) " & " %8.2f (rmse_c) " \\" _n // Write time RMSE in milliseconds using the RSS/15 calculation from the analysis stage.
file write fh "\bottomrule" _n // Finish the summary table body with its bottom rule.
file close fh // Close the file so all summary-table rows are written to disk.

foreach sample in individual pooled { // Write an appendix bin table for each of the two samples.
    use "${PROCESSED}/`sample'-bins.dta", clear // Load the bin values used in this sample's figures.
    file open fh using "${TABLES}/tab-`sample'-bins.tex", write replace // Open this sample's table fragment, replacing the prior version.
    forvalues j=1/15 { // Write one table row for each of the 15 bins, in dataset order.
        if "`sample'" == "individual" { // Use the individual table's extra column for the number of retained time observations.
            file write fh %2.0f (bin[`j']) " & \(" %6.3f (dj[`j']) "\) & " %4.0f (n[`j']) " & " /// Write the bin index, inverse-probability coordinate in math mode, and choice count.
                %4.0f (n_rt[`j']) " & " %5.3f (pj[`j']) " & " %5.3f (ej[`j']) " & " /// Continue with the retained-time count and observed and predicted choice frequencies.
                %6.4f (gj[`j']) " & " %7.1f (tj[`j']) " & " %7.1f (tfit[`j']) " \\" _n // Finish the row with the time-function mean, observed time mean, and fitted time mean.
        } // Finish the individual row format.
        else { // Use the shorter row format for the class data, which have no separate retained-time count.
            file write fh %2.0f (bin[`j']) " & \(" %6.4f (dj[`j']) "\) & " %5.0fc (n[`j']) " & " /// Write the class bin index, supplied difference in math mode, and total choice count.
                %5.3f (pj[`j']) " & " %5.3f (ej[`j']) " & " %6.4f (gj[`j']) " & " /// Continue with observed and predicted right-choice shares and the time-function value.
                %7.1f (tj[`j']) " & " %7.1f (tfit[`j']) " \\" _n // Finish the class row with observed and fitted time means in milliseconds.
        } // Finish the class row format.
    } // Move to the next bin until all 15 rows have been written.
    file write fh "\bottomrule" _n // Add the bottom rule after the final bin.
    file close fh // Close this sample's table file before moving on.
} // Finish exporting the two appendix table bodies.

* Round prose values separately; the estimation and bin tables retain detail.
local llp_prose_text : display %12.2fc ll_c // Format the class likelihood to two decimals for its appearance in prose.
local llp_prose_math = subinstr(strtrim("`llp_prose_text'"), ",", "{,}", .) // Protect its thousands separator in LaTeX math mode, as in the summary table.
file open fh using "${TABLES}/numeric-values.tex", write replace // Open the file of numeric macros used by the report's prose.
file write fh "\newcommand{\AlphaIndividual}{" %5.3f (alpha_i) "}" _n // Define the individual alpha macro with three decimal places.
file write fh "\newcommand{\AlphaPooled}{" %5.3f (alpha_c) "}" _n // Define the class alpha macro with three decimal places.
file write fh "\newcommand{\LLIndividual}{" %8.2f (ll_i) "}" _n // Define the individual likelihood macro with two decimal places.
file write fh "\newcommand{\LLPooled}{" "`llp_prose_math'" "}" _n // Define the class likelihood macro from the formatted math-safe text.
file write fh "\newcommand{\TimeScaleIndividual}{" %5.0f (A_i) "}" _n // Round the individual time scale to a whole number for running prose.
file write fh "\newcommand{\OffsetIndividual}{" %4.0f (t0_i) "}" _n // Round individual nondecision time to a whole millisecond.
file write fh "\newcommand{\TimeScalePooled}{" %5.0f (A_c) "}" _n // Round the class time scale to a whole number for running prose.
file write fh "\newcommand{\OffsetPooled}{" %4.0f (t0_c) "}" _n // Round class nondecision time to a whole millisecond.
file write fh "\newcommand{\RSquaredIndividual}{" %5.3f (r2_i) "}" _n // Define individual R-squared to three decimal places.
file write fh "\newcommand{\RSquaredPooled}{" %5.3f (r2_c) "}" _n // Define class R-squared to three decimal places.
file write fh "\newcommand{\RTMean}{" %6.0fc (rt_mean) "}" _n // Write the full-sample time mean in whole milliseconds with a thousands separator.
file write fh "\newcommand{\RTSD}{" %6.0fc (rt_sd) "}" _n // Write the sample SD of response time in the same whole-millisecond format.
file write fh "\newcommand{\RTLower}{" %8.2f (rt_lower) "}" _n // Keep two decimal places for the lower trimming cutoff.
file write fh "\newcommand{\RTUpper}{" %8.2f (rt_upper) "}" _n // Keep two decimal places for the upper trimming cutoff.
file write fh "\newcommand{\RTKept}{" %3.0f (rt_n) "}" _n // Write the exact number of individual times retained after trimming.
file write fh "\newcommand{\RTRemoved}{" %3.0f (rt_removed) "}" _n // Write the exact number of individual times excluded by trimming.
file write fh "\newcommand{\MaxDelta}{" %4.2f (delta_max) "}" _n // Round the maximum absolute bid difference to two decimals for prose.
file write fh "\newcommand{\BinWidth}{" %5.3f (bin_width) "}" _n // Round the common individual bin width to three decimals for prose.
file write fh "\newcommand{\ChoiceRMSEIndividual}{" %4.1f (100*choice_rmse_i) "}" _n // Write individual choice RMSE in percentage points with one decimal.
file write fh "\newcommand{\ChoiceRMSEPooled}{" %4.1f (100*choice_rmse_c) "}" _n // Write class choice RMSE in percentage points with one decimal.
file write fh "\newcommand{\RTRMSEIndividual}{" %4.0f (rmse_i) "}" _n // Round individual time RMSE to whole milliseconds for prose.
file write fh "\newcommand{\RTRMSEPooled}{" %4.0f (rmse_c) "}" _n // Round class time RMSE to whole milliseconds for prose.
file write fh "\newcommand{\CurveGapIndividual}{" %4.1f (max_curve_gap_i) "}" _n // Write the largest individual dot-to-curve time gap with one decimal.
file write fh "\newcommand{\AltAIndividual}{" %5.0f (A_all_i) "}" _n // Write the alternative individual time-scale estimate for the appendix.
file write fh "\newcommand{\AltOffsetIndividual}{" %4.0f (t0_all_i) "}" _n // Write the alternative nondecision-time estimate for the appendix.
file write fh "\newcommand{\AltRSquaredIndividual}{" %5.3f (r2_all_i) "}" _n // Write the alternative regression's R-squared to three decimals.
file close fh // Close the numeric-macro file after writing all 26 values.
