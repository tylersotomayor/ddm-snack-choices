* Inputs: estimated scalars and bin datasets. Outputs: four vector PDF figures.
version 19.0 // Use Stata 19's graph syntax and defaults.
foreach sample in individual pooled { // Draw a choice plot and a time plot for each sample.
    if "`sample'" == "individual" { // Set the individual graph's parameter names and horizontal range.
        local suffix i // Read the individual estimates, whose scalar names end in _i.
        local bound = scalar(delta_max) // Plot the individual curve from -M to M, matching the bin range.
        local tick "-4(2)4" // Label the individual horizontal axis every two dollars from -4 to 4.
    } // Finish the individual graph settings.
    else { // Use the following settings for the pooled class plots.
        local suffix c // Read the class estimates, whose scalar names end in _c.
        local bound = 10 // Use a symmetric range covering all supplied class differences.
        local tick "-10(2)10" // Label the class horizontal axis every two dollars from -10 to 10.
    } // Finish selecting settings for this sample.
    use "${PROCESSED}/`sample'-bins.dta", clear // Load the sample's observed and predicted bin means.
    local alpha = scalar(alpha_`suffix') // Copy fitted choice sensitivity into a local macro for the curve expression.
    local scale = scalar(A_`suffix') // Copy the fitted time scale for the response-time curve.
    local offset = scalar(t0_`suffix') // Copy the fitted nondecision time for the response-time curve.
    twoway /// Combine the theoretical choice curve with two sets of bin points.
        (function y=invlogit(2*`alpha'*x), range(-`bound' `bound') n(1001) lcolor(black) lwidth(medthin)) /// Draw P(d) at 1,001 evenly spaced points as a thin black curve.
        (scatter ej dj, msymbol(O) msize(small) mcolor(black)) /// Use small filled circles for the model's mean choice predictions.
        (scatter pj dj, msymbol(Oh) msize(medlarge) mcolor(black) mlwidth(medthin)), /// Use larger open circles for observed right-choice shares so overlaps remain visible.
        xtitle("Value difference, right minus left ($)", margin(t=3)) /// Label the bid difference in dollars; the top margin keeps the title clear of the tick labels in Latin Modern.
        ytitle("Frequency of choosing right") /// Label the vertical axis as a right-choice frequency.
        xlabel(`tick', labsize(small) nogrid) ylabel(0(.2)1, labsize(small) angle(horizontal) nogrid) /// Show the selected dollar ticks and frequencies from zero to one, without grid lines.
        legend(order(3 "Observed frequency" 2 "DDM bin prediction" 1 "DDM curve") /// Order the legend as observed points, predicted points, and theoretical curve.
            rows(1) size(small) region(lstyle(none)) position(6)) /// Put the legend below the plot on one row, with small text and no border.
        graphregion(color(white)) plotregion(color(white)) /// Use white backgrounds for both the plot and surrounding graph.
        xsize(6.5) ysize(4.0) name(`sample'_choice, replace) // Set a 6.5-by-4-inch canvas and give this sample's choice graph a distinct name.
    graph export "${FIGURES}/fig-`sample'-choice.pdf", as(pdf) replace // Write the choice figure as a vector PDF, replacing the prior run's file.
    twoway /// Combine the theoretical time curve with predicted and observed time means.
        (function y=`scale'*cond(abs(x)<1e-8, `alpha', tanh(`alpha'*x)/x)+`offset', /// Draw A*G(d)+t0, using G(0)=alpha within a tiny neighborhood of zero to avoid division by zero.
            range(-`bound' `bound') n(1001) lcolor(black) lwidth(medthin)) /// Use the same range, point density, and black line style as in the choice plot.
        (scatter tfit dj, msymbol(O) msize(small) mcolor(black)) /// Mark fitted bin times with small filled circles.
        (scatter tj dj, msymbol(Oh) msize(medlarge) mcolor(black) mlwidth(medthin)), /// Mark observed bin times with larger open circles.
        xtitle("Value difference, right minus left ($)", margin(t=3)) /// Keep the horizontal-axis label and its margin consistent with the choice figures.
        ytitle("Mean response time (ms)") /// Label mean response time in milliseconds.
        xlabel(`tick', labsize(small) nogrid) ylabel(, labsize(small) angle(horizontal) nogrid) /// Use the sample's dollar ticks and horizontal time labels, without grid lines.
        legend(order(3 "Observed mean" 2 "DDM bin prediction" 1 "DDM curve") /// List observed means first, followed by fitted bin means and the theoretical curve.
            rows(1) size(small) region(lstyle(none)) position(6)) /// Place the time legend below the plot with the same sizing as the choice legend.
        graphregion(color(white)) plotregion(color(white)) /// Keep both graph backgrounds white.
        xsize(6.5) ysize(4.0) name(`sample'_rt, replace) // Use the same canvas dimensions and a separate name for the response-time graph.
    graph export "${FIGURES}/fig-`sample'-response-time.pdf", as(pdf) replace // Write the time figure as a vector PDF, replacing the prior run's file.
} // Finish both plots for this sample, then continue to the next sample.
