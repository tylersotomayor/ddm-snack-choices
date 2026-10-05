* Inputs: clean trial and pooled datasets. Outputs: fits and 15-bin summaries.
* Lecture 5 normalization: P=invlogit(2*alpha*delta); G=tanh(alpha*delta)/delta.
* Both time regressions give equal weight to the 15 bin means.
* Scalar suffix i denotes the individual; c denotes the pooled class.
version 19.0 // Use the same Stata version as in data preparation.
use "${PROCESSED}/individual-trials.dta", clear // Load all individual choices, including trials whose times were excluded.
logit chose_r delta, noconstant iterate(200) tolerance(1e-12) // Fit log odds = 2*alpha*delta, with no intercept and a tight convergence tolerance.
assert e(converged) == 1 // Do not use a choice estimate if the optimizer failed to converge.
scalar alpha_i = _b[delta] / 2 // Divide the logit slope by two to recover the lecture's alpha.
scalar alpha_se_i = _se[delta] / 2 // Apply the same rescaling to the slope's standard error.
scalar ll_i = e(ll) // Save the maximized log likelihood before another regression replaces it.
scalar n_i = e(N) // Record the number of trials entering the choice likelihood.
assert alpha_i > 0 // Check that right-choice probability increases with the right snack's bid advantage.
gen double pfit = invlogit(2 * alpha_i * delta) // Evaluate each trial's probability of choosing right at the fitted alpha.
gen double gfit = cond(delta == 0, alpha_i, tanh(alpha_i * delta) / delta) // Use G(0)=alpha at the equal-bid trial to avoid division by zero.
assert !missing(pfit, gfit) // Both model functions must be defined on every trial.
gen double score = 2 * delta * (chose_r - pfit) // Calculate each trial's derivative of log likelihood with respect to alpha.
quietly summarize score // Sum those derivative contributions across the choice sample.
scalar score_i = r(sum) // Store the score at the unrounded maximum-likelihood estimate.
assert abs(score_i) < 1e-5 // Check that the estimate nearly solves the first-order condition.
gen double curvature = -4 * delta^2 * pfit * (1 - pfit) // Calculate each trial's second derivative with respect to alpha.
quietly summarize curvature // Add the curvature contributions for the full log likelihood.
scalar hessian_i = r(sum) // Keep the total second derivative for inspection.
assert hessian_i < 0 // Negative curvature confirms a maximum at a stationary point.
gen double loglik = ln(invlogit(2 * alpha_i * sigma * delta)) // Evaluate log P(sigma*delta), the assignment's trial-level likelihood expression.
quietly summarize loglik // Add the log probabilities over trials.
assert abs(r(sum) - ll_i) < 1e-7 // Confirm that this likelihood formula agrees with Stata's logit likelihood.
gen double rt_used = rt if rt_keep // Set excluded times to missing so they do not enter the bin means.
gen double g_used = gfit if rt_keep // Average G over exactly the same retained trials as response time.
save "${PROCESSED}/individual-trials.dta", replace // Add fitted probabilities and time-function values to the saved trial data.
export delimited using "${PROCESSED}/individual-trials.csv", replace // Update the readable trial file with those fitted quantities.
collapse (count) n=chose_r n_rt=rt_used (sum) n_right=chose_r /// Count all choices, retained times, and right choices within each bin.
    (mean) pj=chose_r ej=pfit mean_delta=delta tj=rt_used gj=g_used g_all=gfit, by(bin) // Form bin means; g_all also keeps the time-function mean over untrimmed trials.
assert _N == 15 // All 15 individual intervals must survive aggregation.
assert n > 0 & n_rt > 0 // Every bin needs choices and at least one retained response time.
gen double lower = -delta_max + (bin - 1) * bin_width // Recover each bin's lower boundary from its index.
gen double upper = -delta_max + bin * bin_width // Recover each bin's upper boundary using the same width.
gen double dj = ln(ej / (1 - ej)) / (2 * alpha_i) // Invert the mean predicted probability to obtain the coordinate required in 1(d).
assert dj >= lower - 1e-10 & dj <= upper + 1e-10 // The inverse coordinate lies within its bin, allowing tiny numerical error.
assert abs(invlogit(2 * alpha_i * dj) - ej) < 1e-12 // Substituting d_j back into P must recover the bin's predicted probability.
regress tj gj // Fit T_j=t0+A*g_j by OLS, giving all 15 means equal weight.
scalar A_i = _b[gj] // The slope estimates the individual time scale.
scalar t0_i = _b[_cons] // The intercept estimates time outside evidence accumulation, in milliseconds.
scalar r2_i = e(r2) // Save the fraction of variation in bin-mean times explained by this regression.
scalar rmse_i = sqrt(e(rss) / e(N)) // Use RSS/15 for descriptive RMSE, rather than Stata's residual variance with 13 degrees of freedom.
assert A_i > 0 // Check that the fitted time scale is positive.
predict double tfit // Obtain fitted bin times from the regression just estimated.
gen double residual = tj - tfit // Positive residuals mean decisions took longer than predicted.
gen double choice_resid = pj - ej // Measure observed minus predicted right-choice frequency in each bin.
gen double weighted_choice_sq = n * choice_resid^2 // Weight each squared choice error by the number of trials in its bin.
quietly summarize weighted_choice_sq // Sum the weighted squared choice errors.
scalar choice_rmse_i = sqrt(r(sum) / n_i) // Normalize by all choice trials and take the square root.
gen double choice_sq = choice_resid^2 // Also retain squared errors without trial-count weights.
quietly summarize choice_sq // Average those squared errors over the 15 bins.
scalar choice_bin_rmse_i = sqrt(r(mean)) // Keep the unweighted choice RMSE as a diagnostic in the log.
gen double rt_curve = A_i * cond(dj == 0, alpha_i, tanh(alpha_i * dj) / dj) + t0_i // Evaluate the smooth time curve at the plotted coordinates.
gen double dot_curve_gap = abs(tfit - rt_curve) // Compare A*mean(G)+t0 with A*G(d_j)+t0; they need not coincide for individual bins.
quietly summarize dot_curve_gap // Find the largest absolute separation between the fitted dots and curve.
scalar max_curve_gap_i = r(max) // Save that separation for the appendix's explanation of Figure 2.
* The appendix also fits the same time means with G averaged over all trials.
quietly regress tj g_all // Change only the averaging sample for G; retain the trimmed time means and equal bin weights.
scalar A_all_i = _b[g_all] // Store the alternative slope without replacing the main A estimate.
scalar t0_all_i = _b[_cons] // Store the alternative intercept for the appendix.
scalar r2_all_i = e(r2) // Save the alternative regression's fit for comparison.
save "${PROCESSED}/individual-bins.dta", replace // Keep the main fitted dots and both averaging conventions in the bin dataset.
export delimited using "${PROCESSED}/individual-bins.csv", replace // Export every individual plotted coordinate and prediction in readable form.

use "${PROCESSED}/pooled-input.dta", clear // Start the class analysis from the 15 supplied summaries.
* Expand each bin into two rows, weighted by counts of R and L choices.
preserve // Keep a copy of the bin data while constructing the likelihood representation.
expand 2 // Make one row for right choices and one for left choices in each bin.
bysort bin: gen byte chose_r = _n == 1 // Label the first copy as right and the second as left.
gen long weight = cond(chose_r, nj, mj) // Attach the supplied right or left count to each copy.
logit chose_r dj [fw=weight], noconstant iterate(200) tolerance(1e-12) // Frequency weights fit the grouped likelihood, treating each bin mean as the trial difference.
assert e(converged) == 1 // Require the class choice fit to converge before using its estimates.
scalar alpha_c = _b[dj] / 2 // Convert the class logit slope to the same alpha normalization as the individual fit.
scalar alpha_se_c = _se[dj] / 2 // Rescale the class slope's standard error accordingly.
scalar ll_c = e(ll) // Save the likelihood without binomial coefficients, which do not affect its maximizer.
assert alpha_c > 0 // Check that class right-choice probabilities increase with the bid advantage.
restore // Return to one row per class bin while keeping the estimated scalars.
gen double ej = invlogit(2 * alpha_c * dj) // Predict the right-choice frequency at each supplied mean difference.
gen double gj = cond(dj == 0, alpha_c, tanh(alpha_c * dj) / dj) // Evaluate the class time function, using its limit if a mean difference is zero.
gen double loglik = nj * ln(ej) + mj * ln(1 - ej) // Calculate the grouped likelihood directly from the two choice counts.
quietly summarize loglik // Sum the 15 bin contributions to the class likelihood.
assert abs(r(sum) - ll_c) < 1e-7 // Check that the grouped formula matches the frequency-weighted logit result.
gen double score = 2 * dj * (nj - n * ej) // Evaluate each bin's contribution to the derivative with respect to class alpha.
quietly summarize score // Add the derivative contributions across bins.
scalar score_c = r(sum) // Save the score at the fitted class parameter.
assert abs(score_c) < 1e-5 // Check that the class estimate nearly satisfies the first-order condition.
regress tj gj // Fit the 15 supplied time means with equal OLS weights; individual times are unavailable for trimming.
scalar A_c = _b[gj] // Store the class time-scale estimate from the slope.
scalar t0_c = _b[_cons] // Store the class nondecision-time estimate from the intercept.
scalar r2_c = e(r2) // Record the fit across the 15 class time means.
scalar rmse_c = sqrt(e(rss) / e(N)) // Compute the class time RMSE with the same RSS/15 convention as the individual fit.
assert A_c > 0 // Check the positive time-scale condition in the class fit.
predict double tfit // Generate the fitted class response-time means.
gen double residual = tj - tfit // Retain each class bin's time discrepancy for inspection.
gen double choice_resid = pj - ej // Subtract predicted from observed class right-choice shares.
gen double weighted_choice_sq = n * choice_resid^2 // Weight squared choice discrepancies by their class trial counts.
quietly summarize weighted_choice_sq // Sum those count-weighted errors over bins.
scalar choice_rmse_c = sqrt(r(sum) / pool_n) // Express choice RMSE as a frequency; the tables convert it to percentage points.
gen double choice_sq = choice_resid^2 // Keep the unweighted squared choice discrepancies too.
quietly summarize choice_sq // Average them across the 15 supplied bins.
scalar choice_bin_rmse_c = sqrt(r(mean)) // Store unweighted choice RMSE as an additional log diagnostic.
save "${PROCESSED}/pooled-bins.dta", replace // Save the class coordinates, observed means, and fitted values for plotting.
export delimited using "${PROCESSED}/pooled-bins.csv", replace // Export the same class bin data for inspection without Stata.
scalar list // Print the saved counts, parameter estimates, and fit checks in the run log.
