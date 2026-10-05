* Inputs: assignment CSVs in data/raw. Outputs: individual and pooled .dta files.
* item1 is the left good; item2 is the right good. RT is in milliseconds.
* Choice estimation keeps all 435 trials. Trimming is applied only to RT analysis.
* The replication copy contains only ts3686's original rows; the full course file works too.
version 19.0 // Use the Stata 19 command definitions used for the report.
import delimited using "${RAW}/indv-Data.csv", varnames(1) case(lower) clear // Read the header as variable names and convert those names to lowercase.
foreach v in item1id item2id chosenid columbiaid { // Clean surrounding spaces in the four identifier fields.
    replace `v' = strtrim(`v') // Make equality checks insensitive to padding in the CSV.
} // Finish cleaning identifiers before selecting the sample.
keep if columbiaid == "ts3686" // Retain the participant whose choices the report analyzes.
assert _N == 435 // The supplied individual sample has 435 trials, including repeated pairs.
assert !missing(item1v, item2v, rt) // Each trial needs both bids and a response time.
assert rt > 0 // Response times must be strictly positive.
assert chosenid == item1id | chosenid == item2id // The choice must be one of the two displayed snacks.
gen long trial = _n // Number the retained trials in their original order.
gen str60 pair = cond(item1id < item2id, item1id + "|" + item2id, item2id + "|" + item1id) // Alphabetize each pair's IDs so left/right reversals count as the same pair.
* Repeated pairs are distinct trials in the supplied file; do not deduplicate.
egen byte pair_tag = tag(pair) // Mark one record per unordered pair without deleting repeated trials.
quietly count if pair_tag // Count the marked pairs rather than the number of trials.
scalar n_pairs = r(N) // Save the distinct-pair count for inspection in the log.
isid trial // Check that the trial index uniquely identifies each retained record.
gen double delta = item2v - item1v // A positive difference means the right snack has the higher bid.
gen byte chose_r = chosenid == item2id // Code a right choice as one and a left choice as zero.
gen byte sigma = 2 * chose_r - 1 // Convert choices to +1 for right and -1 for left in the likelihood formula.
gen double absdelta = abs(delta) // Measure the bid gap regardless of which side has the higher bid.
quietly summarize absdelta // Find the largest absolute bid difference in this sample.
scalar delta_max = r(max) // Save M, the extent of the symmetric bin range.
scalar bin_width = 2 * delta_max / 15 // Divide the interval from -M to M into 15 equal parts.
* Intervals are left-closed/right-open, except the final interval includes +M.
gen byte bin = min(15, floor((delta + delta_max) / bin_width) + 1) // Use left-closed bins and keep the endpoint +M in bin 15.
assert inrange(bin, 1, 15) // Every trial must belong to one of the 15 intervals.
quietly summarize rt // Compute the mean and sample SD once, using all 435 response times.
scalar rt_mean = r(mean) // Keep the full-sample mean before excluding any time observations.
scalar rt_sd = r(sd) // Stata's sample SD uses N-1 in its variance denominator.
scalar rt_lower = rt_mean - 2 * rt_sd // Record the lower two-SD cutoff for the report.
scalar rt_upper = rt_mean + 2 * rt_sd // Record the upper two-SD cutoff for the report.
gen byte rt_keep = abs(rt - rt_mean) < 2 * rt_sd // Flag times strictly within both cutoffs without dropping choice records.
quietly count if rt_keep // Count the trials available for the response-time means.
scalar rt_n = r(N) // Save the number retained by this single trimming step.
scalar rt_removed = _N - rt_n // Count the excluded times relative to the original 435 trials.
quietly count if delta == 0 // Check how often the two bids are exactly equal.
scalar n_ties = r(N) // Save the tie count; these trials need the limit of the time function at zero.
quietly count if sigma * delta > 0 // Count choices of the snack with the higher bid.
scalar n_consistent = r(N) // Retain that count as a descriptive check on choice coding.
quietly count if sigma * delta < 0 // Count choices of the lower-bid snack, leaving ties separate.
scalar n_reversals = r(N) // Save the number of choices against the bid ranking.
save "${PROCESSED}/individual-trials.dta", replace // Save all trials, their bins, and the time-retention flag for estimation.
export delimited using "${PROCESSED}/individual-trials.csv", replace // Provide a readable copy of the same trial-level data.

import delimited using "${RAW}/popData.csv", varnames(1) clear // Load the supplied class means and choice counts.
assert _N == 15 // The class file must have one row for each of its 15 bins.
assert !missing(dj, nj, mj, tj) // Require a bid difference, both choice counts, and a time mean in every bin.
assert nj >= 0 & mj >= 0 & tj > 0 // Counts cannot be negative and mean times must be positive.
gen byte bin = _n // Preserve the supplied ordering of class bins.
gen long n = nj + mj // Add right and left choices to obtain each bin's total count.
gen double pj = nj / n // Compute the observed fraction choosing right.
assert n > 0 // Every bin needs at least one choice to define its choice frequency.
quietly summarize n // Add the bin counts to recover the class sample size.
scalar pool_n = r(sum) // Keep the total for count-weighted choice fit and the summary table.
save "${PROCESSED}/pooled-input.dta", replace // Save the class inputs; individual times are unavailable for trimming.
