* Inputs: project directories. Outputs: reproducible environment and paths.
version 19.0 // Keep this stage under the same Stata version as the master file.
set more off // Do not pause if a command produces a long results listing.
set type double // Store newly created numeric variables in double precision by default.
set scheme s2mono // Use Stata's built-in monochrome graph scheme.
graph set window fontface "Latin Modern Roman" // Use the report's plot font in the graph window.
graph set print fontface "Latin Modern Roman" // Use the same font in the exported figures.
foreach path in "data/intermediate" "data/processed" "output/figures" "output/tables" "output/pdf" { // Create the project's output directories one at a time.
    capture mkdir `path' // Leave an existing directory in place and continue.
} // Finish creating directories before checking the inputs.
confirm file "${RAW}/indv-Data.csv" // Stop here if the individual-trial file is missing.
confirm file "${RAW}/popData.csv" // Stop here if the supplied class summaries are missing.
