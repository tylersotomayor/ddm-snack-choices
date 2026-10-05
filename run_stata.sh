#!/bin/bash
# The first line selects Bash, which is included with macOS.
# Inputs: the numbered do-files and the two CSVs in data/raw.
# Outputs: four figures, table fragments, processed data, and logs.
set -euo pipefail # Stop on failed shell commands, unset variables, or failed pipeline stages.

cd -- "$(dirname -- "$0")" # Resolve every relative path from this script's folder, even when its path contains spaces.

STATA_BIN="${STATA_BIN:-/Applications/StataNow/StataBE.app/Contents/MacOS/StataBE}" # Allow another macOS Stata installation through the STATA_BIN environment variable.

if [ ! -x "$STATA_BIN" ]; then # Check the selected executable before starting the analysis.
  printf 'Stata was not found at: %s\nSet STATA_BIN to your macOS Stata executable.\n' "$STATA_BIN" >&2 # Show the missing path and how to choose another installation.
  exit 1 # Return a failure status without claiming that any outputs were rebuilt.
fi # Continue only when Stata can be launched.

mkdir -p logs data/intermediate data/processed output/figures output/tables # Make the parent directories needed by a fresh copy of the package.
rm -f logs/stata-master.log logs/stata-console.log 00_master.log # Remove old run logs so an earlier completion message cannot pass today's check.

stata_status=0 # Start with a successful shell status; replace it if Stata reports a failure.
"$STATA_BIN" -e do code/00_master.do || stata_status=$? # Run the master do-file in macOS batch mode and retain any failing exit status.

if [ -f 00_master.log ]; then # Stata may also write a batch console log in the working folder.
  mv 00_master.log logs/stata-console.log # Keep that console record with the main results log.
fi # An absent secondary log does not invalidate a completed main log.

if [ "$stata_status" -ne 0 ] || ! grep -qx 'PROJECT 1 PIPELINE COMPLETED SUCCESSFULLY' logs/stata-master.log; then # Check both the process status and the final printed result; some Stata errors leave a zero process status.
  printf 'The analysis stopped before completion. See logs/stata-master.log and logs/stata-console.log.\n' >&2 # Point to the command and error message that stopped the run.
  exit 1 # Report failure to the caller rather than accepting partial outputs.
fi # A completed log means all estimation and export stages ran.

for sample in individual pooled; do # Check both samples' figure exports.
  for kind in choice response-time; do # Each sample should have one choice plot and one time plot.
    figure="output/figures/fig-${sample}-${kind}.pdf" # Construct the filename used by 04_figures.do.
    if [ ! -s "$figure" ]; then # Require the expected PDF to exist and contain data.
      printf 'Missing or empty figure: %s\n' "$figure" >&2 # Identify the export that needs attention.
      exit 1 # Do not call the run successful if a required figure is absent.
    fi # Continue once this figure exists and is nonempty.
  done # Finish checking both plot types for the current sample.
done # Finish checking all four plots.
printf 'Done. All four plots are in output/figures/.\n' # Tell the reader where to find the regenerated figures.
