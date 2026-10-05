* Inputs: the two assignment CSVs in data/raw.
* Outputs: processed data, four figure PDFs, table fragments, and a run log.
* Run from the project folder through ./run_stata.sh.
version 19.0 // Use the Stata 19 command definitions used for the report.
clear all // Clear data and stored estimation results before this run.
set more off // Let the run finish without pausing at a full results window.
capture log close _all // Close any open logs; an absent log is harmless here.
capture mkdir "logs" // A fresh copy of the package may not have a log directory yet.
log using "logs/stata-master.log", text replace // Record commands and results in a readable text file.
global RAW "data/raw" // Read source files from this directory without modifying them.
global PROCESSED "data/processed" // Keep every derived dataset in this directory.
global FIGURES "output/figures" // Export the four plots here.
global TABLES "output/tables" // Write the report's table bodies and numeric macros here.

do "code/01_setup.do" // Set numerical storage, graph styling, and output directories.
do "code/02_import_clean.do" // Select the individual sample and prepare both datasets.
do "code/03_analysis.do" // Estimate choice sensitivity, form bins, and fit response times.
do "code/04_figures.do" // Draw all four figures from those estimates and bin means.
do "code/05_tables.do" // Export the estimates and bin values used in the report's tables.
display as result "PROJECT 1 PIPELINE COMPLETED SUCCESSFULLY" // The shell script checks for this exact completion line.
log close // Flush and close the results log before Stata exits.
exit, clear // Close the batch session without a prompt about the data in memory.
