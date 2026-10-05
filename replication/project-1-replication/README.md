# Project 1 replication files

**Snack Choices and Decision Times under the Drift-Diffusion Model**  
Tyler Sotomayor · ECON GU4850 · October 5, 2026

This folder reproduces the report's four figures from the individual choice
records and supplied class summaries. The six numbered Stata do-files also
write the table bodies and numeric values used in the report, and a single
Python script reproduces the same analysis without Stata. Each executable
Stata line has a comment explaining its purpose, including the separate
lines of multiline plotting commands.

## Run

The package uses **Stata 19 on macOS** and no community-written commands.
It was tested with StataNow/StataBE. LaTeX is not needed. Python is needed
only for the alternative path described below.

The plots use Latin Modern Roman, the typeface of the report's body text, so
figure labels match the surrounding type. It is not a stock macOS font. The
OpenType files (`lmroman10-regular.otf` and its companions) ship with TeX Live
and MacTeX under `texmf-dist/fonts/opentype/public/lm/` and are also
distributed by GUST. Install them with Font Book before running. Without them,
Stata substitutes another font; the numbers and geometry are unaffected.

Extract the ZIP and open Terminal in the `project-1-replication` folder. Run:

```sh
./run_stata.sh # Rebuild the analysis, all four figures, and the table fragments.
```

The script looks for Stata at
`/Applications/StataNow/StataBE.app/Contents/MacOS/StataBE`. If your macOS
installation is elsewhere, supply the full executable path when running it:

```sh
STATA_BIN="/path/to/your/Stata/executable" ./run_stata.sh # Use your installed Stata binary.
```

A completed run prints `Done. All four plots are in output/figures/.` Existing
derived data, figures, tables, and logs are replaced when the script runs again;
the two input files stay unchanged. If the run stops, the last command in
`logs/stata-master.log` identifies where it stopped. Stata's separate batch
record is saved as `logs/stata-console.log` when available.

## Python alternative

`code/replicate.py` reproduces the same analysis without Stata. It needs
Python 3.8 or later and matplotlib. Every number is computed with the standard
library alone: the choice parameter solves the likelihood's first-order
condition by bisection, and the response-time regressions use closed-form OLS.
matplotlib is used only to draw. Run:

```sh
python3 code/replicate.py # Rebuild the datasets, table fragments, and figures under output/python/.
```

It writes `output/python/data/`, `output/python/tables/`, and
`output/python/figures/`, then compares its table fragments byte for byte with
the Stata copies in `output/tables/` and, when `data/processed/` exists from a
Stata run, its datasets with Stata's to a relative tolerance of 1e-9. It exits
with status 1 on any mismatch. The four fragments and 26 numeric macros match
the Stata versions byte for byte. The figures show the same points, curves,
labels, and legend as the Stata figures, but a different drawing engine
produces them, so they are not pixel-identical. They use Latin Modern Roman
when it is installed and fall back to Times New Roman or DejaVu Serif. The
report itself uses the Stata figures.

## Figures and code

| Report figure | PDF in `output/figures/` |
|---|---|
| 1: Individual choice frequencies | `fig-individual-choice.pdf` |
| 2: Individual mean response times | `fig-individual-response-time.pdf` |
| 3: Pooled choice frequencies | `fig-pooled-choice.pdf` |
| 4: Pooled mean response times | `fig-pooled-response-time.pdf` |

Copies of the finished figures and table fragments are included. The script
recomputes them from the inputs; it does not read the existing copies.

`code/00_master.do` is the entry point. It calls the remaining files in order:

| Do-file | Work done |
|---|---|
| `01_setup.do` | Set numerical storage, create directories, and set the plot style. |
| `02_import_clean.do` | Select the individual trials, assign bins, and flag response-time outliers. |
| `03_analysis.do` | Estimate choice sensitivity, form bin means, and fit the time regressions. |
| `04_figures.do` | Draw the observed points, fitted points, and theoretical curves. |
| `05_tables.do` | Write the summary table, both bin tables, and 26 numeric macros. |
| `replicate.py` | Python replication of all five stages, with a self-check against the Stata outputs. |

The master file closes Stata after the batch run. Running a later do-file alone
will not work because it uses datasets or scalars created by an earlier stage.

## Data

`data/raw/indv-Data.csv` contains the **435 original records for `ts3686`**
selected from the course's individual-trial file. The header, selected record
bytes, and their order are preserved. No bids or times have been rounded, and
repeated pairs have not been removed. The original file also contains other
students' records, which are unnecessary for these figures and are omitted here.

| Original column | Meaning |
|---|---|
| `item1ID`, `item2ID` | IDs of the left and right snacks, respectively. |
| `item1V`, `item2V` | Their bid values, in dollars. |
| `chosenID` | ID of the chosen snack. |
| `RT` | Response time, in milliseconds. |
| `columbiaID` | Participant identifier; every included row is `ts3686`. |

`data/raw/popData.csv` is an unchanged copy of the supplied class file. Its
15 rows summarize 30,450 choices.

| Column | Meaning |
|---|---|
| `dj` | Supplied mean bid difference, right minus left, in dollars. |
| `nj`, `mj` | Counts of right and left choices. |
| `tj` | Supplied mean response time, in milliseconds. |

All derived datasets are written to `data/processed/`. The two `*-bins.csv`
files contain the coordinates behind the figures: `dj` is the horizontal
coordinate, `pj` the observed right-choice share, `ej` its model prediction,
`tj` the observed time mean, and `tfit` the fitted time mean. `gj` is the time
function used in the regression. The individual file also records the total
choice count `n` and retained-time count `n_rt` in each bin.

## Calculation notes

The bid difference is right minus left. The choice model is
`P(delta) = invlogit(2 * alpha * delta)`, so alpha is half the slope from a
logit regression with no intercept. The comments in `03_analysis.do` connect
that regression to the likelihood, score, and curvature in Question 1.
The scalar suffixes `_i` and `_c` identify the individual and class estimates.

For the individual figures, the code divides `[-M, M]` into 15 equal intervals,
where `M` is the largest absolute bid difference. Intervals include their left
endpoint; the final interval also includes `M`. The plotted coordinate is
the inverse of the bin's mean predicted choice probability. This puts each
filled choice dot on the theoretical curve, as in Question 1(d).

Response-time trimming uses the original 435-trial mean and sample standard
deviation once. Times at or beyond two standard deviations are excluded from
the time calculations, leaving 418 observations. All 435 choices still enter
the likelihood and choice means. Both the time mean and the mean of
`G(delta) = tanh(alpha * delta) / delta` use the same retained trials. At zero,
the code uses the limit `G(0) = alpha`.

Both time regressions fit `T_j = t0 + A * g_j` with equal weight on the 15 bins.
The individual fitted time dots use the average of G within each bin, whereas
the smooth curve evaluates G at the plotted coordinate. Their small separation
is expected. The appendix's alternative averaging convention is calculated
separately and does not change the plotted dots.

For the class data, the supplied mean difference is treated as each trial's
difference within its bin. Right and left counts enter as frequency weights
in the choice likelihood. The time file has only bin means, so there is no
individual-level trimming step. Under this approximation, the supplied `dj`
also equals the inverse-probability coordinate, and the fitted time dots lie
on the smooth curve.

The formulas follow Michael Woodford's *Lecture 5: The Drift-Diffusion Model*
(September 23, 2026, slides 10–11 and 24–26) and Andrew Souther's recitation
notes from September 18 (pp. 4–5) and September 25 (pp. 3–4), for ECON GU4850.
The course notes are not needed to execute the code.

`checksums.sha256` lists the code, the Python script, data, README, runner, and
table fragments.
Run `shasum -a 256 -c checksums.sha256` in this folder to confirm that an
extracted or regenerated copy matches. The four figure PDFs are deliberately
not listed: Stata writes creation metadata into each export, so their bytes
change from run to run even when the drawing does not. Compare regenerated
figures by viewing them, or rasterize them with `pdftoppm -r 150` and compare
the images.

The package was checked from a fresh extracted copy with the bundled outputs
removed. All four regenerated plots matched the report figures when rendered
at 150 dpi, and the processed CSVs and table fragments matched byte for byte.
