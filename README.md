# Snack Choices and Decision Times under the Drift-Diffusion Model

Project Report 1 for ECON GU4850, Cognitive Mechanisms and Economic Behavior,
Columbia University, Fall 2026. The report fits the drift-diffusion model to
435 of my own snack choices and to 30,450 choices pooled across the class,
using bid differences as the value signal, and then tests the model's
response-time prediction with the choice parameter held fixed.

Report: `output/pdf/project-1-report.pdf`.

The standalone replication package is in `replication/project-1-replication/`,
with a distribution ZIP at `output/replication/project-1-replication.zip`.
Its README gives the one-command run instructions. All six numbered do-files
and the shell runner have line-by-line comments. The package includes only
the 435 individual records used in this report and the supplied class means;
its code is identical to the analysis code in this project.

The report answers Questions 1(a)–1(d), 2, 3(a), and 3(b) using the rows for
UNI `ts3686` and the supplied pooled data. Its LaTeX preamble, frontmatter,
heading styles, bibliography, captions, and directory structure are adapted
from the existing economics working-paper template. Only the template's
reusable source and build files were copied; unrelated research projects and
sample outputs were not imported. Authorship is adapted to this individual
course assignment.

## Rebuild

From this directory:

```sh
./run_stata.sh
python3 code/verify_results.py
./compile.sh
./check_refs.sh --strict
```

StataBE is expected at the path in `run_stata.sh`. The LaTeX build requires
`pdflatex` and `biber`. The figures use Latin Modern Roman, which Stata can
only embed if the font is installed in macOS; its OpenType files ship with
TeX Live. `main.tex` sets `\pdfinclusioncopyfonts=1` so pdfTeX keeps the
fonts embedded in the figure PDFs rather than substituting the body font,
which would drop the ff ligature in the axis titles.

`python3 code/replicate.py` rebuilds the same datasets, table fragments, and
figures without Stata under `output/python/` and checks them against the Stata
outputs; the table fragments match byte for byte. The report itself uses only
the Stata figures.

Final empirical figures and table fragments come only
from the numbered Stata do-files. The Python script independently audits
the results using bisection and closed-form OLS; it does not generate final
empirical outputs.

## Data and conventions

- `data/raw/indv-Data.csv` holds the 435 original records for UNI `ts3686`,
  extracted byte for byte from the class file. The full class file carries
  other students' identifiers and is not distributed. `data/raw/popData.csv`
  is the supplied class summary, unchanged.
- The individual sample keeps all 435 records for `ts3686`, including repeated
  goods pairs. `item1` is interpreted as left and `item2` as right.
- Choice probabilities use `invlogit(2*alpha*delta)` with no intercept.
- The 15 individual intervals are left-closed and right-open; the final one
  includes the maximum. Their common width is `2*max(abs(delta))/15`.
- Response-time trimming uses the original sample mean and sample standard
  deviation once. Observations at or beyond two standard deviations are
  excluded only from response-time calculations.
- Individual time means and time-function means use the same retained trials.
  The choice estimate and inverse-probability coordinates remain unchanged.
  The appendix also reports the alternative convention using all trials for
  the time-function mean.
- The pooled file is used as supplied. Its 15 mean differences are treated as
  exact trial differences for the approximate likelihood and time calculation.
- Both response-time regressions are unweighted OLS over the 15 bins.
- `data/processed/*-bins.csv` gives every plotted coordinate and prediction.
  The LaTeX numeric macros are regenerated from the estimated Stata scalars.

The source notes in `literature/notes/` identify the lecture and recitation
passages used for the report.
