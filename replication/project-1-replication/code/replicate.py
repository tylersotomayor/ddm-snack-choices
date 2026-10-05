#!/usr/bin/env python3
"""Python replication of the Project 1 analysis, with no Stata required.

Reproduces what the six numbered do-files produce:

  output/python/data/     the trial-level and 15-bin datasets
  output/python/tables/   the four LaTeX fragments and 26 numeric macros,
                          byte-identical to the Stata versions
  output/python/figures/  the four report figures, drawn with matplotlib

Every number comes from the standard library: the choice parameter solves
the likelihood's first-order condition by bisection, and the response-time
regressions use closed-form OLS. matplotlib is imported only to draw.

Run from the package root:

    python3 code/replicate.py

If Stata outputs exist in output/tables/ and data/processed/, the script
compares its results to them and exits with status 1 on any mismatch.
"""

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "output" / "python"
SUBJECT = "ts3686"
NBINS = 15


# ----------------------------------------------------------------- model ----
def choice_probability(alpha, delta):
    """P(right) = e^{a d} / (e^{a d} + e^{-a d}) = 1 / (1 + e^{-2 a d})."""
    return 1.0 / (1.0 + math.exp(-2.0 * alpha * delta))


def time_function(alpha, delta):
    """G(d) = tanh(a d) / d, with its limit a at d = 0."""
    return alpha if delta == 0 else math.tanh(alpha * delta) / delta


def solve_alpha(groups):
    """Maximum-likelihood alpha for grouped data (delta, right count, count).

    The log likelihood is strictly concave in alpha, so the maximum is the
    unique zero of the score. Bisection brackets it and halves 200 times,
    which pins the root to the last bit of a double.
    """

    def score(alpha):
        return math.fsum(2.0 * d * (r - n * choice_probability(alpha, d)) for d, r, n in groups)

    low, high = 0.0, 1.0
    if score(low) <= 0:
        sys.exit("the score at alpha = 0 is not positive; no positive maximum")
    while score(high) > 0:
        high *= 2.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if score(middle) > 0:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def ols(x, y):
    """Slope, intercept, R^2, RMSE with RSS/n, and fitted values."""
    n = len(x)
    xm, ym = math.fsum(x) / n, math.fsum(y) / n
    sxx = math.fsum((a - xm) ** 2 for a in x)
    sxy = math.fsum((a - xm) * (b - ym) for a, b in zip(x, y))
    slope = sxy / sxx
    intercept = ym - slope * xm
    fit = [intercept + slope * a for a in x]
    rss = math.fsum((b - f) ** 2 for b, f in zip(y, fit))
    tss = math.fsum((b - ym) ** 2 for b in y)
    return slope, intercept, 1.0 - rss / tss, math.sqrt(rss / n), fit


def mean(values):
    values = list(values)
    return math.fsum(values) / len(values)


# ------------------------------------------------------------------ data ----
def read_csv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        return [{k.strip(): v.strip() for k, v in row.items()} for row in reader]


def individual_analysis():
    rows = [r for r in read_csv(RAW / "indv-Data.csv") if r["columbiaID"] == SUBJECT]
    if len(rows) != 435:
        sys.exit(f"expected 435 trials for {SUBJECT}, found {len(rows)}")
    trials = []
    seen_pairs = set()
    for i, r in enumerate(rows, start=1):
        if r["chosenID"] not in (r["item1ID"], r["item2ID"]):
            sys.exit(f"trial {i}: chosen item is not one of the pair")
        v1, v2, rt = float(r["item1V"]), float(r["item2V"]), float(r["RT"])
        if rt <= 0:
            sys.exit(f"trial {i}: response time must be positive")
        pair = "|".join(sorted((r["item1ID"], r["item2ID"])))
        chose_r = int(r["chosenID"] == r["item2ID"])
        trials.append({
            "item1id": r["item1ID"], "item2id": r["item2ID"], "item1v": v1, "item2v": v2,
            "chosenid": r["chosenID"], "rt": rt, "columbiaid": r["columbiaID"], "trial": i,
            "pair": pair, "pair_tag": int(pair not in seen_pairs), "delta": v2 - v1,
            "chose_r": chose_r, "sigma": 2 * chose_r - 1, "absdelta": abs(v2 - v1),
        })
        seen_pairs.add(pair)

    # Bins: 15 equal intervals on [-M, M], left-closed, with +M in bin 15.
    delta_max = max(t["absdelta"] for t in trials)
    bin_width = 2.0 * delta_max / NBINS
    for t in trials:
        t["bin"] = min(NBINS, math.floor((t["delta"] + delta_max) / bin_width) + 1)

    # Response-time trimming: one pass, full-sample mean and sample SD.
    rts = [t["rt"] for t in trials]
    rt_mean = mean(rts)
    rt_sd = math.sqrt(math.fsum((x - rt_mean) ** 2 for x in rts) / (len(rts) - 1))
    for t in trials:
        t["rt_keep"] = int(abs(t["rt"] - rt_mean) < 2.0 * rt_sd)

    # Choice sensitivity by maximum likelihood over all 435 trials.
    alpha = solve_alpha([(t["delta"], t["chose_r"], 1) for t in trials])
    for t in trials:
        p = choice_probability(alpha, t["delta"])
        t["pfit"] = p
        t["gfit"] = time_function(alpha, t["delta"])
        t["score"] = 2.0 * t["delta"] * (t["chose_r"] - p)
        t["curvature"] = -4.0 * t["delta"] ** 2 * p * (1.0 - p)
        t["loglik"] = math.log(choice_probability(alpha, t["sigma"] * t["delta"]))
        t["rt_used"] = t["rt"] if t["rt_keep"] else None
        t["g_used"] = t["gfit"] if t["rt_keep"] else None
    loglik = math.fsum(t["loglik"] for t in trials)

    bins = []
    for j in range(1, NBINS + 1):
        members = [t for t in trials if t["bin"] == j]
        kept = [t for t in members if t["rt_keep"]]
        if not members or not kept:
            sys.exit(f"bin {j} has no trials or no retained response times")
        ej = mean(t["pfit"] for t in members)
        bins.append({
            "bin": j, "n": len(members), "n_rt": len(kept),
            "n_right": sum(t["chose_r"] for t in members),
            "pj": mean(t["chose_r"] for t in members), "ej": ej,
            "mean_delta": mean(t["delta"] for t in members),
            "tj": mean(t["rt"] for t in kept), "gj": mean(t["gfit"] for t in kept),
            "g_all": mean(t["gfit"] for t in members),
            "lower": -delta_max + (j - 1) * bin_width, "upper": -delta_max + j * bin_width,
            "dj": math.log(ej / (1.0 - ej)) / (2.0 * alpha),
        })
    A, t0, r2, rmse, fit = ols([b["gj"] for b in bins], [b["tj"] for b in bins])
    for b, f in zip(bins, fit):
        b["tfit"] = f
        b["residual"] = b["tj"] - f
        b["choice_resid"] = b["pj"] - b["ej"]
        b["weighted_choice_sq"] = b["n"] * b["choice_resid"] ** 2
        b["choice_sq"] = b["choice_resid"] ** 2
        b["rt_curve"] = A * time_function(alpha, b["dj"]) + t0
        b["dot_curve_gap"] = abs(b["tfit"] - b["rt_curve"])
    A_all, t0_all, r2_all, _, _ = ols([b["g_all"] for b in bins], [b["tj"] for b in bins])
    return {
        "trials": trials, "bins": bins, "alpha": alpha, "ll": loglik, "n": len(trials),
        "delta_max": delta_max, "bin_width": bin_width,
        "rt_mean": rt_mean, "rt_sd": rt_sd, "rt_lower": rt_mean - 2 * rt_sd,
        "rt_upper": rt_mean + 2 * rt_sd, "rt_n": sum(t["rt_keep"] for t in trials),
        "A": A, "t0": t0, "r2": r2, "rmse": rmse,
        "choice_rmse": math.sqrt(math.fsum(b["weighted_choice_sq"] for b in bins) / len(trials)),
        "max_curve_gap": max(b["dot_curve_gap"] for b in bins),
        "A_all": A_all, "t0_all": t0_all, "r2_all": r2_all,
    }


def pooled_analysis():
    rows = read_csv(RAW / "popData.csv")
    if len(rows) != NBINS:
        sys.exit(f"expected {NBINS} class bins, found {len(rows)}")
    bins = []
    for j, r in enumerate(rows, start=1):
        nj, mj = int(r["nj"]), int(r["mj"])
        bins.append({"dj": float(r["dj"]), "nj": nj, "mj": mj, "tj": float(r["tj"]),
                     "bin": j, "n": nj + mj, "pj": nj / (nj + mj)})
    pool_n = sum(b["n"] for b in bins)
    alpha = solve_alpha([(b["dj"], b["nj"], b["n"]) for b in bins])
    for b in bins:
        b["ej"] = choice_probability(alpha, b["dj"])
        b["gj"] = time_function(alpha, b["dj"])
        b["loglik"] = b["nj"] * math.log(b["ej"]) + b["mj"] * math.log(1.0 - b["ej"])
        b["score"] = 2.0 * b["dj"] * (b["nj"] - b["n"] * b["ej"])
    A, t0, r2, rmse, fit = ols([b["gj"] for b in bins], [b["tj"] for b in bins])
    for b, f in zip(bins, fit):
        b["tfit"] = f
        b["residual"] = b["tj"] - f
        b["choice_resid"] = b["pj"] - b["ej"]
        b["weighted_choice_sq"] = b["n"] * b["choice_resid"] ** 2
        b["choice_sq"] = b["choice_resid"] ** 2
    return {
        "bins": bins, "alpha": alpha, "ll": math.fsum(b["loglik"] for b in bins), "n": pool_n,
        "A": A, "t0": t0, "r2": r2, "rmse": rmse,
        "choice_rmse": math.sqrt(math.fsum(b["weighted_choice_sq"] for b in bins) / pool_n),
    }


# --------------------------------------------------------------- outputs ----
TRIAL_COLUMNS = ["item1id", "item2id", "item1v", "item2v", "chosenid", "rt", "columbiaid",
                 "trial", "pair", "pair_tag", "delta", "chose_r", "sigma", "absdelta", "bin",
                 "rt_keep", "pfit", "gfit", "score", "curvature", "loglik", "rt_used", "g_used"]
IBIN_COLUMNS = ["bin", "n", "n_rt", "n_right", "pj", "ej", "mean_delta", "tj", "gj", "g_all",
                "lower", "upper", "dj", "tfit", "residual", "choice_resid",
                "weighted_choice_sq", "choice_sq", "rt_curve", "dot_curve_gap"]
PBIN_COLUMNS = ["dj", "nj", "mj", "tj", "bin", "n", "pj", "ej", "gj", "loglik", "score",
                "tfit", "residual", "choice_resid", "weighted_choice_sq", "choice_sq"]


def write_csv(path, rows, columns):
    def cell(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return repr(v)
        return str(v)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([cell(row[c]) for c in columns])


def f(value, width, decimals, comma=False):
    """Stata's %w.df and %w.dfc: right-justified, widened when needed."""
    return format(value, f"{width}{',' if comma else ''}.{decimals}f")


def braced(text):
    """Protect a thousands separator inside LaTeX math mode, as 05_tables.do does."""
    return text.strip().replace(",", "{,}")


def write_tables(ind, pool):
    tables = OUT / "tables"
    lines = [
        f"Choice trials & {f(ind['n'], 9, 0, True)} & {f(pool['n'], 9, 0, True)} \\\\",
        f"Choice sensitivity & {f(ind['alpha'], 8, 6)} & {f(pool['alpha'], 8, 6)} \\\\",
        f"Log likelihood & \\({f(ind['ll'], 12, 3)}\\) & \\({braced(f(pool['ll'], 12, 3, True))}\\) \\\\",
        f"Choice RMSE (percentage points) & {f(100 * ind['choice_rmse'], 6, 2)} & {f(100 * pool['choice_rmse'], 6, 2)} \\\\",
        "\\midrule",
        "Response-time bins & 15 & 15 \\\\",
        f"Time scale, $\\widehat{{A}}$ & {f(ind['A'], 10, 3, True)} & {f(pool['A'], 10, 3, True)} \\\\",
        f"Nondecision time, $\\widehat{{t}}_0$ (ms) & {f(ind['t0'], 10, 3, True)} & {f(pool['t0'], 10, 3, True)} \\\\",
        f"\\(R^2\\) of bin means & {f(ind['r2'], 6, 4)} & {f(pool['r2'], 6, 4)} \\\\",
        f"Response-time RMSE (ms) & {f(ind['rmse'], 8, 2)} & {f(pool['rmse'], 8, 2)} \\\\",
        "\\bottomrule",
    ]
    (tables / "tab-estimates.tex").write_text("\n".join(lines) + "\n")

    rows = [f"{f(b['bin'], 2, 0)} & \\({f(b['dj'], 6, 3)}\\) & {f(b['n'], 4, 0)} & "
            f"{f(b['n_rt'], 4, 0)} & {f(b['pj'], 5, 3)} & {f(b['ej'], 5, 3)} & "
            f"{f(b['gj'], 6, 4)} & {f(b['tj'], 7, 1)} & {f(b['tfit'], 7, 1)} \\\\"
            for b in ind["bins"]]
    (tables / "tab-individual-bins.tex").write_text("\n".join(rows + ["\\bottomrule"]) + "\n")

    rows = [f"{f(b['bin'], 2, 0)} & \\({f(b['dj'], 6, 4)}\\) & {f(b['n'], 5, 0, True)} & "
            f"{f(b['pj'], 5, 3)} & {f(b['ej'], 5, 3)} & {f(b['gj'], 6, 4)} & "
            f"{f(b['tj'], 7, 1)} & {f(b['tfit'], 7, 1)} \\\\"
            for b in pool["bins"]]
    (tables / "tab-pooled-bins.tex").write_text("\n".join(rows + ["\\bottomrule"]) + "\n")

    macros = [
        ("AlphaIndividual", f(ind["alpha"], 5, 3)), ("AlphaPooled", f(pool["alpha"], 5, 3)),
        ("LLIndividual", f(ind["ll"], 8, 2)), ("LLPooled", braced(f(pool["ll"], 12, 2, True))),
        ("TimeScaleIndividual", f(ind["A"], 5, 0)), ("OffsetIndividual", f(ind["t0"], 4, 0)),
        ("TimeScalePooled", f(pool["A"], 5, 0)), ("OffsetPooled", f(pool["t0"], 4, 0)),
        ("RSquaredIndividual", f(ind["r2"], 5, 3)), ("RSquaredPooled", f(pool["r2"], 5, 3)),
        ("RTMean", f(ind["rt_mean"], 6, 0, True)), ("RTSD", f(ind["rt_sd"], 6, 0, True)),
        ("RTLower", f(ind["rt_lower"], 8, 2)), ("RTUpper", f(ind["rt_upper"], 8, 2)),
        ("RTKept", f(ind["rt_n"], 3, 0)), ("RTRemoved", f(ind["n"] - ind["rt_n"], 3, 0)),
        ("MaxDelta", f(ind["delta_max"], 4, 2)), ("BinWidth", f(ind["bin_width"], 5, 3)),
        ("ChoiceRMSEIndividual", f(100 * ind["choice_rmse"], 4, 1)),
        ("ChoiceRMSEPooled", f(100 * pool["choice_rmse"], 4, 1)),
        ("RTRMSEIndividual", f(ind["rmse"], 4, 0)), ("RTRMSEPooled", f(pool["rmse"], 4, 0)),
        ("CurveGapIndividual", f(ind["max_curve_gap"], 4, 1)),
        ("AltAIndividual", f(ind["A_all"], 5, 0)), ("AltOffsetIndividual", f(ind["t0_all"], 4, 0)),
        ("AltRSquaredIndividual", f(ind["r2_all"], 5, 3)),
    ]
    (tables / "numeric-values.tex").write_text(
        "".join(f"\\newcommand{{\\{name}}}{{{value}}}\n" for name, value in macros))


def draw_figures(ind, pool):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Latin Modern Roman", "Times New Roman", "DejaVu Serif"],
        "pdf.fonttype": 3, "axes.unicode_minus": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.labelsize": 9,
        "legend.fontsize": 8, "legend.frameon": False,
    })
    metadata = {"Creator": "replicate.py", "Producer": "matplotlib",
                "CreationDate": None, "ModDate": None}
    figures = OUT / "figures"

    def frame(xticks, xlabel, ylabel):
        fig, ax = plt.subplots(figsize=(6.5, 4.0))
        fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.25)
        ax.set_xticks(xticks)
        ax.set_xlabel(xlabel, labelpad=6)
        ax.set_ylabel(ylabel)
        ax.tick_params(length=3, width=0.6)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_linewidth(0.6)
        return fig, ax

    def finish(fig, ax, handles, labels, path):
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.17),
                  ncol=3, handletextpad=0.4, columnspacing=2.0)
        fig.savefig(path, format="pdf", metadata=metadata)
        plt.close(fig)

    curve_style = dict(color="black", linewidth=0.8)
    dot_style = dict(linestyle="none", marker="o", markersize=3.6, color="black")
    circle_style = dict(linestyle="none", marker="o", markersize=6.5, markerfacecolor="none",
                        markeredgecolor="black", markeredgewidth=0.8)
    xlabel = "Value difference, right minus left ($)"

    for name, res, bound, step in (("individual", ind, ind["delta_max"], 2),
                                   ("pooled", pool, 10.0, 2)):
        alpha, bins = res["alpha"], res["bins"]
        grid = [-bound + 2.0 * bound * i / 1000 for i in range(1001)]
        ticks = list(range(-int(bound), int(bound) + 1, step))
        dj = [b["dj"] for b in bins]

        fig, ax = frame(ticks, xlabel, "Frequency of choosing right")
        curve, = ax.plot(grid, [choice_probability(alpha, d) for d in grid], **curve_style)
        dots, = ax.plot(dj, [b["ej"] for b in bins], **dot_style)
        circles, = ax.plot(dj, [b["pj"] for b in bins], **circle_style)
        ax.yaxis.set_major_locator(FixedLocator([0, 0.2, 0.4, 0.6, 0.8, 1.0]))
        ax.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: "0" if v == 0 else "1" if v == 1 else f"{v:.1f}".lstrip("0")))
        ax.set_ylim(-0.03, 1.03)
        finish(fig, ax, [circles, dots, curve],
               ["Observed frequency", "DDM bin prediction", "DDM curve"],
               figures / f"fig-{name}-choice.pdf")

        fig, ax = frame(ticks, xlabel, "Mean response time (ms)")
        curve, = ax.plot(grid, [res["A"] * time_function(alpha, d) + res["t0"] for d in grid],
                         **curve_style)
        dots, = ax.plot(dj, [b["tfit"] for b in bins], **dot_style)
        circles, = ax.plot(dj, [b["tj"] for b in bins], **circle_style)
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        finish(fig, ax, [circles, dots, curve],
               ["Observed mean", "DDM bin prediction", "DDM curve"],
               figures / f"fig-{name}-response-time.pdf")


# ------------------------------------------------------------ comparison ----
def compare_with_stata():
    """Check the Python outputs against Stata's, when those are present."""
    problems = 0
    checked = 0
    for name in ("numeric-values.tex", "tab-estimates.tex",
                 "tab-individual-bins.tex", "tab-pooled-bins.tex"):
        stata = ROOT / "output" / "tables" / name
        if stata.exists():
            checked += 1
            same = stata.read_bytes() == (OUT / "tables" / name).read_bytes()
            problems += not same
            print(f"  {name}: {'byte-identical to Stata' if same else 'DIFFERS from Stata'}")
    for name, skip in (("individual-bins.csv", ()), ("pooled-bins.csv", ()),
                       ("individual-trials.csv", ("pair_tag",))):
        stata = ROOT / "data" / "processed" / name
        if not stata.exists():
            continue
        checked += 1
        a, b = read_csv(stata), read_csv(OUT / "data" / name)
        worst = 0.0
        for ra, rb in zip(a, b):
            for key in ra:
                if key in skip:
                    continue
                try:
                    x, y = float(ra[key] or "nan"), float(rb[key] or "nan")
                except ValueError:
                    worst = max(worst, 0.0 if ra[key] == rb[key] else math.inf)
                    continue
                if math.isnan(x) and math.isnan(y):
                    continue
                worst = max(worst, abs(x - y) / max(1.0, abs(x)))
        ok = len(a) == len(b) and worst < 1e-9
        problems += not ok
        print(f"  {name}: {'agrees with Stata' if ok else 'DIFFERS from Stata'} "
              f"(largest relative difference {worst:.1e})")
    if not checked:
        print("  no Stata outputs found to compare against")
    return problems


def main():
    for sub in ("data", "tables", "figures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    ind = individual_analysis()
    pool = pooled_analysis()
    write_csv(OUT / "data" / "individual-trials.csv", ind["trials"], TRIAL_COLUMNS)
    write_csv(OUT / "data" / "individual-bins.csv", ind["bins"], IBIN_COLUMNS)
    write_csv(OUT / "data" / "pooled-bins.csv", pool["bins"], PBIN_COLUMNS)
    write_tables(ind, pool)
    draw_figures(ind, pool)
    print(f"Individual: alpha={ind['alpha']:.6f}  LL={ind['ll']:.3f}  A={ind['A']:.3f}  "
          f"t0={ind['t0']:.3f}  R2={ind['r2']:.4f}  ({ind['rt_n']} of {ind['n']} times retained)")
    print(f"Pooled:     alpha={pool['alpha']:.6f}  LL={pool['ll']:.3f}  A={pool['A']:.3f}  "
          f"t0={pool['t0']:.3f}  R2={pool['r2']:.4f}")
    print(f"Outputs written under {OUT.relative_to(ROOT)}/")
    print("Comparison with Stata outputs:")
    problems = compare_with_stata()
    if problems:
        sys.exit(f"{problems} comparison(s) failed")
    print("Done.")


if __name__ == "__main__":
    main()
