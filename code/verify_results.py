"""Independent numerical audit, using only Python's standard library.

This script reads the original CSVs and checks the final Stata-generated
results. It does not author empirical figures, tables, or replacement results.
Run from the project root: python3 code/verify_results.py.
"""

from bisect import bisect_right
import csv
from math import exp, isclose, log, sqrt, tanh
from pathlib import Path
from statistics import mean, stdev


def read_csv(path):
    with Path(path).open(newline="") as stream:
        return [
            {key.strip(): value.strip() for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]


def probability(alpha, delta):
    return 1 / (1 + exp(-2 * alpha * delta))


def time_function(alpha, delta):
    return tanh(alpha * delta) / delta if delta else alpha


def close(actual, expected, tol=1e-7):
    assert isclose(float(actual), expected, rel_tol=tol, abs_tol=tol), (actual, expected)


def estimate_alpha(observations):
    # observations: (delta, right count, total count). Solve the score by bisection.
    def score(alpha):
        return sum(2 * d * (right - n * probability(alpha, d)) for d, right, n in observations)

    lower, upper = 0.0, 1.0
    assert score(lower) > 0
    while score(upper) > 0:
        upper *= 2
    for _ in range(100):
        middle = (lower + upper) / 2
        if score(middle) > 0:
            lower = middle
        else:
            upper = middle
    alpha = (lower + upper) / 2
    close(score(alpha), 0, 1e-8)
    return alpha


def ols(x, y):
    xm, ym = mean(x), mean(y)
    slope = sum((a - xm) * (b - ym) for a, b in zip(x, y)) / sum((a - xm) ** 2 for a in x)
    intercept = ym - slope * xm
    fit = [intercept + slope * a for a in x]
    rss = sum((a - b) ** 2 for a, b in zip(y, fit))
    r_squared = 1 - rss / sum((a - ym) ** 2 for a in y)
    return slope, intercept, r_squared, sqrt(rss / len(y)), fit


raw = [r for r in read_csv("data/raw/indv-Data.csv") if r["columbiaID"] == "ts3686"]
assert len(raw) == 435
assert len({tuple(r.values()) for r in raw}) == 435
assert len({tuple(sorted((r["item1ID"], r["item2ID"]))) for r in raw}) == 368
assert len({r[k] for r in raw for k in ("item1ID", "item2ID")}) == 30
delta = [float(r["item2V"]) - float(r["item1V"]) for r in raw]
choice = [int(r["chosenID"] == r["item2ID"]) for r in raw]
rt = [float(r["RT"]) for r in raw]
alpha = estimate_alpha(list(zip(delta, choice, [1] * len(raw))))
cutoff_mean, cutoff_sd = mean(rt), stdev(rt)
keep = [abs(t - cutoff_mean) < 2 * cutoff_sd for t in rt]
assert sum(keep) == 418
assert sum(d == 0 for d in delta) == 1
extent = max(map(abs, delta))
width = 2 * extent / 15
edges = [-extent + j * width for j in range(16)]
groups = [min(14, bisect_right(edges, d) - 1) for d in delta]
bins = read_csv("data/processed/individual-bins.csv")
assert len(bins) == 15
assert sum(int(r["n"]) for r in bins) == 435
assert sum(int(r["n_rt"]) for r in bins) == 418
gs, ts = [], []
for j, result in enumerate(bins):
    indexes = [i for i, group in enumerate(groups) if group == j]
    retained = [i for i in indexes if keep[i]]
    expected = mean(probability(alpha, delta[i]) for i in indexes)
    g = mean(time_function(alpha, delta[i]) for i in retained)
    t = mean(rt[i] for i in retained)
    close(result["n"], len(indexes))
    close(result["n_rt"], len(retained))
    close(result["pj"], mean(choice[i] for i in indexes))
    close(result["ej"], expected)
    close(result["dj"], log(expected / (1 - expected)) / (2 * alpha))
    close(result["gj"], g)
    close(result["tj"], t)
    gs.append(g)
    ts.append(t)
A, t0, r2, rmse, fit = ols(gs, ts)
for r, value in zip(bins, fit):
    close(r["tfit"], value)
print(f"Individual audit passed: alpha={alpha:.9f}, A={A:.6f}, t0={t0:.6f}, R2={r2:.9f}, RMSE={rmse:.6f}")

pooled = read_csv("data/raw/popData.csv")
triples = [(float(r["dj"]), int(r["nj"]), int(r["nj"]) + int(r["mj"])) for r in pooled]
assert sum(n for _, _, n in triples) == 30450
alpha = estimate_alpha(triples)
bins = read_csv("data/processed/pooled-bins.csv")
gs, ts = [], []
for source, result, (d, right, n) in zip(pooled, bins, triples):
    g = time_function(alpha, d)
    close(result["pj"], right / n)
    close(result["ej"], probability(alpha, d))
    close(result["gj"], g)
    gs.append(g)
    ts.append(float(source["tj"]))
A, t0, r2, rmse, fit = ols(gs, ts)
for result, value in zip(bins, fit):
    close(result["tfit"], value)
print(f"Pooled audit passed: alpha={alpha:.9f}, A={A:.6f}, t0={t0:.6f}, R2={r2:.9f}, RMSE={rmse:.6f}")
print("All trial counts, bin means, inverse coordinates, likelihood fits, and response-time predictions agree with Stata.")
