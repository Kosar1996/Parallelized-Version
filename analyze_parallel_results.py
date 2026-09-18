"""
Post-processing for the 1/2/4/8/16-processor correctness + benchmark runs.

Reads out_parallel_1timestep_N<P>.mat for P in {1,2,4,8,16}, compares each
against the N=1 baseline for correctness, and builds the speedup curve
(wall-clock time and speedup vs. processor count).

Usage: run from this folder once all 5 out_parallel_1timestep_N*.mat files
have been pulled back from the cluster:
    python3 analyze_parallel_results.py
"""
import scipy.io as sio
import numpy as np
import os

PROC_COUNTS = [1, 2, 4, 8, 16]
HERE = os.path.dirname(os.path.abspath(__file__))

def load_result(n):
    path = os.path.join(HERE, f"out_parallel_1timestep_N{n}.mat")
    if not os.path.exists(path):
        return None
    d = sio.loadmat(path, struct_as_record=False, squeeze_me=True)
    out = d["out"]
    wallclock = float(d["wallClockSeconds"])
    return out, wallclock

def get_field(out, name):
    return getattr(out, name)

results = {}
for n in PROC_COUNTS:
    r = load_result(n)
    if r is None:
        print(f"N={n}: out_parallel_1timestep_N{n}.mat not found -- skipping")
        continue
    results[n] = r

if 1 not in results:
    print("No N=1 baseline found -- cannot verify correctness or compute speedup. Stopping.")
    raise SystemExit(1)

out1, t1 = results[1]
stopStep1 = int(get_field(out1, "stopStep"))
finalT1 = float(np.atleast_1d(get_field(out1, "t"))[-1])
p2DMax1 = get_field(out1, "p2DMaxHist")
p2DMax1 = float(np.atleast_1d(p2DMax1)[-1]) if p2DMax1 is not None else None

print("=" * 70)
print("CORRECTNESS VERIFICATION (each N vs. N=1 baseline)")
print("=" * 70)
print(f"{'N':>4} {'stopStep':>10} {'final t match':>15} {'max|P2D| match':>16} {'max|P2D| diff':>15}")

correctness_rows = []
for n in PROC_COUNTS:
    if n not in results:
        continue
    out_n, t_n = results[n]
    stopStep_n = int(get_field(out_n, "stopStep"))
    finalT_n = float(np.atleast_1d(get_field(out_n, "t"))[-1])
    p2DMax_n = get_field(out_n, "p2DMaxHist")
    p2DMax_n = float(np.atleast_1d(p2DMax_n)[-1]) if p2DMax_n is not None else None

    t_match = "IDENTICAL" if finalT_n == finalT1 else f"DIFFERS ({finalT_n} vs {finalT1})"
    if p2DMax1 is not None and p2DMax_n is not None:
        p_diff = abs(p2DMax_n - p2DMax1)
        p_match = "IDENTICAL" if p_diff == 0 else ("within tol" if p_diff < 1e-6 else "DIFFERS")
    else:
        p_diff = float("nan")
        p_match = "N/A"

    print(f"{n:>4} {stopStep_n:>10} {t_match:>15} {p_match:>16} {p_diff:>15.6e}")
    correctness_rows.append((n, stopStep_n, finalT_n, p2DMax_n, p_diff))

print()
print("=" * 70)
print("SCALING BENCHMARK")
print("=" * 70)
print(f"{'N':>4} {'wall-clock (s)':>16} {'speedup vs N=1':>16}")
benchmark_rows = []
for n in PROC_COUNTS:
    if n not in results:
        continue
    _, t_n = results[n]
    speedup = t1 / t_n
    print(f"{n:>4} {t_n:>16.3f} {speedup:>16.3f}")
    benchmark_rows.append((n, t_n, speedup))

# Save a CSV for the report / plotting
import csv
csv_path = os.path.join(HERE, "parallel_benchmark_results.csv")
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["N", "wallClockSeconds", "speedup", "stopStep", "finalT", "maxP2D", "maxP2D_diff_vs_N1"])
    bench_by_n = {n: (t, s) for n, t, s in benchmark_rows}
    corr_by_n = {n: (ss, ft, p, pd) for n, ss, ft, p, pd in correctness_rows}
    for n in PROC_COUNTS:
        if n not in results:
            continue
        t_n, s_n = bench_by_n[n]
        ss, ft, p, pd = corr_by_n[n]
        w.writerow([n, t_n, s_n, ss, ft, p, pd])

print(f"\nSaved: {csv_path}")
