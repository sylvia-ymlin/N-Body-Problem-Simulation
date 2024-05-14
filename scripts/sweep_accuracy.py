"""
Sweep theta values and measure pos_maxdiff against v1 (exact) reference.
Usage: python3 scripts/sweep_accuracy.py
"""
import json
import os
import subprocess
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from benchmark_suite import PYTHON_BIN
from benchmark_full_suite import ensure_input

BINARY = "./build/nbody_simulate"
N = 2000
NSTEPS = 200
DT = 1e-5
DATA_FILE = f"data/inputs/accuracy_{N}.gal"
OUT_DIR = "data/outputs"
THETAS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
TOLERANCE = 1e-3


def run(version, theta):
    out = f"{OUT_DIR}/result_v{version}.gal"
    cmd = [BINARY, str(version), str(N), DATA_FILE, str(NSTEPS), str(DT), "1", str(theta), "0"]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0:
        print(f"  ERROR v{version} theta={theta}: {r.stderr.decode()[:200]}")
        return None
    return out


def read_gal(path):
    size = os.path.getsize(path)
    n = size // 40
    return np.fromfile(path, dtype=np.float64).reshape(n, 5)


def pos_maxdiff(ref_path, test_path):
    d1 = read_gal(ref_path)
    d2 = read_gal(test_path)
    idx1 = np.argsort(d1[:, 2])
    idx2 = np.argsort(d2[:, 2])
    d1s = d1[idx1]
    d2s = d2[idx2]
    dist = np.sqrt((d1s[:, 0] - d2s[:, 0])**2 + (d1s[:, 1] - d2s[:, 1])**2)
    return float(np.max(dist)), float(np.mean(dist))


def main():
    ensure_input(DATA_FILE, N, "disk")
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"Generating reference solution (v1, theta=0, N={N}, nsteps={NSTEPS}, dt={DT})...")
    ref = run(1, 0.0)
    if ref is None:
        sys.exit(1)

    results = []
    print(f"\nTheta sweep (v4, N={N}, comparing against v1):")
    print(f"{'theta':>8} | {'max_diff':>12} | {'mean_diff':>12} | {'pass':>6}")
    print("-" * 48)

    for theta in THETAS:
        test = run(4, theta)
        if test is None:
            continue
        mx, mn = pos_maxdiff(ref, test)
        passed = mx < TOLERANCE
        print(f"{theta:>8.2f} | {mx:>12.2e} | {mn:>12.2e} | {'✓' if passed else '✗':>6}")
        results.append({"theta": theta, "max_diff": mx, "mean_diff": mn, "pass": passed})

    passing = [r for r in results if r["pass"]]
    if passing:
        best = min(passing, key=lambda r: r["theta"])
        print(f"\nOptimal theta (max_diff < {TOLERANCE}): {best['theta']} (max_diff={best['max_diff']:.2e})")
    else:
        print(f"\nNo theta passes tolerance {TOLERANCE}")

    out_path = "data/metrics/sweep_accuracy.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=4)
    print(f"Results saved to {out_path}")


if __name__ == "__main__":
    main()
