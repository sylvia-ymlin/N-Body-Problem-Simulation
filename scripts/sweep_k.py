import json
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from benchmark_suite import run_bench
from benchmark_full_suite import ensure_input

def sweep_k(binary="./build/nbody_simulate",
            output_path="data/metrics/sweep_k.json",
            theta=0.5):
    N = 100000
    data_file = "data/inputs/bench_strong_100000.gal"
    ensure_input(data_file, N, "disk")

    ks = [1, 2, 4, 8, 16, 32, 64]
    results = []

    print(f"K sweep (v4, N={N}, serial, theta={theta})...")
    for k in ks:
        print(f"  k={k}...")
        res = run_bench(binary, 4, N, data_file, 1, 1, theta, k)
        if res:
            res.update({"k": k, "N": N, "theta": theta})
            results.append(res)
            print(f"    -> {res['runtime']:.3f}s")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)
    print(f"\nResults saved to {output_path}")

    best = min(results, key=lambda r: r["runtime"])
    print(f"Optimal k: {best['k']} -> {best['runtime']:.3f}s")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./build/nbody_simulate")
    parser.add_argument("--output", default="data/metrics/sweep_k.json")
    parser.add_argument("--theta", type=float, default=0.5)
    args = parser.parse_args()
    sweep_k(args.binary, args.output, args.theta)
