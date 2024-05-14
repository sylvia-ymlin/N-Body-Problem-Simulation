import json
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
from benchmark_suite import run_bench
from benchmark_full_suite import ensure_input

def sweep_theta(binary="./build/nbody_simulate",
                output_path="data/metrics/sweep_theta.json"):
    N = 100000
    data_file = "data/inputs/bench_strong_100000.gal"
    ensure_input(data_file, N, "disk")

    thetas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    results = []

    print(f"Theta sweep (v4, N={N}, serial)...")
    for theta in thetas:
        print(f"  theta={theta}...")
        res = run_bench(binary, 4, N, data_file, 1, 1, theta, 0)
        if res:
            res.update({"theta": theta, "N": N})
            results.append(res)
            print(f"    -> {res['runtime']:.3f}s")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)
    print(f"\nResults saved to {output_path}")

    best = min(results, key=lambda r: r["runtime"])
    print(f"Optimal theta: {best['theta']} -> {best['runtime']:.3f}s")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./build/nbody_simulate")
    parser.add_argument("--output", default="data/metrics/sweep_theta.json")
    args = parser.parse_args()
    sweep_theta(args.binary, args.output)
