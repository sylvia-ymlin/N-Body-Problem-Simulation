import subprocess
import time
import os
import json
import sys
import argparse

PYTHON_BIN = sys.executable or "python3"

def run_bench(binary, version, N, filename, nsteps, n_threads, theta=0.5, k=0, extra_env=None):
    cmd = [
        binary,
        str(version),
        str(N),
        filename,
        str(nsteps),
        "0.001", # dt
        str(n_threads),
        str(theta),
        str(k)
    ]
    start = time.time()
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    end = time.time()
    
    if result.returncode != 0:
        return None
    
    lines = result.stdout.splitlines()
    # Fallback to wall time if the binary output format changes.
    sim_time = end - start
    force_time = 0.0
    tree_time = 0.0
    kmeans_time = 0.0
    schedule = None
    for line in lines:
        if "Simulation Complete:" in line:
            sim_time = float(line.split(":")[1].replace("s", "").strip())
        # Hotspot line is printed by v5 once per process.
        if "[HP] Schedule:" in line:
            parts = line.split("|")
            schedule = parts[0].split(":", 1)[1].strip()
            kmeans_time = float(parts[1].split(":")[1].replace("s", "").strip())
            tree_time = float(parts[2].split(":")[1].replace("s", "").strip())
            force_time = float(parts[3].split(":")[1].replace("s", "").strip())
            
    return {
        "runtime": sim_time,
        "force_time": force_time,
        "tree_time": tree_time,
        "kmeans_time": kmeans_time,
        "schedule": schedule,
    }

def parse_threads(raw):
    return [int(part) for part in raw.split(",") if part.strip()]

def run_comprehensive_suite(binary="./build/nbody_simulate",
                            output_path="data/metrics/scaling_results.json",
                            strong_scaling_n=100000,
                            threads_list=None,
                            nsteps=1,
                            theta=0.5):
    # 1. Strong Scaling: N=100,000, Threads=1, 2, 4, 8, 16
    threads_list = threads_list or [1, 2, 4, 8, 16]
    strong_results = []
    
    print(f"Starting Strong Scaling Benchmarks (N={strong_scaling_n})...")
    data_file = f"data/inputs/bench_strong_{strong_scaling_n}.gal"
    if not os.path.exists(data_file):
        subprocess.run([PYTHON_BIN, "scripts/generate_data.py", str(strong_scaling_n), data_file, "disk"], capture_output=True, check=True)
        
    for t in threads_list:
        print(f"  Threads={t}...")
        res = run_bench(binary, 5, strong_scaling_n, data_file, nsteps, t, theta)
        if res:
            res.update({"threads": t, "N": strong_scaling_n})
            strong_results.append(res)

    # 2. Weak Scaling: N_base=100,000, N = N_base * threads
    weak_results = []
    print("\nStarting Weak Scaling Benchmarks (N=100k * threads)...")
    for t in threads_list:
        n_weak = strong_scaling_n * t
        print(f"  Threads={t}, N={n_weak}...")
        weak_file = f"data/inputs/bench_weak_{n_weak}.gal"
        if not os.path.exists(weak_file):
            subprocess.run([PYTHON_BIN, "scripts/generate_data.py", str(n_weak), weak_file, "disk"], capture_output=True, check=True)
        
        res = run_bench(binary, 5, n_weak, weak_file, nsteps, t, theta)
        if res:
            res.update({"threads": t, "N": n_weak})
            weak_results.append(res)

    final_results = {
        "strong": strong_results,
        "weak": weak_results
    }
    
    with open(output_path, "w") as f:
        json.dump(final_results, f, indent=4)
    
    print(f"\nResults saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./build/nbody_simulate")
    parser.add_argument("--output", default="data/metrics/scaling_results.json")
    parser.add_argument("--threads", default="1,2,4,8,16")
    parser.add_argument("--n-base", type=int, default=100000)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--theta", type=float, default=0.5)
    args = parser.parse_args()

    if not os.path.exists("data"):
        os.makedirs("data")

    run_comprehensive_suite(binary=args.binary,
                            output_path=args.output,
                            strong_scaling_n=args.n_base,
                            threads_list=parse_threads(args.threads),
                            nsteps=args.steps,
                            theta=args.theta)
