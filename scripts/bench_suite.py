import subprocess
import time
import os
import json
import sys

def run_bench(version, N, filename, nsteps, n_threads, theta=0.5, k=0):
    cmd = [
        "./build/nbody_simulate",
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
    result = subprocess.run(cmd, capture_output=True, text=True)
    end = time.time()
    
    if result.returncode != 0:
        return None
    
    lines = result.stdout.splitlines()
    sim_time = end - start
    force_time = 0.0
    for line in lines:
        if "Simulation Complete:" in line:
            sim_time = float(line.split(":")[1].replace("s", "").strip())
        if "[HP] Tree Build:" in line:
            parts = line.split("|")
            force_time = float(parts[1].split(":")[1].replace("s", "").strip())
            
    return {
        "runtime": sim_time,
        "force_time": force_time
    }

def run_comprehensive_suite():
    # 1. Strong Scaling: N=100,000, Threads=1, 2, 4, 8, 16
    strong_scaling_n = 100000
    threads_list = [1, 2, 4, 8, 16]
    strong_results = []
    
    print("Starting Strong Scaling Benchmarks (N=100k)...")
    data_file = f"data/bench_strong_{strong_scaling_n}.gal"
    if not os.path.exists(data_file):
        subprocess.run(["python3", "scripts/generate_data.py", str(strong_scaling_n), data_file, "disk"], capture_output=True)
        
    for t in threads_list:
        print(f"  Threads={t}...")
        res = run_bench(5, strong_scaling_n, data_file, 1, t)
        if res:
            res.update({"threads": t, "N": strong_scaling_n})
            strong_results.append(res)

    # 2. Weak Scaling: N_base=100,000, N = N_base * threads
    weak_results = []
    print("\nStarting Weak Scaling Benchmarks (N=100k * threads)...")
    for t in threads_list:
        n_weak = strong_scaling_n * t
        print(f"  Threads={t}, N={n_weak}...")
        weak_file = f"data/bench_weak_{n_weak}.gal"
        if not os.path.exists(weak_file):
            subprocess.run(["python3", "scripts/generate_data.py", str(n_weak), weak_file, "disk"], capture_output=True)
        
        res = run_bench(5, n_weak, weak_file, 1, t)
        if res:
            res.update({"threads": t, "N": n_weak})
            weak_results.append(res)

    final_results = {
        "strong": strong_results,
        "weak": weak_results
    }
    
    with open("xeon_results.json", "w") as f:
        json.dump(final_results, f, indent=4)
    
    print("\nResults saved to xeon_results.json")

if __name__ == "__main__":
    if not os.path.exists("data"): os.makedirs("data")
    run_comprehensive_suite()
