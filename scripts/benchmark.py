import subprocess
import time
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

# Configuration
N_VALUES = [1000, 2000, 5000, 10000, 20000, 50000]
NSTEPS = 20
THETA = 0.5
NUM_RUNS = 3
DATA_FILE = "bench_data.bin"

# Pre-recorded data from user to avoid re-running long benchmarks
LATEST_RESULTS = {
    "v1_naive_serial": [0.0362, 0.1264, 0.7495, 2.9942, 11.1612, 67.0591],
    "v1_naive_parallel": [0.0167, 0.0407, 0.2250, 0.8707, 3.2780, 21.3802],
    "v2_barnes_hut": [0.0305, 0.0681, 0.1846, 0.4329, 0.9636, 3.0551],
    "v3_arena": [0.0278, 0.0642, 0.1724, 0.4108, 0.9106, 2.9420],
    "v4_morton": [0.0238, 0.0526, 0.1821, 0.3810, 0.7653, 2.5314],
    "v5_parallel_4t": [0.0217, 0.0355, 0.0801, 0.1553, 0.2810, 0.8144],
    "v5_parallel_8t": [0.0231, 0.0368, 0.0754, 0.1421, 0.2393, 0.7011]
}

# Version definitions
VERSIONS = {
    "v1_naive_serial": "env OMP_NUM_THREADS=1 ./build/v1_naive {N} {file} {steps}",
    "v1_naive_parallel": "env OMP_NUM_THREADS=4 ./build/v1_naive {N} {file} {steps}",
    "v2_barnes_hut": "./build/v2_barnes_hut {N} {file} {steps} {theta}",
    "v3_arena": "./build/v3_arena {N} {file} {steps} {theta}",
    "v4_morton": "./build/v4_morton {N} {file} {steps} {theta}",
    "v5_parallel_4t": "./build/v5_parallel {N} {file} {steps} 0.001 4 {theta} 8",
    "v5_parallel_8t": "./build/v5_parallel {N} {file} {steps} 0.001 8 {theta} 8"
}

def run_command(cmd):
    start = time.time()
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        end = time.time()
        return end - start
    except Exception as e:
        print(f"Error executing command: {e}")
        return None

def run_benchmark():
    results = {name: [] for name in VERSIONS}
    for N in N_VALUES:
        print(f"\n--- Benchmarking N = {N} ---")
        subprocess.run(["python3", "scripts/generate_data.py", str(N), DATA_FILE], check=True, capture_output=True)
        for name, cmd_template in VERSIONS.items():
            cmd = cmd_template.format(N=N, file=DATA_FILE, steps=NSTEPS, theta=THETA)
            durations = []
            for i in range(NUM_RUNS):
                dur = run_command(cmd)
                if dur is not None: durations.append(dur)
            min_duration = min(durations) if durations else None
            results[name].append(min_duration)
            if min_duration: print(f"{name:18}: {min_duration:.4f}s")
    if os.path.exists(DATA_FILE): os.remove(DATA_FILE)
    return results

def plot_all_in_one(results):
    os.makedirs("docs", exist_ok=True)
    plt.style.use('seaborn-v0_8-muted')
    
    # Create a two-panel figure: Absolute Time and Speedup
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), dpi=120)
    
    colors = {
        "v1_naive_serial": "#d8dee9", "v1_naive_parallel": "#4c566a",
        "v2_barnes_hut": "#88c0d0", "v3_arena": "#81a1c1", "v4_morton": "#5e81ac",
        "v5_parallel_4t": "#d08770", "v5_parallel_8t": "#bf616a"
    }
    
    # Panel 1: Execution Time (Log-Log)
    for name, times in results.items():
        valid_indices = [i for i, t in enumerate(times) if t is not None]
        if not valid_indices: continue
        x = [N_VALUES[i] for i in valid_indices]
        y = [times[i] for i in valid_indices]
        
        # Style naive versions differently to reduce clutter
        is_naive = "naive" in name
        alpha = 0.4 if is_naive else 0.9
        lw = 1.5 if is_naive else 3.0
        ls = "--" if is_naive else "-"
        
        ax1.plot(x, y, label=name.replace("_", " ").title(), 
                 color=colors.get(name), linestyle=ls, marker='o', 
                 markersize=6, linewidth=lw, alpha=alpha, zorder=2 if is_naive else 5)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Number of Particles (N)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Execution Time (seconds)', fontsize=12, fontweight='bold')
    ax1.set_title('Absolute Performance (Log-Log Scale)', fontsize=14, fontweight='bold')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend(loc='upper left', frameon=True, fontsize=9)

    # Panel 2: Speedup vs v1_naive_serial (Linear)
    baseline = results.get("v1_naive_serial")
    if baseline:
        # Plot the baseline (v1_naive_serial) explicitly at y=1
        ax2.axhline(y=1, color=colors["v1_naive_serial"], linestyle="--", alpha=0.6, label="V1 Naive Serial (Baseline)")
        
        for name, times in results.items():
            if name == "v1_naive_serial": continue
            valid_indices = [i for i, t in enumerate(times) if t is not None and baseline[i] is not None]
            if not valid_indices: continue
            
            x = [N_VALUES[i] for i in valid_indices]
            speedup = [baseline[i] / times[i] for i in valid_indices]
            
            # v1_naive_parallel now with markers and will be visible on Log scale
            is_naive_para = name == "v1_naive_parallel"
            ax2.plot(x, speedup, label=name.replace("_", " ").title(), 
                     color=colors.get(name), marker='s', 
                     markersize=6, linewidth=2.5, alpha=0.9,
                     linestyle='-' if not is_naive_para else '--')

        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel('Number of Particles (N)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Speedup (x-times faster than Naive Serial)', fontsize=12, fontweight='bold')
        ax2.set_title('Relative Speedup Evolution (Log Scale)', fontsize=14, fontweight='bold')
        
        # Set nice looking y-ticks for speedup
        from matplotlib.ticker import ScalarFormatter
        ax2.yaxis.set_major_formatter(ScalarFormatter())
        ax2.set_yticks([1, 2, 5, 10, 20, 50, 100])
        
        ax2.grid(True, which="both", alpha=0.3)
        ax2.legend(loc='upper left', frameon=True, fontsize=9)
        
        # Highlight the 100x line if reached
        ax2.axhline(y=100, color='r', linestyle=':', alpha=0.5)

    fig.suptitle('N-Body Simulation Performance Analysis', 
                 fontsize=18, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('docs/performance_comparison.png', dpi=300, bbox_inches='tight')
    print("\nEnhanced dual-panel plot saved to docs/performance_comparison.png")

if __name__ == "__main__":
    if "--plot-only" in sys.argv:
        print("Using pre-recorded data for visualization...")
        plot_all_in_one(LATEST_RESULTS)
    else:
        results = run_benchmark()
        plot_all_in_one(results)
