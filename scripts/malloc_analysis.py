#!/usr/bin/env python3
"""
Memory Allocator Performance Analysis
Test different particle scales to validate malloc overhead claims
"""

import subprocess
import time
import json
import matplotlib.pyplot as plt
import numpy as np
import os

def generate_test_data(n_values):
    """Generate test data for different particle counts"""
    for n in n_values:
        filename = f"data/inputs/test_{n}.bin"
        if not os.path.exists(filename):
            print(f"Generating {n} particles...")
            subprocess.run([
                "python3", "scripts/generate_data.py", 
                str(n), filename
            ], check=True, capture_output=True)

def run_benchmark(version, n, steps=10, runs=3):
    """Run benchmark for a specific version and particle count"""
    times = []
    input_file = f"data/inputs/test_{n}.bin"
    
    for i in range(runs):
        cmd = [
            f"./build/{version}",
            str(n), input_file, str(steps), 
            "0.001", "1", "0.5", "8"
        ]
        
        start_time = time.time()
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            end_time = time.time()
            
            if result.returncode == 0:
                # Extract time from output if available
                output = result.stdout
                if "took" in output:
                    # Try to parse time from the output
                    try:
                        time_str = output.split("took")[1].split("seconds")[0].strip()
                        elapsed = float(time_str)
                    except:
                        elapsed = end_time - start_time
                else:
                    elapsed = end_time - start_time
                
                times.append(elapsed)
                print(f"  {version} N={n} Run {i+1}: {elapsed:.3f}s")
            else:
                print(f"  {version} N={n} failed: {result.stderr}")
                return None
                
        except Exception as e:
            print(f"  {version} N={n} error: {e}")
            return None
    
    return np.mean(times) if times else None

def main():
    # Test different particle scales
    n_values = [1000, 2000, 5000, 10000, 20000]
    versions = ["v2_barnes_hut", "v3_arena"]
    steps = 10
    runs = 3
    
    print("Generating test data...")
    generate_test_data(n_values)
    
    results = {v: [] for v in versions}
    
    print("\nRunning benchmarks...")
    for n in n_values:
        print(f"\nTesting N={n}:")
        for version in versions:
            avg_time = run_benchmark(version, n, steps, runs)
            if avg_time is not None:
                results[version].append(avg_time)
            else:
                # If failed, try with smaller steps
                avg_time = run_benchmark(version, n, 5, runs)
                if avg_time is not None:
                    results[version].append(avg_time)
    
    # Save results
    with open("malloc_analysis_results.json", "w") as f:
        json.dump({"n_values": n_values, "results": results}, f, indent=2)
    
    # Create visualization
    create_visualization(n_values, results)
    
    print("\nAnalysis complete! Results saved to malloc_analysis_results.json")
    print("Visualization saved to malloc_analysis_plot.png")

def create_visualization(n_values, results):
    """Create visualization of the results"""
    plt.style.use('seaborn-v0_8')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot absolute times
    for version, times in results.items():
        if times and len(times) == len(n_values):
            label = "v2 (malloc)" if "v2" in version else "v3 (arena)"
            ax1.plot(n_values, times, 'o-', label=label, linewidth=2, markersize=8)
    
    ax1.set_xlabel('Number of Particles (N)')
    ax1.set_ylabel('Time per 10 Steps (seconds)')
    ax1.set_title('Absolute Performance: malloc vs Arena Allocator')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    # Plot speedup ratio
    if len(results["v2_barnes_hut"]) == len(results["v3_arena"]) == len(n_values):
        speedup_ratio = [
            results["v2_barnes_hut"][i] / results["v3_arena"][i] 
            for i in range(len(n_values))
        ]
        
        ax2.plot(n_values, speedup_ratio, 's-', color='red', 
                label='Speedup (v2/v3)', linewidth=2, markersize=8)
        ax2.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7)
        
        ax2.set_xlabel('Number of Particles (N)')
        ax2.set_ylabel('Speedup Ratio (v2_time / v3_time)')
        ax2.set_title('Arena Allocator Speedup Over System malloc')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig('malloc_analysis_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print summary
    print("\n=== Performance Summary ===")
    for i, n in enumerate(n_values):
        if (i < len(results["v2_barnes_hut"]) and 
            i < len(results["v3_arena"])):
            v2_time = results["v2_barnes_hut"][i]
            v3_time = results["v3_arena"][i]
            speedup = v2_time / v3_time
            print(f"N={n}: v2={v2_time:.3f}s, v3={v3_time:.3f}s, speedup={speedup:.2f}x")

if __name__ == "__main__":
    main()