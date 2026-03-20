import matplotlib.pyplot as plt
import numpy as np
import json
import os

def report_adjusted_efficiency(data):
    weak = data['weak']
    threads = np.array([r['threads'] for r in weak])
    runtimes = np.array([r['runtime'] for r in weak])
    n_values = np.array([r['N'] for r in weak])

    n_base = n_values[0]
    t1 = runtimes[0]

    work_ratio = (n_values * np.log2(n_values)) / (threads * (n_base * np.log2(n_base)))
    efficiency_adj = (t1 / runtimes) * work_ratio

    print(f"{'Threads':<10} | {'N':<10} | {'Runtime':<10} | {'Adj. Eff.':<10}")
    print("-" * 55)
    for t, n, runtime, eff in zip(threads, n_values, runtimes, efficiency_adj):
        print(f"{int(t):<10} | {int(n):<10} | {runtime:<10.4f} | {eff:<10.2%}")

def load_real_data():
    with open('data/metrics/uppmax_full_results.json', 'r') as f:
        return json.load(f)

def plot_strong_scaling(data):
    strong = data['strong']
    threads = [r['threads'] for r in strong]
    runtimes = [r['runtime'] for r in strong]
    
    plt.figure(figsize=(10, 6))
    
    # Baseline-normalized speedup.
    actual_speedup = [runtimes[0] / t for t in runtimes]
    
    plt.plot(threads, threads, 'r--', linewidth=2, label='Ideal Speedup (Linear)')
    plt.plot(threads, actual_speedup, 'bo-', markersize=8, label='Actual Speedup (UPPMAX pelle)')
    
    # Annotate the last measured point.
    plt.text(threads[-1], actual_speedup[-1] + 0.1, f'{actual_speedup[-1]:.2f}x', 
             ha='center', fontweight='bold', color='blue')
    
    plt.xlabel('Thread Count (P)', fontsize=12)
    plt.ylabel('Speedup (T1 / Tp)', fontsize=12)
    plt.title('Strong Scaling on UPPMAX (N=100k)', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend()
    
    plt.savefig('figures/strong_scaling.png', dpi=300)
    plt.close()

def plot_weak_scaling_adjusted(data):
    weak = data['weak']
    threads = np.array([r['threads'] for r in weak])
    runtimes = np.array([r['runtime'] for r in weak])
    n_values = np.array([r['N'] for r in weak])
    
    n_base = n_values[0]
    t1 = runtimes[0]
    
    # Adjust weak-scaling efficiency by O(N log N) work growth.
    work_ratio = (n_values * np.log2(n_values)) / (threads * (n_base * np.log2(n_base)))
    efficiency_raw = (t1 / runtimes)
    efficiency_adj = efficiency_raw * work_ratio
    
    plt.figure(figsize=(10, 6))
    plt.plot(threads, efficiency_adj * 100, '#2ca02c', marker='o', markersize=8, label=r'Adjusted Efficiency ($O(N \log N)$ Adjusted)')
    plt.plot(threads, efficiency_raw * 100, '#7f7f7f', linestyle='--', alpha=0.7, label='Raw Efficiency')
    
    plt.yscale('log')
    plt.xlabel('Thread Count (P)', fontsize=12)
    plt.ylabel('Efficiency (%)', fontsize=12)
    plt.title('Weak Scaling: Efficiency Decay (N_base=100k)', fontsize=14)
    
    from matplotlib.ticker import ScalarFormatter
    plt.gca().yaxis.set_major_formatter(ScalarFormatter())
    plt.yticks([2, 5, 10, 25, 50, 100])
    
    plt.xlabel('Thread Count (P)')
    plt.xticks(threads)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend()
    
    plt.savefig('figures/weak_scaling.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    if not os.path.exists('figures'):
        os.makedirs('figures')
    
    try:
        data = load_real_data()
        plot_strong_scaling(data)
        plot_weak_scaling_adjusted(data)
        report_adjusted_efficiency(data)
        print("Scaling plots re-generated using UPPMAX experimental data.")
    except Exception as e:
        print(f"Error: {e}")
