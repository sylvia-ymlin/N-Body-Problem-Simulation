import matplotlib.pyplot as plt
import numpy as np
import json
import os

def report_adjusted_efficiency(weak):
    threads = np.array([r['threads'] for r in weak], dtype=float)
    runtimes = np.array([r['median_after_warmup'] for r in weak], dtype=float)
    n_values = np.array([r['N'] for r in weak], dtype=float)

    n_base = n_values[0]
    t1 = runtimes[0]

    work_ratio = (n_values * np.log2(n_values)) / (threads * (n_base * np.log2(n_base)))
    efficiency_adj = (t1 / runtimes) * work_ratio

    print(f"{'Threads':<10} | {'N':<10} | {'Runtime':<10} | {'Adj. Eff.':<10}")
    print("-" * 55)
    for t, n, runtime, eff in zip(threads, n_values, runtimes, efficiency_adj):
        print(f"{int(t):<10} | {int(n):<10} | {runtime:<10.4f} | {eff:<10.2%}")

def load_scaling_data():
    with open('data/metrics/autodl_strong_scaling.json', 'r') as f:
        strong = json.load(f)['results']
    with open('data/metrics/autodl_weak_scaling.json', 'r') as f:
        weak = json.load(f)['results']
    return strong, weak

def load_locality_data():
    with open('data/metrics/autodl_locality_ablation.json', 'r') as f:
        return json.load(f)

def plot_locality_ablation(data):
    runs = data['results']
    labels = []
    runtimes = []
    colors = []
    for r in runs:
        if r['k'] == 0:
            labels.append('Morton')
            colors.append('#1f77b4')
        else:
            labels.append(f'K-Means (k={r["k"]})')
            colors.append('#ff7f0e')
        runtimes.append(r['median_after_warmup'])

    plt.figure(figsize=(9, 5))
    bars = plt.bar(labels, runtimes, color=colors, width=0.55)
    for bar, val in zip(bars, runtimes):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                 f'{val:.2f}s', ha='center', va='bottom', fontweight='bold')

    best_kmeans = min(r['median_after_warmup'] for r in runs if r['k'] != 0)
    ratio = best_kmeans / runtimes[0]
    plt.title(f'Locality Strategy Ablation (v4, N=100k, serial)\nBest K-Means is {ratio:.2f}× slower than Morton', fontsize=13)
    plt.ylabel('Wall Time (s)', fontsize=12)
    plt.ylim(0, max(runtimes) * 1.2)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('figures/locality_ablation.png', dpi=300)
    plt.close()

def plot_strong_scaling(strong):
    threads = [r['threads'] for r in strong]
    runtimes = [r['median_after_warmup'] for r in strong]
    
    plt.figure(figsize=(10, 6))
    
    # Baseline-normalized speedup.
    actual_speedup = [runtimes[0] / t for t in runtimes]
    
    plt.plot(threads, threads, 'r--', linewidth=2, label='Ideal Speedup (Linear)')
    plt.plot(threads, actual_speedup, 'bo-', markersize=8, label='Actual Speedup (AutoDL AMD EPYC)')
    
    # Annotate the last measured point.
    plt.text(threads[-1], actual_speedup[-1] + 0.1, f'{actual_speedup[-1]:.2f}x', 
             ha='center', fontweight='bold', color='blue')
    
    plt.xlabel('Thread Count (P)', fontsize=12)
    plt.ylabel('Speedup (T1 / Tp)', fontsize=12)
    plt.title('Strong Scaling on AutoDL (N=100k)', fontsize=14)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend()
    
    plt.savefig('figures/strong_scaling.png', dpi=300)
    plt.close()

def plot_weak_scaling_adjusted(weak):
    threads = np.array([r['threads'] for r in weak], dtype=float)
    runtimes = np.array([r['median_after_warmup'] for r in weak], dtype=float)
    n_values = np.array([r['N'] for r in weak], dtype=float)
    
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
    plt.title('Weak Scaling on AutoDL: Efficiency Decay (N_base=100k)', fontsize=14)
    
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
        strong, weak = load_scaling_data()
        plot_strong_scaling(strong)
        plot_weak_scaling_adjusted(weak)
        report_adjusted_efficiency(weak)

        locality = load_locality_data()
        plot_locality_ablation(locality)

        print("Plots re-generated using AutoDL experimental data.")
    except Exception as e:
        print(f"Error: {e}")
