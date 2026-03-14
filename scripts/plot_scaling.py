import matplotlib.pyplot as plt
import numpy as np
import json
import os

def load_real_data():
    with open('xeon_results_real.json', 'r') as f:
        return json.load(f)

def plot_strong_scaling(data):
    strong = data['strong']
    threads = [r['threads'] for r in strong]
    runtimes = [r['runtime'] for r in strong]
    
    plt.figure(figsize=(10, 6))
    
    # Ideal speedup
    actual_speedup = [runtimes[0] / t for t in runtimes]
    
    plt.plot(threads, threads, 'r--', linewidth=2, label='Ideal Speedup (Linear)')
    plt.plot(threads, actual_speedup, 'bo-', markersize=8, label='Actual Speedup (Xeon 6430)')
    
    # Annotate Point of Collapse
    plt.text(threads[-1], actual_speedup[-1] + 0.1, f'{actual_speedup[-1]:.2f}x', 
             ha='center', fontweight='bold', color='blue')
    
    plt.xlabel('Thread Count (P)', fontsize=12)
    plt.ylabel('Speedup (T1 / Tp)', fontsize=12)
    plt.title('Strong Scaling: Performance Collapse at High P (N=100k)', fontsize=14)
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
    
    work_ratio = (n_values * np.log2(n_values)) / (threads * (n_base * np.log2(n_base)))
    efficiency_raw = (t1 / runtimes)
    efficiency_adj = efficiency_raw * work_ratio
    
    plt.figure(figsize=(10, 6))
    plt.plot(threads, efficiency_adj * 100, '#2ca02c', marker='o', markersize=8, label='Adjusted Efficiency ($O(N \log N)$ Adjusted)')
    plt.plot(threads, efficiency_raw * 100, '#7f7f7f', linestyle='--', alpha=0.7, label='Raw Efficiency')
    
    plt.yscale('log') # Log scale to visualize the decay
    plt.xlabel('Thread Count (P)', fontsize=12)
    plt.ylabel('Efficiency (%)', fontsize=12)
    plt.title('Weak Scaling: Efficiency Decay (N_base=100k)', fontsize=14)
    
    # Formatting the log axis to look natural
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
        print("Scaling plots re-generated using REAL Xeon experimental data.")
    except Exception as e:
        print(f"Error: {e}")
