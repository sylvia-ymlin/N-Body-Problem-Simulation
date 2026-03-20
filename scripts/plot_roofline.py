import numpy as np
import matplotlib.pyplot as plt

def generate_roofline():
    # Hardware Characteristics (Estimated for Xeon Gold 6430 / M3)
    # Bandwidth in GB/s
    peak_bw = 60 
    # Peak Performance in GFLOPS
    peak_perf = 250 

    # Operational Intensity range (FLOP/Byte)
    oi_range = np.logspace(-2, 1, 100)
    
    # Roofline performance calculation
    roofline_perf = np.minimum(peak_perf, peak_bw * oi_range)

    plt.figure(figsize=(10, 7))
    plt.loglog(oi_range, roofline_perf, 'r-', label='Xeon Gold 6430 Roofline', linewidth=3)

    # Data Points (Calculated from benchmarks)
    # Naive v1: OI ~ 0.18, GFLOPS ~ 1.5
    # Parallel v5: OI ~ 0.28, GFLOPS ~ 12.0 (Estimated for 4T Peak)
    versions = {
        'v1 Naive': (0.18, 1.5, 'bo'),
        'v5 Parallel (Optimal 4T)': (0.28, 12.0, 'go'),
        'v5 Parallel (Stalled 16T)': (0.28, 6.0, 'ro')
    }

    for label, (oi, gflops, style) in versions.items():
        plt.plot(oi, gflops, style, label=f'{label} (N=1M)', markersize=10)
        plt.annotate(label, (oi, gflops), textcoords="offset points", xytext=(0,10), ha='center')

    # Formatting
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.xlabel('Operational Intensity (FLOPs/Byte)', fontsize=12)
    plt.ylabel('Performance (GFLOPS)', fontsize=12)
    plt.title('Roofline Analysis: N-Body Simulation Engine', fontsize=14)
    plt.legend()
    
    # Shade Memory-Bound and Compute-Bound regions
    knee_point = peak_perf / peak_bw
    plt.axvline(x=knee_point, color='k', linestyle='--', alpha=0.3)
    plt.text(knee_point * 1.1, 1e-1, f'Knee Point: {knee_point:.1f} FLOP/B', rotation=90)
    
    plt.fill_between(oi_range, 1e-2, roofline_perf, where=(oi_range < knee_point), color='blue', alpha=0.1)
    plt.text(1e-2, 1e1, 'Memory Bound Region', color='blue', alpha=0.6, fontsize=12)
    plt.text(knee_point * 2, 1e1, 'Compute Bound Region', color='red', alpha=0.6, fontsize=12)

    plt.ylim(1e-1, 1e3)
    plt.xlim(1e-2, 1e1)
    
    plt.tight_layout()
    plt.savefig('figures/roofline_analysis.png', dpi=300)
    print("Roofline plot saved to figures/roofline_analysis.png")

if __name__ == "__main__":
    generate_roofline()
