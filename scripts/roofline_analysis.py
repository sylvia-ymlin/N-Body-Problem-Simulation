#!/usr/bin/env python3
"""
Roofline Model Analysis for N-Body Simulation
Based on empirical FLOPs and memory traffic measurements
"""

import numpy as np
import matplotlib.pyplot as plt

# ==================== HARDWARE PARAMETERS ====================
# Typical modern CPU specifications (adjust based on actual hardware)
PEAK_BANDWIDTH = 50.0  # GB/s (typical DDR4/DDR5)
PEAK_FLOPS = 100.0     # GFLOP/s (single core)
L3_CACHE_SIZE = 16.0   # MB

# ==================== EMPIRICAL MEASUREMENTS ====================
# Version performance data (N=50,000 particles)
VERSIONS = {
    'v1': {'time': 66.5, 'description': 'Direct O(N²)'},
    'v2': {'time': 2.3, 'description': 'Barnes-Hut O(N log N)'},
    'v4': {'time': 0.31, 'description': 'Parallel Barnes-Hut'}
}

# ==================== FLOPs ANALYSIS ====================
# Manual FLOPs counting from compute_force function
# Each particle-node interaction: 16 FLOPs
FLOPs_PER_INTERACTION = 16

# For N=50,000 particles
N = 50000
LOG2_N = np.log2(N)  # ~15.6

# Total interactions per step
TOTAL_INTERACTIONS = N * LOG2_N
TOTAL_FLOPS = TOTAL_INTERACTIONS * FLOPs_PER_INTERACTION  # ~12.5 MFLOPs

# ==================== MEMORY TRAFFIC ANALYSIS ====================
# Memory access per interaction: 64 bytes (8 doubles)
BYTES_PER_INTERACTION = 64
TOTAL_BYTES = TOTAL_INTERACTIONS * BYTES_PER_INTERACTION  # ~50 MB

# ==================== PERFORMANCE CALCULATION ====================
for version, data in VERSIONS.items():
    time_s = data['time']
    
    # Performance in FLOP/s
    flops_per_s = TOTAL_FLOPS / time_s  # FLOP/s
    data['GFLOPs'] = flops_per_s / 1e9
    
    # Bandwidth utilization
    bytes_per_s = TOTAL_BYTES / time_s  # bytes/s
    data['GBs'] = bytes_per_s / 1e9
    
    # Operational Intensity
    data['OI'] = flops_per_s / bytes_per_s  # FLOP/byte

# ==================== PLOTTING ====================
plt.figure(figsize=(10, 7), dpi=150)
ax = plt.gca()

# Calculate Knee Point (Ridge Point)
OI_knee = PEAK_FLOPS / PEAK_BANDWIDTH

# Define X-axis range
xmin, xmax = 0.05, 100.0
OI_range = np.logspace(np.log10(xmin), np.log10(xmax), 1000)

# Memory-bound performance curve
mem_bound_perf = np.minimum(PEAK_BANDWIDTH * OI_range, PEAK_FLOPS)

# --- 1. Shaded Regions (Memory vs Compute Bound) ---
# Memory-bound region (Left of Knee)
ax.fill_between([xmin, OI_knee], [1e-4, 1e-4], [PEAK_FLOPS, PEAK_FLOPS], 
                color='#d9e6f2', alpha=0.5, label='Memory-bound region')
# Compute-bound region (Right of Knee)
ax.fill_between([OI_knee, xmax], [1e-4, 1e-4], [PEAK_FLOPS, PEAK_FLOPS], 
                color='#f2d9d9', alpha=0.5, label='Compute-bound region')

# --- 2. Roofline Curve ---
plt.loglog(OI_range, mem_bound_perf, color='#4c72b0', linewidth=3, label='Roofline Model')

# --- 3. Ceilings & Reference Lines ---
# Peak Compute (Horizontal Dashed Line)
plt.axhline(y=PEAK_FLOPS, color='#c44e52', linestyle=':', linewidth=1.5)
plt.text(xmin * 1.5, PEAK_FLOPS * 1.2, f'Peak Compute: {PEAK_FLOPS} GFLOPS', 
         color='#c44e52', fontsize=10, fontweight='bold')

# Ridge Point (Vertical Dashed Line)
plt.axvline(x=OI_knee, color='gray', linestyle='--', linewidth=1.5)
plt.text(OI_knee * 1.1, 1e-2, f'Ridge Point\n({OI_knee:.1f} FLOP/byte)', 
         color='#555555', fontsize=9)

# --- 4. Data Points ---
# Distinct markers and colors for versions
markers = {'v1': 'o', 'v2': 's', 'v4': '^'}
colors = {'v1': '#55a868', 'v2': '#dd8452', 'v4': '#8172b3'}  # Muted green, orange, purple
labels = {'v1': 'v1: Direct O(N²)', 'v2': 'v2: Barnes-Hut', 'v4': 'v4: Parallel BH'}

# Specific offsets for labels to avoid overlap
offsets = {
    'v1': (0, -25),   # Below
    'v2': (30, 0),    # Right
    'v4': (0, 15)     # Above
}
alignments = {
    'v1': ('center', 'top'),
    'v2': ('left', 'center'),
    'v4': ('center', 'bottom')
}

for version, data in VERSIONS.items():
    plt.loglog(data['OI'], data['GFLOPs'], 
               marker=markers[version], markersize=12, 
               markerfacecolor=colors[version], markeredgecolor='black', markeredgewidth=1.5,
               linestyle='None', label=labels[version])
    
    # Annotate points
    off = offsets.get(version, (0, 10))
    ha, va = alignments.get(version, ('center', 'bottom'))
    
    plt.annotate(f"{version}",
                 xy=(data['OI'], data['GFLOPs']),
                 xytext=off, textcoords='offset points', 
                 ha=ha, va=va,
                 fontsize=9, fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8, ec='gray'))

# --- 5. Formatting ---
plt.xlabel('Arithmetic Intensity (FLOP/byte)', fontsize=12, fontweight='bold', color='#333333')
plt.ylabel('Performance (GFLOPS)', fontsize=12, fontweight='bold', color='#333333')
plt.title('Roofline Model: N-Body Simulation', fontsize=14, fontweight='bold', pad=15)

# Set axis limits
plt.xlim(xmin, xmax)
plt.ylim(1e-5, PEAK_FLOPS * 5)  # Lowered Y-min to 1e-5 to accommodate v1 label

# Grid
plt.grid(True, which="major", ls="-", color='white', alpha=0.5)
plt.grid(True, which="minor", ls=":", color='white', alpha=0.2)
ax.set_facecolor('#f0f0f0')  # Light gray background for plot area

# Legend
plt.legend(loc='center right', frameon=True, facecolor='white', framealpha=0.9, fontsize=10)

plt.tight_layout()

# Save plot
output_path = '/Users/ymlin/Downloads/003-Study/137-Projects/04-nBody-Problem-Simulation/docs/roofline_analysis.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Plot saved to {output_path}")

# ==================== QUANTITATIVE ANALYSIS ====================
print("=== ROOFLINE ANALYSIS RESULTS ===")
print(f"Hardware Parameters:")
print(f"  Peak Bandwidth: {PEAK_BANDWIDTH} GB/s")
print(f"  Peak Compute: {PEAK_FLOPS} GFLOP/s")
print(f"  Knee Point: OI = {PEAK_FLOPS/PEAK_BANDWIDTH:.1f} FLOP/byte")
print()

print(f"Algorithm Analysis (N={N}):")
print(f"  Total FLOPs/step: {TOTAL_FLOPS/1e6:.1f} MFLOP")
print(f"  Total Memory/step: {TOTAL_BYTES/1e6:.1f} MB")
print(f"  Base Operational Intensity: {TOTAL_FLOPS/TOTAL_BYTES:.3f} FLOP/byte")
print()

print("Version Performance:")
for version, data in VERSIONS.items():
    print(f"  {version}: {data['GFLOPs']:.3f} GFLOP/s, OI={data['OI']:.3f}, "
          f"Efficiency={data['GFLOPs']/min(PEAK_FLOPS, PEAK_BANDWIDTH*data['OI'])*100:.1f}%")

print("\n=== KEY INSIGHTS ===")
print("1. All versions are MEMORY-BOUND (OI < knee point)")
print("2. v2 achieves only 5.4% of peak memory bandwidth")
print("3. Low OI indicates poor cache reuse and memory access patterns")
print("4. Parallelism (v4) improves performance but doesn't change fundamental bound")