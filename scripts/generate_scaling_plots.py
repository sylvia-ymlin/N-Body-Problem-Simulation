#!/usr/bin/env python3
"""
Generate professional scaling visualization plots for n-body simulation report
Using the specified color palette - Resized to smaller dimensions
"""

import matplotlib.pyplot as plt
import numpy as np

# Define the specified color palette
NBody_Palette = {
    "gray_light": "#D8DEE9",
    "gray_dark":  "#4C566A", 
    "blue_light": "#88C0D0",
    "blue_mid":   "#81A1C1",
    "blue_deep":  "#5E81AC",
    "orange":     "#D08770",
    "red":        "#BF616A"
}

# Set matplotlib style for academic papers
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.figsize': (5, 3.5),
    'figure.dpi': 300,
    'axes.edgecolor': NBody_Palette["gray_dark"],
    'axes.labelcolor': NBody_Palette["gray_dark"],
    'text.color': NBody_Palette["gray_dark"],
    'xtick.color': NBody_Palette["gray_dark"],
    'ytick.color': NBody_Palette["gray_dark"]
})

# Strong scaling data (fixed N=100,000)
cores_strong = np.array([1, 2, 4, 8, 16])
time_strong = np.array([285, 148, 82, 72, 68])
speedup_strong = time_strong[0] / time_strong

# Weak scaling data (fixed work per core)
cores_weak = np.array([1, 2, 4, 8, 16])
particles_weak = np.array([12500, 25000, 50000, 100000, 200000])
time_weak = np.array([36, 38, 41, 45, 52])
efficiency_weak = time_weak[0] / time_weak * 100

# Create Strong Scaling Speedup chart
fig, ax = plt.subplots(figsize=(5, 3.5))

ax.plot(cores_strong, cores_strong, '--', color=NBody_Palette["gray_dark"], 
         label='Ideal Scaling', linewidth=1.5, alpha=0.7)
ax.plot(cores_strong, speedup_strong, 'o-', color=NBody_Palette["blue_deep"], 
         linewidth=1.5, markersize=5, label='Measured Speedup')
ax.set_xlabel('Number of Cores', fontweight='bold')
ax.set_ylabel('Speedup (vs 1 core)', fontweight='bold')
ax.set_title('Strong Scaling Performance\n(Fixed N=100,000 particles)', fontweight='bold')
ax.grid(True, color=NBody_Palette["gray_light"], alpha=0.8)
ax.legend()
ax.set_facecolor(NBody_Palette["gray_light"])

plt.tight_layout()
plt.savefig('strong_scaling_speedup.png', dpi=300, bbox_inches='tight', 
           facecolor='white', edgecolor='none')
plt.savefig('strong_scaling_speedup.pdf', bbox_inches='tight')
plt.close()

# Create Weak Scaling Efficiency chart with particle count labels
fig, ax = plt.subplots(figsize=(5, 3.5))

bars = ax.plot(cores_weak, efficiency_weak, 'D-', color=NBody_Palette["blue_mid"], 
        linewidth=1.5, markersize=5, label='Scaling Efficiency')
ax.axhline(y=80, color=NBody_Palette["red"], linestyle='--', alpha=0.8, 
          label='80% Efficiency Target')

# Add particle count annotations
for i, (core, eff, particles) in enumerate(zip(cores_weak, efficiency_weak, particles_weak)):
    ax.annotate(f'N={particles//1000}K', 
               xy=(core, eff), 
               xytext=(3, 3), 
               textcoords='offset points',
               fontsize=8,
               bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

ax.set_xlabel('Number of Cores', fontweight='bold')
ax.set_ylabel('Scaling Efficiency (%)', fontweight='bold')
ax.set_title('Weak Scaling Efficiency\n(Fixed Work per Core)', fontweight='bold')
ax.grid(True, color=NBody_Palette["gray_light"], alpha=0.8)
ax.legend()
ax.set_facecolor(NBody_Palette["gray_light"])

plt.tight_layout()
plt.savefig('weak_scaling_efficiency.png', dpi=300, bbox_inches='tight', 
           facecolor='white', edgecolor='none')
plt.savefig('weak_scaling_efficiency.pdf', bbox_inches='tight')
plt.close()

print("Resized scaling plots generated (5x3.5 inches):")
print("- strong_scaling_speedup.png: Strong scaling speedup")
print("- weak_scaling_efficiency.png: Weak scaling with particle counts")

# Generate markdown for report
print("\nMarkdown for REPORT.md:")
print("""
### Parallel Scaling Analysis (v4 Engine)

#### Strong Scaling Performance
![Strong Scaling](../strong_scaling_speedup.png)

*Strong scaling demonstrates excellent performance up to 4 cores (88% efficiency), with near-linear speedup. Beyond 4 cores, memory bandwidth saturation limits further scaling, which is typical for memory-bound applications.*

#### Weak Scaling and Large-Scale Performance  
![Weak Scaling](../weak_scaling_efficiency.png)

*Weak scaling maintains 80% efficiency at 8 cores, confirming effective workload distribution and minimal parallel overhead in the hybrid parallel strategy. Particle counts scale proportionally with core count (12.5K particles per core).*
""")