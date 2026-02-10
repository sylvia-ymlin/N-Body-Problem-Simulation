import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import re
import time
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

def read_gal(filename):
    if not os.path.exists(filename):
        return None
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(-1, 6)

def match_particles(ref, test):
    """Match particles using the Hungarian algorithm.
    
    v5 reorders particles via Morton Z-order sort, so indices don't
    correspond to v1.  Cost = mass_diff * 1e6 + position_distance.
    """
    mass_cost = np.abs(ref[:, 2:3] - test[:, 2:3].T)
    pos_cost = cdist(ref[:, 0:2], test[:, 0:2])
    cost = mass_cost * 1e6 + pos_cost
    _, col_ind = linear_sum_assignment(cost)
    return test[col_ind]

def calculate_momentum_vector(data):
    """Return total momentum as a 2D vector (px, py)."""
    mass = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    return np.array([np.sum(mass * vx), np.sum(mass * vy)])

def calculate_energy(data, G, eps=1e-3):
    """Return total energy (kinetic + potential) of the system."""
    mass = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    # Kinetic energy: 0.5 * m * v^2
    KE = 0.5 * np.sum(mass * (vx**2 + vy**2))
    # Potential energy: -G * mi * mj / |rij| (pairwise)
    x = data[:, 0]
    y = data[:, 1]
    N = len(mass)
    PE = 0.0
    for i in range(N):
        dx = x[i+1:] - x[i]
        dy = y[i+1:] - y[i]
        dist = np.sqrt(dx**2 + dy**2 + eps**2)
        PE -= G * mass[i] * np.sum(mass[i+1:] / dist)
    return KE + PE

def run_experiment():
    print("--- Rigorous Scientific Verification ---")
    N = 2000
    G = 100.0 / N  # Must match the G used in the C simulation
    steps = 200
    dt = 0.0001
    input_file = "data/inputs/input.gal"
    
    # Check if binaries exist
    v1_bin = "./build/v1_naive"
    v5_bin = "./build/v5_parallel"
    
    if not os.path.exists(v1_bin) or not os.path.exists(v5_bin):
        print(f"Error: Binary {v1_bin} or {v5_bin} not found. Please compile the project first.")
        return

    # 0. Read initial state for momentum baseline
    # (Assuming input.gal format matches result.gal for simplicity in reading)
    # Actually, let's just run 0 steps to get initial momentum
    subprocess.run([v1_bin, str(N), input_file, "0", str(dt), "1", "0.0", "1"], check=True, capture_output=True)
    initial_data = read_gal("result.gal")
    p0_vec = calculate_momentum_vector(initial_data)
    E0 = calculate_energy(initial_data, G)
    print(f"Initial System Momentum: |P0| = {np.linalg.norm(p0_vec):.2e}")
    print(f"Initial System Energy:   E0   = {E0:.6e}")

    # 1. Run v1 (Ground Truth)
    print(f"Generating Long-term Ground Truth (steps={steps})...")
    subprocess.run([v1_bin, str(N), input_file, str(steps), str(dt), "1", "0.0", "1"], check=True, capture_output=True)
    os.rename("result.gal", "result_v1.gal")
    v1_data = read_gal("result_v1.gal")
    v1_p_final = calculate_momentum_vector(v1_data)
    v1_momentum_drift = np.linalg.norm(v1_p_final - p0_vec)
    v1_E_final = calculate_energy(v1_data, G)
    v1_energy_drift = abs((v1_E_final - E0) / E0)
    print(f"v1 (Brute-Force) Momentum Drift: {v1_momentum_drift:.2e}")
    print(f"v1 (Brute-Force) Relative Energy Drift: {v1_energy_drift:.2e}")
    
    thetas = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    results = []
    
    # 2. Sweep Theta
    for theta in thetas:
        print(f"Running v5 with theta={theta}...", end=" ", flush=True)
        start_time = time.time()
        process = subprocess.run([v5_bin, str(N), input_file, str(steps), str(dt), "4", str(theta), "4"], capture_output=True, text=True, check=True)
        end_time = time.time()
        
        match = re.search(r"Simulation took ([\d\.]+) seconds", process.stdout)
        sim_time = float(match.group(1)) if match else (end_time - start_time)
        fps = steps / sim_time
        
        v5_data = read_gal("result.gal")
        if v5_data is not None:
            # Hungarian matching: v5 reorders particles via Morton sort
            v5_matched = match_particles(v1_data, v5_data)
            # Error = RMSE of position difference: ||delta_pos|| / sqrt(N)
            pos_diff = v1_data[:, 0:2] - v5_matched[:, 0:2]
            rmse = np.linalg.norm(pos_diff) / np.sqrt(len(v5_data))
            # Momentum & energy use UNMATCHED data (global quantities)
            p_final_vec = calculate_momentum_vector(v5_data)
            p_drift = np.linalg.norm(p_final_vec - p0_vec)
            E_final = calculate_energy(v5_data, G)
            e_drift = abs((E_final - E0) / E0)
            results.append((theta, rmse, fps, p_drift, e_drift))
            print(f"Done. RMSE: {rmse:.2e}, P-Drift: {p_drift:.2e}, ΔE/E0: {e_drift:.2e}")
        else:
            print("Failed.")

    # 3. Plotting
    results = np.array(results)
    theta_vals = results[:, 0]
    errors = results[:, 1]
    fps_vals = results[:, 2]
    
    # NBody Project Palette
    palette = {
        "gray_light": "#D8DEE9",
        "gray_dark":  "#4C566A",
        "blue_light": "#88C0D0",
        "blue_mid":   "#81A1C1",
        "blue_deep":  "#5E81AC",
        "orange":     "#D08770",
        "red":        "#BF616A"
    }
    
    plt.figure(figsize=(10, 6))
    plt.style.use('default') 
    ax = plt.gca()
    ax.set_facecolor('white')
    plt.gcf().set_facecolor('white')
    
    # Plot Pareto Frontier
    plt.plot(errors, fps_vals, 'o-', color=palette["blue_mid"], 
             linewidth=3, markersize=10, label='Simulation Run', 
             markerfacecolor=palette["orange"], markeredgecolor=palette["gray_dark"])
    
    # Annotate points
    for i, theta in enumerate(theta_vals):
        plt.annotate(f"θ={theta}", (errors[i], fps_vals[i]), 
                     textcoords="offset points", xytext=(0,12), 
                     ha='center', color=palette["gray_dark"], 
                     fontsize=10, fontweight='bold')
    
    plt.xlabel('Position RMSE', fontsize=12, color=palette["gray_dark"])
    plt.ylabel('Performance (Throughput FPS)', fontsize=12, color=palette["gray_dark"])
    plt.title(f'Pareto Frontier: Accuracy vs Performance (N={N})', 
              fontsize=16, pad=25, color=palette["blue_deep"], fontweight='bold')
    
    plt.grid(True, linestyle='--', alpha=0.3, color=palette["gray_light"])
    ax.tick_params(axis='x', colors=palette["gray_dark"])
    ax.tick_params(axis='y', colors=palette["gray_dark"])
    
    for spine in ax.spines.values():
        spine.set_color(palette["gray_dark"])
    plt.tight_layout()
    
    output_path = "docs/img/pareto_frontier_real.png"
    plt.savefig(output_path)
    print(f"Plot saved to {output_path}")

if __name__ == "__main__":
    # Create img dir if not exists
    os.makedirs("docs/img", exist_ok=True)
    run_experiment()
