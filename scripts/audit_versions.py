import subprocess
import os
import struct
import math
from scipy.optimize import linear_sum_assignment
import numpy as np

def read_gal(filename, N):
    data = []
    try:
        with open(filename, "rb") as f:
            for _ in range(N):
                rx = struct.unpack("d", f.read(8))[0]
                ry = struct.unpack("d", f.read(8))[0]
                m = struct.unpack("d", f.read(8))[0]
                vx = struct.unpack("d", f.read(8))[0]
                vy = struct.unpack("d", f.read(8))[0]
                data.append({'x': rx, 'y': ry, 'm': m, 'vx': vx, 'vy': vy})
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return []
    return data

def match_particles(ref_data, test_data):
    if not ref_data or not test_data:
        return []
    N = len(ref_data)
    cost_matrix = np.zeros((N, N))
    
    # Pre-compute cost matrix (squared distance)
    # Using only positions for matching
    ref_pos = np.array([[p['x'], p['y']] for p in ref_data])
    test_pos = np.array([[p['x'], p['y']] for p in test_data])
    
    # Broadcast to compute pairwise differences
    # dist[i,j] = sum((ref[i] - test[j])^2)
    # This is expensive O(N^2), but N=500 is fine (250k).
    for i in range(N):
        diff = test_pos - ref_pos[i]
        dist_sq = np.sum(diff**2, axis=1)
        cost_matrix[i, :] = dist_sq
        
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    matched_data = [None] * N
    for i, j in zip(row_ind, col_ind):
        matched_data[i] = test_data[j]
        
    return matched_data

def compute_rmse(ref_data, test_data):
    if not ref_data or not test_data:
        return float('inf')
    
    mse = 0
    N = len(ref_data)
    for i in range(N):
        dx = ref_data[i]['x'] - test_data[i]['x']
        dy = ref_data[i]['y'] - test_data[i]['y']
        mse += dx*dx + dy*dy
    return math.sqrt(mse / N)

def run_audit():
    N = 500
    input_file = "data/inputs/input_500.bin"
    steps = 50
    dt = 0.001
    
    binaries = [
        ("v1_naive", 1), # name, n_threads expected arg (dummy)
        ("v2_barnes_hut", 1),
        ("v3_arena", 1),
        ("v4_morton", 1),
        ("v5_parallel", 4)
    ]
    
    print("--- Running Audit ---")
    
    # Run v1
    print("Running v1_naive (Reference)...")
    subprocess.run(["./build/v1_naive", str(N), input_file, str(steps), str(dt), "1", "0.0", "1"], check=True)
    os.rename("result.gal", "result_v1.gal")
    ref_data = read_gal("result_v1.gal", N)
    
    results = []
    
    for bin_name, threads in binaries:
        if bin_name == "v1_naive":
            continue
            
        print(f"Running {bin_name}...")
        try:
            # Theta = 0.0 to force consistency
            subprocess.run([f"./build/{bin_name}", str(N), input_file, str(steps), str(dt), str(threads), "0.0", "1"], check=True)
            output_name = f"result_{bin_name}.gal"
            if os.path.exists("result.gal"):
                os.rename("result.gal", output_name)
            
            test_data = read_gal(output_name, N)
            
            # v4 and v5 sort particles. v2 and v3 likely do not (unless implemented in barnes_hut? no).
            # v2/v3 use pointer tree, no sorting of array.
            # v4/v5 use Z_order_sort.
            
            # Blindly apply matching if RMSE is high?
            # Or just always apply matching.
            matched_test = match_particles(ref_data, test_data)
            
            rmse = compute_rmse(ref_data, matched_test)
            results.append((bin_name, rmse))
            print(f"  RMSE: {rmse:.2e}")
            
        except subprocess.CalledProcessError as e:
            print(f"  FAILED: {e}")
            results.append((bin_name, "FAILED"))
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append((bin_name, "ERROR"))

    print("\n--- Summary ---")
    print(f"{'Version':<15} | {'RMSE (vs v1)':<15}")
    print("-" * 35)
    for name, res in results:
        if isinstance(res, float):
            print(f"{name:<15} | {res:.2e}")
        else:
            print(f"{name:<15} | {res}")

if __name__ == "__main__":
    run_audit()
