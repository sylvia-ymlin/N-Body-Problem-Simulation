import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import time

def read_gal(filename, N):
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(N, 5)

def match_particles(ref, test):
    """Match particles using the Hungarian algorithm.
    
    v5 reorders particles via Morton Z-order sort, so indices don't
    correspond to v1.  We build a cost matrix where mass difference
    (heavily weighted) is the primary criterion and Euclidean position
    distance is the tiebreaker, then solve the optimal assignment.
    """
    # Primary: mass difference (mass is invariant across simulations)
    mass_cost = np.abs(ref[:, 2:3] - test[:, 2:3].T)  # N x N
    # Secondary: position distance (tiebreak same-mass particles)
    pos_cost = cdist(ref[:, 0:2], test[:, 0:2])
    # Weight mass heavily so same-mass particles are matched first
    cost = mass_cost * 1e6 + pos_cost
    _, col_ind = linear_sum_assignment(cost)
    return test[col_ind]


def compute_rmse(ref_data, test_data):
    N = ref_data.shape[0]
    pos_diff = ref_data[:, 0:2] - test_data[:, 0:2]
    rmse = np.linalg.norm(pos_diff) / np.sqrt(N)
    return rmse

def run_validation():
    print("--- Scientific Accuracy Validation ---")
    N = 500
    original_input = "data/inputs/input_500.bin"
    steps = 1
    dt = 0.001

    print(f"--- Generating Sorted Baseline (Using v4 Internal Sort) ---")
    # Run v4 just to generate sorted.gal (dummy run)
    # sort_mode=2 (Sort + Dump)
    subprocess.run(["./build/v4_morton", str(N), original_input, "0", str(dt), "1", "0.0", "1", "2"], check=True)
    
    if not os.path.exists("sorted.gal"):
        print("Error: sorted.gal not generated.")
        return
        
    sorted_input = "sorted.gal"
    print(f"Generated {sorted_input}")

    # Run v1 (Exact O(N^2)) on SORTED input
    print("Running v1 (Exact O(N^2)) on Internal Sorted Input...")
    start_v1 = time.time()
    subprocess.run(["./build/v1_naive", str(N), sorted_input, str(steps), str(dt), "1", "0.0", "1"], check=True)
    end_v1 = time.time()
    os.rename("result.gal", "result_v1.gal")

    v1_data = read_gal("result_v1.gal", N)
    if v1_data.size == 0: return

    # Check v4 Accuracy with varying theta
    targets = [
        ("theta=0.0", 0.0, 1e-12), # Single step should be machine precision
        ("theta=0.5", 0.5, 0.05),
        ("theta=1.0", 1.0, 0.05)
    ]
    
    overall_pass = True
    
    for label, theta, threshold in targets:
        print(f"\n--- {label} ---")
        start_bh = time.time()
        # Use v4 with EnableSort=0 (Use the pre-dumped sorted file as is)
        # So v4 uses input order (which is Sorted).
        subprocess.run(["./build/v4_morton", str(N), sorted_input, str(steps), str(dt), "1", str(theta), "1", "0"], check=True)
        end_bh = time.time()
        
        os.rename("result.gal", "result_bh.gal")
        bh_data = read_gal("result_bh.gal", N)
        
        rmse = compute_rmse(v1_data, bh_data)
        print(f"Position RMSE: {rmse:.2e}")
        
        if rmse < threshold:
            print(f"RESULT: [PASS] Deviation {rmse:.2e} within threshold {threshold}")
        else:
            print(f"RESULT: [FAIL] Deviation {rmse:.2e} exceeds threshold {threshold}")
            overall_pass = False

    print("\n" + "="*40)
    if overall_pass:
        print("Overall: ALL PASSED")
    else:
        print("Overall: SOME FAILED")

if __name__ == "__main__":
    if not os.path.exists("./build/v1_naive"):
        print("Error: Binaries not found. Please run 'cmake --build build' first.")
    else:
        run_validation()
