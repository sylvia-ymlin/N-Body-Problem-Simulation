import sys
import os
import numpy as np

def read_gal(filename):
    # io_write_result writes 5 doubles per particle: px, py, m, vx, vy.
    file_size = os.path.getsize(filename)
    N = file_size // 40
    data = np.fromfile(filename, dtype=np.float64).reshape(N, 5)
    return data # px, py, m, vx, vy

def sort_indices(data, force_lexsort=False):
    if force_lexsort:
        return np.lexsort((data[:, 2], data[:, 0], data[:, 1])), "lexsort"

    masses = data[:, 2]
    if len(np.unique(masses)) == len(masses):
        return np.argsort(masses), "mass"

    return np.lexsort((masses, data[:, 0], data[:, 1])), "lexsort"

def calculate_error(file1, file2, sample_k=None, force_lexsort=False):
    d1 = read_gal(file1)
    d2 = read_gal(file2)
    
    if len(d1) != len(d2):
        print(f"Error: Particle count mismatch ({len(d1)} vs {len(d2)})")
        return None
    
    # Align particles first because some kernels reorder arrays (e.g., Morton sort).
    idx1, method1 = sort_indices(d1, force_lexsort)
    idx2, method2 = sort_indices(d2, force_lexsort)
    method = method1 if method1 == method2 else "mixed"

    if len(d1) > 50000 or (sample_k is not None):
        print(f"Info: Using stochastic sampling (K={sample_k or 1000}) for large N={len(d1)}")
        d1 = d1[idx1]
        d2 = d2[idx2]
        
        K = sample_k if sample_k else 1000
        # Sampling keeps checks fast for large N while still catching drift.
        indices = np.random.choice(len(d1), min(K, len(d1)), replace=False)
        d1_s = d1[indices]
        d2_s = d2[indices]
        dist = np.sqrt((d1_s[:,0] - d2_s[:,0])**2 + (d1_s[:,1] - d2_s[:,1])**2)
    else:
        d1_sorted = d1[idx1]
        d2_sorted = d2[idx2]
        dist = np.sqrt((d1_sorted[:,0] - d2_sorted[:,0])**2 + (d1_sorted[:,1] - d2_sorted[:,1])**2)
    
    max_err = np.max(dist)
    mean_err = np.mean(dist)
    return max_err, mean_err, method

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 accuracy.py <ref.gal> <test.gal> [sample_k] [tolerance]")
        sys.exit(1)
        
    sample_k = int(sys.argv[3]) if len(sys.argv) > 3 else None
    tolerance = float(sys.argv[4]) if len(sys.argv) > 4 else 5e-2
    res = calculate_error(sys.argv[1], sys.argv[2], sample_k)
    if res:
        print(f"Alignment Strategy: {res[2]}")
        print(f"Max Position Error: {res[0]:.2e}")
        print(f"Mean Position Error: {res[1]:.2e}")
        if res[0] < tolerance: 
            print("✅ Precision Check Passed")
        else:
            print("⚠️ Significant Divergence Detected")
