import sys
import os
import numpy as np

def read_gal(filename):
    # io_write_result in io.c writes:
    # 5 doubles per particle: px, py, m, vx, vy
    file_size = os.path.getsize(filename)
    N = file_size // 40
    data = np.fromfile(filename, dtype=np.float64).reshape(N, 5)
    return data # px, py, m, vx, vy

def calculate_error(file1, file2, sample_k=None):
    d1 = read_gal(file1)
    d2 = read_gal(file2)
    
    if len(d1) != len(d2):
        print(f"Error: Particle count mismatch ({len(d1)} vs {len(d2)})")
        return None
    
    # Sort both datasets by (mass, p_x, p_y) to align reordered particles
    # We use a tolerance for p_x, p_y in sorting if they are slightly different, 
    # but for bit-perfect theta=0, lexsort is fine.
    # For approximation comparisons, we might need a more robust matcher.
    
    # If using sampling, we don't need to sort the whole thing if we are comparing 
    # specific particles, but since d1/d2 might be reordered, we still need 
    # a way to find the "same" particle.
    
    # For now, let's stick to full sort for N < 50k, and implement sampling matcher for large N.
    if len(d1) > 50000 or (sample_k is not None):
        print(f"Info: Using stochastic sampling (K={sample_k or 1000}) for large N={len(d1)}")
        # Match by mass (heuristic, might not be unique but works for our generators)
        # Better: use the original index if we added it, but we didn't.
        # Let's sort both and then sample the sorted list.
        idx1 = np.lexsort((d1[:,2], d1[:,0], d1[:,1]))
        idx2 = np.lexsort((d2[:,2], d2[:,0], d2[:,1]))
        d1 = d1[idx1]
        d2 = d2[idx2]
        
        K = sample_k if sample_k else 1000
        indices = np.random.choice(len(d1), min(K, len(d1)), replace=False)
        d1_s = d1[indices]
        d2_s = d2[indices]
        dist = np.sqrt((d1_s[:,0] - d2_s[:,0])**2 + (d1_s[:,1] - d2_s[:,1])**2)
    else:
        idx1 = np.lexsort((d1[:,2], d1[:,0], d1[:,1]))
        idx2 = np.lexsort((d2[:,2], d2[:,0], d2[:,1]))
        d1_sorted = d1[idx1]
        d2_sorted = d2[idx2]
        dist = np.sqrt((d1_sorted[:,0] - d2_sorted[:,0])**2 + (d1_sorted[:,1] - d2_sorted[:,1])**2)
    
    max_err = np.max(dist)
    mean_err = np.mean(dist)
    return max_err, mean_err

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 accuracy_gate.py <ref.gal> <test.gal> [sample_k]")
        sys.exit(1)
        
    sample_k = int(sys.argv[3]) if len(sys.argv) > 3 else None
    res = calculate_error(sys.argv[1], sys.argv[2], sample_k)
    if res:
        print(f"Max Position Error (Sorted): {res[0]:.2e}")
        print(f"Mean Position Error (Sorted): {res[1]:.2e}")
        # Barnes-Hut with theta=0.5 can have small divergence from O(N^2)
        if res[0] < 5e-2: 
            print("✅ Precision Check Passed")
        else:
            print("⚠️ Significant Divergence Detected")
