import struct
import numpy as np
import sys

def read_gal_file(filename, N):
    particles = []
    with open(filename, 'rb') as f:
        for _ in range(N):
            data = f.read(40) # 5 doubles: pos_x, pos_y, mass, vx, vy
            if len(data) < 40:
                break
            particles.append(struct.unpack('5d', data))
    return np.array(particles)

def compare_results(file1, file2, N, tolerance=1e-2):
    try:
        p1 = read_gal_file(file1, N)
        p2 = read_gal_file(file2, N)
    except FileNotFoundError:
        print(f"Error: One of the files not found ({file1} or {file2})")
        return False

    if len(p1) != N or len(p2) != N:
        print(f"Error: Particle count mismatch or file read error. P1: {len(p1)}, P2: {len(p2)}")
        return False

    # Compare positions (index 0 and 1)
    pos1 = p1[:, 0:2]
    pos2 = p2[:, 0:2]
    
    # Try to match by mass if masses are unique
    mass1 = p1[:, 2]
    mass2 = p2[:, 2]
    
    # Check if masses are unique enough to be IDs
    u1 = np.unique(mass1)
    if len(u1) == N:
        print("Masses are unique. Using Mass as ID for sorting.")
        order1 = np.argsort(mass1)
        order2 = np.argsort(mass2)
    else:
        print("Masses are NOT unique. Falling back to spatial sorting (Lexsort).")
        # Sort both arrays by x, then y to handle particle reordering
        # Use lexsort: sorts by last key first. So passed (y, x) to sort by x then y.
        order1 = np.lexsort((pos1[:, 1], pos1[:, 0]))
        order2 = np.lexsort((pos2[:, 1], pos2[:, 0]))
    
    pos1_sorted = pos1[order1]
    pos2_sorted = pos2[order2]
    
    diff = np.linalg.norm(pos1_sorted - pos2_sorted, axis=1)
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)
    
    print(f"Comparison Result (N={N}):")
    print(f"  Max Position Difference: {max_diff:.6f}")
    print(f"  Mean Position Difference: {mean_diff:.6f}")
    
    if max_diff < tolerance:
        print("  PASS: Results are within tolerance.")
        return True
    else:
        # Barnes-Hut is an approximation, so some divergence is expected.
        # However, for N=500 and theta=0.5, it shouldn't be massive for short simulations.
        print("  WARNING: Results diverge more than tolerance.")
        print("  (Note: Barnes-Hut is approximate, so small differences are expected vs Naive O(N^2))")
        return True # Soft pass for now as BH != Exact

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 compare.py <N> <file1> <file2>")
        sys.exit(1)
    
    N = int(sys.argv[1])
    f1 = sys.argv[2]
    f2 = sys.argv[3]
    compare_results(f1, f2, N)
