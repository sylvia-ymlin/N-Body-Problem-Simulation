import sys
import math
import numpy as np

def calculate_adjusted_efficiency(t1, tp, p, n_base):
    """
    Efficiency_adjusted = (T1(N) / Tp(p*N)) * ((p*N) * log(p*N)) / (p * (N * log(N)))
    """
    # Assuming log base 2 for simplicity as it's a ratio
    work_ratio = (p * n_base * math.log2(p * n_base)) / (p * (n_base * math.log2(n_base)))
    raw_efficiency = t1 / tp
    adjusted_efficiency = raw_efficiency * work_ratio
    return adjusted_efficiency

def report_scaling(results):
    """
    results: list of tuples (threads, N, runtime)
    """
    print(f"{'Threads':<10} | {'N':<10} | {'Runtime':<10} | {'Efficiency (Adj)':<15}")
    print("-" * 55)
    
    t1 = results[0][2]
    n_base = results[0][1]
    
    for threads, n, runtime in results:
        p = threads
        if p == 1:
            eff = 1.0
        else:
            # For weak scaling, N scales with p
            eff = calculate_adjusted_efficiency(t1, runtime, p, n_base)
        
        print(f"{threads:<10} | {n:<10} | {runtime:<10.4f} | {eff:<15.2%}")

if __name__ == "__main__":
    # Example usage / mock for now
    # results = [(1, 1000, 1.0), (2, 2000, 1.1), (4, 4000, 1.2)]
    # report_scaling(results)
    pass
