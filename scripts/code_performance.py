import subprocess
import time
import os
import matplotlib.pyplot as plt

def measure_execution_time(command, env=None):
    start_time = time.time()
    subprocess.run(command, shell=True, env=env)
    end_time = time.time()
    return end_time - start_time

# Compile the C program
# subprocess.run("gcc -fopenmp galsim.c -o galsim", shell=True)

# Measure execution time in serial and in parallel with different numbers of threads and different N
N_values = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
num_threads_list = [1, 2, 4, 8, 16]

plt.figure(figsize=(10, 6))

num_runs = 5  # number of times to run the command for each N

for num_threads in num_threads_list:
    execution_times = []
    for N in N_values:
        min_execution_time = float('inf')  # initialize to infinity
        for _ in range(num_runs):

            command = f"./galsim {N} ./input_data/ellipse_N_{N:05d}.gal 200 1e-5 {num_threads}"

            execution_time = measure_execution_time(command)
            min_execution_time = min(min_execution_time, execution_time)  # update if current run is faster
        execution_times.append(min_execution_time)
        print(f"Minimum execution time with {num_threads} threads and N={N}: {min_execution_time} seconds")
    
    # Plot the results for this thread count
    plt.plot(N_values, execution_times, marker='o', label=f'{num_threads} threads')

plt.xlabel('N')
plt.ylabel('Execution time (seconds)')
plt.title('Execution time vs N')
plt.grid(True)
plt.legend()
plt.savefig('execution_time_vs_N.png', dpi=300)
plt.show()