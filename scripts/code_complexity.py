import math
from statistics import mean
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
import subprocess

# Define the form of the function you want to fit
def func(x, a, b, c):
    return a * x**2 + b * x + c

# Values of N to test
N_values = [10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]

# List to store execution times
times = []

for N in N_values:
    print(f"Running for N={N}")
    # Generate the filename based on N
    filename = f"./input_data/ellipse_N_{N:05d}.gal"

    # Command to run the code
    cmd = ["./galsim", str(N), filename, "200", "1e-5", "8"]

    # List to store the execution times for each run
    run_times = []

    for _ in range(10):
        # Measure the execution time
        start_time = time.time()
        subprocess.run(cmd)
        end_time = time.time()

        # Store the execution time
        run_times.append((end_time - start_time))

    # Store the minimum execution time
    times.append(mean(run_times)/N)


# Fit the data using curve_fit
popt, pcov = curve_fit(func, N_values, times)

# Generate y-values based on the fit
y_fit = [func(x, *popt) for x in N_values]

# Plot the execution times
plt.figure(figsize=(10, 6))

# Plot the scatter plot
plt.scatter(N_values, times, label='Execution time')

# Plot the fitted curve
plt.plot(N_values, y_fit, 'r-', label='Fitted Curve')

plt.xlabel('N values')
plt.ylabel('Execution time')
# Add a legend
plt.legend()

# Save the figure
plt.savefig('code_complexity.png')