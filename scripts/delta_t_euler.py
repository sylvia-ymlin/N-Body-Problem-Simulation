import subprocess
import matplotlib.pyplot as plt
import numpy as np

# Change the executable inputs to run with different nsteps and delta_t
command_template = "./galsim 3 input_data/sun_and_planets_N_3.gal {} {} 0"
compare_command = "./compare_gal_files/compare_gal_files 3 result.gal ref_output_data/sun_and_planets_3.gal"

# nsteps and delta_t values
nsteps_values = [20*i for i in range(1, 100)]
delta_t_values = [1e-3 / i for i in range(1, 100)]
nsteps_values.append(30)
delta_t_values.append(1e-3 / 1.5)

pos_maxdiff_values = []
# Run the executable for each nsteps and delta_t
for nsteps, delta_t in zip(nsteps_values, delta_t_values):
    command = command_template.format(nsteps, delta_t)
    subprocess.run(command, shell=True)  # Run the command

    compare_result = subprocess.run(compare_command, shell=True, capture_output=True, text=True)
    output_lines = compare_result.stdout.strip().split('\n')  # Split the output into lines
    pos_maxdiff_line = [line for line in output_lines if line.startswith('pos_maxdiff')][0]  # Find the line with pos_maxdiff
    pos_maxdiff_value = float(pos_maxdiff_line.split('=')[1].strip())  # Extract the value after '='
    pos_maxdiff_values.append(pos_maxdiff_value)

print(pos_maxdiff_values)


# # Fit a linear function to the data
# log_delta_t_values = np.log(delta_t_values)
# log_pos_maxdiff_values = np.log(pos_maxdiff_values)
# slope, intercept = np.polyfit(log_delta_t_values, log_pos_maxdiff_values, 1)
# fit_line = np.exp(intercept) * np.power(delta_t_values, slope)
# label='Linear Fit'

# Fit a quadratic function to the data
delta_t_values = np.array(delta_t_values)
coefficients = np.polyfit(delta_t_values, pos_maxdiff_values, 2)
# Generate a larger set of x-values for a smoother curve
x_smooth = np.linspace(min(delta_t_values), max(delta_t_values), 500)

# Evaluate the polynomial at the points x_smooth
fit_line = coefficients[0]*np.power(x_smooth, 2) + coefficients[1]*x_smooth + coefficients[2]
label='Quadratic Fit'

# Plot the results
plt.figure(figsize=(10, 6))
# plt.xscale('log')
# plt.yscale('log')
plt.gca().invert_xaxis()
plt.scatter(delta_t_values, pos_maxdiff_values)
plt.plot(x_smooth, fit_line, color='red', linestyle='--', linewidth=1, label=label)
plt.xlabel(r'$\triangle$ t', fontsize=18)
plt.ylabel('pos_maxdiff', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.grid(True)  # Add grid

# Change the name here if you modify N, else it would overwrite the previous plot
plt.savefig('delta_t_euler.png', dpi=300, bbox_inches='tight')
plt.legend()
plt.show()
