import matplotlib.pyplot as plt

# Initialize variables
data = {}
max_input_size = 400
step_size = 2
input_sizes_full = list(range(10, max_input_size + step_size, step_size))  # Full range of input sizes

# Read the file
with open('l2_error_convergence.txt', 'r') as file:
    lines = file.readlines()

    # Assume each algorithm's name is followed by its times
    for i in range(0, len(lines), 2):
        name = lines[i].strip()  # Get the name of the algorithm
        times = list(map(float, lines[i + 1].strip().split(',')))  # Convert times to float
        input_sizes = input_sizes_full[:len(times)]  # Adjust input_sizes based on the number of timings
        data[name] = (input_sizes, times)  # Map name to (input_sizes, times)

# Check if data is parsed correctly
print(data)

# Create a plot
plt.figure(figsize=(10, 6))

# Plot each algorithm's time against its corresponding input sizes
for name, (input_sizes, times) in data.items():
    plt.plot(input_sizes, times, label=name)

# Set x-axis to logarithmic scale
# plt.xscale('log')
# plt.yscale('log')

# Label the axes
plt.xlabel('N (mesh_fineness)')
plt.ylabel('Time')
plt.title('Algorithm time comparison')

# Display legend
plt.legend()

# Show the plot
plt.show()
