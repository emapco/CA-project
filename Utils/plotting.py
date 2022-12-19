import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import csv

def plot_last(log_file: str):
    """Showing the plot of the last tensor grid.
    Args:
        log_file (str): The log file generated from the CA.
    """
    with open(log_file,'r') as f:
        reader = csv.reader(f)
        # We do not need the n_state in this function
        next(reader)
        # Read the second line (which contains the shape of the tensor)
        shape_line = next(reader)
        # Convert the shape line to a list of integers
        shape = [int(x) for x in shape_line[0:3]]
        size = shape[0] * shape[1] * shape[2]
        data_line = None
        for line in reader:
            data_line = line
        # Convert the data line to a NumPy array
        data = np.array(data_line[:size])
        # Reshape the data into a tensor with the specified shape
        tensor = np.reshape(data, shape)
    # Get the x, y, and z coordinates for each element in the tensor
    x, y, z = np.indices((shape[0], shape[1], shape[2]))
    # Flatten the coordinates into 1D arrays
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    # Create a figure and a 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Plot the tensor, using the flattened tensor as the size of the points
    ax.scatter(x, y, z, c=tensor, cmap='viridis')
    plt.show()

def plot_step(result_file: str):
    """Show the plot with the n_states in one plot, be able to provide interesting information for t_end 
    Args:
        result_file (str): The Result file generated by the c++ function
    """
    data  = []
    n_steps = 0
    with open(result_file,'r') as f:
        reader = csv.reader(f)
        # We need the number of states in the CA
        n_states = int(next(reader)[0])
        for i in range(n_states):
            sublist = []
            data.append(sublist)
        next(reader)
        for line in reader:
            n_steps += 1
            for i in range(n_states):
                data[i].append(line[i])
    fig, ax = plt.subplots()
    x = []
    for i in range(n_steps):
        x.append(i)
    # Plot the lines
    for i in range(n_states):
        ax.plot(x,data[i], label=f'Line {i}')
    ax.legend(loc='upper left')
    plt.show()

def plot_last_galaxy(log_file: str):
    """Showing the plot of the last tensor grid, especially for the galaxy cell.
    Args:
        log_file (str): The log file generated from the Galaxy Cell.
    """
    with open(log_file,'r') as f:
        reader = csv.reader(f)
        # We do not need the n_state in this function
        next(reader)
        # Read the second line (which contains the shape of the tensor)
        shape_line = next(reader)
        # Convert the shape line to a list of integers
        shape = [int(x) for x in shape_line[0:3]]
        size = shape[0] * shape[1] * shape[2]
        data_line = None
        for line in reader:
            data_line = line
        # Convert the data line to a NumPy array
        data = np.array(data_line[:size])
        # Reshape the data into a tensor with the specified shape
        tensor = np.reshape(data, shape)
    # Get the x, y, and z coordinates for each element in the tensor
    x, y, z = np.indices((shape[0], shape[1], shape[2]))
    # Flatten the coordinates into 1D arrays
    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    # Flatten the tensor into a 1D array
    sizes = tensor.flatten()
    # Create a figure and a 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Plot the tensor, using the flattened tensor as the size of the points
    ax.scatter(x, y, z, c=tensor, cmap='viridis', size=sizes*100)
    plt.show()
