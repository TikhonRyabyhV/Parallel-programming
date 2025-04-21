import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import struct

# Path to the binary file
print("Enter path for input file: ")
INPUT_FILE = input()

#INPUT_FILE = "./results/conveq.bin"

def read_binary_file(filename):
    with open(filename, "rb") as f:
        # Read K and M (integers)
        K = struct.unpack("i", f.read(4))[0]
        M = struct.unpack("i", f.read(4))[0]
        # Read tau, h, c (doubles)
        tau = struct.unpack("d", f.read(8))[0]
        h = struct.unpack("d", f.read(8))[0]
        c = struct.unpack("d", f.read(8))[0]
        # Read the u grid (K+1) x (M+1) doubles
        u = np.zeros((K + 1, M + 1))
        for k in range(K + 1):
            for m in range(M + 1):
                u[k, m] = struct.unpack("d", f.read(8))[0]
    return K, M, tau, h, c, u

def plot_surface(K, M, tau, h, u):
    # Create coordinate arrays
    t = np.linspace(0, K * tau, K + 1)
    x = np.linspace(0, M * h, M + 1)
    T, X = np.meshgrid(t, x)  # Note: swapped to match u[m, k] indexing

    # Create the 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface (transpose u to match T, X grid)
    surf = ax.plot_surface(T, X, u.T, cmap='viridis', edgecolor='none')
    
    # Add labels and title
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_zlabel('u(x, t)')
    ax.set_title('3D Surface Plot of u(x, t)')
    
    # Add a color bar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    
    if INPUT_FILE == "./results/conveq.bin":
        plt.savefig("./pictures/3D_plot.png")
    elif INPUT_FILE == "./results/parconveq.bin":
        plt.savefig("./pictures/3D_plot_par.png")
    else:
        print("Incorrect path!")


def plot_heatmap(K, M, tau, h, u):
    # Create coordinate arrays
    t = np.linspace(0, K * tau, K + 1)
    x = np.linspace(0, M * h, M + 1)
    
    # Create the 2D heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot the heatmap (transpose u to align with x, t axes)
    cax = ax.imshow(u.T, origin='lower', cmap='viridis', 
                    extent=[0, K * tau, 0, M * h], aspect='auto')
    
    # Add labels and title
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_title('2D Heatmap of u(x, t)')
    
    # Add a color bar
    fig.colorbar(cax, ax=ax, label='u(x, t)')
    
    if INPUT_FILE == "./results/conveq.bin":
        plt.savefig("./pictures/2D_plot.png")
    elif INPUT_FILE == "./results/parconveq.bin":
        plt.savefig("./pictures/2D_plot_par.png")
    else:
        print("Incorrect path!")


def main():
    try:
        # Read the binary file
        K, M, tau, h, c, u = read_binary_file(INPUT_FILE)
        print(f"Read file: K={K}, M={M}, tau={tau}, h={h}, c={c}")
        
        # Plot the 3D surface
        plot_surface(K, M, tau, h, u)
        
        # Plot the 2D heatmap
        plot_heatmap(K, M, tau, h, u)
        
    except FileNotFoundError:
        print(f"Error: File {INPUT_FILE} not found.")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()
