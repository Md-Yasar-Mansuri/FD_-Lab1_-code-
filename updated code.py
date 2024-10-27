# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 10:16:01 2024

@author: MD YASAR
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 18:29:47 2024

@author: MD YASAR
"""

import numpy as np
import matplotlib.pyplot as plt

def FinDiff(n, method, source_term):
    l = 1
    k = 1
    qn = 1
    h = 1 / n
    cmatrix = np.zeros((n + 1, n + 1))
    bvector = np.zeros((n + 1))
    x_points = np.linspace(0, 1, n + 1)

    # Select source terms
    if source_term == 1:
        r_values = np.ones_like(x_points)  # Constant r for first source term
        u_analytical = -((x_points**2)/2) + 2 * x_points
    elif source_term == 2:
        r_values = np.sin(np.pi * x_points / l)  # Sinusoidal r for second source term
        u_analytical = (1 / k) * (np.sin(np.pi * x_points / l) * (l / np.pi) ** 2 + (qn + l / np.pi) * x_points)
    elif source_term == 3:
        r_values = np.sin(np.pi * x_points / (2 * l))  # Sinusoidal r with different frequency for third source term
        u_analytical = (1 / k) * (np.sin(np.pi * x_points / (2 * l)) * (4 * l ** 2 / np.pi ** 2) + x_points)
    else:
        print("Invalid source term selected. Choose 1, 2, or 3.")
        return None

    # Construct coefficient matrix
    for i in range(1, n):
        cmatrix[i, i] = -2     # Adding value diagonally r_c = -2
        cmatrix[i, i - 1] = 1  # r-_j = 1
        cmatrix[i, i + 1] = 1  # r+_j = 1
    for i in range(1, n + 1):
        bvector[i] = -(r_values[i] * (h ** 2) / k)

    cmatrix[0, 0] = 1  # Dirichlet boundary condition

    if method == "backward":
        cmatrix[-1, -1] = 1  # Neumann backward condition
        cmatrix[-1, -2] = -1
        bvector[-1] = qn * h
    elif method == "centered":
        cmatrix[-1, -1] = -2
        cmatrix[-1, -2] = 2
        bvector[-1] = bvector[-1] - (2 * h * qn)

    uvector = np.linalg.solve(cmatrix, bvector)
    error.append(np.max(np.abs(u_analytical - uvector)))
    return uvector, u_analytical

# Define intervals and source term to test
n_intervals = [5, 10, 50, 100, 500, 1000]  # Different values for n
source_terms = [1, 2, 3]  # Select source term 1, 2, or 3

# Plotting for each source term
for source_term in source_terms:
    plt.figure()
    error = []
    for n in n_intervals:
        FinDiff(n, method="centered", source_term=source_term)
    plt.loglog(n_intervals, error, label="NM Centered", linestyle='--', marker='o')
    plt.ylabel('log(error)')
    plt.xlabel('log(1/h)')
    plt.title(f' showing the error for Neumann backward and centered for different step size on source term {source_term}')
    plt.legend()

    error = []
    for n in n_intervals:
        FinDiff(n, method="backward", source_term=source_term)
    plt.loglog(n_intervals, error, label="NM Backward")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.show()

# Test example with n = 10 for centered and backward methods


# Assuming FinDiff function and other necessary code is already defined

# Define parameters
n_values = [4,10,15]  # List of different n values

for n in n_values:
    x_points = np.linspace(0, 1, n + 1)

# Compute solutions
    u_centered, u_analytical = FinDiff(n, method="centered", source_term=1)
    u_backward, _ = FinDiff(n, method="backward", source_term=1)
    
# finding maximum errro between analytical and numarical solutions 
 
    max_error_centered = np.max(np.abs(u_analytical - u_centered))
    max_error_backward = np.max(np.abs(u_analytical - u_backward))
    
    print(f'Maximum error for Neumann (Centered) method with N = {n}: {max_error_centered}')
    print(f'Maximum error for Neumann (Backward) method with N = {n}: {max_error_backward}')


# Plot 1: Analytical vs Numerical (Centered)
    plt.figure()
    plt.plot(x_points, u_analytical, label="Analytical", linestyle='--', color='black')  # Analytical in black
    plt.plot(x_points, u_centered, label="Neumann(Centered)", marker='o', color='green')  # Centered in blue
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Analytical vs Centered Method for  the First Source Term')
    plt.legend()
    plt.show()

# Plot 2: Analytical vs Numerical (Backward)
    plt.figure()
    plt.plot(x_points, u_analytical, label="Analytical", linestyle='--', color='black')  # Analytical in black
    plt.plot(x_points, u_backward, label="Neumann (Backward)", marker='x', color='red')  # Backward in red
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Analytical vs Backward Method for the  First Source Term')
    plt.legend()
    plt.show()
