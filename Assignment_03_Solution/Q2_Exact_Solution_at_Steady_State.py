"""
Question1 of Assignment2
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
"""
# Importing relevant libraries and modules
import numpy as np
import math as m
import matplotlib.pyplot as plt


# Function to calculate successive terms of final solution
def sum_series(h, l, x, y, ta, tb):
    sums = 0.0
    for k in range(1, 151):
        c = (k*m.pi)/l
        t1 = (1.0 - ((-1.0) ** k)) / k
        t2 = tb*m.sinh(c*(h-y)) + ta*m.sinh(c*y)
        t3 = m.sin(c*x)
        t4 = m.sinh(c*h)
        sums += (t1 * t2 * t3) / t4
    ans = (2*sums)/m.pi
    return ans


# Defining required quantities
num_x = 30
num_y = 40
del_x = 0.01
del_y = 0.01
T_B = 40.0
T_A = 10.0
H = 0.4
L = 0.3
x_mesh = np.zeros(num_x + 1)
y_mesh = np.zeros(num_y + 1)
num_sol_numerical = np.zeros((num_x + 1, num_y + 1))

# Calculating grid points
for i in range(0, num_x + 1):
    x_mesh[i] = i * del_x
for i in range(0, num_y + 1):
    y_mesh[i] = i * del_y

# Calculating the exact solution
for o in range(0, num_x + 1):
    for l in range(0, num_y + 1):
        T_val = sum_series(H, L, x_mesh[o], y_mesh[l], T_A, T_B)
        num_sol_numerical[o][l] = T_val

# Plotting the steady state temperature distribution
d = np.transpose(num_sol_numerical)
c = plt.imshow(d, cmap='hot', vmin=d.min(), vmax=d.max(), extent=[0, 30, 0, 40], interpolation='bicubic',
               origin='lower', aspect='auto')
plt.colorbar()
plt.clim(0, 40)
plt.title('Exact Analytical Solution using first 150 Terms of the Infinite Sum', fontweight="bold")
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.show()
