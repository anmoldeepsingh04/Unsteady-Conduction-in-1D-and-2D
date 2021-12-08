"""
Question1 of Assignment2
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
"""
# Importing relevant libraries and modules
import numpy as np
import math as m
import matplotlib.pyplot as plt

# Defining required quantities
num_x = 30
num_y = 40
num_sol_numerical = np.zeros((num_x + 1, num_y + 1))
x_mesh = np.zeros(num_x + 1)
y_mesh = np.zeros(num_y + 1)


# Function to calculate successive terms of final solution
def sum_series(h, l, x, y, ta, tb, n):
    sums = 0.0
    for k in range(1, int(n)):
        constant = (k * m.pi) / l
        t1 = (1.0 - ((-1.0) ** k)) / k
        t2 = tb * m.sinh(constant * (h - y)) + ta * m.sinh(constant * y)
        t3 = m.sin(constant * x)
        t4 = m.sinh(constant * h)
        sums += (t1 * t2 * t3) / t4
    ans = (2 * sums) / m.pi
    return ans


# Calculating grid points
for i in range(0, num_x + 1):
    x_mesh[i] = i * 0.01
for i in range(0, num_y + 1):
    y_mesh[i] = i * 0.01

reps = np.zeros(150)  # To control the total number of terms in the infinite sum_val

for f in range(1, len(reps)):
    reps[f - 1] = f

# Calculating the exact solution
for p in range(0, len(reps)):
    for o in range(0, num_x + 1):
        for l in range(0, num_y + 1):
            T_val = sum_series(0.4, 0.3, x_mesh[o], y_mesh[l], 10.0, 40.0, p+1)
            num_sol_numerical[o][l] = T_val

    # if p % 4 == 0:  # To control the number of plots being displayed
    #     d = np.transpose(num_sol_numerical)
    #     c = plt.imshow(d, cmap='hot', vmin=d.min(), vmax=d.max(), extent=[0, 30, 0, 40], interpolation='bicubic',
    #                    origin='lower', aspect='auto')
    #     plt.colorbar()
    #     plt.clim(0, 40)
    #     plt.title('Exact Analytical Solution with {0} terms'.format(p), fontweight="bold")
    #     plt.xlabel('$Spatial Dimension$', fontsize=14)
    #     plt.ylabel('$Temperature$', fontsize=14)
    #     plt.get_current_fig_manager().window.state('zoomed')
    #     plt.show(block=False)
    #     plt.pause(0.5)
    #     plt.close()

# Plotting the steady state temperature distribution
d = np.transpose(num_sol_numerical)
c = plt.imshow(d, cmap='hot', vmin=d.min(), vmax=d.max(), extent=[0, 30, 0, 40], interpolation='bicubic',
               origin='lower', aspect='auto')
plt.colorbar()
plt.clim(0, 40)
plt.title('Exact Analytical Solution with {0} terms'.format(p), fontweight="bold")
plt.xlabel('$Spatial Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.get_current_fig_manager().window.state('zoomed')
plt.show()
