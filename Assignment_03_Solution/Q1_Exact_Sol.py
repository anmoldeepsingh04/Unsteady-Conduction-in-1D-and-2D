"""
Question2.1 of Assignment3
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
"""

# Importing relevant libraries and modules
import numpy as np
import matplotlib.pyplot as plt
import math as m

# Defining required quantities
num1 = 1000
del_x1 = 1.0 / num1
Ti = 100.0
Ts = 300.0
time_iter = 500

# Defining required arrays
x_mesh1 = np.zeros(num1)
exact_sol_at_1 = np.zeros(num1)
exact_sol_at_2 = np.zeros(num1)
exact_sol_at_3 = np.zeros(num1)
exact_sol_at_4 = np.zeros(num1)
exact_sol_at_5 = np.zeros(num1)

# Calculating grid points
for i in range(0, num1):
    x_mesh1[i] = (i + 0.5) * del_x1

# Calculating the numerical solution
for i in range(0, num1):
    sum_val = 0.0
    t = 0.1
    for j in range(1, 10):
        sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
    exact_sol_at_1[i] = Ts + 2 * (Ti - Ts) * sum_val
for i in range(0, num1):
    sum_val = 0.0
    t = 0.2
    for j in range(1, 10):
        sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
    exact_sol_at_2[i] = Ts + 2 * (Ti - Ts) * sum_val
for i in range(0, num1):
    sum_val = 0.0
    t = 0.3
    for j in range(1, 10):
        sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
    exact_sol_at_3[i] = Ts + 2 * (Ti - Ts) * sum_val
for i in range(0, num1):
    sum_val = 0.0
    t = 0.4
    for j in range(1, 10):
        sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
    exact_sol_at_4[i] = Ts + 2 * (Ti - Ts) * sum_val
for i in range(0, num1):
    sum_val = 0.0
    t = 0.5
    for j in range(1, 10):
        sum_val += m.exp(-1 * ((j * m.pi) ** 2) * t) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
    exact_sol_at_5[i] = Ts + 2 * (Ti - Ts) * sum_val

# Plotting the final numerically calculated solution
plt.plot(x_mesh1, exact_sol_at_1, 'r', label='Exact Solution @ t=0.1hr')
plt.plot(x_mesh1, exact_sol_at_2, 'g', label='Exact Solution @ t=0.2hr')
plt.plot(x_mesh1, exact_sol_at_3, 'y', label='Exact Solution @ t=0.3hr')
plt.plot(x_mesh1, exact_sol_at_4, 'b', label='Exact Solution @ t=0.4hr')
plt.plot(x_mesh1, exact_sol_at_5, 'k', label='Exact Solution @ t=0.5hr')
plt.title('Exact solution')
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.legend(fontsize='small', shadow=True, loc='lower right')
plt.grid()
plt.show()