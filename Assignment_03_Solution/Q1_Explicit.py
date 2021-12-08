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
num = 20
num1 = 1000
del_x = 1.0 / num
del_x1 = 1.0 / num1
del_t = 0.001
alpha = 1.0
beta = (alpha * del_t) / (del_x ** 2)
L = 1.0
T = 0.5
Ti = 100.0
Ts = 300.0
time_iter = 500

# Defining required arrays
x_mesh = np.zeros(num)
num_sol_old = np.zeros(num)
num_sol_new = np.zeros(num)
num_sol_new_at_100 = np.zeros(num)
num_sol_new_at_200 = np.zeros(num)
num_sol_new_at_300 = np.zeros(num)
num_sol_new_at_400 = np.zeros(num)
num_sol_new_at_500 = np.zeros(num)
x_mesh1 = np.zeros(num1)
exact_sol_at_1 = np.zeros(num1)
exact_sol_at_2 = np.zeros(num1)
exact_sol_at_3 = np.zeros(num1)
exact_sol_at_4 = np.zeros(num1)
exact_sol_at_5 = np.zeros(num1)

# Calculating grid points and setting initial values
for i in range(0, num):
    x_mesh[i] = (i + 0.5) * del_x
    num_sol_old[i] = Ti

# Calculating the numerical solution
for j in range(0, time_iter):
    if j == 0:
        T_A = Ti
        T_B = Ti
        for i in range(0, num):
            if i == 0:
                num_sol_new[i] = 2 * beta * T_A + (1 - 3 * beta) * num_sol_old[i] + beta * num_sol_old[i + 1]

            if 1 <= i <= num - 2:
                num_sol_new[i] = beta * num_sol_old[i - 1] + (1 - 2 * beta) * num_sol_old[i] + beta * num_sol_old[i + 1]

            if i == num - 1:
                num_sol_new[i] = beta * num_sol_old[i - 1] + (1 - 3 * beta) * num_sol_old[i] + 2 * beta * T_B

    else:
        for i in range(0, num):
            T_A = Ts
            T_B = Ts
            if i == 0:
                num_sol_new[i] = 2 * beta * T_A + (1 - 3 * beta) * num_sol_old[i] + beta * num_sol_old[i + 1]

            if 1 <= i <= num - 2:
                num_sol_new[i] = beta * num_sol_old[i - 1] + (1 - 2 * beta) * num_sol_old[i] + beta * num_sol_old[i + 1]

            if i == num - 1:
                num_sol_new[i] = beta * num_sol_old[i - 1] + (1 - 3 * beta) * num_sol_old[i] + 2 * beta * T_B

    for i in range(0, num):
        num_sol_old[i] = num_sol_new[i]
        if j == time_iter * 0.2 - 1:
            num_sol_new_at_100[i] = num_sol_new[i]
        if j == time_iter * 0.4 - 1:
            num_sol_new_at_200[i] = num_sol_new[i]
        if j == time_iter * 0.6 - 1:
            num_sol_new_at_300[i] = num_sol_new[i]
        if j == time_iter * 0.8 - 1:
            num_sol_new_at_400[i] = num_sol_new[i]
        if j == time_iter * 1.0 - 1:
            num_sol_new_at_500[i] = num_sol_new[i]

# Calculating the exact solution
for a in range(0, 1):
    for i in range(0, num1):
        x_mesh1[i] = (i + 0.5) * del_x1

    for i in range(0, num1):
        sum_val = 0.0
        t = T * (1 / 5)
        for j in range(1, 10):
            sum_val += m.exp(-1 * ((j * m.pi) ** 2) * (t)) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
        exact_sol_at_1[i] = Ts + 2 * (Ti - Ts) * sum_val
    for i in range(0, num1):
        sum_val = 0.0
        t = T * (2 / 5)
        for j in range(1, 10):
            sum_val += m.exp(-1 * ((j * m.pi) ** 2) * (t)) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
        exact_sol_at_2[i] = Ts + 2 * (Ti - Ts) * sum_val
    for i in range(0, num1):
        sum_val = 0.0
        t = T * (3 / 5)
        for j in range(1, 10):
            sum_val += m.exp(-1 * ((j * m.pi) ** 2) * (t)) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
        exact_sol_at_3[i] = Ts + 2 * (Ti - Ts) * sum_val
    for i in range(0, num1):
        sum_val = 0.0
        t = T * (4 / 5)
        for j in range(1, 10):
            sum_val += m.exp(-1 * ((j * m.pi) ** 2) * (t)) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
        exact_sol_at_4[i] = Ts + 2 * (Ti - Ts) * sum_val
    for i in range(0, num1):
        sum_val = 0.0
        t = T * (5 / 5)
        for j in range(1, 10):
            sum_val += m.exp(-1 * ((j * m.pi) ** 2) * (t)) * ((1 - (-1) ** j) / (j * m.pi)) * m.sin(j * m.pi * x_mesh1[i])
        exact_sol_at_5[i] = Ts + 2 * (Ti - Ts) * sum_val

# Extending arrays to boundaries for plotting
for i in range(0, 1):
    num_sol_new_at_100 = np.insert(num_sol_new_at_100, 0, Ts)
    num_sol_new_at_100 = np.insert(num_sol_new_at_100, len(num_sol_new_at_100), Ts)
    num_sol_new_at_200 = np.insert(num_sol_new_at_200, 0, Ts)
    num_sol_new_at_200 = np.insert(num_sol_new_at_200, len(num_sol_new_at_200), Ts)
    num_sol_new_at_300 = np.insert(num_sol_new_at_300, 0, Ts)
    num_sol_new_at_300 = np.insert(num_sol_new_at_300, len(num_sol_new_at_300), Ts)
    num_sol_new_at_400 = np.insert(num_sol_new_at_400, 0, Ts)
    num_sol_new_at_400 = np.insert(num_sol_new_at_400, len(num_sol_new_at_400), Ts)
    num_sol_new_at_500 = np.insert(num_sol_new_at_500, 0, Ts)
    num_sol_new_at_500 = np.insert(num_sol_new_at_500, len(num_sol_new_at_500), Ts)
    exact_sol_at_1 = np.insert(exact_sol_at_1, 0, Ts)
    exact_sol_at_1 = np.insert(exact_sol_at_1, len(exact_sol_at_1), Ts)
    exact_sol_at_2 = np.insert(exact_sol_at_2, 0, Ts)
    exact_sol_at_2 = np.insert(exact_sol_at_2, len(exact_sol_at_2), Ts)
    exact_sol_at_3 = np.insert(exact_sol_at_3, 0, Ts)
    exact_sol_at_3 = np.insert(exact_sol_at_3, len(exact_sol_at_3), Ts)
    exact_sol_at_4 = np.insert(exact_sol_at_4, 0, Ts)
    exact_sol_at_4 = np.insert(exact_sol_at_4, len(exact_sol_at_4), Ts)
    exact_sol_at_5 = np.insert(exact_sol_at_5, 0, Ts)
    exact_sol_at_5 = np.insert(exact_sol_at_5, len(exact_sol_at_5), Ts)
    x_mesh = np.insert(x_mesh, 0, 0.0)
    x_mesh = np.insert(x_mesh, len(x_mesh), L)
    x_mesh1 = np.insert(x_mesh1, 0, 0.0)
    x_mesh1 = np.insert(x_mesh1, len(x_mesh1), L)

# Plotting the final numerically calculated solution
plt.plot(x_mesh, num_sol_new_at_100, 'r*', label='Numerical Solution @ t=0.1hr')
plt.plot(x_mesh1, exact_sol_at_1, 'r', label='Exact Solution @ t=0.1hr')
plt.plot(x_mesh, num_sol_new_at_200, 'gd', label='Numerical Solution @ t=0.2hr')
plt.plot(x_mesh1, exact_sol_at_2, 'g', label='Exact Solution @ t=0.2hr')
plt.plot(x_mesh, num_sol_new_at_300, 'ys', label='Numerical Solution @ t=0.3hr')
plt.plot(x_mesh1, exact_sol_at_3, 'y', label='Exact Solution @ t=0.3hr')
plt.plot(x_mesh, num_sol_new_at_400, 'kx', label='Numerical Solution @ t=0.4hr')
plt.plot(x_mesh1, exact_sol_at_4, 'k', label='Exact Solution @ t=0.4hr')
plt.plot(x_mesh, num_sol_new_at_500, 'b^', label='Numerical Solution @ t=0.5hr')
plt.plot(x_mesh1, exact_sol_at_5, 'b', label='Exact Solution @ t=0.5hr')
plt.title('Computational solution and Exact solution for %d grid points using TDMA and Fully Explicit Method' % num)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.legend(fontsize='small', shadow=True, loc='lower right')
plt.grid()
plt.show()
