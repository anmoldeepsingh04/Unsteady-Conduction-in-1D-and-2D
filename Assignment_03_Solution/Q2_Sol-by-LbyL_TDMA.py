"""
Question2.1 of Assignment3
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
Note: The following code automatically displays the temperature distribution at different time steps.
Do not close the plots in between the simulation as it will lead to termination with error
Kindly read the relevant comments for more information
"""

# Importing relevant libraries and modules
import timeit  # To calculate the execution time
import numpy as np
import math as m
import matplotlib.pyplot as plt

start = timeit.default_timer()  # Starting the timer as code starts from here


# defining the TDMA Algorithm
def tdma(num, low, up, dia, func, sol):
    d1 = np.zeros(num, float)
    rhs1 = np.zeros(num, float)
    d1[0] = dia[0]
    rhs1[0] = func[0]
    for h in range(1, len(dia)):
        d1[h] = dia[h] - (low[h] * up[h - 1] / d1[h - 1])
        rhs1[h] = func[h] - (low[h] * rhs1[h - 1] / d1[h - 1])
    sol[num - 1] = rhs1[num - 1] / d1[num - 1]
    for c in range(len(dia) - 2, -1, -1):
        sol[c] = (rhs1[c] - up[c] * sol[c + 1]) / d1[c]


# Function to calculate successive terms of final solution
def sum_series(h, l, x, y, ta, tb):
    sums = 0.0
    for k in range(1, 100):
        constant = (k * m.pi) / l
        t1 = (1.0 - ((-1.0) ** k)) / k
        t2 = tb * m.sinh(constant * (h - y)) + ta * m.sinh(constant * y)
        t3 = m.sin(constant * x)
        t4 = m.sinh(constant * h)
        sums += (t1 * t2 * t3) / t4
    ans = (2 * sums) / m.pi
    return ans


# Defining required quantities
num_x = 30
num_y = 40
del_x = 0.01
del_y = 0.01
del_t = 1.0
alpha = 0.00011234
k = 380.0
beta = (alpha * del_t) / (del_x ** 2)
T_B = 40.0
T_A = 10.0
H = 0.4
L = 0.3
converged = False
error = 1.0e-5

# Defining required arrays
x_mesh = np.zeros(num_x + 1)
y_mesh = np.zeros(num_y + 1)
num_sol_old = np.zeros((num_x - 1) * (num_y - 1))
num_sol_new = np.zeros((num_x - 1) * (num_y - 1))
num_sol_numerical = np.zeros((num_x + 1, num_y + 1))
main_dia = np.zeros(num_y - 1)
lower_dia = np.zeros(num_y - 1)
upper_dia = np.zeros(num_y - 1)
rhs = np.zeros(num_y - 1)
u_E = np.zeros(num_y - 1)
u_W = np.zeros(num_y - 1)
Monitor_point1 = np.array([])
Monitor_point2 = np.array([])
Monitor_point3 = np.array([])
Monitor_point4 = np.array([])
Monitor_point5 = np.array([])
Monitor_point6 = np.array([])
Monitor_point7 = np.array([])
Monitor_point8 = np.array([])
Monitor_point9 = np.array([])
Monitor_point10 = np.array([])

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

# Converting the 2D Exact solution array to 1D array to compare with numerical solution
num_check = np.zeros((num_x - 1) * (num_y - 1))
flat = num_sol_numerical.flatten()
for h in range(0, (num_x - 1) * (num_y - 1)):
    num_check[h] = flat[h + 38]

# Defining a counter to count the number of iterations
iter_count = 0

# Calculating the numerical solution
while not converged:
    iter_count = iter_count + 1
    for k in range(0, num_x - 1):
        temp_sol = np.zeros(num_y - 1)
        # Solving for left boundary point

        # Inserting values in diagonals
        if k == 0:
            for f in range(0, len(u_W)):
                u_W[f] = 0.0
            for m in range(0, len(u_E)):
                u_E[m] = num_sol_old[m + (k + 1) * (num_y - 1)]
            for i in range(1, num_y - 2):
                lower_dia[i] = -beta
                main_dia[i] = 1.0 + 4.0 * beta
                upper_dia[i] = -beta
                rhs[i] = beta * (u_W[i] + u_E[i]) + num_sol_old[i + k * (num_y - 1)]

            # Defining the boundary conditions
            def apply_bc_dia():
                lower_dia[0] = 0.0
                upper_dia[0] = -beta
                main_dia[0] = 1.0 + 4.0 * beta
                rhs[0] = beta * (T_B + u_W[0] + u_E[0]) + num_sol_old[0]
                lower_dia[num_y - 2] = -beta
                upper_dia[num_y - 2] = 0.0
                main_dia[num_y - 2] = 1.0 + 4.0 * beta
                rhs[num_y - 2] = beta * (T_A + u_W[num_y - 2] + u_E[num_y - 2]) + num_sol_old[num_y - 2]


            # Applying the boundary conditions
            apply_bc_dia()

            # Solving the system
            tdma(num_y - 1, lower_dia, upper_dia, main_dia, rhs, temp_sol)
            for g in range(0, len(temp_sol)):
                num_sol_new[g] = temp_sol[g]

        # Solving for right boundary point

        # Inserting values in diagonals
        elif k == num_x - 2:
            for f in range(0, len(u_W)):
                u_W[f] = num_sol_new[f + (num_y - 1) * k]
            for m in range(0, len(u_E)):
                u_E[m] = 0.0
            for i in range(1, num_y - 2):
                lower_dia[i] = -beta
                main_dia[i] = 1.0 + 4.0 * beta
                upper_dia[i] = -beta
                rhs[i] = beta * (u_W[i] + u_E[i]) + num_sol_old[i + (num_x - 1) * (num_y - 1) - num_y]

            # Defining the boundary conditions
            def apply_bc_dia():
                lower_dia[0] = 0.0
                upper_dia[0] = -beta
                main_dia[0] = 1.0 + 4.0 * beta
                rhs[0] = beta * (T_B + u_W[0] + u_E[0]) + num_sol_old[(num_x - 1) * (num_y - 1) - num_y]
                lower_dia[num_y - 2] = -beta
                upper_dia[num_y - 2] = 0.0
                main_dia[num_y - 2] = 1.0 + 4.0 * beta
                rhs[num_y - 2] = beta * (T_A + u_W[num_y - 2] + u_E[num_y - 2]) + num_sol_old[
                    (num_x - 1) * (num_y - 1) - 1]


            # Applying the boundary conditions
            apply_bc_dia()

            # Solving the system
            tdma(num_y - 1, lower_dia, upper_dia, main_dia, rhs, temp_sol)

            for g in range(0, len(temp_sol)):
                num_sol_new[g + (num_x - 1) * (num_y - 1) - num_y] = temp_sol[g]

        # Solving for internal points

        # Inserting values in diagonals
        else:
            for f in range(0, len(u_W)):
                u_W[f] = num_sol_new[f + (k - 1) * (num_y - 1)]
            for m in range(0, len(u_E)):
                u_E[m] = num_sol_old[m + (k + 1) * (num_y - 1)]
            for i in range(1, num_y - 2):
                lower_dia[i] = -beta
                main_dia[i] = 1.0 + 4.0 * beta
                upper_dia[i] = -beta
                rhs[i] = beta * (u_W[i] + u_E[i]) + num_sol_old[i + k * (num_y - 1)]

            # Defining the boundary conditions
            def apply_bc_dia():
                lower_dia[0] = 0.0
                upper_dia[0] = -beta
                main_dia[0] = 1.0 + 4.0 * beta
                rhs[0] = beta * (T_B + u_W[0] + u_E[0]) + num_sol_old[k * (num_y - 1)]
                lower_dia[num_y - 2] = -beta
                upper_dia[num_y - 2] = 0.0
                main_dia[num_y - 2] = 1.0 + 4.0 * beta
                rhs[num_y - 2] = beta * (T_A + u_W[num_y - 2] + u_E[num_y - 2]) + num_sol_old[
                    k * (num_y - 1) + (num_y - 2)]

            # Applying the boundary conditions
            apply_bc_dia()

            # Solving the system
            tdma(num_y - 1, lower_dia, upper_dia, main_dia, rhs, temp_sol)
            for g in range(0, len(temp_sol)):
                num_sol_new[g + k * (num_y - 1)] = temp_sol[g]

    # Transferring the values of present solution matrix into previous solution matrix
    for p in range(0, len(num_sol_new)):
        num_sol_old[p] = num_sol_new[p]

    # Checking if the solution is converged by employing L(infinity) norm at each grid point
    '''In case iteration count exceeds 15000, stopping as it provides a reasonably close answer'''
    for p in range(0, len(num_sol_new)):
        if abs(num_sol_new[p] - num_check[p]) < error or iter_count == 15000:
            converged = True
    print(iter_count)

    # Inserting values for different monitor points on the entire plate
    Monitor_point1 = np.append(Monitor_point1, [num_sol_new[157]], axis=0)
    Monitor_point2 = np.append(Monitor_point2, [num_sol_new[586]], axis=0)
    Monitor_point3 = np.append(Monitor_point3, [num_sol_new[976]], axis=0)
    Monitor_point4 = np.append(Monitor_point4, [num_sol_new[392]], axis=0)
    Monitor_point5 = np.append(Monitor_point5, [num_sol_new[742]], axis=0)
    Monitor_point6 = np.append(Monitor_point6, [num_sol_new[1054]], axis=0)
    Monitor_point7 = np.append(Monitor_point7, [num_sol_new[253]], axis=0)
    Monitor_point8 = np.append(Monitor_point8, [num_sol_new[387]], axis=0)
    Monitor_point9 = np.append(Monitor_point9, [num_sol_new[621]], axis=0)
    Monitor_point10 = np.append(Monitor_point10, [num_sol_new[894]], axis=0)

    # For plotting the solution at intermediate time steps
    if iter_count % 200 == 0:  # Controls the number of time iterations that you want to see the plot for
        for l in range(0, 1):
            # Plotting the required results
            temp_arr = np.zeros(39)
            x = np.append(temp_arr, num_sol_new)
            y = np.append(x, temp_arr)
            real = np.zeros(len(y))
            for i in range(0, len(y)):
                real[i] = y[i]
            f_arr = real.reshape(31, 39)
            arr = np.transpose(f_arr)
            temp_10 = np.zeros(31)
            temp_40 = np.zeros(31)
            for i in range(0, 31):
                if i == 0 or i == 31:
                    temp_10[i] = 0.0
                    temp_40[i] = 0.0
                temp_10[i] = 10.0
                temp_40[i] = 40.0
            temp_arr_1 = np.r_[[temp_40], arr]
            Final_sol = np.r_[temp_arr_1, [temp_10]]
            x_mesh1 = np.zeros(len(y))
            print(len(real))
            c = plt.imshow(Final_sol, cmap='hot', vmin=arr.min(), vmax=arr.max(), extent=[0, 30, 0, 40], interpolation='bicubic', origin='lower', aspect='auto')
            plt.colorbar()
            plt.clim(0, 40)
            plt.title('Numerical Solution at time step {0}'.format(iter_count), fontweight="bold")
            plt.xlabel('$x$', fontsize=14)
            plt.ylabel('$Temperature$', fontsize=14)
            plt.get_current_fig_manager().window.state('zoomed')
            plt.show(block=False)
# Do not close the plots in between the simulation as it will lead to termination with error
            plt.pause(1.5)  # Waiting time before next plot.
            plt.close()
            '''plt.savefig("plt{val}".format(val=j), dpi = 200 )'''  # For saving the plot

# Plotting the final numerically calculated solution
for l in range(0, 1):
    temp_arr = np.zeros(39)
    x = np.append(temp_arr, num_sol_new)
    y = np.append(x, temp_arr)
    real = np.zeros(len(y))
    for i in range(0, len(y)):
        real[i] = y[i]
    arr = real.reshape(31, 39)
    arr = np.transpose(arr)
    x_mesh1 = np.zeros(len(y))
    temp_10 = np.zeros(31)
    temp_40 = np.zeros(31)
    for i in range(0, 31):
        if i == 0 or i == 31:
            temp_10[i] = 0.0
            temp_40[i] = 0.0
        temp_10[i] = 10.0
        temp_40[i] = 40.0
    temp_arr_1 = np.r_[[temp_40], arr]
    Final_sol = np.r_[temp_arr_1, [temp_10]]

    x_mesh1 = np.zeros(len(y))
    c = plt.imshow(Final_sol, cmap='hot', vmin=arr.min(), vmax=arr.max(), extent=[0, 30, 0, 40], interpolation='bicubic', origin='lower', aspect='auto')
    plt.colorbar()
    plt.clim(0, 40)
    plt.title('Numerical Solution at time step {0}'.format(iter_count), fontweight="bold")
    plt.xlabel('$Spatial Dimension$', fontsize=14)
    plt.ylabel('$Temperature$', fontsize=14)
    plt.get_current_fig_manager().window.state('zoomed')
    plt.show()

# Plotting the evolution of temperature at the monitor points
time_step = np.zeros(iter_count)
for i in range(0, len(time_step)):
    time_step[i] = i
plt.plot(time_step, Monitor_point1, label='Monitor Point 1')
plt.plot(time_step, Monitor_point2, label='Monitor Point 2')
plt.plot(time_step, Monitor_point3, label='Monitor Point 3')
plt.plot(time_step, Monitor_point4, label='Monitor Point 4')
plt.plot(time_step, Monitor_point5, label='Monitor Point 5')
plt.plot(time_step, Monitor_point6, label='Monitor Point 6')
plt.plot(time_step, Monitor_point7, label='Monitor Point 7')
plt.plot(time_step, Monitor_point8, label='Monitor Point 8')
plt.plot(time_step, Monitor_point9, label='Monitor Point 9')
plt.plot(time_step, Monitor_point10, label='Monitor Point 10')
plt.title('Monitor Point Values at time step {0}'.format(iter_count), fontweight="bold")
plt.xlabel('$Temporal Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.legend(fontsize='small', shadow=True, loc='lower right')
plt.grid()
plt.get_current_fig_manager().window.state('zoomed')
plt.show()

# Plotting the final exact solution at steady state
d = np.transpose(num_sol_numerical)
c = plt.imshow(d, cmap='hot', vmin=d.min(), vmax=d.max(), extent=[0, 30, 0, 40], interpolation='bicubic',
               origin='lower', aspect='auto')
plt.colorbar()
plt.clim(0, 40)
plt.title('Exact Analytical Solution', fontweight="bold")
plt.xlabel('$Spatial Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.get_current_fig_manager().window.state('zoomed')
plt.show()

# Program execution ends here. Displaying the time of execution and number of iterations
stop = timeit.default_timer()
execution_time = stop - start

print("Program Executed in " + str(execution_time), "seconds")
print("Program Executed took ", iter_count, "iterations")
