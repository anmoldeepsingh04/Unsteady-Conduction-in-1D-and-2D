"""
Question2.1 of Assignment3
Completed by: Anmoldeep Singh 180030002
Dated: 30-04-2021
"""

# Importing relevant libraries and modules
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
import math as m


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


def gauss_seidel(A, b, tolerance, x):
    iter1 = 0  # Counter to count the number of iterations
    while True:
        iter1 = iter1 + 1
        x_old = x.copy()
        for h in range(A.shape[0]):
            x[h] = (b[h] - np.dot(A[h, :h], x[:h]) - np.dot(A[h, (h + 1):], x_old[(h + 1):])) / A[h, h]
        # Stopping criteria, evaluating L(inf) norm at each point
        if max(abs(x_old)) == 0:
            continue
        L_norm_Inf = max(abs((x - x_old))) / max(abs(x_old))
        if L_norm_Inf < tolerance:
            break
    print("Total number of iterations: " + str(iter1))
    return x


# Defining required quantities
num_x = 30
num_y = 40
del_x = 0.01
del_y = 0.01
del_t = 1.0
t = num_x + num_y
alpha = 0.00011234
k = 380.0
beta = (alpha * del_t) / (del_x ** 2)
T_B = 40.0
T_A = 10.0
H = 0.4
L = 0.3
tol = 1.0e-5

# Defining required arrays
num_sol_calc = np.zeros((num_x, num_y))
num_sol_calc1 = np.zeros((num_x + 1, num_y + 1))
coeff_matrix = np.zeros(((num_x - 1) * (num_y - 1), (num_x - 1) * (num_y - 1)))
b = np.zeros(((num_x - 1) * (num_y - 1), 1))
x_mesh = np.zeros(num_x + 1)
y_mesh = np.zeros(num_y + 1)
init_guess_temp = 0.0000125
temp_values = init_guess_temp*np.ones(((num_x - 1) * (num_y - 1), 1))

# Calculating grid points
for i in range(0, num_x + 1):
    x_mesh[i] = i * del_x
for i in range(0, num_y + 1):
    y_mesh[i] = i * del_y

# Calculating the exact solution
for o in range(0, num_x):
    for l in range(0, num_y):
        T_val = sum_series(H, L, x_mesh[o], y_mesh[l], T_A, T_B)
        num_sol_calc[o][l] = T_val
for o in range(0, num_x + 1):
    for l in range(0, num_y + 1):
        T_val = sum_series(H, L, x_mesh[o], y_mesh[l], T_A, T_B)
        num_sol_calc1[o][l] = T_val

c = b.reshape((num_x - 1), (num_y - 1))
# Inserting values in RHS matrix
for i in range(0, num_x - 1):
    c[i][0] = 40.0
    c[i][num_y - 2] = 10.0
b = b.reshape((num_x - 1) * (num_y - 1), 1)

# Inserting values in coefficient matrix
coeff_matrix[1][2] = -beta
coeff_matrix[(num_x - 1) * (num_y - 1) - 2][(num_x - 1) * (num_y - 1) - 2] = -beta

for i in range(0, (num_x - 1) * (num_y - 1)):
    coeff_matrix[i][i] = 1.0 + 4 * beta
for i in range(0, (num_x - 1) * (num_y - 1) - 1):
    coeff_matrix[i][i + 1] = -beta
    try:
        coeff_matrix[i][i + 4] = -beta
    except:
        pass
for i in range(1, (num_x - 1) * (num_y - 1) - 1):
    coeff_matrix[i][i - 1] = -beta
    try:
        coeff_matrix[i + 3][i - 1] = -beta
    except:
        pass
for i in range(1, (num_x - 1) * (num_y - 1) - 1):
    if (2 * i + 1) % t == 0 or (2 * i) % t == 0:
        coeff_matrix[i][i + 1] = 0.0
        coeff_matrix[i + 1][i] = 0.0

# Computing the numerical solution
gauss_seidel(coeff_matrix, b, tol, temp_values)
temp_values = solve(coeff_matrix, b)
temp_values = np.transpose(num_sol_calc)

# Plotting the final exact solution at steady state
c = plt.imshow(temp_values, cmap='hot', vmin=temp_values.min(), vmax=temp_values.max(), extent=[0, 30, 0, 40],
               interpolation='bicubic', origin='lower', aspect='auto')
plt.colorbar()
plt.title('Numerical Solution using Gauss-Seidel Point-by-Point Iterative Method', fontweight="bold")
plt.xlabel('$Spatial Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.get_current_fig_manager().window.state('zoomed')
plt.show()

# Plotting exact solution
c = plt.imshow(np.transpose(num_sol_calc1), cmap='hot', vmin=num_sol_calc.min(), vmax=num_sol_calc.max(),
               extent=[0, 30, 0, 40], interpolation='bicubic', origin='lower', aspect='auto')
plt.colorbar()
plt.title('Exact Analytical Solution', fontweight="bold")
plt.xlabel('$Spatial Dimension$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.get_current_fig_manager().window.state('zoomed')
plt.show()
