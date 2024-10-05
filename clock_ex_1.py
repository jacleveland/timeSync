import numpy as np
from scipy.linalg import svd
import matplotlib.pyplot as plt
"""
    Example One:
    
    Flow of time synchronization:
    1 --> 2 < --- > 3

    Tells us that the border B = {1} and the non-border values Y = V\B = {2, 3}

    Graph Laplacian:

    L(G) =   1  2  3
        1   [[0, 0, 0], 
        2    [-1, 2, -1], 
        3    [0, -1, 1]]

    A = I - h*alpha*L(G)


    b = [1, 0] * (h*alpha*x_b)

    Starting times: x_1 = 40023.054, x_2 = 40023.043, x_3 = 40023.067
"""

# h = dt (ish), alpha between (0, 1/h*d) s.t. d is largest diag value
# here, d = 2 and let h = 1/100 --> alpha in (0, 50)
params = [.001, 0]

# Graph Laplacian (non-boundary values)
m = np.array([[2, -1], [-1,1]])

# Times (in seconds)
x1 = 40023.054
x2 = 40023.045
x3 = 40023.067
clocks = [x1, x2, x3]
xv = np.array([[x2-x1], [x3-x1]])
print(f"Initial times: clock A: {x1:.3f} clock B: {x2:.3f} clock C: {x3:.3f}")

# one euler step
def euler_step(xv, m):
    A = np.identity(2) - params[0]*params[1]*m
    xv = A @ xv

    return xv

"""
    Alternative methods
"""
def rk2_step(xv, m):
    t = params[0]*params[1]
    xv = xv - t * m @ (xv - t * 0.5 * m @ xv)

    return xv

def rk4_step(xv, m):
    t = -1*params[0] * params[1]
    k1 = m @ xv
    k2 = m @ (xv + 0.5*t*k1)
    k3 = m @ (xv + 0.5*t*k2)
    k4 = m @ (xv + t*k3)
    xv = xv + (t/6) * (k1 + 2*k2 + 2*k3 + k4)

    return xv

def trap_step(xv, m):
    mfwd = np.identity(2) - params[0]*params[1]*0.5*m
    mbwd = np.identity(2) + params[0]*params[1]*0.5*m

    # Using svd
    #u, s, v  = svd(mbwd)
    #xv_prod = mfwd @ xv
    #xv_prod = v.T @ ((u.T @ xv_prod) / s)
    # return xv_prod

    # Using built-in solver
    xv = np.linalg.solve(mbwd, mfwd @ xv)
    return xv 

# count number of steps needed to reach desired convergence
def itr_count(A, xv, xb, x2, x3, store, store_indx, method=euler_step, dif=1e-6):
    i = 0

    # iterate through method while boundary time not reached
    while (i < 1000):
        xv_old = np.array([[xv[0, 0]], [xv[1, 0]]])
        xv = method(xv, A)
        # want to shrink errors to zero --> "take your current errors and shrink them a little bit"
        x2 -= (xv_old[0, 0] - xv[0, 0])
        x3 -= (xv_old[1, 0] - xv[1, 0])
        if (x2 > 1e12 or x3 > 1e12):
            i = 1000

        xb += 0.001
        x2 += 0.001
        x3 += 0.001

        i+=1

        if (np.abs(xb - x2) < dif and np.abs(xb - x3) < dif):
            break
    
   
    # output results
    store[store_indx] = i
    


# Plotting for different values of alpha for a fixed dt
a_range = np.arange(0, 1400, 5)
euler_step_count = np.zeros_like(a_range)
rk2_step_count = np.zeros_like(a_range)
rk4_step_count = np.zeros_like(a_range)
trap_step_count = np.zeros_like(a_range)

for i in range(len(a_range)):
    params[1] = a_range[i]
    itr_count(m, xv, x1, x2, x3, euler_step_count, i, method=euler_step, dif=1e-5)
    itr_count(m, xv, x1, x2, x3, rk2_step_count, i, method=rk2_step, dif=1e-5)
    itr_count(m, xv, x1, x2, x3, rk4_step_count, i, method=rk4_step, dif=1e-5)
    itr_count(m, xv, x1, x2, x3, trap_step_count, i, method=trap_step, dif=1e-5)


# Graph
plt.rcParams['font.size'] = 24
plt.plot(a_range, euler_step_count)
plt.plot(a_range, rk2_step_count, color='mediumpurple')
plt.plot(a_range, rk4_step_count, color='maroon')
plt.plot(a_range, trap_step_count, color='darkgreen')
plt.xlabel("Alpha Value")
plt.ylabel("Number of Iterations")
plt.legend(["Euler", "RK2", "RK4", "Trap"])
plt.title("Steps to Converge vs Size of Alpha for Timestep h = 0.001")
plt.show()

# Optimal Parameters
print("RK2:", np.min(rk2_step_count), a_range[np.argmin(rk2_step_count)], 
      "EULER:", np.min(euler_step_count), a_range[np.argmin(euler_step_count)], 
      "RK4:", np.min(rk4_step_count), a_range[np.argmin(rk4_step_count)], 
      "TRAP:", np.min(trap_step_count), a_range[np.argmin(trap_step_count)])