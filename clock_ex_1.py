import numpy as np
from scipy.linalg import null_space, svd
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

# Graph Laplacian
m = np.array([[2, -1], [-1,1]])
eval = np.linalg.eigvals(m)
#print("Eigenvalues of L(G): ", eval)

# For Sheaf Laplacian, alpha = 2 / h(lambda max + lambda min) acc to problem 6 in paper
params[1] = 2 / (params[0] * (np.max(eval) + np.min(eval)))
#params[1] = 67
#print("Ideal Alpha Value:", params[1])
# A = 1 - h * alpha * m
A = np.identity(2) - params[0]*params[1]*m
#print(A)
#print("Eigenvalues of A: ", np.linalg.eig(A)[0])
#print("Nullspace of M", null_space(np.array([[0, 0, 0], [-1, 2, -1], [0, -1, 1]])).flatten())

# Times (in seconds)
x1 = 40023.054
x2 = 40023.045
x3 = 40023.067
clocks = [x1, x2, x3]
xv = np.array([[x2-x1], [x3-x1]])
print(f"Initial times: clock A: {x1:.3f} clock B: {x2:.3f} clock C: {x3:.3f}")
#print(xv[0, 0], xv[1, 0], x2, x3)

# one euler step
def euler_step(xv, m):
    #A = np.identity(2) - params[0]*params[1]*m
    xv = A @ xv
    # check time difference
    #print("err value:", xv[0, 0], xv[1, 0], end=" ")
    return xv

"""
    Alternative methods
"""
def rk2_step(xv, m):
    m1 = np.identity(2) - (params[1]*params[0]*0.5*m)
    #print ("M1 evals", np.linalg.eigvals(m1))
    xk = m1 @ xv
    m2 = np.identity(2) - (params[1]*params[0]*m)
    #print ("M2 evals", np.linalg.eigvals(m2))
    xv = m2 @ xk
    return xv

def rk4_step(xv, m):
    k1 = -1*params[0]*params[1]*m @ xv
    k2 = -1*params[0]*params[1]*0.5*m @ (xv + k1)
    k3 = -1*params[0]*params[1]*0.5*m @ (xv + k2)
    k4 = -1*params[0]*params[1]*m @ (xv + k3)
    xv = xv + (1/6) * (k1 + 2*k2 + 2*k3 + k4)
    return xv

def trap_step(xv, m):
    mfwd = np.identity(2) - params[0]*params[1]*0.5*m
    mbwd = np.identity(2) + params[0]*params[1]*0.5*m

    # TEST -- Using svd
    u, s, v  = svd(mbwd)
    xv_prod = mfwd @ xv
    #for i in range(len(v)):
    #    if (s[i] > 1e-6):
    #        v[i, i] *= 1/s[i]

    xv_prod = v.T @ ((u.T @ xv_prod) / s)

    # was there supposed to be an issue with this or
    xv = np.linalg.solve(mbwd, mfwd @ xv)

    # I guess np.linalg.solve is more accurate than using svd
    #if ((np.abs(xv_prod - xv) > np.full_like(xv, 1e-6)).any()):
    #    print("false", xv_prod, xv)
    return xv_prod

# count number of steps needed to reach desired convergence
def itr_count(A, xv, xb, x2, x3, store, store_indx, method=euler_step, dif=1e-6):
    i = 0
    # For cumulative difference plotting
    #temp = [0]*23
    #temp2 = [0]*7
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
        #print(x2, x3)
        #if (i < 23):
            #temp[i] = np.abs(xv_old[0, 0] - xv[0, 0])
        #    temp2[i] =  np.abs(x3 - x1)
            #temp[i] = np.abs(xv[0, 0] - xv_old[0, 0])
            #temp2[i] = np.abs(xv[1, 0] - xv_old[1, 0])

        i+=1

        if (np.abs(xb - x2) < dif and np.abs(xb - x3) < dif):
            break

    # output results
    #store[store_indx] = i
    print(f"Clocks converged after {i} iterations or by time: {x2 : .3f}")
    
    #plt.rcParams['font.size'] = 16
    #plt.plot(np.arange(0, 23, 1), temp)
    #plt.plot(np.arange(0, 27, 1), temp2, color='maroon')
    #plt.legend(["Clock B's Time", "Clock C's Time"])
    #plt.xlabel("Number of Iterations")
    #plt.ylabel("Clock Error")
    #plt.title("Clock B's Change in Error vs. Iteration Count for RK Methods")
    #plt.show()
    

"""
# Plotting for different values of alpha for a fixed dt
# Though it may be more valulable to do the inverse? Or just the constants in general
a_range = np.arange(850, 980, 5)
euler_step_count = np.zeros_like(a_range)
rk2_step_count = np.zeros_like(a_range)
rk4_step_count = np.zeros_like(a_range)
trap_step_count = np.zeros_like(a_range)

for i in range(len(a_range)):
    params[1] = a_range[i]
    #itr_count(m, xv, x1, x2, x3, euler_step_count, i, method=euler_step, dif=1e-5)
    itr_count(m, xv, x1, x2, x3, rk2_step_count, i, method=rk2_step, dif=1e-5)
    itr_count(m, xv, x1, x2, x3, rk4_step_count, i, method=rk4_step, dif=1e-5)
    #itr_count(m, xv, x1, x2, x3, trap_step_count, i, method=trap_step, dif=1e-5)

#print(np.min())

plt.rcParams['font.size'] = 16
plt.plot(a_range, euler_step_count)
plt.plot(a_range, rk2_step_count, color='mediumpurple')
plt.plot(a_range, rk4_step_count, color='maroon')
plt.plot(a_range, trap_step_count, color='darkgreen')
plt.xlabel("Alpha Value")
plt.ylabel("Number of Iterations")
plt.legend(["Euler", "RK2", "RK4", "Trap"])
plt.title("Steps to Converge vs Size of Alpha for Timestep h = 0.001")
plt.show()
"""


print("Solving the heat equation (and incrementing our clock values by 0.001s) until all \ntimes are within 1e-4s of each other...")
# Different choice for alpha needed to min iterations
params[1] = 66
print(f"... using {'Euler Method:'}", end=" ")
itr_count(m, xv, x1, x2, x3, None, None, method=euler_step, dif=1e-4)
params[1] = 980
print(f"... using {'Runge-Kutta 2:'}", end=" ")
itr_count(m, xv, x1, x2, x3, None, None, method=rk2_step, dif=1e-4)
params[1] = 820
print(f"... using {'Runge-Kutta 4:'}", end=" ")
itr_count(m, xv, x1, x2, x3, None, None, method=rk4_step, dif=1e-4)
params[1] = 2150
print(f"... using {'Trapezoidal Method:'}", end=" ")
itr_count(m, xv, x1, x2, x3, None, None, method=trap_step, dif=1e-4)

#plt.rcParams['font.size'] = 16
#plt.plot(np.arange(0, 23, 1), strk2)
#plt.plot(np.arange(0, 23, 1), strk4, color='maroon')
#plt.legend(["RK 2", "RK 4"])
#plt.xlabel("Number of Iterations")
#plt.ylabel("Clock Error")
#plt.title("Convergence of Runge-Kutta Methods for Example 1")
#plt.show()




#### OLD
"""
# WORKING CODE FOR FULL TIMES (As in it works)

# h = dt (ish), alpha between (0, 1/h*d) s.t. d is largest diag value
# here, d = 2 and let h = 1/100 --> alpha in (0, 50)
params = [.01, 0]

# Graph Laplacian
m = np.array([[2, -1], [-1,1]])
eval = np.linalg.eigvals(m)
print("Eigenvalues of L(G): ", eval)

# For Sheaf Laplacian, alpha = 2 / h(lambda max + lambda min) acc to problem 6 in paper
params[1] = 2 / (params[0] * (np.sum(eval)))
print("Ideal Alpha Value:", params[1])
# A = 1 - h * alpha * m
A = [[1, 0], [0, 1]] - params[0]*params[1]*m
print("Eigenvalues of A: ", np.linalg.eig(A)[0])

# Times (in seconds)
x1 = 40023.054
x2 = 40023.043
x3 = 40023.067
xb = np.array([[1*params[0]*params[1]], [0]]) * x1
xv = np.array([[x2], [x3]])
print("Input times: 1:", x1, "2:", x2, "3:", x3)


# one euler step
def euler_step(xb, xv, A):
    xv = A @ xv + xb 
    # check time difference
    print("dif:", np.abs(xv[0] - x1), np.abs(xv[1] - x1))
    return xv


# test if boundary condition satsified
def conv(dif, xb, xv):
    for i in range(len(xv)):
        if (np.abs(xv[i] - xb) > dif):
            return False
    
    return True

# count number of steps needed to reach desired convergence
def itr_count(A, B, xv, xb, method=euler_step, dif=1e-6):
    i = 0
    # iterate through method while boundary time not reached
    while (i < 100 or not conv(dif, xb, xv)):
        xv = method(B, xv, A)
        i+=1

        if (conv(dif, xb, xv)):
            break

    # output results
    print("i:", i, end=" ")
    for j in range(len(xv)):
        print(j+2, ":", xv[j], end=" ")
    print('')

itr_count(A, xb, xv, x1, method=euler_step, dif=1e-3)
"""