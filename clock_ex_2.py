import numpy as np
import matplotlib.pyplot as plt
import random

"""
    Example 2: Stochastic Block Model
    Idea: Let us look at the affect that clustering of clocks in a network may have on
    time synchronization by generating a network that has 3 distinct clusters.

    For: 
    p - probability that two nodes in the same cluster are connected together
    q - probability that two nodes in different clusters are connected together = 1 - p
    m - nodes per cluster
    n - total nodes = 3 * n

    In this example, we will set n = 150 and vary p between [0.5, 0.95] at intervals of +/- 0.05. 
    After generating the Laplacian matrix and a series of time stamps with random noise added to them,
    we will then apply the iterative methods to the network and record the average iterations it takes for convergence
    to be reached across 10 trials. Similar to example 1, we will estimate an "optimial" minimum number of iterations through varying alpha
    during each trial across all four methods.

    Again, we will let h = 0.001

"""

"""
    Generating the network
"""

m = 50
n = 3*m

# assign probability for connection i --> j
def add(i, j, m, p, q):
    if (i <= m and j <= m):
        return (np.random.rand() < p)
    elif (i > m and j > m and i <= 2*m and j <= 2*m):
        return (np.random.rand() < p)
    elif (i > 2*m and j > 2*m):
        return (np.random.rand() < p)
    return (np.random.rand() < q)

# generate the network
def gennetwork(p, q, n, m):
    adj = np.zeros((n, n), dtype=float)
    x = []
    y = []
    for i in range(n):
        for j in range(n):
                if (not (i == j)) and add(i, j, m, p, q):
                    adj[i, j] = 1
                    x.append(i)
                    y.append(j)

    return adj, x, y


"""
    Numerical Methods
"""

# h, alpha
params = [0.001, 0.001]

def euler_step(xv, L, params):
    A = np.identity(n) - params[0]*params[1]*L
    xv = A @ xv
    return xv

def rk2_step(xv, L, params):
    t = params[0]*params[1]
    xv = xv - t * L @ (xv - t * 0.5 * L @ xv)
    return xv

def rk4_step(xv, L, params):
    t = -1*params[0] * params[1]
    k1 = L @ xv
    k2 = L @ (xv + 0.5*t*k1)
    k3 = L @ (xv + 0.5*t*k2)
    k4 = L @ (xv + t*k3)
    xv = xv + (t/6) * (k1 + 2*k2 + 2*k3 + k4)
    return xv

def trap_step(xv, L, params):
    mfwd = np.identity(n) - params[0]*params[1]*0.5*L
    mbwd = np.identity(n) + params[0]*params[1]*0.5*L
    xv = np.linalg.solve(mbwd, mfwd @ xv)
    return xv

"""
    Generate starting times
"""

x0 = np.full((n, 1), 40023.054, dtype=float)
# Range of Error
s = 0.2

for i in range(n):
    x0[i, 0] += random.uniform(-s, s)


"""
    Plot Time Differences between Nodes
"""

def timedif(xv, s, t):
    xv_dif = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            xv_dif[i, j] = np.abs(xv[i, 0] - xv[j, 0])
    
    title = "Time Difference Matrix at " + str(t) + " Iterations"
    plt.title(title)
    plt.imshow(xv_dif, vmin=0, vmax=s)

    # Grid showing different communities
    ax = plt.gca()
    ax.set_xticks(np.arange(-.5, 90, 30), minor=True)
    ax.set_yticks(np.arange(-.5, 90, 30), minor=True)
    ax.grid(which='minor', color='w', linestyle='-', linewidth=2)
    ax.tick_params(which='minor', bottom=False, left=False)

    plt.colorbar()
    plt.show()


"""
    Calculate iterations until convergence is reached
"""

# Threshold for convergence met
def converge(xv, dif):
    tot = 0
    xtest = random.choice(xv)
    for i in range(n):
        tot += np.abs(xv[i]-xtest)
    return (tot < dif)

# Threshold for divergence met
def diverge(xv, div):
    for i in range(n):
        if (xv[i] >= div):
            return True
    return False

def itr_count(L, xv, params, method=euler_step, dif=1e-5, div = 1e8):
    i = 0

    # iterate through method while boundary time not reached
    while (i < 1000):
        if (diverge(xv, div)):
            i = 1000
            break

        xv = method(xv, L, params)
        # Increase all clock times
        for j in range(n):
            xv[j, 0] += 0.001

        i+=1

        if (converge(xv, dif)):
            break
    
    return i


"""
    Calculate Results
"""

# Number of networks to test per p value
num_trials = 10

# Range of p values to test
p_val = np.arange(0.5, 1.0, 0.05)
i_val_euler = np.zeros_like(p_val)
i_val_rk2 = np.zeros_like(p_val)
i_val_rk4 = np.zeros_like(p_val)
i_val_trap = np.zeros_like(p_val)

# For all p values in the range
for j in range(len(p_val)):
    # Get p, q
    p = p_val[j]
    q = 1.0 - p

    # Initialize sums to 0
    i_sum_euler = 0
    i_sum_rk2 = 0
    i_sum_rk4 = 0
    i_sum_trap = 0

    # Loop for trials for that p value 
    for t in range(num_trials):
        # Generate random network and find its laplacian
        adj, x, y = gennetwork(p, q, n, m)
        L = np.diag(np.sum(adj, axis=1)) - adj

        # Initialize test iteration counts
        a_range = np.arange(15, 70, 1)
        euler_step_count = np.zeros_like(a_range)
        rk2_step_count = np.zeros_like(a_range)
        rk4_step_count = np.zeros_like(a_range)
        trap_step_count = np.zeros_like(a_range)

        # Test for min iterations
        for i in range(len(a_range)):
            params[1] = a_range[i]
            
            # Solver - Euler
            xv = np.copy(x0)
            euler_step_count[i] = itr_count(L, xv, params, method=euler_step)
            
            # Solver - RK2
            xv = np.copy(x0)
            rk2_step_count[i] = itr_count(L, xv, params, method=rk2_step)
            
            # Solver - RK4
            xv = np.copy(x0)
            rk4_step_count[i] = itr_count(L, xv, params, method=rk4_step)

            # Solver - Trapezoidal
            xv = np.copy(x0)
            trap_step_count[i] = itr_count(L, xv, params, method=trap_step)

        # increment trial values by min iterations
        i_sum_euler += np.min(euler_step_count)
        i_sum_rk2 += np.min(rk2_step_count)
        i_sum_rk4 += np.min(rk4_step_count)
        i_sum_trap += np.min(trap_step_count)

    # average across all trials and store results
    i_val_euler[j] = i_sum_euler / num_trials
    i_val_rk2[j] = i_sum_rk2 / num_trials
    i_val_rk4[j] = i_sum_rk4 / num_trials
    i_val_trap[j] = i_sum_trap / num_trials   


"""
    Display Results
"""

plt.title("Average Steps to Converge vs p")
plt.xlabel("p")
plt.ylabel("Number of Iterations")
plt.plot(p_val, i_val_euler, '.-')
plt.plot(p_val, i_val_rk2, '.-')
plt.plot(p_val, i_val_rk4, '.-')
plt.plot(p_val, i_val_trap, '.-')
plt.legend(["Euler", "RK2", "RK4", "Trap"])
plt.xticks(p_val)
plt.yticks(np.arange(0, 105, 10))
plt.show()
print(i_val_euler)
print(i_val_rk2)
print(i_val_rk4)
print(i_val_trap)