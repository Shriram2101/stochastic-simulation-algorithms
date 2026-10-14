# LOTKA MODEL EXACT STOCHASTIC SIMULATION ALGORITHM FOR SINGLE RUN

# Importing libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from itertools import zip_longest

# Declaring variables
Y1 = [1000]    # Number of Y1 molecules
Y2 = [1000]    # Number of Y2 molecules

# Declaring rate constants
c1 = c3 = 10
c2 = 0.01

# Initializing start time & end time
t = [0]   # INITIALIZING START TIME AS A LIST
tend = 10

# IMPLEMENTING GILLESPIE ALGORITHM
while t[-1] <= tend:
    h = [c1 * Y1[-1], c2 * Y1[-1] * Y2[-1], c3 * Y2[-1]]
    h_sum = sum(h)

    if h_sum == 0:
        break

    # Choose first random number
    r1 = np.random.uniform()
    tau = (1 / h_sum) * math.log(1 / r1)
    t.append(t[-1] + tau)                           # Update time

    # Choose second random number
    r2 = np.random.uniform()

    # Selection of reactions
    if r2 * h_sum <= h[0]:
        Y1.append(Y1[-1] + 1)
        Y2.append(Y2[-1])

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1]):
        Y1.append(Y1[-1] - 1)
        Y2.append(Y2[-1] + 1)

    elif (r2 * h_sum) > (h[0]) and r2 * h_sum <= (h[0] + h[1] + h[2]):
        Y1.append(Y1[-1])
        Y2.append(Y2[-1] - 1)

# Plotting the graph using matplotlib
plt.xlim(0, 10)
plt.ylim(0, 3000)
plt.plot(t, Y1, label='Prey')
plt.plot(t, Y2, label='Predator')
plt.plot(Y1, Y2)
plt.xticks(np.arange(0, 11, 1))
plt.yticks(np.arange(0, 3001, 500))
plt.legend()
plt.xlabel("Time")
plt.ylabel("Number of Preys and Predator")
plt.show()


