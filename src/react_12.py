# EXPONENTIAL DECAY STOCHASTIC SIMULATION ALGORITHM FOR SINGLE RUN

# IMPORTING NECESSARY LIBRARIES
import math
import numpy as np
import matplotlib.pyplot as plt

# DECLARING MOLECULAR CONCENTRATION AND RATE CONSTANT
X = [10000]
c = 0.5

# INITIALIZING START & END TIME
t = [0]
tend = 10

# IMPLEMENTING GILLESPIE ALGORITHM
while t[-1] <= tend:
    h = c * X[-1]

    if h == 0:
        break

    # CHOOSING THE RANDOM NUMBER
    r1 = np.random.uniform()
    tau = (1 / h) * math.log(1 / r1)

    t.append(t[-1] + tau)              # UPDATE TIME VALUE
    X.append(X[-1] - 1)                # UPDATE MOLECULAR CONCENTRATION

# PLOTTING GRAPH USING MATPLOTLIB
plt.plot(t, X)
plt.xlim(0, 10)
plt.ylim(0, 10000)
plt.xlabel("Time")
plt.ylabel("Number of molecules")
plt.show()
