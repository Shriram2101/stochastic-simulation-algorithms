import math
import numpy as np
import matplotlib.pyplot as plt
from itertools import zip_longest

count = 1

for i in range(count):
    X = [5000]
    c = 0.5

    t = [0]
    tend = 10
    while t[-1] <= tend:
        h = c * X[-1]

        if h == 0:
            break

        r1 = np.random.uniform()
        tau = ((1 / h) * math.log(1 / r1))
        t.append(t[-1] + tau)
        X.append(X[-1] - 1)

plt.plot(t, X)
plt.xlim(0, 10)
plt.ylim(0, 5000)
plt.xlabel("Time")
plt.ylabel("Number of X molecules")
plt.xticks(np.arange(0, 11, 1))
plt.yticks(np.arange(0, 5001, 500))
plt.show()





