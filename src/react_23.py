import math
import numpy as np
import matplotlib.pyplot as plt

from itertools import zip_longest

Y = [12000]

c1 = 5
c2 = 0.00125

t = [0]
tend = 5

while t[-1] <= tend:
    h = [c1 * Y[-1], c2 * (Y[-1] * (Y[-1] - 1))]
    h_sum = sum(h)

    if h_sum == 0:
        break

    r1 = np.random.uniform()
    tau = (1 / h_sum) * math.log(1 / r1)
    t.append(t[-1] + tau)

    r2 = np.random.uniform()

    if r2 * h_sum <= h[0]:
        Y.append(Y[-1] + 2)

    elif (r2 * h_sum > h[0]) and (r2 * h_sum) <= (h[0] + h[1]):
        Y.append(Y[-1] - 2)

'''with open('react_23.txt', 'w') as f:
    f.write('Time \t'.expandtabs(50))
    f.write('Concentration')
    f.write('\n')
    for (a, b) in zip_longest(t, Y):
        f.write('\n')
        f.write("{}\t{}".expandtabs(50).format(str(a), str(b)))
    f.close()'''

plt.plot(t, Y)
plt.xlim(0, 5)
plt.ylim(0, 12000)
plt.xlabel("Time")
plt.ylabel("Number of molecules")
plt.show()
