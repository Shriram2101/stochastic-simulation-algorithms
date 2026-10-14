import math
import numpy as np
import matplotlib.pyplot as plt

from itertools import zip_longest

Y1 = [1000]
Y2 = [2000]

c1 = 5000
c2 = 50
c3 = 0.00005
c4 = 5

t = [0]
tend = 14

while t[-1] <= tend:
    h = [c1, c2 * Y1[-1], c3 * Y2[-1] * Y1[-1] * ((Y1[-1] - 1) / 2), c4 * Y1[-1]]
    h_sum = sum(h)

    if h_sum == 0:
        break

    r1 = np.random.uniform()
    tau = (1 / h_sum) * math.log(1 / r1)
    t.append(t[-1] + tau)

    r2 = np.random.uniform()

    if r2 * h_sum <= h[0]:
        Y1.append(Y1[-1] + 1)
        Y2.append(Y2[-1])

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1]):
        Y1.append(Y1[-1] - 1)
        Y2.append(Y2[-1] + 1)

    elif (r2 * h_sum) > (h[0]) and r2 * h_sum <= (h[0] + h[1] + h[2]):
        Y1.append(Y1[-1] + 1)
        Y2.append(Y2[-1] - 1)

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1] + h[2] + h[3]):
        Y1.append(Y1[-1] - 1)
        Y2.append(Y2[-1])

plt.xlim(0, 8000)
plt.ylim(0, 8000)
#plt.plot(t, Y1, label='Y1 Plot')
#plt.plot(t, Y2, label='Y2 Plot')
#plt.legend()
plt.plot(Y1, Y2)
plt.xticks(np.arange(0, 8001, 1000))
plt.yticks(np.arange(0, 8001, 1000))
plt.xlabel("Number of Y1 molecules")
plt.ylabel("Number of Y2 molecules")
plt.show()

with open('Brusselator.txt', 'w') as f:
    f.write('Time \t'.expandtabs(50))
    f.write('Concentration of Y1 molecules \t'.expandtabs(50))
    f.write('Concentration of Y2 molecules \t'.expandtabs(50))
    f.write('\n')
    for (A, B, C) in zip_longest(t, Y1, Y2):
        f.write('\n')
        f.write("{}\t{}\t{}".expandtabs(50).format(str(A), str(B), str(C)))
    f.close()