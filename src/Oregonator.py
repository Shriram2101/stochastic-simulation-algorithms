import math
import numpy as np
import matplotlib.pyplot as plt

from itertools import zip_longest

Y1 = [500]
Y2 = [1000]
Y3 = [2000]

rho1 = [2000]
rho2 = [50000]

c1 = 2 #rho1[-1] / Y2[-1]
c2 = 0.1 #rho2[-1] / Y1[-1] * Y2[-1]
c3 = 104.0 #(rho1[-1] + rho2[-1]) / Y1[-1]
c4 = 0.016 #(2 * rho1[-1]) / Y1[-1] ** 2
c5 = 26 #(rho1[-1] + rho2[-1]) / Y3[-1]

t = [0]
tend = 3

while t[-1] <= tend:
    h = [c1 * Y2[-1], c2 * Y1[-1] * Y2[-1], c3 * Y1[-1], c4 * Y1[-1] * ((Y1[-1] - 1) / 2), c5 * Y3[-1]]
    h_sum = sum(h)

    if h_sum == 0:
        break

    r1 = np.random.uniform()
    tau = (1 / h_sum) * math.log(1 / r1)
    t.append(t[-1] + tau)

    r2 = np.random.uniform()

    if r2 * h_sum <= h[0]:
        Y1.append(Y1[-1] + 1)
        Y2.append(Y2[-1] - 1)
        Y3.append(Y3[-1])

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1]):
        Y1.append(Y1[-1] - 1)
        Y2.append(Y2[-1] - 1)
        Y3.append(Y3[-1])

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1] + h[2]):
        Y1.append(Y1[-1] + 1)
        Y2.append(Y2[-1])
        Y3.append(Y3[-1] + 1)

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1] + h[2] + h[3]):
        Y1.append(Y1[-1] - 2)
        Y2.append(Y2[-1])
        Y3.append(Y3[-1])

    elif (r2 * h_sum) > h[0] and r2 * h_sum <= (h[0] + h[1] + h[2] + h[3] + h[4]):
        Y1.append(Y1[-1])
        Y2.append(Y2[-1] + 1)
        Y3.append(Y3[-1] - 1)

with open('Oreganator.txt', 'w') as f:
    f.write('Time \t'.expandtabs(50))
    f.write('Concentration of Y1 molecules \t'.expandtabs(50))
    f.write('Concentration of Y2 molecules \t'.expandtabs(50))
    f.write('Concentration of Y3 molecules \t'.expandtabs(50))
    f.write('\n')
    for (A, B, C, D) in zip_longest(t, Y1, Y2, Y3):
        f.write('\n')
        f.write("{}\t{}\t{}\t{}".expandtabs(50).format(str(A), str(B), str(C), str(D)))
    f.close()

#plt.plot(t, Y1, label='Y1 Plot')
#plt.plot(t, Y2, label='Y2 Plot')
#plt.plot(t, Y3, label='Y3 Plot')
#plt.legend()
plt.plot(Y2, Y3)
plt.xlim(0, 7000)
plt.ylim(0, 10000)
plt.xticks(np.arange(0, 7001, 1000))
plt.yticks(np.arange(0, 10001, 1000))
plt.xlabel("Number of Y1 molecules")
plt.ylabel("Number of Y3 molecules")
plt.show()


















