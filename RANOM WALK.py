import numpy as np
import matplotlib.pyplot as plt

end_step = 1000
position = [0]

for _ in range(end_step):
    r = np.random.uniform(0, 1)
    if r < 0.5:
        position.append(position[-1] - 1)
    else:
        position.append(position[-1] + 1)

plt.plot(position)
plt.xlabel('Steps')
plt.ylabel('Position')
plt.show()
