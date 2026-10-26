import numpy as np
import matplotlib.pyplot as plt

end_time = 1000
x_position = [0]
y_position = [0]

for _ in range(end_time):
    rx = np.random.uniform(0, 1)
    ry = np.random.uniform(0, 1)

    if rx < 0.5:
        x_position.append(x_position[-1] + 1)
    else:
        x_position.append(x_position[-1] - 1)

    if ry < 0.5:
        y_position.append(y_position[-1] + 1)
    else:
        y_position.append(y_position[-1] - 1)
# plt.figure(facecolor="#3aafa9")
plt.plot(x_position, y_position, color="#17252a")

'''ax = plt.gca()
ax.set_facecolor("#3aafa9")
ax.tick_params(axis="x", colors="#17252a")
ax.tick_params(axis="y", colors="#17252a")
ax.xaxis.label.set_color("#17252a")
ax.yaxis.label.set_color("#17252a")
ax.spines['left'].set_color("#17252a")
ax.spines['right'].set_color("#3aafa9")
ax.spines['top'].set_color("#3aafa9")
ax.spines['bottom'].set_color("#17252a")'''
plt.xlabel('Time')
plt.ylabel('Position')
plt.show()
