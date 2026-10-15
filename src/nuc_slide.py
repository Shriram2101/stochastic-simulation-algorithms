import random
import math
import time
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

clock_time = time.time()

lattice_size = 5000
nucleosome = 30
bp = 147
k_slide = 1.0

initial_config = 1
final_config = 2000
occupancy = []
initial = []
initial_separation = []
final_separation = []

while initial_config <= final_config:
    if nucleosome * bp > lattice_size:
        print("Not enough lattice size")
    else:
        lattice = ['0'] * lattice_size
        segment_size = lattice_size // nucleosome
        for i in range(nucleosome):
            position = random.randint(i * segment_size, (i + 1) * segment_size - bp)
            for j in range(bp):
                lattice[position + j] = '1'
        tri_indices = [i + 1 for i, x in enumerate(lattice) if x == '1']
        initial.append(list(lattice))

        nucleosome_initial_sep = []

        nucleosome_initial_start_position = None
        for indices in range(0, len(tri_indices), bp):
            nucleosome_initial_start_position = tri_indices[indices]
            break
        gap_initial_left = nucleosome_initial_start_position - 1
        nucleosome_initial_sep.append(gap_initial_left)

        for i in range(0, len(tri_indices), bp):
            current_group = tri_indices[i:i + bp]
            next_group = tri_indices[i + bp:i + 2 * bp]
            if current_group and next_group:
                separation = (next_group[0] - current_group[-1]) - 1
                nucleosome_initial_sep.append(separation)

        for indices in range(0, len(tri_indices), bp):
            nucleosome_initial_start_position = tri_indices[indices]
        gap_right = (lattice_size - nucleosome_initial_start_position) - (bp - 1)
        nucleosome_initial_sep.append(gap_right)

        initial_separation.append(list(nucleosome_initial_sep))

        transpose_initial_list = list(zip(*initial))

        with open('den88_sep10_initial_dist_0sep.txt', 'w') as z:
            for i, column in enumerate(transpose_initial_list):
                count_list = column.count('1')
                probability = count_list / final_config
                z.write(f"{i + 1}\t {count_list}\t {probability}\n")

        t = [0]
        end_time = 1000
        new_t = [0]
        time_step = 20

        while t[-1] <= end_time:
            nucleosome_separation = []

            nucleosome_start_position = None
            for indices in range(0, len(tri_indices), bp):
                nucleosome_start_position = tri_indices[indices]
                break
            gap_left = nucleosome_start_position - 1
            nucleosome_separation.append(gap_left)

            for i in range(0, len(tri_indices), bp):
                current_group = tri_indices[i:i + bp]
                next_group = tri_indices[i + bp:i + 2 * bp]
                if current_group and next_group:
                    separation = (next_group[0] - current_group[-1]) - 1
                    nucleosome_separation.append(separation)

            for indices in range(0, len(tri_indices), bp):
                nucleosome_start_position = tri_indices[indices]
            gap_right = (lattice_size - nucleosome_start_position) - (bp - 1)
            nucleosome_separation.append(gap_right)

            left_space = nucleosome_separation[:-1]
            right_space = nucleosome_separation[1:]

            r1 = np.random.uniform()
            rates = np.zeros(2 * nucleosome)

            for element in range(nucleosome):
                rates[element * 2] = k_slide if left_space[element] >= 11 else 0
                rates[(element * 2) + 1] = k_slide if right_space[element] >= 11 else 0

            total_rate = np.sum(rates)

            cumulative_rate = []
            rate_sum = 0
            for rate_choice in rates:
                rate_sum += rate_choice
                cumulative_rate.append(rate_sum)

            r2 = np.random.uniform()
            tau = (1 / rate_sum) * math.log(1 / r2)

            event_choice = None
            for i, loop_rate in enumerate(cumulative_rate):
                if r1 > 0 and r1 * total_rate <= loop_rate:
                    event_choice = i
                    break

            if event_choice is not None:
                nucleosome_index = event_choice // 2
                direction = event_choice % 2

                if direction == 0 and left_space[nucleosome_index] >= 10:
                    tri_positions = tri_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp]
                    new_positions = [pos - 10 for pos in tri_positions]
                    for i in range(bp):
                        lattice[new_positions[i] - 1] = '1'
                        lattice[tri_positions[i] - 1] = '0'

                    tri_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp] = new_positions

                elif direction == 1 and right_space[nucleosome_index] >= 10:
                    tri_positions = tri_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp]
                    new_positions = [pos + 10 for pos in tri_positions]
                    for i in range(bp):
                        lattice[new_positions[i] - 1] = '1'
                    for j in range(0, 10):
                        lattice[tri_positions[j] - 1] = '0'

                    tri_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp] = new_positions

            t.append(t[-1] + tau)

            if t[-1] > new_t[-1]:
                t_med = t[-1] - new_t[-1]
                t_med_p = t_med / time_step
                itr = int(t_med_p)
                for j in range(itr + 1):
                    if new_t[-1] != end_time:
                        new_t.append(new_t[-1] + time_step)
                    else:
                        break

        occupancy.append(list(lattice))

        nucleosome_final_sep = []

        nucleosome_final_start_position = None
        for indices in range(0, len(tri_indices), bp):
            nucleosome_final_start_position = tri_indices[indices]
            break
        gap_final_left = nucleosome_final_start_position - 1
        nucleosome_final_sep.append(gap_final_left)

        for i in range(0, len(tri_indices), bp):
            current_group = tri_indices[i:i + bp]
            next_group = tri_indices[i + bp:i + 2 * bp]
            if current_group and next_group:
                separation = (next_group[0] - current_group[-1]) - 1
                nucleosome_final_sep.append(separation)

        for indices in range(0, len(tri_indices), bp):
            nucleosome_final_start_position = tri_indices[indices]
        gap_right = (lattice_size - nucleosome_final_start_position) - (bp - 1)
        nucleosome_final_sep.append(gap_right)

        final_separation.append(list(nucleosome_final_sep))

        print(initial_config, nucleosome_final_sep)
        initial_config += 1

    transposed_list = list(zip(*occupancy))

    with open('den88_sep10_final_dist_0sep.txt', 'w') as z:
        for i, column in enumerate(transposed_list):
            count_list = column.count('1')
            probability = count_list / final_config
            z.write(f"{i + 1}\t {count_list}\t {probability}\n")

end_clock_time = time.time()
elapsed_clock_time = end_clock_time - clock_time
print("ELAPSED CLOCK TIME:\t", elapsed_clock_time)

all_initial_separation = np.concatenate(initial_separation)
all_final_separation = np.concatenate(final_separation)

log_initial_separation = np.log10(all_initial_separation)
log_final_separation = np.log10(all_final_separation)

az = sns.kdeplot(all_initial_separation)
plt.title('Distribution of Nucleosome Separations Initially')
plt.xlabel('Separation (bp)')
plt.ylabel('Frequency')
plt.show()

ay = sns.kdeplot(all_final_separation)
plt.title('Distribution of Nucleosome Separations Finally')
plt.xlabel('Separation (bp)')
plt.ylabel('Frequency')
plt.show()

aw = sns.kdeplot(log_initial_separation)
plt.title('Distribution of Nucleosome Separations Initially')
plt.xlabel('Separation (bp)')
plt.ylabel('Frequency')
plt.show()

ax = sns.kdeplot(log_final_separation)
plt.title('Distribution of Nucleosome Separations Finally')
plt.xlabel('Separation (bp)')
plt.ylabel('Frequency')
plt.show()
