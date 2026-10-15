import random
import math
import time
import numpy as np

clock_time = time.time()

lattice_size = 2000
bp = 147
k_add = 12
k_slide = 0.0024
k_fall = math.exp(-7)

initial_config = 1
final_config = 1000
nuc_count = []
occupancy = []

while initial_config <= final_config:
    lattice = ['0'] * lattice_size
    position = random.randint(0, lattice_size - bp)
    for i in range(bp):
        lattice[position + i] = '1'
    nuc_indices = [i + 1 for i, x in enumerate(lattice) if x == '1']
    vacant_indices = [i + 1 for i, x in enumerate(lattice) if x == '0']

    t = [0]
    end_time = 60

    while t[-1] <= end_time:

        nucleosome = len(nuc_indices) // bp

        nucleosome_separation = []

        nucleosome_start_position = None
        for indices in range(0, len(nuc_indices), bp):
            nucleosome_start_position = nuc_indices[indices]
            break
        gap_left = nucleosome_start_position - 1
        nucleosome_separation.append(gap_left)

        for i in range(0, len(nuc_indices), bp):
            current_group = nuc_indices[i:i + bp]
            next_group = nuc_indices[i + bp:i + 2 * bp]
            if current_group and next_group:
                separation = (next_group[0] - current_group[-1]) - 1
                nucleosome_separation.append(separation)

        for indices in range(0, len(nuc_indices), bp):
            nucleosome_start_position = nuc_indices[indices]
        gap_right = (lattice_size - nucleosome_start_position) - (bp - 1)
        nucleosome_separation.append(gap_right)

        left_space = nucleosome_separation[:-1]
        right_space = nucleosome_separation[1:]

        r1 = np.random.uniform()
        rates = np.zeros(4 * nucleosome)

        for element in range(nucleosome):
            rates[element * 4] = k_add
            rates[(element * 4) + 1] = k_slide if left_space[element] >= 10 else 0
            rates[(element * 4) + 2] = k_slide if right_space[element] >= 10 else 0
            rates[(element * 4) + 3] = k_fall

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
            nucleosome_index = event_choice // 4
            direction = event_choice % 4

            if direction == 0:
                vacant_indices = [i for i, x in enumerate(lattice) if x == '0']

                vacant_segments = []
                current_segment = []

                for index in vacant_indices:
                    if not current_segment or index == current_segment[-1] + 1:
                        current_segment.append(index)
                    else:
                        if len(current_segment) >= bp:
                            vacant_segments.append(current_segment)
                        current_segment = [index]

                if len(current_segment) >= bp:
                    vacant_segments.append(current_segment)

                if vacant_segments:
                    selected_segment = random.choice(vacant_segments)
                    start_position = random.choice(selected_segment[:len(selected_segment) - bp + 1])

                    for i in range(bp):
                        lattice[start_position + i] = '1'

                    nuc_indices.extend([start_position + i + 1 for i in range(bp)])
                    nuc_indices.sort()

            elif direction == 1 and left_space[nucleosome_index] >= 10:
                tri_positions = nuc_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp]
                new_positions = [pos - 10 for pos in tri_positions]
                for i in range(bp):
                    lattice[new_positions[i] - 1] = '1'
                    lattice[tri_positions[i] - 1] = '0'

                nuc_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp] = new_positions

            elif direction == 2 and right_space[nucleosome_index] >= 10:
                tri_positions = nuc_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp]
                new_positions = [pos + 10 for pos in tri_positions]
                for i in range(bp):
                    lattice[new_positions[i] - 1] = '1'
                for j in range(0, 10):
                    lattice[tri_positions[j] - 1] = '0'

                nuc_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp] = new_positions

            elif direction == 3:
                tri_positions = nuc_indices[nucleosome_index * bp: (nucleosome_index + 1) * bp]
                for i in range(bp):
                    lattice[tri_positions[i] - 1] = '0'
                nucleosome -= 1
                tri_indices = nuc_indices[:nucleosome_index * bp] + nuc_indices[(nucleosome_index + 1) * bp:]

        t.append(t[-1] + tau)

    occupancy.append(list(lattice))

    nuc_list = int(len(nuc_indices) / 147)
    nuc_count.append(nuc_list)

    print(initial_config, nuc_list)
    initial_config += 1

nuc_avg = sum(nuc_count) / final_config

transposed_list = list(zip(*occupancy))

with open('Dynamics_60.txt', 'w') as z:
    for i, column in enumerate(transposed_list):
        count_list = column.count('1')
        probability = count_list / final_config
        z.write(f"{i + 1}\t {count_list}\t {probability}\n")

end_clock_time = time.time()
elapsed_clock_time = end_clock_time - clock_time
print("ELAPSED CLOCK TIME:\t", elapsed_clock_time)
print("AVERAGE NUCLEOSOME:", nuc_avg)
