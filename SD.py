"""STOCHASTIC GENE EXPRESSION CODE FOR MODEL INVOLVING TF COMPLEXES, mRNA, PROTEINS,
miRNA, mRNA-miRNA COMPLEXES. THE CODE IMPLEMENTS EXACT GILLESPIE ALGORITHM, AND CALCULATES,
MEAN, VARIANCE AND STANDARD DEVIATION AND ALSO GENERATES THE FILES CORRESPONDINGLY FOR 5,000 ITERATIONS."""

# IMPORTING ESSENTIAL LIBRARIES IN PYTHON
import math
import random
import numpy as np
import time
from itertools import zip_longest

# SET THE START TIME FOR CLOCK AND CPU
clock_time = time.time()                             # CLOCK TIME
process_time = time.process_time()                   # CPU TIME

# INITIALIZING THE MINIMUM AND MAXIMUM SIMULATION ITERATIONS
loop = 1
count = 5000

# INITIALIZING THE END POINT FOR THE SIMULATION
tend = 20000

"""INITIALIZING LISTS TO STORE THE CUMULATIVE SUM OF RESULTS OF TF COMPLEXES, mRNA, PROTEINS,
miRNA, mRNA-miRNA COMPLEXES"""
cumulative_t = [0] * tend
cumulative_n = [0] * tend
cumulative_m = [0] * tend
cumulative_p = [0] * tend
cumulative_mi = [0] * tend
cumulative_mi_com = [0] * tend

"""INITIALIZING LISTS TO STORE THE CUMULATIVE SQUARED RESULTS OF TF-COMPLEXES, mRNA, PROTEINS,
miRNA, mRNA-miRNA COMPLEXES"""
cumulative_squared_n = [0] * tend
cumulative_squared_m = [0] * tend
cumulative_squared_p = [0] * tend
cumulative_squared_mi = [0] * tend
cumulative_squared_mi_com = [0] * tend

# IMPLEMENTING GILLESPIE ALGORITHM
while loop <= count:

    # DECLARING VARIABLES
    k_on = 0.0027               # TF BINDING RATE
    k_off = 0.0023              # TF UNBINDING RATE

    # THE INITIAL POPULATION FOR EACH OF THESE IS ASSUMED TO BE AT ZERO
    n = [0]                     # INITIAL NUMBER OF TF-PROMOTER COMPLEX
    m = [0]                     # INITIAL NUMBER OF mRNA
    p = [0]                     # INITIAL NUMBER OF PROTEIN MOLECULES
    mi = [0]                    # INITIAL NUMBER OF miRNA
    mi_com = [0]                # INITIAL NUMBER OF mRNA AND miRNA complexes

    N_tf = [10]                 # NUMBER OF TF
    N_p = [1]                   # NUMBER OF PROMOTERS. VARIED AS (1,2,3,4)

    r_bound = 0.33              # TF BOUND TRANSCRIPTION RATE
    r_free = 0.00               # TF FREE TRANSCRIPTION RATE
    gamma = 0.011               # mRNA DECAY RATE
    p_form = 0.5                # PROTEIN FORMATION RATE
    p_degrade = 0.01            # PROTEIN DEGRADATION RATE

    gamma_mi = 0.0004           # mRNA-miRNA COMPLEX DECAY RATE
    gs = 0.0003                 # miRNA DECAY RATE
    ks = 0.05                   # miRNA PRODUCTION FROM CERTAIN SOURCE
    k_plus = 0.005              # mRNA-miRNA ASSOCIATION RATE
    alpha = 0.5                 # CATALYTICITY PARAMETER

    t = [0]                     # THE STARTING TIME FOR THE SIMULATION
    time_step = 100              # INCREMENTING THE TIME STEP BY 10

    new_t = [0]                 # ARRAY TO STORE THE VALUE OF TIME_STEP
    new_n = [0]                 # ARRAY TO STORE THE VALUE OF TF COMPLEX CORRESPONDING TO TIME_STEP
    new_m = [0]                 # ARRAY TO STORE THE VALUE OF mRNA CORRESPONDING TO TIME_STEP
    new_p = [0]                 # ARRAY TO STORE THE VALUE OF PROTEINS CORRESPONDING TO TIME_STEP

    new_mi = [0]                # ARRAY TO STORE THE VALUE OF miRNA CORRESPONDING TO TIME_STEP
    new_mi_com = [0]            # ARRAY TO STORE THE VALUE OF mRNA-miRNA COMPLEXES CORRESPONDING TO TIME_STEP

    probability_miRNA = [0]

    random.seed(a=loop, version=2)  # INITIALIZING A SEED VALUE FOR THE PURPOSE OF GENERATION OF RANDOM VALUES

    # START OF GILLESPIE ALGORITHM BY DECLARING PROPENSITIES FOR EACH REACTIONS
    while t[-1] <= tend:
        h = [k_on * ((N_tf[-1] - n[-1]) * (N_p[-1] - n[-1])),
             k_off * n[-1],
             r_bound * n[-1],
             r_free * (N_p[-1] - n[-1]),
             gamma * m[-1],
             p_form * m[-1],
             p_degrade * p[-1],
             ks,
             gs * mi[-1],
             k_plus * m[-1] * mi[-1],
             alpha * gamma_mi * mi_com[-1],
             (alpha * (1 - gamma_mi)) * mi_com[-1]]
        h_sum = sum(h)

        # CHOOSING THE FIRST RANDOM NUMBER
        r1 = np.random.uniform()
        tau = (1 / h_sum) * math.log(1 / r1)
        t.append(t[-1] + tau)                       # UPDATE TIME

        # CHOOSING SECOND RANDOM NUMBER
        r2 = np.random.uniform()

        # SELECTION OF REACTIONS
        if r2 * h_sum <= h[0]:
            n.append(n[-1] + 1)
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1]:
            n.append(n[-1] - 1)
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2]:
            n.append(n[-1])
            m.append(m[-1] + 1)
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3]:
            n.append(n[-1])
            m.append(m[-1] + 1)
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4]:
            n.append(n[-1])
            m.append(m[-1] - 1)
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5]:
            n.append(n[-1])
            m.append(m[-1] - 1)
            p.append(p[-1] + 1)
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6]:
            n.append(n[-1])
            m.append(m[-1])
            p.append(p[-1] - 1)
            mi.append(mi[-1])
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7]:
            n.append(n[-1])
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1] + 1)
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + h[8]:
            n.append(n[-1])
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1] - 1)
            mi_com.append(mi_com[-1])

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + h[8] + h[9]:
            n.append(n[-1])
            m.append(m[-1] - 1)
            p.append(p[-1])
            mi.append(mi[-1] - 1)
            mi_com.append(mi_com[-1] + 1)

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + h[8] + h[9] + h[10]:
            n.append(n[-1])
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1])
            mi_com.append(mi_com[-1] - 1)

        elif r2 * h_sum <= h[0] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + h[8] + h[9] + h[10] + h[11]:
            n.append(n[-1])
            m.append(m[-1])
            p.append(p[-1])
            mi.append(mi[-1] + 1)
            mi_com.append(mi_com[-1] - 1)

        # ALGORITHM TO IMPLEMENT TIME STEP OF 10s
        if t[-1] > new_t[-1]:
            t_med = t[-1] - new_t[-1]
            t_med_p = t_med / time_step
            itr = int(t_med_p)
            for j in range(itr+1):
                if new_t[-1] != tend:
                    new_t.append(new_t[-1] + time_step)
                    new_n.append(n[-1])
                    new_m.append(m[-1])
                    new_p.append(p[-1])
                    new_mi.append(mi[-1])
                    new_mi_com.append(mi_com[-1])
                else:
                    break

    # AFTER EACH ITERATION ADD RESULTS TO THE CUMULATIVE SUMS FOR CALCULATING MEAN
    cumulative_t = [x + y for x, y in zip(cumulative_t, new_t)]
    cumulative_n = [x + y for x, y in zip(cumulative_n, new_n)]
    cumulative_m = [x + y for x, y in zip(cumulative_m, new_m)]
    cumulative_p = [x + y for x, y in zip(cumulative_p, new_p)]
    cumulative_mi = [x + y for x, y in zip(cumulative_mi, new_mi)]
    cumulative_mi_com = [x + y for x, y in zip(cumulative_mi_com, new_mi_com)]

    # CALCULATION OF MEAN
    mean_t = [x / count for x in cumulative_t]
    mean_n = [x / count for x in cumulative_n]
    mean_m = [x / count for x in cumulative_m]
    mean_p = [x / count for x in cumulative_p]
    mean_mi = [x / count for x in cumulative_mi]
    mean_mi_com = [x / count for x in cumulative_mi_com]

    # AFTER EACH ITERATION ADD RESULTS TO THE CUMULATIVE SQUARED FOR CALCULATION OF STANDARD DEVIATION
    cumulative_squared_n = [pow((x - y), 2) for x, y in zip(new_n, mean_n)]
    cumulative_squared_m = [pow((x - y), 2) for x, y in zip(new_m, mean_m)]
    cumulative_squared_p = [pow((x - y), 2) for x, y in zip(new_p, mean_p)]
    cumulative_squared_mi = [pow((x - y), 2) for x, y in zip(new_mi, mean_mi)]
    cumulative_squared_mi_com = [pow((x - y), 2) for x, y in zip(new_mi_com, mean_mi_com)]

    # CALCULATION OF STANDARD DEVIATION
    sd_n = [math.sqrt(x / count) for x in cumulative_squared_n]
    sd_m = [math.sqrt(x / count) for x in cumulative_squared_m]
    sd_p = [math.sqrt(x / count) for x in cumulative_squared_p]
    sd_mi = [math.sqrt(x / count) for x in cumulative_squared_mi]
    sd_mi_com = [math.sqrt(x / count) for x in cumulative_squared_mi_com]

    # GENERATE DIFFERENT FILES FOR EACH AND EVERY RUN.
    '''filename = f"file{loop}.txt"
    with open(filename, 'w') as z:
        for (a, b, c, d, e, f) in zip_longest(new_t, new_n, new_m, new_p, new_mi, new_mi_com):
            z.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(str(a), str(b), str(c), str(d), str(e), str(f)))'''

    # PRINTING THE MEAN DATA AS AN OUTPUT FILE
    with open("mean_bound_1.txt", 'w') as z:
        for (a, b, c, d, e, f) in zip_longest(mean_t, mean_n, mean_m, mean_p, mean_mi, mean_mi_com):
            z.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(str(a), str(b), str(c), str(d), str(e), str(f)))

    # PRINTING THE STANDARD DEVIATION DATA AS AN OUTPUT FILE
    with open("sd_bound_1.txt", 'w') as z:
        for (a, b, c, d, e, f) in zip_longest(mean_t, sd_n, sd_m, sd_p, sd_mi, sd_mi_com):
            z.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(str(a), str(b), str(c), str(d), str(e), str(f)))

    # PRINT THE ITERATIONS WHICH ARE COMPLETED AND INCREMENT THE ITERATIONS ONE-BY-ONE
    print(loop)
    loop += 1

# GETTING THE OUTPUT OF CLOCK AS WELL AS CPU TIME
end_clock_time = time.time()
end_process_time = time.process_time()
elapsed_clock_time = end_clock_time - clock_time
elapsed_process_time = end_process_time - process_time
print("ELAPSED CLOCK TIME:\t", elapsed_clock_time)
print("ELAPSED PROCESS TIME\t", elapsed_process_time)
