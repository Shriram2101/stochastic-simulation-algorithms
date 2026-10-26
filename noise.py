# OPEN THE FILES GENERATED FROM MEAN AND STANDARD DEVIATION CODE
with open("mean_bound_1.txt", 'r') as noise_mean,\
     open("sd_bound_1.txt", 'r') as noise_std:

    # READ THE DATA FILES
    data_mean = noise_mean.readlines()
    data_std = noise_std.readlines()

# INITIALIZING AN ARRAY TO STORE THE DATA VALUES
result_std = []

# LOOP THROUGH EACH DATA LINES
for line_mean, line_std in zip(data_mean, data_std):
    # SPLIT THE LINES INTO VALUES
    values_mean = line_mean.split()
    values_std = line_std.split()

    # PERFORM DIVISION SO AS TO GET CV
    result_r2 = [float(x_std) / float(x_mean) for x_std, x_mean in zip(values_std, values_mean)]
    result_std.append(result_r2)

# GENERATE ANOTHER TEXT FILE FOR CV
with open("CV_bound_1.txt", 'w') as f:
    for row in result_std:
        f.write(" ".join(map(str, row)) + "\n")
