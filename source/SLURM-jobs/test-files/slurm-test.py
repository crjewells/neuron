import os
import csv
import random
import numpy as np

id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
group = int(np.ceil(id/100)) # GROUP NUMBER: All have the same experimental variable
trial = id % 100         # TRIAL NUMBER: 1-100 for each experimental variable

sim_var = np.random.rand() + group

filename = "test/group" + str(group) + "-t" + str(trial) + ".csv"

with open(filename, 'w', newline='') as f:
     
    # using csv.writer method from CSV package
    write = csv.writer(f, delimiter=',')
    write.writerow([sim_var])