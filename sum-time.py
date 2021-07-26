import csv
import numpy
import os
import string
import sys

analysis = sys.argv[1].rstrip("/")
n_chains = 4

total_time = 0.0

for chain_i in range(n_chains):
	time_filename = "%s.%s.time" % (analysis, string.ascii_lowercase[chain_i])
	time_path = os.path.join(analysis, time_filename)

	time_file = open(time_path)
	time_line = time_file.readline()
	user_time, sys_time, wall_time, cpu_usage, mem = time_line.split(maxsplit = 4)

	total_time += float(user_time[:-4]) + float(sys_time[:-6])

print(total_time / 3600.0)
