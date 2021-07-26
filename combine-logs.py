import csv
import os
import string
import sys
import subprocess

analysis = sys.argv[1].rstrip("/")

n_burnin = 64
n_chains = 4
subsample_rate = 8

combined_trees_fn = "%s.combined.trees" % (analysis)
combined_trees_file = open(combined_trees_fn, "w")

combined_log_fn = "%s.combined.log" % (analysis)
combined_log_file = open(combined_log_fn, "w")
combined_log_writer = csv.writer(combined_log_file, dialect = csv.excel_tab)

trees_combined_i = 0
log_combined_i = 0
for chain_i in range(n_chains):
	chain_code = string.ascii_lowercase[chain_i]

	chain_name = analysis + "." + chain_code
	chain_trees_fn = chain_name + ".trees"
	chain_log_fn = chain_name + ".log"

	chain_trees_path = os.path.join(analysis, chain_trees_fn)
	chain_log_path = os.path.join(analysis, chain_log_fn)

	chain_trees_file = open(chain_trees_path)

	trees_sample_i = 0
	is_header = True
	for line in chain_trees_file:
		if is_header:
			if line.startswith("tree STATE_"):
				is_header = False
			elif (chain_i == 0):
				combined_trees_file.write(line)

		if not is_header:
			if line.startswith("tree"):
				if (trees_sample_i > n_burnin) and (trees_sample_i % subsample_rate == 0):
					newick_string = line.strip().split()[-1]
					output = "tree STATE_%d = %s\n" % (trees_combined_i, newick_string)
					combined_trees_file.write(output)
					trees_combined_i += 1
				trees_sample_i += 1

	chain_trees_file.close()

	chain_log_file = open(chain_log_path)

	l = chain_log_file.readline()
	while l[0] == "#":
		l = chain_log_file.readline()

	if chain_i == 0:
		combined_log_file.write(l)

	chain_log_reader = csv.reader(chain_log_file, dialect = csv.excel_tab)
	log_sample_i = 0
	for row in chain_log_reader:
		if (log_sample_i > n_burnin) and (log_sample_i % subsample_rate == 0):
			row[0] = log_combined_i
			combined_log_writer.writerow(row)
			log_combined_i += 1
		log_sample_i += 1

	chain_log_file.close()

combined_trees_file.write("End;\n")
combined_trees_file.close()
combined_log_file.close()
