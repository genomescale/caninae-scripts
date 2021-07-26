import dendropy
import numpy
import string
import csv
import sys

from dendropy.calculate import statistics

trees_fn = sys.argv[1]

trees = dendropy.TreeList.get(path = trees_fn, schema = "nexus", rooting = "force-rooted")
n_trees = len(trees)
n_steps = 1024
ltt_height = 32.0
ultrametric_threshold = 1e-6 # corresponds to 1 year in the past

ltt_matrix = numpy.ones((n_trees, n_steps + 1), dtype = numpy.int32)

for tree_i, tree in enumerate(trees):
	root_node_distances = tree.calc_node_root_distances()
	root_height = max(root_node_distances)

	lineage_changes = []
	for node in tree.postorder_node_iter():
		node_age = root_height - node.root_distance
		node_quantile = int(n_steps * node_age / ltt_height)

		if not node.is_leaf():
			ltt_matrix[tree_i][:node_quantile + 1] += 1
		elif node_age > ultrametric_threshold:
			ltt_matrix[tree_i][:node_quantile + 1] -= 1

ltt_fn = trees_fn[:trees_fn.rfind(".")] + ".ltt.csv"
ltt_file = open(ltt_fn, "w")
ltt_writer = csv.writer(ltt_file)
ltt_header = ["mya", "mean", "hpd_low", "hpd_high"]
ltt_writer.writerow(ltt_header)

for step_i, step_row in enumerate(ltt_matrix.transpose()):
	step_mya = -(ltt_height * step_i / n_steps)
	step_mean = numpy.mean(step_row)
	step_hpd_low, step_hpd_high = statistics.empirical_hpd(step_row)
	output_row = [step_mya, step_mean, step_hpd_low, step_hpd_high]
	ltt_writer.writerow(output_row)

ltt_file.close()
