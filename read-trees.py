import collections
import csv
import dendropy
import numpy
import os
import string
import sys

from dendropy.calculate import statistics

taxon_order = [
	"Atelocynus_microtis",
	"Canis_aureus",
	"Canis_latrans",
	"Canis_lupus",
	"Canis_simensis",
	"Cerdocyon_thous",
	"Chrysocyon_brachyurus",
	"Cuon_alpinus",
	"Lupulella_adusta",
	"Lupulella_mesomelas",
	"Lycalopex_culpaeus",
	"Lycalopex_fulvipes",
	"Lycalopex_griseus",
	"Lycalopex_gymnocercus",
	"Lycalopex_sechurae",
	"Lycalopex_vetulus",
	"Lycaon_pictus",
	"Nyctereutes_procyonoides",
	"Otocyon_megalotis",
	"Speothos_venaticus",
	"Urocyon_cinereoargenteus",
	"Urocyon_littoralis",
	"Vulpes_cana",
	"Vulpes_chama",
	"Vulpes_corsac",
	"Vulpes_ferrilata",
	"Vulpes_lagopus",
	"Vulpes_macrotis",
	"Vulpes_rueppellii",
	"Vulpes_vulpes",
	"Vulpes_zerda",
]

node_label_order = list(string.digits + string.ascii_uppercase)

common_namespace = dendropy.TaxonNamespace()
common_namespace.new_taxa(taxon_order)

trees_path = sys.argv[1]
log_path = sys.argv[2]
burnin = 0

output_base = trees_path[:trees_path.rfind(".")]
if output_base.endswith(".combined"):
	output_base = trees_path[:trees_path.rfind(".combined")]

clock_path    = output_base + ".clock.csv"
nexus_path    = output_base + ".mcc.nexus"
newick_path   = output_base + ".mcc.tree"
nodes_path    = output_base + ".nodes.csv"
branches_path = output_base + ".branches.csv"

trees = dendropy.TreeList.get(path = trees_path,
	schema = "nexus",
	rooting = "force-rooted",
	preserve_underscores = True,
	taxon_namespace = common_namespace,
	tree_offset = burnin
)
n_samples = len(trees)

log_file = open(log_path)
log_reader = csv.reader(log_file, dialect = csv.excel_tab)
log_header = next(log_reader)

if output_base != "morphology-only":
	molecular_rate_i = log_header.index("strictClockRate.c:molecular")
if output_base != "molecular-only":
	morphology_rate_i = log_header.index("strictClockRate.c:morphology")

clock_file = open(clock_path, "w")
clock_writer = csv.writer(clock_file)
clock_header = ["Rate", "Bitmask", "Age", "SubsPerSite", "ChildLeafCount"]

clock_writer.writerow(clock_header)

node_counts = {}
node_ages = {}

branch_counts = {}
branch_lengths = {}
branch_popsizes = {}

for tree_i, tree in enumerate(trees):
	child_leaf_counts = collections.Counter()

	log_row = next(log_reader)

	if output_base != "morphology-only":
		tree_molecular_rate = float(log_row[molecular_rate_i])
	if output_base != "molecular-only":
		tree_morphology_rate = float(log_row[morphology_rate_i])

	tree.encode_bipartitions()
	tree.calc_node_ages(ultrametricity_precision = 1e-3)

	for branch in tree.postorder_edge_iter():
		if branch.tail_node == None:
			parent_bitmask = -1
			length = 0.0
		else:
			parent_bitmask = branch.tail_node.bipartition.leafset_bitmask
			length = branch.length

		child_node = branch.head_node
		child_bitmask = child_node.bipartition.leafset_bitmask

		if child_node.is_leaf():
			child_leaf_counts[parent_bitmask] += 1
			n_child_leaves = 0
		else:
			n_child_leaves = child_leaf_counts[child_bitmask]

		if output_base != "morphology-only":
			molecular_subs_per_site = child_node.age * tree_molecular_rate
			molecular_row = ["molecular", child_bitmask, child_node.age, molecular_subs_per_site, n_child_leaves]
			clock_writer.writerow(molecular_row)
		if output_base != "molecular-only":
			morphology_subs_per_site = child_node.age * tree_morphology_rate
			morphology_row = ["morphology", child_bitmask, child_node.age, morphology_subs_per_site, n_child_leaves]
			clock_writer.writerow(morphology_row)

		popsize_annotation = child_node.annotations.get_value("dmv")
		if popsize_annotation == None:
			popsize = 0.0
		else:
			popsize = float(popsize_annotation[0])

		if child_bitmask in node_counts:
			node_counts[child_bitmask] += 1
			node_ages[child_bitmask].append(child_node.age)
			if parent_bitmask in branch_counts[child_bitmask]:
				branch_counts[child_bitmask][parent_bitmask] += 1
				branch_lengths[child_bitmask][parent_bitmask].append(length)
				branch_popsizes[child_bitmask][parent_bitmask].append(popsize)
			else:
				branch_counts[child_bitmask][parent_bitmask] = 1
				branch_lengths[child_bitmask][parent_bitmask] = [length]
				branch_popsizes[child_bitmask][parent_bitmask] = [popsize]
		else:
			node_counts[child_bitmask] = 1
			node_ages[child_bitmask] = [child_node.age]
			branch_counts[child_bitmask] = {parent_bitmask: 1}
			branch_lengths[child_bitmask] = {parent_bitmask: [length]}
			branch_popsizes[child_bitmask] = {parent_bitmask: [popsize]}

clock_file.close()

mcc_node_labels = {}
mcc_branch_labels = {}

mcc_tree = trees.maximum_product_of_split_support_tree()
for node in mcc_tree.postorder_node_iter():
	bitmask = node.bipartition.leafset_bitmask
	node.age = numpy.mean(node_ages[bitmask])
	node.annotations.drop(name = "dmv")

mcc_tree.set_edge_lengths_from_node_ages()

label_i = 1
for node in mcc_tree.ageorder_node_iter(descending = True):
	bitmask = node.bipartition.leafset_bitmask
	node.annotations.add_new("id", bitmask)

	if node.is_leaf():
		label = node.taxon.label
	else:
		label = node_label_order[label_i]
		node.annotations.add_new("tag", label)
		label_i += 1

	mcc_node_labels[bitmask] = label

for branch in mcc_tree.postorder_edge_iter():
	if branch.tail_node == None:
		parent_bitmask = -1
		length = 0.0
	else:
		parent_bitmask = branch.tail_node.bipartition.leafset_bitmask
		length = branch.length

	child_node = branch.head_node
	child_bitmask = child_node.bipartition.leafset_bitmask

	mcc_branch_labels[child_bitmask] = (parent_bitmask, mcc_node_labels[child_bitmask])

mcc_tree.write(path = nexus_path, schema = "nexus", unquoted_underscores = True, suppress_rooting = True)
mcc_tree.write(path = newick_path, schema = "newick", unquoted_underscores = True, suppress_annotations = False, annotations_as_nhx = True, suppress_rooting = True)

nodes_file = open(nodes_path, "w")
nodes_writer = csv.writer(nodes_file)
nodes_header = ["bitmask", "label", "frequency", "probability", "age_mean", "age_sd", "age_low", "age_high", "log_age_mean", "log_age_sd"]
nodes_writer.writerow(nodes_header)

for bitmask in sorted(node_counts):
	if bitmask in mcc_node_labels:
		label = mcc_node_labels[bitmask]
	else:
		label = ""

	frequency = node_counts[bitmask]
	probability = float(frequency) / float(n_samples)

	ages = node_ages[bitmask]
	age_mean = numpy.mean(ages)

	if 0.0 in ages:
		log_age_mean = "NA"
	else:
		log_ages = numpy.log(ages)
		log_age_mean = numpy.mean(log_ages)

	if frequency > 10:
		age_sd = numpy.std(ages)
		age_low, age_high = statistics.empirical_hpd(ages)

		if 0.0 in ages:
			log_age_sd = "NA"
		else:
			log_age_sd = numpy.std(log_ages)
	else:
		age_sd = "NA"
		age_low = min(ages)
		age_high = max(ages)

		log_age_sd = "NA"

	output_row = [bitmask, label, frequency, probability, age_mean, age_sd, age_low, age_high, log_age_mean, log_age_sd]
	nodes_writer.writerow(output_row)

nodes_file.close()

branches_file = open(branches_path, "w")
branches_writer = csv.writer(branches_file)
branches_header = ["child_bitmask", "parent_bitmask", "child_label", "parent_label", "frequency", "probability", "length_mean", "length_sd", "length_low", "length_high", "popsize_mean", "popsize_sd", "popsize_low", "popsize_high"]
branches_writer.writerow(branches_header)

for child_bitmask in sorted(branch_counts):
	for parent_bitmask in sorted(branch_counts[child_bitmask]):
		if child_bitmask in mcc_branch_labels and mcc_branch_labels[child_bitmask][0] == parent_bitmask:
			child_label = mcc_branch_labels[child_bitmask][1]
		else:
			child_label = ""
			parent_label = ""

		if child_label != "" and parent_bitmask in mcc_branch_labels:
			parent_label = mcc_branch_labels[parent_bitmask][1]

		frequency = branch_counts[child_bitmask][parent_bitmask]
		probability = float(frequency) / float(n_samples)

		lengths = branch_lengths[child_bitmask][parent_bitmask]
		popsizes = branch_popsizes[child_bitmask][parent_bitmask]

		length_mean = numpy.mean(lengths)
		popsize_mean = numpy.mean(popsizes)

		if frequency > 10:
			length_sd = numpy.std(lengths)
			length_low, length_high = statistics.empirical_hpd(lengths)

			popsize_sd = numpy.std(popsizes)
			popsize_low, popsize_high = statistics.empirical_hpd(popsizes)
		else:
			length_sd = "NA"
			length_low = min(lengths)
			length_high = max(lengths)

			popsize_sd = "NA"
			popsize_low = min(popsizes)
			popsize_high = max(popsizes)

		output_row = [child_bitmask, parent_bitmask, child_label, parent_label, frequency, probability, length_mean, length_sd, length_low, length_high, popsize_mean, popsize_sd, popsize_low, popsize_high]
		branches_writer.writerow(output_row)

branches_file.close()
