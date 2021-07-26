import csv
import dendropy
import os
import sys

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

mcc_path = sys.argv[1]
nodes_path = sys.argv[2] # for node labels
branches_path = sys.argv[3] # for branch length differences

mcc_split = os.path.split(mcc_path)[1].split(".")
branches_split = os.path.split(branches_path)[1].split(".")
output_path = mcc_split[0] + "-" + branches_split[0] + ".tree"

common_namespace = dendropy.TaxonNamespace()
common_namespace.new_taxa(taxon_order)

mcc_tree = dendropy.Tree.get(path = mcc_path,
	schema = "newick",
	preserve_underscores = True,
	rooting = "force-rooted",
	taxon_namespace = common_namespace,
)

mcc_tree.calc_node_ages()

node_labels = {}

nodes_file = open(nodes_path)
nodes_reader = csv.reader(nodes_file)
for row_i, row in enumerate(nodes_reader):
	if row_i > 0:
		bitmask = int(row[0])
		label = row[1]

		if label != "":
			node_labels[bitmask] = label

nodes_file.close()

branch_lengths = {}
branch_probabilities = {}

branches_file = open(branches_path)
branches_reader = csv.reader(branches_file)
for row_i, row in enumerate(branches_reader):
	if row_i > 0:
		child_bitmask = int(row[0])
		parent_bitmask = int(row[1])
		probability = float(row[5])
		length = float(row[6])

		if child_bitmask in branch_lengths:
			branch_lengths[child_bitmask][parent_bitmask] = length
			branch_probabilities[child_bitmask][parent_bitmask] = probability
		else:
			branch_lengths[child_bitmask] = {parent_bitmask: length}
			branch_probabilities[child_bitmask] = {parent_bitmask: probability}

branches_file.close()

diff_path = mcc_split[0] + "-" + branches_split[0] + ".csv"
diff_file = open(diff_path, "w")
diff_writer = csv.writer(diff_file)
diff_header = ["id", "start", "end", "ldiff"]
diff_writer.writerow(diff_header)

for branch in mcc_tree.postorder_edge_iter():
	child_node = branch.head_node
	child_bitmask = int(child_node.annotations.get_value("id"))

	if not child_node.is_leaf():
		if child_bitmask in node_labels:
			child_node.annotations["tag"] = node_labels[child_bitmask]
		else:
			child_node.annotations.drop(name = "tag")

	parent_node = branch.tail_node
	if parent_node != None:
		parent_bitmask = int(parent_node.annotations.get_value("id"))

		if (child_bitmask in branch_lengths) and (parent_bitmask in branch_lengths[child_bitmask]):
			concat_length = branch_lengths[child_bitmask][parent_bitmask]
			sb2_length = branch.length
			length_difference = concat_length - sb2_length

			probability = branch_probabilities[child_bitmask][parent_bitmask]

			parent_height = parent_node.age
			child_height = child_node.age

			if length_difference > 0.0:
				diff_row = [child_bitmask, -parent_height, -parent_height - length_difference, length_difference]
			else:
				diff_row = [child_bitmask, -child_height, -child_height + length_difference, length_difference]

			diff_writer.writerow(diff_row)
		else:
			length_difference = 0.0
			probability = 0.0

		child_node.annotations.add_new("length_difference", length_difference)
		child_node.annotations.add_new("probability", probability)

diff_file.close()

mcc_tree.write(path = output_path, schema = "newick", annotations_as_nhx = True, suppress_annotations = False, suppress_rooting = True)
