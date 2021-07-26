import dendropy
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

tree = dendropy.Tree.get(path = input_path,
	schema = "nexus",
	rooting = "force-rooted",
	preserve_underscores = True,
)

for leaf_node in tree.leaf_nodes():
	leaf_node.annotations.drop()
	leaf_node.annotations.add_new("sa", "FALSE")

for internal_node in tree.internal_nodes():
	internal_node.annotations.drop()
	is_sampled_ancestor = False
	left_child, right_child = internal_node.child_nodes()
	if left_child.is_leaf() or right_child.is_leaf():
		if left_child.edge_length != right_child.edge_length:
			if left_child.edge_length == 0:
				sa_leaf = left_child
			elif right_child.edge_length == 0:
				sa_leaf = right_child
			else:
				sa_leaf = None

			if sa_leaf != None:
				is_sampled_ancestor = True

				sa_taxon = sa_leaf.taxon
				internal_node.remove_child(sa_leaf)
				internal_node.taxon = sa_taxon

	if is_sampled_ancestor:
		internal_node.annotations.add_new("sa", "TRUE")
	else:
		internal_node.annotations.add_new("sa", "FALSE")

tree.write(path = output_path, schema = "newick", unquoted_underscores = True, suppress_annotations = False, annotations_as_nhx = True, suppress_rooting = True)
