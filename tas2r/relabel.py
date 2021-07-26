import dendropy

t = dendropy.Tree.get(path = "pastajob.tre", schema = "newick")

for taxon in t.taxon_namespace:
	taxon.label = taxon.label.replace(" ", "_")

	if taxon.label.startswith("NA_"):
		taxon.label = taxon.label[3:]
	elif "," in taxon.label:
		taxon.label = taxon.label[:taxon.label.rfind("_", 0, taxon.label.rfind(","))]

t.write(path = "tas2r.tree", schema = "newick")
