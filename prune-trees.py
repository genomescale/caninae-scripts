import dendropy
import sys

molecular_taxa = [
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

fossil_taxa = []

analysis = sys.argv[1].rstrip("/")

combined_trees_fn = "%s.combined.trees" % (analysis)
pruned_trees_fn = "%s-pruned.combined.trees" % (analysis)

common_namespace = dendropy.TaxonNamespace()
combined_trees = dendropy.TreeList.get(path = combined_trees_fn, schema = "nexus", preserve_underscores = True, rooting = "force-rooted", taxon_namespace = common_namespace)

for taxon in common_namespace:
	if taxon.label.replace(" ", "_") not in molecular_taxa:
		fossil_taxa.append(taxon)

for tree_i, complete_tree in enumerate(combined_trees):
	complete_tree.prune_taxa(fossil_taxa)

combined_trees.purge_taxon_namespace()
combined_trees.write(path = pruned_trees_fn, schema = "nexus", unquoted_underscores = True, suppress_rooting = True)
