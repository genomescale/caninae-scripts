import collections
import subprocess
import csv
from fasta import *

# replace these paths with the right ones for your system

prank_path = "/path/to/prank"
trimal_path = "/path/to/trimal"

locus_symbol_map = [
	("APOB", "apolipoprotein B"),
	("BDNF", "brain-derived neurotrophic factor"),
	("BRCA1", "breast and ovarian cancer susceptibility protein"),
	("CHST12", "carbohydrate sulfotransferase 12"),
	("CMKOR1", "chemokine orphan receptor 1"),
	("GHR", "GHR (GHR) gene, exon 10"),
	("GHR", "growth hormone receptor gene, exon 10"),
	("RAG1", "recombination activating protein"),
	("TMEM20", "transmembrane protein"),
	("VANGL2", "vang-like 2"),
	("VWF", "von Willebrand factor"),
	("Ch14", "12.149"),
	("Ch21", "7.92"),
	("Ch24", "1.75"),
	("FGFR3", "growth factor receptor"),
	("CHRNA1", "cholinergic receptor"),
	("CHRNA1", "nicotinic alpha polypeptide 1 precursor"),
	("CYP1A1", "cytochrome P450 gene"),
	("FES", "feline sarcoma protooncogene"),
	("GHR", "growth hormone receptor (GHR) gene, exon 9"),
	("VTN", "VTN"),
	("VTN", "vitronectin"),
	("TRSP", "tRNA-Sec"),
	("RPPH1", "RNase P RNA"),
	("COII", "cytochrome oxidase subunit II"),
]

locus_segments = {
	"APOB": [(240289, 240321, "APOBS1"), (240322, 240355, "APOBS2"), (447809, 447831, "APOBcds")],
	"BRCA1": [(0, 240423, "BRCA1S1"), (240424, 999999, "BRCA1S2")],
	"CHRNA1": [(447855, 447877, "CHRNA1intron"), (885308, 885331, "CHRNA1exons"), (239439, 239446, "CHRNA1exons")],
	"GHR": [(239463, 239470, "GHRex10"), (885379, 885402, "GHRex10"), (240653, 240686, "GHRex09"), (447924, 447946, "GHRcds")],
}

garbage_accessions = [
	"KT426724.1", # Chrysocyon brachyurus T2R1, duplicate of Vulpes T2R1
	"KT426731.1", # Chrysocyon brachyurus T2R9, appears to be a Vulpes T2R9
	"KT426793.1", # Vulpes corsac T2R40, duplicate of Vulpes corsac T2R4
	"KT426799.1", # Vulpes corsac T2R67, duplicate of Chrysocyon T2R67
	"KT426819.1", # Vulpes ferrilata T2R67, duplicate of Vulpes zerda T2R42
	"KT426829.1", # Vulpes vulpes T2R19, many more mutations than in other species
	"KT426846.1", # Vulpes zerda T2R9, many more mutations than in other species
	"KX604064.1", # Vulpes zerda T2R7LIKE2, many more mutations than in other species
]

accession_ranges_file = open("accession-ranges.csv")
accession_ranges_reader = csv.reader(accession_ranges_file)

accession_range_sources = {}
for row in accession_ranges_reader:
	ncbi_prefix = row[0]
	first_accession = int(row[1])
	last_accession = int(row[2])
	source = row[3]

	if ncbi_prefix in accession_range_sources:
		if source in accession_range_sources[ncbi_prefix]:
			accession_range_sources[ncbi_prefix][source].append((first_accession, last_accession))
		else:
			accession_range_sources[ncbi_prefix][source] = [(first_accession, last_accession)]
	else:
		accession_range_sources[ncbi_prefix] = {source: [(first_accession, last_accession)]}

accession_ranges_file.close()

ncbi_sequences = read_fasta("combined.fasta")

overrides = {
	"KT426813.1": ("T2R40", "Chrysocyon_brachyurus"),
	"KT426738.1": ("T2R40", "Vulpes_ferrilata"),
	"KT426830.1": ("T2R62", "Vulpes_ferrilata"),
	"KT426818.1": ("T2R34", "Vulpes_vulpes"),
}

loci = collections.defaultdict(dict)

sequence_sources_file = open("sequence-sources.csv", "w")
sequence_sources_writer = csv.writer(sequence_sources_file)
sequence_sources_header = ["accession", "source", "label", "organism", "locus"]
sequence_sources_writer.writerow(sequence_sources_header)

for label in sorted(ncbi_sequences):
	label_split = label.split()
	accession = label_split[0]
	
	accession_prefix = accession[:2]
	accession_number = None
	source = None
	if accession_prefix == "NA":
		source = "Shang2017"
	else:
		accession_number = int(accession[2:accession.rfind(".")])
		for possible_source in accession_range_sources[accession_prefix]:
			for first_accession, last_accession in accession_range_sources[accession_prefix][possible_source]:
				if accession_number >= first_accession and accession_number <= last_accession:
					assert source == None
					source = possible_source

	assert source != None

	if accession in overrides:
		locus_name, organism = overrides[accession]
		assert organism not in loci[locus_name]
		loci[locus_name][organism] = ncbi_sequences[label]
		sequence_sources_writer.writerow([accession, source, label, organism, locus_name])
	elif accession not in garbage_accessions:
		if label_split[1] == "UNVERIFIED:":
			organism_offset = 2
		else:
			organism_offset = 1

		if "isolate" in label_split:
			description_offset = label_split.index("isolate") + 2
		else:
			description_offset = organism_offset + 2

		organism = "_".join(label_split[organism_offset:description_offset])

		locus_name = None

		# One of the bitter taste receptor genes
		for word in label_split[description_offset:]:
			if word.lower().startswith("t2r") or word.lower().startswith("tas2r"):
				locus_name = word.lower().replace("-", "").replace("tas", "t").replace("partial", "").rstrip("p").upper()
				assert organism not in loci[locus_name]
				loci[locus_name][organism] = ncbi_sequences[label]
				sequence_sources_writer.writerow([accession, source, label, organism, locus_name])

		# a gene with symbols in brackets
		if locus_name == None:
			for word in label_split[description_offset:]:
				if word[0] == "(" and word[-1] == ")":
					locus_name = word[1:-1].upper()
					if locus_name == "CYPIA1":
						locus_name = "CYP1A1" # fucking typos in genbank
					if locus_name in locus_segments:
						accession_number = int(accession[2:-2])
						for min_accession, max_accession, segment_name in locus_segments[locus_name]:
							if accession_number >= min_accession and accession_number <= max_accession:
								assert organism not in loci[segment_name]
								loci[segment_name][organism] = ncbi_sequences[label]
								sequence_sources_writer.writerow([accession, source, label, organism, segment_name])
					else:
						assert organism not in loci[locus_name]
						loci[locus_name][organism] = ncbi_sequences[label]
						sequence_sources_writer.writerow([accession, source, label, organism, locus_name])

		# a gene where the symbol is not necessarily given
		if locus_name == None:
			description = " ".join(label_split[description_offset:])
			accession_number = int(accession[2:-2])
			for symbol, prefix in locus_symbol_map:
				if description.startswith(prefix):
					locus_name = symbol
					if locus_name in locus_segments:
						accession_number = int(accession[2:-2])
						for min_accession, max_accession, segment_name in locus_segments[locus_name]:
							if accession_number >= min_accession and accession_number <= max_accession:
								assert organism not in loci[segment_name]
								loci[segment_name][organism] = ncbi_sequences[label]
								sequence_sources_writer.writerow([accession, source, label, organism, segment_name])
					else:
						assert organism not in loci[locus_name]
						loci[locus_name][organism] = ncbi_sequences[label]
						sequence_sources_writer.writerow([accession, source, label, organism, locus_name])

		if locus_name == None:
			assert "mitochondrion," in label_split

sequence_sources_file.close()

for locus_name in sorted(loci):
	locus = loci[locus_name]
	locus_filename = locus_name + ".fasta"
	prank_filename = locus_name + ".best.fas"
	trimmed_filename = locus_name + ".trimmed.fasta"

	write_fasta(locus_filename, locus)

	prank_command = [prank_path, "-DNA", "-d=" + locus_filename, "-o=" + locus_name]
	subprocess.check_call(prank_command)

	trimal_command = [trimal_path, "-in", prank_filename, "-out", trimmed_filename, "-gappyout"]
	subprocess.check_call(trimal_command)

	# iqtree_command = ["iqtree", "-s", trimmed_filename]
	# subprocess.check_call(iqtree_command)
