import csv
import random
from fasta import *

iupac_codes = {
	"M": ("A", "C"),
	"R": ("A", "G"),
	"W": ("A", "T"),
	"S": ("C", "G"),
	"Y": ("C", "T"),
	"K": ("G", "T"),
}

def phase_sequence(unphased_sequence):
	phased_sequence = ""

	for c in unphased_sequence.upper():
		if c in iupac_codes:
			phased_sequence += iupac_codes[c][random.randint(0, 1)]
		else:
			assert c in "-ACGNT"
			phased_sequence += c

	return phased_sequence

label_mapping_filename = "label-mapping.csv"
whitelist_filename = "whitelist.csv"
msa_stats_filename = "msa-stats.csv"

label_mapping_file = open(label_mapping_filename)
label_mapping_reader = csv.reader(label_mapping_file)

whitelist_file = open(whitelist_filename)
whitelist = whitelist_file.read().strip().split()
whitelist_file.close()

label_mapping = {}

for label, species in label_mapping_reader:
	label_mapping[label] = species

molecular_species = set()

trimmed_locus_lengths = {}
trimmed_missing_counts = {}

for locus_name in whitelist:
	trimmed_filename = locus_name + ".trimmed.fasta"

	trimmed_locus = read_fasta(trimmed_filename)

	filtered_locus = {}

	n_species = 0
	for label in sorted(trimmed_locus):
		if label in label_mapping:
			species = label_mapping[label]
			trimmed_sequence = trimmed_locus[label]

			molecular_species.add(species)

			filtered_locus[species + "_x"] = phase_sequence(trimmed_sequence)

			if n_species == 0:
				trimmed_locus_lengths[locus_name] = len(trimmed_sequence)
				trimmed_missing_counts[locus_name] = trimmed_sequence.count("-")
			else:
				trimmed_missing_counts[locus_name] += trimmed_sequence.count("-")

			n_species += 1

	filtered_locus_filename = locus_name + ".filtered.fasta"

	write_fasta(filtered_locus_filename, filtered_locus)

msa_stats_file = open(msa_stats_filename, "w")
msa_stats_writer = csv.writer(msa_stats_file)

msa_stats_header = ["Locus", "Present Species Count", "Missing Species Count", "Original Length", "Trimmed Length", "Original Missing", "Trimmed Missing"]
msa_stats_writer.writerow(msa_stats_header)

for locus_name in whitelist:
	prank_filename = locus_name + ".best.fas"

	prank_locus = read_fasta(prank_filename)

	prank_locus_length = 0
	prank_locus_missing = 0

	n_species = 0
	for label in sorted(prank_locus):
		if label in label_mapping:
			prank_sequence = prank_locus[label]

			prank_locus_missing += prank_sequence.count("-")

			if n_species == 0:
				prank_locus_length = len(prank_sequence)

			n_species += 1

	msa_stats_writer.writerow([locus_name,
		n_species,
		len(molecular_species) - n_species,
		prank_locus_length,
		trimmed_locus_lengths[locus_name],
		prank_locus_missing,
		trimmed_missing_counts[locus_name]])

label_mapping_file.close()
whitelist_file.close()
msa_stats_file.close()
