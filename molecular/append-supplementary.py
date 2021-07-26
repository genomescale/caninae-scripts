from fasta import *

ncbi_filename = "ncbi.fasta"
supplementary_filename = "726_2017_2422_MOESM3_ESM.txt"
combined_filename = "combined.fasta"

ncbi_sequences = read_fasta(ncbi_filename)
supplementary_sequences = read_fasta(supplementary_filename)

combined_sequences = {}

for label in sorted(ncbi_sequences):
	combined_sequences[label] = ncbi_sequences[label]

for label in sorted(supplementary_sequences):
	binomial = "Lycaon_pictus"
	if label.lower().startswith("huntingdog"):
		new_label = "NA Lycaon pictus " +label[10:]
		combined_sequences[new_label] = supplementary_sequences[label]

write_fasta(combined_filename, combined_sequences)
