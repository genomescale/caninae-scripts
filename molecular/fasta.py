def read_fasta(fasta_path):
	label = None
	sequence = None

	locus_data = {}

	fasta_file = open(fasta_path)

	l = fasta_file.readline()
	while l != "":
		if l[0] == ">":
			if label != None:
				locus_data[label] = sequence

			label = l[1:].strip()
			sequence = ""
		elif label != None:
			sequence += l.strip()

		l = fasta_file.readline()

	fasta_file.close()

	if label != None:
		locus_data[label] = sequence

	return locus_data

def write_fasta(fasta_path, locus_data):
	fasta_file = open(fasta_path, "w")

	for label in sorted(locus_data):
		sequence = locus_data[label]
		fasta_file.write(">%s\n%s\n" % (label, sequence))

	fasta_file.close()
