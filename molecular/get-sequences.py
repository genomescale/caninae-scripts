import csv
from Bio import Entrez

Entrez.email = "your@email"

ranges_filename = "accession-ranges.csv"
ranges_file = open(ranges_filename)
ranges_reader = csv.reader(ranges_file)

output_filename = "ncbi.fasta"
output_file = open(output_filename, "w")

for row in ranges_reader:
	ncbi_prefix = row[0]
	first_accession = int(row[1])
	last_accession = int(row[2])

	locus_ids = []
	for i in range(first_accession, last_accession + 1):
		ncbi_accession_id = "%s%d" % (ncbi_prefix, i)
		locus_ids.append(ncbi_accession_id)

	search_results = Entrez.read(Entrez.epost("nucleotide", id = ",".join(locus_ids)))
	webenv = search_results["WebEnv"]
	query_key = search_results["QueryKey"]
	handle = Entrez.efetch(db = "nucleotide", rettype = "fasta", retmode = "text", webenv = webenv, query_key = query_key)

	output_file.write(handle.read())

output_file.close()
