import csv
import os
import string

def is_all_uppercase(word):
	for c in word:
		if c not in string.ascii_uppercase:
			return False

	return True

dates_folder = "dates"
morphology_folder = "morphology"

dates_filename = "beast-dates.tsv"
morphology_filename = "morphology.nex"
extant_morphology_filename = "morphology-extant.nex"

dates_path = os.path.join(dates_folder, dates_filename)
morphology_path = os.path.join(morphology_folder, morphology_filename)
extant_morphology_path = os.path.join(morphology_folder, extant_morphology_filename)

dates_file = open(dates_path)
dates_reader = csv.reader(dates_file, dialect = csv.excel_tab)

extant_taxa = set() # includes recently extinct, i.e. Dusicyon

for row in dates_reader:
	taxon = row[0]
	date = float(row[1])
	if date == 0:
		extant_taxa.add(taxon)

dates_file.close()

n_taxa = len(extant_taxa)

morphology_file = open(morphology_path)
morphology_nexus = morphology_file.read()

current_word = None
current_word_start = 0

current_block = None
current_block_start = 0

block_boundaries = []

for i, c in enumerate(morphology_nexus):
	if c.isspace() or c == ";":
		if current_word != None:
			if current_block == None and is_all_uppercase(current_word):
				current_block = current_word
				current_block_start = current_word_start
			current_word = None
		if c == ";":
			current_block_end = i + 1
			block_boundaries.append((current_block, current_block_start, current_block_end))
			current_word = None
			current_block = None
	else:
		if current_word == None:
			current_word = c
			current_word_start = i
		else:
			current_word += c

morphology_file.close()

n_blocks = len(block_boundaries)

extant_morphology_file = open(extant_morphology_path, "w")

previous_block_end = 0

for block_i in range(n_blocks):
	block_name, block_start, block_end = block_boundaries[block_i]
	block = morphology_nexus[block_start:block_end]

	extant_morphology_file.write(morphology_nexus[previous_block_end:block_start])

	if block_name == "DIMENSIONS":
		if "NTAX" in block:
			extant_morphology_file.write("DIMENSIONS NTAX=%d;" % (n_taxa))
		else:
			extant_morphology_file.write(block)
	elif block_name == "TAXLABELS" or block_name == "MATRIX":
		block_lines = block.split("\n")
		extant_morphology_file.write(block_lines[0] + "\n")
		for line in block_lines[1:-1]:
			if line.strip().split()[0] in extant_taxa:
				extant_morphology_file.write(line + "\n")
		extant_morphology_file.write(block_lines[-1] + "\n")
	else:
		extant_morphology_file.write(block)

	previous_block_end = block_end

extant_morphology_file.write("\n")

extant_morphology_file.close()
