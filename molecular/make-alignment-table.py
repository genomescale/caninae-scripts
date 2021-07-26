import csv
import os

msa_stats_filename = "msa-stats.csv"

msa_stats_file = open(msa_stats_filename)
msa_stats_reader = csv.reader(msa_stats_file)

concatenated_length_before = 0
concatenated_length_after = 0

concatenated_missing_before = 0
concatenated_missing_after = 0

concatenated_matrix_before = 0
concatenated_matrix_after = 0

concatenated_full_matrix_before = 0
concatenated_full_matrix_after = 0

concatenated_total_missing_before = 0
concatenated_total_missing_after = 0

for row_i, row in enumerate(msa_stats_reader):
	if row_i > 0:
		locus_name = row[0]
		n_species_present = int(row[1])
		n_species_absent = int(row[2])
		prank_locus_length = int(row[3])
		trimmed_locus_length = int(row[4])
		prank_missing = int(row[5])
		trimmed_missing = int(row[6])

		prank_matrix_size = prank_locus_length * n_species_present
		trimmed_matrix_size = trimmed_locus_length * n_species_present

		prank_full_matrix_size = prank_locus_length * (n_species_present + n_species_absent)
		trimmed_full_matrix_size = trimmed_locus_length * (n_species_present + n_species_absent)

		prank_total_missing = prank_missing + (prank_locus_length * n_species_absent)
		trimmed_total_missing = trimmed_missing + (trimmed_locus_length * n_species_absent)

		concatenated_length_before += prank_locus_length
		concatenated_length_after += trimmed_locus_length

		concatenated_missing_before += prank_missing
		concatenated_missing_after += trimmed_missing

		concatenated_matrix_before += prank_matrix_size
		concatenated_matrix_after += trimmed_matrix_size

		concatenated_full_matrix_before += prank_full_matrix_size
		concatenated_full_matrix_after += trimmed_full_matrix_size

		concatenated_total_missing_before += prank_total_missing
		concatenated_total_missing_after += trimmed_total_missing

		prank_present_missing_percent = 100 * prank_missing / prank_matrix_size
		trimmed_present_missing_percent = 100 * trimmed_missing / trimmed_matrix_size

		prank_total_missing_percent = 100 * prank_total_missing / prank_full_matrix_size
		trimmed_total_missing_percent = 100 * trimmed_total_missing / trimmed_full_matrix_size

		latex_row = "%s & %d & %d & %.2f\\%% & %.2f\\%% & %.2f\\%% & %.2f\\%%\\\\" % (locus_name, prank_locus_length, trimmed_locus_length, prank_present_missing_percent, trimmed_present_missing_percent, prank_total_missing_percent, trimmed_total_missing_percent)
		print(latex_row)

concatenated_missing_percent_before = 100 * concatenated_missing_before / concatenated_matrix_before
concatenated_missing_percent_after = 100 * concatenated_missing_after / concatenated_matrix_after

concatenated_total_missing_percent_before = 100 * concatenated_total_missing_before / concatenated_full_matrix_before
concatenated_total_missing_percent_after = 100 * concatenated_total_missing_after / concatenated_full_matrix_after

latex_row = "Concatenated & %d & %d & %.2f\\%% & %.2f\\%% & %.2f\\%% & %.2f\\%%\\\\" % (concatenated_length_before, concatenated_length_after, concatenated_missing_percent_before, concatenated_missing_percent_after, concatenated_total_missing_percent_before, concatenated_total_missing_percent_after)
print("\\hline\n" + latex_row)
