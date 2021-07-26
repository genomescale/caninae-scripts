import csv
import numpy
import os

from dendropy.calculate import statistics

parameter_order = [
	"molecular_clock_rate",
	"morphological_clock_rate",
	"mean_negc",
	"diversification_rate",
	"turnover",
	"sampling_rate",
]

table_labels = [
	"Molecular clock rate ($\\times10^{-3}$)",
	"Morphological clock rate ($\\times10^{-2}$)",
	"Mean population size",
	"Diversification rate ($\\lambda - \\mu$)",
	"Turnover ($\\mu \\div \\lambda$)",
	"Sampling proportion ($\\psi \\div (\\psi + \\mu)$)",
]

def read_log(log_path):
	log_file = open(log_path)
	log_reader = csv.reader(log_file, dialect = csv.excel_tab)

	bitmask_i = -1
	label_i = -1
	probability_i = -1

	sample_values = {parameter: [] for parameter in parameter_order}

	molecular_clock_rate_i = -1
	morphological_clock_rate_i = -1
	mean_negc_i = -1
	diversification_rate_i = -1
	turnover_i = -1
	sampling_rate_i = -1

	for row_i, row in enumerate(log_reader):
		if row_i == 0:
			if "strictClockRate.c:molecular" in row:
				molecular_clock_rate_i = row.index("strictClockRate.c:molecular")

			if "strictClockRate.c:morphology" in row:
				morphological_clock_rate_i = row.index("strictClockRate.c:morphology")
			elif "strictClockRate.c:morphology_extant" in row:
				morphological_clock_rate_i = row.index("strictClockRate.c:morphology_extant")

			if "popMean.Species" in row:
				mean_negc_i = row.index("popMean.Species")

			if "diversificationRateFBD.t:Species" in row:
				diversification_rate_i = row.index("diversificationRateFBD.t:Species")
			else:
				diversification_rate_i = row.index("netDiversificationRate.t:Species")

			if "turnoverFBD.t:Species" in row:
				turnover_i = row.index("turnoverFBD.t:Species")
			else:
				turnover_i = row.index("ExtinctionFraction.t:Species")

			if "samplingProportionFBD.t:Species" in row:
				sampling_rate_i = row.index("samplingProportionFBD.t:Species")
		else:
			if molecular_clock_rate_i != -1:
				molecular_clock_rate = float(row[molecular_clock_rate_i]) * 1000.0
				sample_values["molecular_clock_rate"].append(molecular_clock_rate)

			if morphological_clock_rate_i != -1:
				morphological_clock_rate = float(row[morphological_clock_rate_i]) * 100.0
				sample_values["morphological_clock_rate"].append(morphological_clock_rate)

			if mean_negc_i != -1:
				mean_negc = float(row[mean_negc_i])
				sample_values["mean_negc"].append(mean_negc)

			turnover = float(row[turnover_i])
			sample_values["turnover"].append(turnover)

			diversification_rate = float(row[diversification_rate_i])
			sample_values["diversification_rate"].append(diversification_rate)

			if sampling_rate_i != -1:
				sampling_rate = float(row[sampling_rate_i])
				sample_values["sampling_rate"].append(sampling_rate)

	log_file.close()

	return sample_values

analysis_order = ["concat-total-evidence", "sb2-total-evidence", "sb2-extant-only"]

all_analyses = {}
for analysis in analysis_order:
	log_fn = analysis + ".combined.log"
	sample_values = read_log(log_fn)
	all_analyses[analysis] = sample_values

for parameter_i, parameter in enumerate(parameter_order):
	row = table_labels[parameter_i]
	for analysis in analysis_order:
		row += " & "

		parameter_values = all_analyses[analysis][parameter]
		if len(parameter_values) == 0:
			row += "NA"
		else:
			parameter_mean = numpy.mean(parameter_values)
			parameter_hpd_low, parameter_hpd_high = statistics.empirical_hpd(parameter_values)
			row += "%.2f (%.2f--%.2f)" % (parameter_mean, parameter_hpd_low, parameter_hpd_high)

	row += " \\\\"

	print(row)
