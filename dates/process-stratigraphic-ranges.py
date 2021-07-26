import csv

pixel_scale = 17.593 # there are this many horizontal pixels per million years
pixel_offset = 75.183 # this is the horizontal coordinate in pixels of the beginning of the x-axis
mya_offset = 30.0 # this is the beginning of the x-axis in mya

raw_filename = "raw-stratigraphic-ranges.csv"
processed_filename = "stratigraphic-ranges.csv"
beast_filename = "beast-dates.tsv"

raw_file = open(raw_filename)
processed_file = open(processed_filename, "w")
beast_file = open(beast_filename, "w")

raw_reader = csv.reader(raw_file)
processed_writer = csv.writer(processed_file)
beast_writer = csv.writer(beast_file, dialect = csv.excel_tab)

processed_header = ["species", "start_mya", "end_mya", "midpoint"]
processed_writer.writerow(processed_header)

for row_i, row in enumerate(raw_reader):
	if row_i > 0:
		species = row[0]
		if row[1] == "":
			range_start_mya = 0.0
			range_end_mya = 0.0
		else:
			range_start_pixels = float(row[1])
			range_end_pixels = float(row[2])

			range_start_mya = mya_offset - ((range_start_pixels - pixel_offset) / pixel_scale)
			range_end_mya = mya_offset - ((range_end_pixels - pixel_offset) / pixel_scale)

		midpoint = round(10.0 * (range_start_mya + range_end_mya)) / 20.0 # round midpoint to nearest 0.05 million years

		processed_row = [species, range_start_mya, range_end_mya, midpoint]
		processed_writer.writerow(processed_row)

		if range_end_mya <= 0.05:
			beast_row = [species, 0.0]
		else:
			beast_row = [species, midpoint]

		beast_writer.writerow(beast_row)

raw_file.close()
processed_file.close()
