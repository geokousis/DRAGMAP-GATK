#!/bin/python3
import pybedtools

# Load the BED file
bed_file = "tmp.exons.bed"

# Create a BedTool object
bed_tool = pybedtools.BedTool(bed_file)

# Merge overlapping exons
merged_exons = bed_tool.merge()

# Save the merged exons to a new file
output_file = "merged_exons.bed"
merged_exons.saveas(output_file)

print("BED file Ready")
