#!/usr/bin/env python3

import os
import gzip
import random
from pathlib import Path
from Bio import SeqIO
import pandas as pd

def detect_filetype(ext:str):
	try:
		if ext in ['fa', 'FA', 'fasta']:
			file_type = 'fasta'
		if ext in ['fq', 'FQ', 'fastq']:
			file_type = 'fastq'
	except:
		raise TypeError("File type could not be determined - please ensure file extension reflects either fasta or fastq reads file.")
	
	return file_type

def SeqIO_load(input_file:str):	
	"""
	Usage: 
	records = SeqIO_load("reads.fastq.gz")
	"""	
	
	if input_file.split(".")[-1] in ["gz", "GZ", "gzip"]:
		base_file, gz_ext = os.path.splitext(input_file)
		ext = base_file.split(".")[1]
		with gzip.open(input_file, "rt") as f:
			file_type = detect_filetype(ext)
			records = list(SeqIO.parse(f, f"{file_type}"))
	
	else:
		file_type = detect_filetype(input_file.split(".")[-1])		
		records = list(SeqIO.parse(input_file, f"{file_type}"))

	return records

def get_selections(records, subsample:int):
	
	assert len(records) > subsample, f"Subsample number exceeds total no. reads - please enter subsample number < {len(records)}"

	return [random.randint(0, len(records)) for i in range(subsample)]

def get_sequences(records, selections):
	return [str(records[s].seq) for s in selections]

def subsampled_df(records, selections, input_file):
	if input_file.split(".")[-1] in ["gz", "GZ", "gzip"]:
		base_file, gz_ext = os.path.splitext(input_file)
		file_type = base_file.split(".")[1]
	else:
		file_type = detect_filetype(input_file)
	
	if file_type == 'fasta':
		df = pd.DataFrame({
			'id': [records[s].id for s in selections],
			'seq': [records[s].seq for s in selections]
			})
	
	if file_type == 'fastq':
		df = pd.DataFrame({
			'id': [records[s].id for s in selections],
			'seq': [records[s].seq for s in selections],
			'qual': [records[s].letter_annotations['phred_quality'] for s in selections]
			})

	return df

def export_reads(records, subsample, input_file):
	# Determine output filename and filetype
	if input_file.split(".")[-1] in ["gz", "GZ", "gzip"]:
		base_file, gz_ext = os.path.splitext(input_file)
		file_type = base_file.split(".")[1]
		output_file = str(Path(base_file).parent / Path(base_file).stem) + "_subsampled." + file_type
	else:
		file_type = detect_filetype(input_file)
		output_file = str(Path(input_file).parent / Path(input_file).stem) + "_subsampled." + file_type
	
	# Subsample
	selections = get_selections(records, subsample)
	output_records = [records[s] for s in selections]

	# Write subsampled reads file
	with open(output_file, "wt") as f:
		SeqIO.write(output_records, f, file_type)

	msg = f"""
	{Path(input_file).name} subsampled to {subsample} reads
	Subsampled reads successfully written to: {output_file}
	"""
	print(msg)


if __name__ == "__main__":
	import sys
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description = "Randomly subsample reads file")

	parser.add_argument(
		"--reads_file", 
		action = "store", 
		dest = "reads_file", 
		type = str, 
		help = "Input reads file")
	parser.add_argument(
		"--subsample", 
		action = "store", 
		dest = "subsample", 
		type = int, 
		default = 50,
		help = "Number of reads to subsample (default: 50)")
	parser.add_argument(
		"--export",
		action = "store",
		dest = "export",
		type = bool,
		default = False,
		help = "Export subsampled reads file"
	)

	args = parser.parse_args()
	reads_file = args.reads_file
	subsample = args.subsample
	export = args.export
	
	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)
	
	records = SeqIO_load(reads_file)

	selections = get_selections(records, subsample)
	
	sequences = get_sequences(records, selections)

	if export:
		export_reads(records, subsample, reads_file)

	# Prevent trailing characters from being printed to console
	for seq in sequences:
		print(seq.rstrip())