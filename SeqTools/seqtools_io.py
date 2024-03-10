#!/usr/bin/env python3

import gzip
import time
from pathlib import Path
from Bio import SeqIO

def overwhy(dir, basename, ext):
	"""
	Prevent file overwriting by appending formatted date+timestamp for non-unique filename.
	Example:
		Following non-unique filename query:
			dir: output_dir/, basename: results, ext:txt
		Yields:
			"output_dir/result_20240307_160534.txt"
	"""
	t_str = time.strftime("%Y%m%d_%H%M%S")
	basepath = Path(dir) / f"{basename}.{ext}"
	if not basepath.exists():
		result = str(basepath)
	else:
		result = str(Path(dir) / basename) + "_" + t_str + "." + ext
	return result

def guess_file_format(file):
    if any([file.endswith(suf) or file.endswith(suf + '.gz') for suf in 'fasta fna fas fa fst FA FASTA'.split()]): return 'fasta'
    if any([file.endswith(suf) or file.endswith(suf + '.gz') for suf in 'fastq fq FQ FASTQ'.split()]): return 'fastq'
    return None	
	
def SeqIO_load(input_file:str):	
	"""
	Usage: 
	records = list(SeqIO_load("reads.fastq.gz"))
	"""	
	file_type = guess_file_format(input_file)
	if not file_type:
		raise TypeError("File type could not be determined - please ensure file extension reflects either fasta or fastq reads file.")
	if input_file.split(".")[-1] in ["gz", "GZ", "gzip"]:
		with gzip.open(input_file, "rt") as f:
			records = list(SeqIO.parse(f, file_type))
	else:
		records = list(SeqIO.parse(input_file, file_type))

	return records		