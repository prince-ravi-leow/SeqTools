#!/usr/bin/env python3

import random
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from seqtools_io import SeqIO_load, guess_file_format, overwhy

class SeqSub():
	def __init__(self, reads_file):
		self.reads_file = reads_file
		self.format = guess_file_format(reads_file) 
		if not self.format:
			raise Exception("File format could not be determined - check that input file is type fasta or fastq")
		self.records = SeqIO_load(reads_file)

	def __str__(self):
		if self.format == "seq":
			to_print = "\n".join(self.sequences)
		if self.format == "fasta":
			to_print = ""	
			for s in self.selections:
				to_print += f"{self.records[s].format("fasta")}"
		if self.format == "fastq":
			to_print = ""	
			for s in self.selections:
				to_print += f"{self.records[s].format("fastq")}"
		return to_print
	
	def load(self, reads_file):
		self.records = SeqIO_load(reads_file)
		return self
	
	def subsample(self, n_seq):
		self.n_seq = n_seq
		assert len(self.records) > self.n_seq, f"Subsample number exceeds total no. reads - please enter subsample number < {len(self.records)}"

		self.selections = [random.randint(0, len(self.records)-1) for i in range(self.n_seq)]
		self.sequences = [str(self.records[s].seq) for s in self.selections]

	def get_df(self):
		if self.format == 'fasta':
			self.df = pd.DataFrame({
				'id': [self.records[s].id for s in self.selections],
				'seq': [self.records[s].seq for s in self.selections]
				})
		if self.format == 'fastq':
			self.df = pd.DataFrame({
				'id': [self.records[s].id for s in self.selections],
				'seq': [self.records[s].seq for s in self.selections],
				'qual': [self.records[s].letter_annotations['phred_quality'] for s in self.selections]
				})

		return self.df

	def export(self):
		parent = Path(self.reads_file).parent
		basename = Path(self.reads_file).stem
		query = basename + "_subsampled"
		output_file = overwhy(parent, query, self.format)
		output_records = [self.records[s] for s in self.selections]
		with open(output_file, "wt") as f:
			SeqIO.write(output_records, f, self.format)
		
		msg = f"{Path(self.reads_file).name} subsampled to {self.n_seq} reads\nSubsampled reads successfully written to: {output_file}"
		print(msg)

if __name__ == "__main__":
	import sys
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description = "Randomly subsample reads file")

	parser.add_argument("-i", dest = "reads_file", type = str, help = "Input reads file")
	parser.add_argument("-n", dest = "n_seq", type = int, default = 50, help = "Number of reads to subsample (default: 50)")
	parser.add_argument("--export", action = "store_true", dest = "export", default = False, help = "Export subsampled reads file")
	parser.add_argument("--no-print", action = "store_true", dest = "no_print", default = False, help = "Prevent subsampled being printed to console")

	args = parser.parse_args()
	reads_file = args.reads_file
	n_seq = args.n_seq
	export = args.export
	no_print = args.no_print

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit(1)

	if not export and no_print:
		print("So you don't want to print results to console, and you don't want to export them either... so what DO you want to do? :thinking:")
		sys.exit(1)

	sub = SeqSub(reads_file)
	sub.subsample(n_seq)
	
	if not no_print:
		print(sub)
		
	if export:
		sub.export()