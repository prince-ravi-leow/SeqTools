#!/usr/bin/env python3

import random

letter_dict = {
	"DNA" : ["A", "T", "G", "C"],
	"protein" : [c for c in "ACDEFGHIKLMNPQRSTVWY"]
}

# FASTQ PHRED score
def create_qscore_to_ascii_dict():
	qscore_to_ascii = {}
	offset = 33  # Phred+33 encoding
	for i in range(0, 41):
		qscore_to_ascii[i] = chr(i + offset)
	return qscore_to_ascii

qscore_to_ascii_dict = create_qscore_to_ascii_dict()

# Sequence generator function
def gen_seq(seq_length, molecule_type):
	letters = letter_dict.get(molecule_type)
	sequence = ""
	for i in range(seq_length):
		sequence += random.choice(letters)
	return sequence

def gen_qscore(seq_length, q_min, q_max):
	q_line = ""
	for i in range(seq_length):
		q_line += qscore_to_ascii_dict[random.randint(q_min, q_max)]
	return q_line

class SeqGen():
	"""
	Randomly generate DNA sequences (pure sequence, fasta or fastq format
	"""
	def __init__(self):
		self

	def __str__(self):
		"""
		Print generated sequences in specified (file) format
		"""
		return self.sequence_strings

	def gen_sequences(self, n_seq = 10, len_min = 10, len_max = 100, molecule_type = "DNA"):
		"""
		Generate n_seq number of sequences, with len_min, len_max sequence length range 
		molecule_type == 'DNA' by default - set to 'protein' for amino acid sequences
		"""
		self.molecule_type = molecule_type

		self.seq_lengths = [random.randint(len_min, len_max) for i in range(n_seq)]
		self.sequences = [gen_seq(seq_length, self.molecule_type) for seq_length in self.seq_lengths]
		
		return self.sequences

	def gen_read_file(self, format = "fasta", q_min = 10, q_max = 30):
		"""
		'Sequence-only' mode activated by default; set format to 'fasta' or 'fastq' to output sequences in their respective formats. 

		It is possible to specify min and max Q-scores with the arguments q_min and q_max (Default: 10 - 30) 
		"""
		assert format in ["fasta", "fastq", "seq"], "format must be either 'seq', 'fasta' or 'fastq'"
		
		self.sequence_strings = ""
		self.format = format
		
		if self.molecule_type == "protein" and format == "fastq":
			raise Exception("fastq not valid for protein sequences")

		for i, (seq_length, sequence) in enumerate(zip(self.seq_lengths, self.sequences)):
			if format == "fasta":
				self.sequence_strings += f">sequence_{i+1}_length_{seq_length}\n"
				self.sequence_strings += f"{sequence}\n"
			if format == "fastq":
				self.sequence_strings += f"@sequence_{i+1}_length_{seq_length}\n"
				self.sequence_strings += f"{sequence}\n"
				self.sequence_strings += "+\n"
				self.sequence_strings += f"{gen_qscore(seq_length, q_min, q_max)}\n"
			if format == "seq":
				self.sequence_strings += f"{sequence}\n"

		return print(self)
	
	def export(self, filename = None):
		"""
		Export generated sequences to respective format (plain-text file for sequence-only mode)
		"""
		filename = f"generated_sequences.{self.format}" if not filename else filename
		with open(filename, "wt") as f:
			f.writelines(self.sequence_strings)
		print(f"Generated sequences in {self.format} format written to: {filename}")

if __name__ == "__main__":
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description = "Randomly generate DNA sequences")
	
	parser.add_argument("-n", dest = "n_seq", type = int, default = 10, help = "Number of sequences to generate (default : 10)")
	parser.add_argument("--min", dest = "len_min", type = int, default = 10, help = "Min. sequence length")
	parser.add_argument("--max", dest = "len_max", type = int, default = 100,help = "Max. sequence length")
	parser.add_argument("--q_min", dest = "q_min", type = int, default = 10,help = "Min. PHRED score (default: 10)")
	parser.add_argument("--q_max", dest = "q_max", type = int, default = 30,help = "Max. PHRED score (default: 30)")
	parser.add_argument("--format",dest = "format", type = str, default = "seq",help = "Output sequence file type ('seq', 'fasta', 'fastq')")
	parser.add_argument("--mol_type",dest = "mol_type", type = str, default = "DNA", help = "Molecule type ('DNA' or 'protein')")
	parser.add_argument("--export", action = "store_true", dest = "export",default = False,	help = "Export generated sequences to file")
		
	args = parser.parse_args()
	n_seq = args.n_seq
	len_min = args.len_min
	len_max = args.len_max
	format = args.format
	mol_type = args.mol_type
	q_min = args.q_min
	q_max = args.q_max
	export = args.export

	generator = SeqGen()
	generator.gen_sequences(n_seq, len_min, len_max, mol_type)
	generator.gen_read_file(format, q_min, q_max)
	
	if export:
		generator.export()