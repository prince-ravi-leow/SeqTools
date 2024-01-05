#!/usr/bin/env python3

import random

bases = ["A", "T", "G", "C"]

# Sequence generator function
def gen_seq(seq_length):
	sequence = ""
	for i in range(seq_length):
		sequence += random.choice(bases)
	return sequence

# FASTQ PHRED score
def create_qscore_to_ascii_dict():
	qscore_to_ascii = {}
	offset = 33  # Phred+33 encoding
	for i in range(0, 41):
		qscore_to_ascii[i] = chr(i + offset)
	return qscore_to_ascii

qscore_to_ascii_dict = create_qscore_to_ascii_dict()

def gen_qscore(seq_length, q_min, q_max):
	q_line = ""
	for i in range(seq_length):
		q_line += qscore_to_ascii_dict[random.randint(q_min, q_max)]
	return q_line

# Main function call
def gen_sequences(
		n_seq, 
		len_min, 
		len_max, 
		seq = True, 
		fasta = False, 
		fastq = False,
		q_min = 10, 
		q_max = 30, 
		export = False):
	"""
	*** Randomly generate DNA sequences with user-defined number of sequences (n_seq), within a specified sequence length range (len_min, len_max) ***  

	Output type: str 
	
	'Sequence-only' mode activated by default; set fasta or fastq flags to `True` to output sequences in each respective format. 

	It is possible to specify min and max Q-scores with the arguments q_min and q_max (Default: 10 - 30) 

	Export generated sequences to respective format (plain-text file for sequence-only mode), by setting export flag to `True`
	"""
	
	sequences = ""

	# Generate n_seq number of sequence lengths, bound by len_min and len_max
	seq_lengths = [random.randint(len_min, len_max) for i in range(n_seq)]

	# Iterate over seq_lengths to generate sequences
	for i, seq_length in enumerate(seq_lengths):
		if fasta:
			seq = False # prevent sequence-only mode being invoked
			sequences += f">sequence_{i+1}_length_{seq_length}\n"
			sequences += f"{gen_seq(seq_length)}\n"
		if fastq:
			seq = False # prevent sequence-only mode being invoked
			sequences += f"@sequence_{i+1}_length_{seq_length}\n"
			sequences += f"{gen_seq(seq_length)}\n"
			sequences += "+\n"
			sequences += f"{gen_qscore(seq_length, q_min, q_max)}\n"
		if seq:
			sequences += f"{gen_seq(seq_length)}\n"

	if export:
		if seq:
			file_ext = "txt"
		if fasta:
			file_ext = "fasta"
		if fastq: 
			file_ext = "fastq"
		with open(f"generated_sequences.{file_ext}", "wt") as f:
			f.writelines(sequences)
	
	return sequences

if __name__ == "__main__":
	from argparse import ArgumentParser
	
	parser = ArgumentParser(description = "Generate DNA sequences")
	
	parser.add_argument(
		"-n", 
		action = "store",
		dest = "n_seq", 
		type = int, 
		required = True,
		help = "Number of sequences to generate")
	parser.add_argument(
		"--min", 
		action = "store",
		dest = "len_min", 
		type = int, 
		required = True,
		help = "Min. sequence length")
	parser.add_argument(
		"--max", 
		action = "store",
		dest = "len_max", 
		type = int, 
		required = True,
		help = "Max. sequence length")
	parser.add_argument(
		"--seq", 
		action = "store",
		type = bool,
		dest = "seq", 
		default = True,
		help = "Generate sequences w/out read file elements (default: True)")
	parser.add_argument(
		"--fasta", 
		action = "store",
		type = bool,
		dest = "fasta", 
		default = False,
		help = "Generate fasta files (default: False)")
	parser.add_argument(
		"--fastq", 
		action = "store",
		type = bool,
		dest = "fastq", 
		default = False,
		help = "Generate fastq files (default: False)")
	parser.add_argument(
		"--q_min", 
		action = "store",
		dest = "q_min", 
		type = int, 
		default = 10,
		help = "Min. PHRED score (default: 10)")
	parser.add_argument(
		"--q_max", 
		action = "store",
		dest = "q_max", 
		type = int, 
		default = 30,
		help = "Max. PHRED score (default: 30)")
	parser.add_argument(
		"--export",
		action = "store", 
		dest = "export",
		type = bool,
		default = False,
		help = "Export generated sequences to file"
	)
		
	args = parser.parse_args()
	n_seq = args.n_seq
	len_min = args.len_min
	len_max = args.len_max
	seq = args.seq
	fasta = args.fasta
	fastq = args.fastq
	q_min = args.q_min
	q_max = args.q_max
	export = args.export

	sequences = gen_sequences(
		n_seq, 
		len_min, 
		len_max, 
		seq, 
		fasta, 
		fastq,
		q_min, 
		q_max , 
		export)

	print(sequences)