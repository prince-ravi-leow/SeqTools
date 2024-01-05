#!/usr/bin/env python3

import os
import gzip
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

def SeqIO_load(input_file:str):	
	"""
	Usage: 
	records = SeqIO_load("reads.fastq.gz")
	"""	

	def detect_filetype(ext:str):
		try:
			if ext in ['fa', 'FA', 'fasta']:
				file_type = 'fasta'
			if ext in ['fq', 'FQ', 'fastq']:
				file_type = 'fastq'
		except:
			raise TypeError("File type could not be determined - please ensure file extension reflects either fasta or fastq reads file.")
	
		return file_type	

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

def build_kmers(sequence:str, k:int):
	kmers = []
	n_kmers = len(sequence) - k + 1
	for i in range(n_kmers):
		kmer = sequence[i:i + k]
		kmers.append(kmer)
	return kmers

def kmer_count(k:int, records):
	"""
	Iterates over Biopython SeqRecord object, and counts non-unique k-mers.

	Returns nested dictionary with following format:
	results_dict = {'header' : {'kmer': 'count'}}
	
	Example:
	from Bio import SeqIO
	records = SeqIO.parse("reads.fastq", "fastq")
	records_dict = kmer_count(k=100, records=records)
	"""
	records_dict = {}
	for record in tqdm(records):
		# Only process reads as or longer than kmer length
		if len(record) <= k: 
			continue
		kmer_dict = {}
		read = str(record.seq)
		n_kmers = len(read) - k + 1
		for i in range(n_kmers):
			kmer = str(record.seq)[i:i+k]
			# If kmer not in dict, add and initiate count
			if kmer not in kmer_dict:
				kmer_dict[kmer] = 0
			kmer_dict[kmer] += 1
		# Discard kmers that occur only once
		for kmer in list(kmer_dict.keys()):
			if kmer_dict[kmer] == 1:
				del kmer_dict[kmer]
		# Discard empty sub dicts, add populated sub dicts to main dict, with key = current header/read
		if kmer_dict:
			records_dict[record.id] = kmer_dict
	return records_dict

def length_histplot(input_file:str, init_plot:bool=False, tick_plot:bool=False, l: list = None) -> None:
	"""
	*** Vanilla plot: ***
	read_length_histplot("reads.fastq", init_plot=True)
	
	*** Plot with adjusted x-axis and steps (l = [start, stop, step]): ***
	read_length_histplot("reads.fastq", tick_plot=True, l=[0, 10000, 1000])
	"""     
	# Open reads file, calculate lengths
	records = SeqIO_load(input_file)
	
	lengths = [len(record) for record in records]
	
	# Make general naive hist plot 
	if init_plot:
		ax = sns.histplot(lengths)
		plt.show()
	# Hist plot with user-defined xlim and steps
	if tick_plot and l is not None and len(l) == 3:    
		fig, ax = plt.subplots()
		sns.histplot(lengths, ax=ax)
		ax.set_xlim(l[0], l[1])
		ax.set_xticks(list(range(l[0], l[1]+1, l[2])))
		plt.show()

def get_phred(records, plot = False):
	"""
	Requires a loaded Biopython SeqRecord object (may have to be recast as list if you want to iterate on it more than once.
	
	Example:
	from Bio import SeqIO

	records = SeqIO.parse("reads.fastq", "fastq")
					*** OR ***
	records = list(SeqIO.parse("reads.fastq", "fastq"))

	get_phred(records)
			
	
	Set plot = True to get histogram plot of phred scores:
	get_phred(records, plot=True)
	"""
	
	quality = [record.letter_annotations["phred_quality"] for record in records]
	quality_flat = [item for sublist in quality for item in sublist]
	quality_array = np.array(quality_flat)
	
	if plot:
		counts, bins, _ = plt.hist(quality_array, bins=range(0, 41))
		plt.show()
	
	return quality_array