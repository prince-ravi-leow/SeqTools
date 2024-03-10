#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from seqtools_io import SeqIO_load

def build_kmers(sequence:str, k:int):
	kmers = []
	n_kmers = len(sequence) - k + 1
	for i in range(n_kmers):
		kmer = sequence[i:i + k]
		kmers.append(kmer)
	return kmers

def kmer_count(records:SeqRecord, k:int):
	"""
	Iterates over Biopython SeqRecord object, and counts non-unique (duplicate) k-mers.

	Returns nested dictionary with following format:
	results_dict = {'header' : {'kmer': 'count'}}
	"""
	records_dict = {}
	for record in tqdm(records):
		if len(record) <= k: 
			continue
		kmer_dict = {}
		read = str(record.seq)
		for i in range(len(read) - k + 1):
			kmer = read[i:i+k]
			if kmer not in kmer_dict:
				kmer_dict[kmer] = 0
			kmer_dict[kmer] += 1
		for kmer in kmer_dict.keys():
			if kmer_dict[kmer] == 1:
				del kmer_dict[kmer]
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
	records = SeqIO_load(input_file)
	lengths = [len(record) for record in records]
	
	if init_plot:
		ax = sns.histplot(lengths)
		plt.show()
		
	if tick_plot and l is not None and len(l) == 3:    
		fig, ax = plt.subplots()
		sns.histplot(lengths, ax=ax)
		ax.set_xlim(l[0], l[1])
		ax.set_xticks(list(range(l[0], l[1]+1, l[2])))
		plt.show()

def get_phred(records:SeqRecord, plot:bool = False):
	"""
	Requires a loaded Biopython SeqRecord object (may have to be recast as list if you want to iterate on it more than once.
	
	Example:
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