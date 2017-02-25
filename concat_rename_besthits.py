#!/usr/bin/env python

'''
'''
import sys, os
from Bio import SeqIO

species, DIR = sys.argv[1:]
if DIR[-1] != "/": DIR += "/"
output = open(species, 'w')

for fasta_file in os.listdir(DIR):
	if fasta_file.endswith('.besthit.fasta'):
		count = 0:
		cluster = fasta_file.split('\.')[0]
		with SeqIO.read(DIR+fasta_file, 'fasta') as f:
			for record in f:
				record.id = cluster+'_'+count
				record.description = cluster+'_'+count
				count +=1
				SeqIO.write(record, 'fasta')
