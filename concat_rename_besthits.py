#!/usr/bin/env python

'''
'''
import sys, os
from Bio import SeqIO

species, DIR = sys.argv[1:]
if DIR[-1] != "/": DIR += "/"
output = open(species+'.fasta', 'w')
f_count = 0
t_count = 0
for fasta_file in os.listdir(DIR):
	if fasta_file.endswith('.besthit.fasta'):
		print 'Adding '+fasta_file
		f_count +=1
		counter = 0
		cluster = fasta_file.split('.')[0]
		fasta = SeqIO.parse(DIR+fasta_file, 'fasta')
		for record in fasta:
			t_count +=1
			record.id = cluster+'_'+str(counter)
			record.description = cluster+'_'+str(counter)
			counter +=1
			SeqIO.write(record, output, 'fasta')

print 'Processed %d assemblies with a total of %d transcripts' % (f_count, t_count)