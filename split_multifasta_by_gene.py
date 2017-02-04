#!/usr/bin/env python

'''
Takes a fasta file of target sequences for targeted assembly and creates a directory of
 separate fasta for each target as well as the tab-delimited TargetFile.txt
Assumes that each isoform id is distinguishes by an underscore and a number at the end of the id
requires biopython
'''
import sys, os
from Bio import SeqIO

input, binsize = sys.argv[1:]
if not os.path.isdir('./Targets'):
    os.mkdir('./Targets')
cwd = os.getcwd()
outtext = open('TargetFile.txt', 'w')
handle = open(input, "rU")
for record in SeqIO.parse(handle, "fasta"):
    gene = '_'.join(record.id.split('_')[:-1])
    output = 'Targets/'+gene+'.fasta'
    SeqIO.write(record, output, "fasta")
    outtext.write(gene+'\t'+cwd+'/'+output+'\n')