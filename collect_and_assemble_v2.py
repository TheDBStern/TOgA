#!/usr/bin/env python

import sys, os
from Bio import SeqIO


def run_trinity_first(left_reads, right_reads, target):
	cmd = 'Trinity --max_memory 4G --seqType fq --left %s --right %s --SS_lib_type RF --min_contig_length 100 --CPU 4 --output temp/%s.trinity --full_cleanup' % (left_reads, right_reads, target)
	os.system(cmd)

def run_trinity_next(left_reads, right_reads, target):
	cmd = 'Trinity --max_memory 4G --seqType fq --left %s --right %s --SS_lib_type RF --min_contig_length 200 --CPU 4 --output temp/%s.trinity --full_cleanup' % (left_reads, right_reads, target)
	os.system(cmd)

def run_bowtie2_idx(target_fasta, target_name):
	cmd = 'bowtie2-build -q %s temp/%s' % (target_fasta, target_name)
	os.system(cmd)

def run_bowtie2_first(left_reads, target):
	cmd = 'bowtie2 -p 16 --local --gbar 1 --mp 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x temp/%s -U %s --al temp/mapped_left.fastq -S temp/alignments.sam'  % (target, left_reads)
	os.system(cmd)

def run_bowtie2_next(left_reads, target):
	cmd = 'bowtie2 -p 16 --sensitive-local -x temp/%s -U %s --al temp/mapped_left.fastq -S temp/alignments.sam'  % (target, left_reads)
	os.system(cmd)
	
def best_blast(gene_name, Trinity_fasta, target_fasta):
	make_blast_db = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out temp/blastdb' % (Trinity_fasta)
	os.system(make_blast_db)
	blast_cmd = 'blastn -db temp/blastdb -query %s -out temp/blastout.txt -outfmt 6' % (target_fasta)
	os.system(blast_cmd)
	firstline = open('temp/blastout.txt').readline()
	tophit = firstline.split('\t')[1]
	gene = '_'.join(tophit.split('_')[:-1])
	output = open('temp/%s.besthit.fasta'%gene_name, 'w')
	for seq_record in SeqIO.parse(Trinity_fasta, 'fasta'): #Trinity_fasta
		if gene in seq_record.id:
			SeqIO.write(seq_record, output, 'fasta')

def get_mates(mapped_left, right_idx):
	output = open('temp/right_mates.fastq', 'w')
	left = SeqIO.parse(mapped_left, 'fastq')
	for seq in left:
		right_id = seq.id.replace('/1','/2')
		SeqIO.write(right_idx[right_id], output, 'fastq')
		
	
target_list, left_reads, right_reads, iterations = sys.argv[1:]
print 'Indexing reads...'
right_idx = SeqIO.index(right_reads, 'fastq')
if not os.path.isdir('./temp/'):
	os.mkdir('./temp/')
with open(target_list, 'rU') as t:
	for line in t:
		target_gene = line.split('\t')[0]
		target_fasta = line.split('\t')[1].strip('\n')
		for i in range(0,int(iterations)):
			if i == 0:
				# loose bowtie and trinity
				print '%s Iteration 0'%target_fasta
				run_bowtie2_idx(target_fasta, target_gene)
				print 'Aligning left reads to target'
				run_bowtie2_first(left_reads, target_gene)
				print 'Getting mate pairs'
				get_mates('temp/mapped_left.fastq', right_idx)
				print 'Assembling with Trinity'
				run_trinity_first('temp/mapped_left.fastq', 'temp/right_mates.fastq', target_gene)
			if i == 1:
				# map to all trinity contigs, more strict alignment and get best trinity contig
				print 'Iteration 1'
				run_bowtie2_idx('temp/%s.trinity.Trinity.fasta'%target_gene, target_gene)
				print 'Aligning left reads to new targets'
				run_bowtie2_next(left_reads, target_gene)
				print 'Getting mate pairs'
				get_mates('temp/mapped_left.fastq', right_idx)
				print 'Assembling with Trinity'
				run_trinity_first('temp/mapped_left.fastq', 'temp/right_mates.fastq', target_gene)
				print 'Finding best contigs'
				best_blast(target_gene,'temp/%s.trinity.Trinity.fasta'%target_gene, target_fasta)
			if i > 1:
				# now just use best contigs for mapping
				print 'Iteration %s' % str(i)
				run_bowtie2_idx('temp/%s.besthit.fasta' % (target_gene), target_gene)
				print 'Aligning left reads to new targets'
				run_bowtie2_next(left_reads, target_gene)
				print 'Getting mate pairs'
				get_mates('temp/mapped_left.fastq', right_idx)
				print 'Assembling with Trinity'
				run_trinity_next('temp/mapped_left.fastq', 'temp/right_mates.fastq', target_gene)
				print 'Finding best contigs'
				best_blast(target_gene,'temp/%s.trinity.Trinity.fasta'%target_gene, target_fasta)
		os.rename('temp/%s.trinity.Trinity.fasta' % target_gene, '%s.trinity.Trinity.fasta' % target_gene)
		os.rename('temp/%s.besthit.fasta' % target_gene, '%s.besthit.fasta' % target_gene)
