#!/usr/bin/env python

import sys, os
import shutil
from Bio import SeqIO


def run_trinity_first(left_reads, right_reads, target):
	cmd = "Trinity --max_memory 4G --seqType fq --left %s --right %s --min_contig_length 100 --CPU 4 --output temp/%s.trinity --full_cleanup" % (left_reads, right_reads, target)
	os.system(cmd)

def run_trinity_next(left_reads, right_reads, target):
	cmd = "Trinity --max_memory 4G --seqType fq --left %s --right %s --SS_lib_type RF --min_contig_length 200 --CPU 4 --output temp/%s.trinity --full_cleanup" % (left_reads, right_reads, target)
	os.system(cmd)

def run_bowtie2_idx(target_fasta, target_name):
	cmd = "bowtie2-build %s temp/%s" % (target_fasta, target_name)
	os.system(cmd)

def run_bowtie2_first(left_reads, right_reads, target):
	cmd = "bowtie2 -p 4 --local --gbar 1 --mp 4 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x temp/%s -1 %s -2 %s --al-conc temp/mapped_paired.fastq --al temp/mapped_singles.fastq -S temp/alignments.sam"  % (target, left_reads, right_reads)
	os.system(cmd)

def run_bowtie2_next(left_reads, right_reads, target):
	cmd = "bowtie2 -p 4 --fast-local -x temp/%s -1 %s -2 %s --al-conc temp/mapped_paired.fastq -S temp/alignments.sam"  % (target, left_reads, right_reads)
	os.system(cmd)
	
def best_blast(gene_name, Trinity_fasta, target_fasta):
	make_blast_db = "makeblastdb -in %s -dbtype nucl -parse_seqids -out temp/blastdb" % (Trinity_fasta)
	os.system(make_blast_db)
	blast_cmd = "blastn -db temp/blastdb -query %s -out temp/blastout.txt -outfmt 6" % (target_fasta)
	os.system(blast_cmd)
	firstline = open('temp/blastout.txt').readline()
	tophit = firstline.split('\t')[1]
	gene = '_'.join(tophit.split('_')[:-1])
	output = open('%s.besthit.fasta'%gene_name, 'w')
	for seq_record in SeqIO.parse(Trinity_fasta, "fasta"): #Trinity_fasta
		if gene in seq_record.id:
			print seq_record.id
			SeqIO.write(seq_record, output, 'fasta')
		
	
target_list, left_reads, right_reads, iterations = sys.argv[1:]
targets = open(target_list, 'rU')
if not os.path.isdir('./temp/'):
    os.mkdir('./temp/')
for line in targets:
	gene = line.split('\t')[0]
	fasta = line.split('\t')[1]	
	for i in range(0,int(iterations)):
		if i == 0:
			# loose bowtie and trinity
			run_bowtie2_idx(fasta, gene)
			run_bowtie2_first(left_reads, right_reads, gene)
			cmd = "cat temp/mapped_paired.1.fastq temp/mapped_singles0.fastq > temp/P1wU.fastq"
			os.system(cmd)
			run_trinity_first("temp/P1wU.fastq", "temp/mapped_paired.2.fastq", gene)
		elif i != 0 and i != int(iterations)-1:
			# more stringent bowtie and trinity
			run_bowtie2_idx("temp/%s.trinity.Trinity.fasta" % (gene), gene)
			run_bowtie2_next(left_reads, right_reads, gene)
			run_trinity_next("temp/mapped_paired.1.fastq", "temp/mapped_paired.2.fastq", gene)
		else:
			#final iteration with best blast
			run_bowtie2_idx("temp/%s.trinity.Trinity.fasta" % (gene), gene)
			run_bowtie2_next(left_reads, right_reads, gene)
			run_trinity_next("temp/mapped_paired.1.fastq", "temp/mapped_paired.2.fastq", gene)
			best_blast(gene,"temp/%s.trinity.Trinity.fasta"%gene, fasta)

	os.rename("temp/%s.trinity.Trinity.fasta"%gene, "%s.trinity.Trinity.fasta"%gene)
	os.rename("temp/%s.besthit.fasta"%gene, "%s.besthit.fasta"%gene)


shutil.rmtree('temp')
