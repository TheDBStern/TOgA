# README #

* Quick summary

This is a simple script for iterative assembly of targeted orthologous sequences with RNA-seq data using Bowtie2 to recruit reads and Trinity for assembly
Requires Bowtie2, Trinity and Blast+ in your PATH as well as BioPython
In the first iteration, it maps the left reads to the target using loose bowtie2 parameters. It then retrieves corresponding pairs and performs paired-end assembly with Trinity.
In subsequent iterations, it performs the same steps using the Trinity assemblies as targets
In the final iteration, it finds the best Trinity contig based on blastn bitscore and pulls all 'isoforms' of that 'gene.'


* Version 0.1

### Usage ###

* prepare ortholog fasta file with each sequence in the form: >GeneName_SequenceName
python split_multifasta_by_gene.py <reference fasta>

* Run targeted assembly pipeline
python collect_and_assemble_v2.py -h

* Collect best-hits into single file
python concat_rename_besthits.py <species> <directory with output files>