# README #

* Quick summary

This is a simple script for iterative assembly of targeted orthologous sequences with RNA-seq data using Bowtie2 to recruit reads and Trinity for assembly
Requires Bowtie2, Trinity and Blast+ in your PATH as well as BioPython
Requires paired-end reads in separate fastq files with headers ending in /1 and /2
In the first iteration, it maps the left reads to the target using loose bowtie2 parameters. It then retrieves corresponding pairs and performs paired-end assembly with Trinity.
In subsequent iterations, it performs the same steps using the Trinity assemblies as targets
In the final iteration, it finds the best Trinity contig based on blastn bitscore and pulls all 'isoforms' of that 'gene.'


* Version 0.1

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact