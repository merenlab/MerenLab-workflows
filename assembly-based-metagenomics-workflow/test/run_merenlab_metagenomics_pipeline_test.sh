#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

# copy latest script here
cp ../merenlab-metagenomics-pipeline.snakefile $output_dir
cp ../mock_files_for_merenlab_metagenomics_pipeline/*json $output_dir

# we have to go into the test directory because snakemake requires you run the command from the directory where the snakemake is
cd $output_dir

INFO "create samples.txt"
echo -e "sample\tgroup\tr1\tr2" > samples.txt
echo -e "S01\tG01\tS01_R1.fastq.gz\tS01_R2.fastq.gz" >> samples.txt
echo -e "S02\tG02\tS02_R1.fastq.gz\tS02_R2.fastq.gz" >> samples.txt
echo -e "S03\tG02\tS03_R1.fastq.gz\tS03_R2.fastq.gz" >> samples.txt


INFO "Call snakefile"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --cluster 'clusterize -n {threads} -log {log}' --jobs 4 --latency-wait 100 -np 

INFO "clear all files"
rm *.snakefile config.json cluster.json samples.txt
rm -rf 00_LOGS 01_QC

# go back to the directory where we started
cd -
