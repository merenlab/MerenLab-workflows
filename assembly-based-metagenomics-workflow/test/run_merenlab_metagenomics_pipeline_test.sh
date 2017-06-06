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

# create samples.txt
echo -e "sample\tr1\tr2" > samples.txt
echo -e "S01\tS01-QUALITY_PASSED_R1.fastq.gz\tS01-QUALITY_PASSED_R2.fastq.gz" >> samples.txt

# create empty fastq files
mkdir 01_QC
touch 01_QC/S01-QUALITY_PASSED_R1.fastq.gz 01_QC/S01-QUALITY_PASSED_R2.fastq.gz

# Call snakefile
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --cluster-config cluster.json --cluster 'clusterize -n {cluster.n} -log {cluster.log}' --jobs 4 --latency-wait 100 -np 

# clear all files
rm *.snakefile config.json cluster.json samples.txt
rm -rf 00_LOGS 01_QC

# go back to the directory where we started
cd -
