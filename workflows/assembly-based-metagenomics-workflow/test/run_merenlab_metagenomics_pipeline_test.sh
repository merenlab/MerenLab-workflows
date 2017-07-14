#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

cmd="-np"
# if you want the test to actually run through the pipeline
# then call it like this: bash run_merenlab_metagenomics_pipeline_test.sh sandbox/test-output full
if [ $# -eq 2 ]; then
    if [ $2 == "full" ]; then
        cmd=""
    fi
fi
# copy latest script here
cp ../merenlab-metagenomics-pipeline.snakefile $output_dir
cp -R ../wrappers/ $output_dir
cp ../mock_files_for_merenlab_metagenomics_pipeline/*json $output_dir
cp -R ../mock_files_for_merenlab_metagenomics_pipeline/three_samples_example/ $output_dir
cp ../mock_files_for_merenlab_metagenomics_pipeline/samples.txt $output_dir
cp ../mock_files_for_merenlab_metagenomics_pipeline/samples-no-groups.txt $output_dir

# we have to go into the test directory because snakemake requires you run the command from the directory where the snakemake is
cd $output_dir

INFO "Call snakefile"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd

INFO "Call snakefile with all against all"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config all_against_all='True'

INFO "create samples.txt"
echo -e "sample\tgroup\tr1\tr2" > samples.txt
echo -e "S01\tMYref1\tS01_R1.fastq.gz\tS01_R2.fastq.gz" >> samples.txt
echo -e "S02\totherREF\tS02_R1.fastq.gz\tS02_R2.fastq.gz" >> samples.txt
echo -e "S03\totherREF\tS03_R1.fastq.gz\tS03_R2.fastq.gz" >> samples.txt

INFO "Create a references.txt file"
echo -e "reference\tpath" > references.txt
echo -e "MYref1\tXX1.fa" >> references.txt
echo -e "otherREF\tXX2.fa" >> references.txt

INFO "Creating fake reference fasta files"
touch XX1.fa
touch XX2.fa

INFO "Call snakefile with group list"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          -np \
          --config references_txt='references.txt'

INFO "Call snakefile with group list with all against all"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          -np \
          --config references_txt='references.txt'\
          all_against_all='True'


INFO "create samples.txt with no group column"
echo -e "sample\tr1\tr2" > samples.txt
echo -e "S01\tS01_R1.fastq.gz\tS01_R2.fastq.gz" >> samples.txt
echo -e "S02\tS02_R1.fastq.gz\tS02_R2.fastq.gz" >> samples.txt
echo -e "S03\tS03_R1.fastq.gz\tS03_R2.fastq.gz" >> samples.txt

INFO "Call snakefile with group list"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          --config references_txt='references.txt' -np

INFO "clear all files"
rm *.snakefile config.json samples.txt
rm -rf 00_LOGS 01_QC

# go back to the directory where we started
cd -
