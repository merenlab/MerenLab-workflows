#!/bin/bash
source 00.sh

# Setup #############################
SETUP_WITH_OUTPUT_DIR $1
#####################################

cmd="-pn"
# if you want the test to actually run through the pipeline
# then call it like this: bash run_merenlab_metagenomics_pipeline_test.sh sandbox/test-output full
if [ $# -eq 2 ]; then
    if [ $2 == "full" ]; then
        cmd=""
    fi
fi
# copy latest script here
cp ../merenlab-metagenomics-pipeline.snakefile $output_dir
cp -R ../wrappers/ $output_dir/wrappers/
cp ../mock_files_for_merenlab_metagenomics_pipeline/*json $output_dir
cp -R ../mock_files_for_merenlab_metagenomics_pipeline/three_samples_example/ $output_dir/three_samples_example/
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
          --config all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_ALL_AGAINST_ALL"}'


INFO "create samples.txt"
echo -e "sample\tgroup\tr1\tr2" > samples.txt
echo -e "S01\tG01\tthree_samples_example/S01_R1.fastq.gz\tthree_samples_example/S01_R2.fastq.gz" >> samples.txt
echo -e "S02\tG02\tthree_samples_example/S02_R1.fastq.gz\tthree_samples_example/S02_R2.fastq.gz" >> samples.txt
echo -e "S03\tG02\tthree_samples_example/S03_R1.fastq.gz\tthree_samples_example/S03_R2.fastq.gz" >> samples.txt

INFO "decompress mock reference files"
gzip -d three_samples_example/*.fa.gz 

INFO "Create a references.txt file"
echo -e "reference\tpath" > references.txt
echo -e "G01\tthree_samples_example/G01-contigs.fa" >> references.txt
echo -e "G02\tthree_samples_example/G02-contigs.fa" >> references.txt

INFO "Creating fake reference fasta files"
touch XX1.fa
touch XX2.fa

INFO "Call snakefile with group list"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt' \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE"}'


INFO "Call snakefile with group list with all against all"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt'\
          all_against_all=True \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE_all_against_all"}'



INFO "create samples.txt with no group column"
echo -e "sample\tr1\tr2" > samples.txt
echo -e "S01\tthree_samples_example/S01_R1.fastq.gz\tthree_samples_example/S01_R2.fastq.gz" >> samples.txt
echo -e "S02\tthree_samples_example/S02_R1.fastq.gz\tthree_samples_example/S02_R2.fastq.gz" >> samples.txt
echo -e "S03\tthree_samples_example/S03_R1.fastq.gz\tthree_samples_example/S03_R2.fastq.gz" >> samples.txt

INFO "Call snakefile with no group list in reference mode"
INFO "This one shouldn't do anything and just say 'Nothing to be done.'"
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \
          $cmd \
          --config references_txt='references.txt' \
          output_dirs='{"MERGE_DIR": "06_MERGED_REFERENCE_MODE_all_against_all"}'


# go back to the directory where we started
cd -

INFO "clear all files"
rm -rf $output_dir
