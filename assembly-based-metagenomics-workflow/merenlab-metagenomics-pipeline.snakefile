'''
    This is a snakemake for the metagenomics workflow in the Meren Lab using
    anvi'o.

    It includes the following steps:
    Quality filtering
    Assembly using megahit
    Mapping of metagenomes to assemblies using bowtie2
    generating anvio contigs database (including running hmm profile)
    generating anvio profile database

    The following files must exist in the working directory:
    config.json - this file contains essential configuration information for
    the pipeline. Here is an example of the format of the file:

    {
        "samples_txt": "samples.txt",
        "remove_human_contamination": "no",
        "memory_portion_usage_for_assembly": "0.4",
        "MIN_CONTIG_LENGTH_FOR_ASSEMBLY": "1000",
        "MIN_CONTIG_SIZE_FOR_PROFILE_DB": "2500",
        "CLUSTER_CONTIGS": "--cluster-contigs"
    }

    samples.txt - 
        TAB-delimited file to describe where samples are. The
        header line should be "sample", "r1", and "r2". Each
        row should list the sample name in the first column,
        and full path for r1 and r2.



    An example run of this workflow on the barhal server:
    $ snakemake --snakefile merenlab-metagenomics-pipeline.snakefile \ 
                --cluster-config cluster.json --cluster 'clusterize  \
                -n {threads} -log {cluster.log}' --jobs 4 --latency-wait 100 -p 

    Note on rule order: whenever the order of rule execution was ambiguous
        mypreferred approach was to use the rule dependencies. See:
        http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-dependencies
    '''
import os
import anvio
import anvio.utils as u


__author__ = "Alon Shaiber"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"

# Setting the names of all directories
LOGS_DIR = "00_LOGS"
QC_DIR = "01_QC"
ASSEMBLY_DIR = "02_ASSEMBLY"
CONTIGS_DIR = "03_CONTIGS"
MAPPING_DIR = "04_MAPPING"
PROFILE_DIR = "05_ANVIO_PROFILE"

#If it doesn't already exist then create a 00_LOGS folder
os.makedirs(LOGS_DIR, exist_ok=True)

# The config file contains many essential configurations for the workflow
configfile: "config.json"

# loading the samples.txt file
samples_txt_file = config["samples_txt"]

# getting the names of samples from samples.txt
fastq_files = u.get_TAB_delimited_file_as_dictionary(samples_txt_file)
SAMPLES = set(fastq_files.keys())


rule all:
    '''
        The final product of the workflow is an anvi'o profile directory
        for each sample
    '''
    input: expand("{DIR}/{sample}", DIR=PROFILE_DIR, sample=SAMPLES)


rule gen_configs:
    '''
        Generating a config file for each sample. Notice that this step
        is ran only once and generates the config files for all samples
    '''
    version: 1.0
    input: samples_txt_file
    output: expand("{DIR}/{sample}.ini", DIR = QC_DIR, sample = SAMPLES)
    params: dir=QC_DIR
    shell: "iu-gen-configs {input} -o {params.dir}"


rule qc:
    ''' Run QC using iu-filter-quality-minoche '''
    version: 1.0
    input: QC_DIR + "/{sample}.ini"
    output: 
        r1= QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq", 
        r2= QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq"
    threads: 4
    shell: "iu-filter-quality-minoche {input} --ignore-deflines"


rule gzip_fastas:
    ''' Compressing the quality controlled fastq files'''
    version: 1.0
    input: QC_DIR + "/{sample}-QUALITY_PASSED_{R}.fastq"
    output: QC_DIR + "/{sample}-QUALITY_PASSED_{R}.fastq.gz"
    shell: "gzip {input}"

rule megahit:
    ''' 
        Assembling fastq files using megahit.
        Notice that megahit requires a directory to be specified as 
        output. If the directory already exists then megahit will not
        run. To avoid this, the output for this rule is defined as the 
        directory (and not the assembly fasta file), because if the 
        fasta file was defined as the output of the rule, then snakemake
        would automaticaly creates the directory.
        All files created by megahit are stored in a temporary folder,
        and only the fasta file is kept for later analysis.
    '''
    version: 1.0
    input:
        r1= QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq.gz", 
        r2= QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq.gz"
    params:
        # the minimum length for contig (smaller contigs will be discarded)
        MIN_CONTIG_LENGTH_FOR_ASSEMBLY = config["MIN_CONTIG_LENGTH_FOR_ASSEMBLY"],
        # portion of total memory to use by megahit
        memory_portion_usage_for_assembly = config["memory_portion_usage_for_assembly"]
    # output folder for megahit is temporary
    output: temp(ASSEMBLY_DIR + "/{sample}_TEMP")
    threads: 11
    shell: "megahit -1 {input.r1} -2 {input.r2} --min-contig-len {params.MIN_CONTIG_LENGTH_FOR_ASSEMBLY} -m {params.memory_portion_usage_for_assembly} -o {output} -t {threads} "

rule reformat_fasta:
    '''
        Reformating the headers of the contigs fasta files in order to
        give contigs meaningful names; so that if the sample name is
        'MYSAMPLE01', the contigs would look like this:
        > MYSAMPLE01_000000000001
        > MYSAMPLE01_000000000002
    '''
    version: 1.0
    input:
        ASSEMBLY_DIR + "/{sample}_TEMP"
    output:
        contig = protected(ASSEMBLY_DIR + "/{sample}/{sample}-contigs.fa"),
        report = ASSEMBLY_DIR + "/{sample}/{sample}-reformat-report.txt"
    shell: "anvi-script-reformat-fasta {input}/final.contigs.fa -o {output.contig} -r {output.report} --simplify-names --prefix {wildcards.sample}"

if config["remove_human_contamination"] == "yes":
    # These rules will only run if the user asked for removal of Human contamination
    rule remove_human_dna_using_centrifuge:
        """ this is just a placeholder for now """
        version: 1.0
        input: ASSEMBLY_DIR + "/{sample}/{sample}-contigs.fa"
        output: ASSEMBLY_DIR + "/{sample}/{sample}-contigs-filtered.fa"
        shell: "touch {output}"

rule gen_contigs_db:
    """ Generates a contigs database using anvi-gen-contigs-database """
    version: 1.0
    # depending on whether human contamination using centrifuge was done
    # or not, the input to this rule will be the raw assembly or the 
    # filtered.
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: CONTIGS_DIR + "/{sample}-contigs.db"
    threads: 5
    shell: "anvi-gen-contigs-database -f {input} -o {output}"

rule anvi_run_hmms:
    """ Run anvi-run-hmms"""
    version: 1.0
    input: rules.gen_contigs_db.output
    # Since this rule doesn't create a new file, then snakemake will
    # touch this file at the end to show that this rule executed.
    output: touch("anvi_run_hmms-{sample}.done")
    threads: 20
    shell: "anvi-run-hmms -c {input} -T {threads}"

rule bowtie_build:
    """ Run bowtie-build on the contigs fasta"""
    # TODO: consider runnig this as a shadow rule
    version: 1.0
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    # I touch this file because the files created have different suffix
    output: touch("%s/{sample}/{sample}-contigs" % MAPPING_DIR) 
    threads: 4
    shell: "bowtie2-build {input} {output}"

rule bowtie:
    """ Run mapping with bowtie2,  sort and convert to bam with samtools"""
    version: 1.0
    input:
        build_output = rules.bowtie_build.output,
        r1 = rules.megahit.input.r1,
        r2 = rules.megahit.input.r2,
    # setting the output as temp, since we only want to keep the bam file.
    output: temp("%s/{sample}/{sample}.sam" % MAPPING_DIR)
    params: dir = MAPPING_DIR + "/{sample}"
    threads: 10
    shell: "bowtie2 --threads {threads} -x {input.build_output} -1 {input.r1} -2 {input.r2} --no-unal -S {output}"

rule samtools_view:
    """ sort sam file with samtools and create a RAW.bam file"""
    version: 1.0
    input: rules.bowtie.output
    # output as temp. we only keep the final bam file
    output: temp("%s/{sample}/{sample}-RAW.bam" % MAPPING_DIR)
    threads: 4
    shell: "samtools view -F 4 -bS {input} > {output}"

rule anvi_init_bam:
    """
        run anvi-init-bam on RAW bam file to create a bam file ready for
        anvi-profile.
    """
    version: 1.0 # later we can decide if we want the version to use the version of anvi'o
    input: rules.samtools_view.output
    output:
        bam = "%s/{sample}/{sample}.bam" % MAPPING_DIR,
        bai = "%s/{sample}/{sample}.bam.bai" % MAPPING_DIR
    threads: 4
    shell: "anvi-init-bam {input} -o {output.bam}"

rule anvi_profile:
    """ run anvi-profile on the bam file"""
    version: 1.0
    input:
        bam = rules.anvi_init_bam.output.bam,
        contigs = rules.gen_contigs_db.output,
        # this is here just so snakemake would run the hmms before running this rule
        hmms = rules.anvi_run_hmms.output 
    output: "%s/{sample}" % PROFILE_DIR
    params:
        # minimal length of contig to include in the profiling
        MIN_CONTIG_SIZE_FOR_PROFILE_DB = config["MIN_CONTIG_SIZE_FOR_PROFILE_DB"],
        # see --cluster-contigs in the help manu of anvi-profile
        CLUSTER_CONTIGS = config["CLUSTER_CONTIGS"],
    threads: 5
    shell: "anvi-profile -i {input.bam} -c {input.contigs} -o {output} -M {params.MIN_CONTIG_SIZE_FOR_PROFILE_DB} -T {threads} --cluster-contigs"

