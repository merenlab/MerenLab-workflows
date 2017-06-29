#If it doesn't already exist then create a 00_LOGS folder
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

LOGS_DIR = "00_LOGS"
QC_DIR = "01_QC"
ASSEMBLY_DIR = "02_ASSEMBLY"
CONTIGS_DIR = "03_CONTIGS"
MAPPING_DIR = "04_MAPPING"
PROFILE_DIR = "05_ANVIO_PROFILE"
os.makedirs(LOGS_DIR, exist_ok=True)

configfile: "config.json"

# loading the samples.txt file 
samples_txt_file = config["samples_txt"]
fastq_files = u.get_TAB_delimited_file_as_dictionary(samples_txt_file)
SAMPLES = set(fastq_files.keys())


rule all:
    input: expand("{DIR}/{sample}", DIR=PROFILE_DIR, sample=SAMPLES)

rule gen_configs:
    version: 1.0
    input: samples_txt_file
    output: expand("{DIR}/{sample}.ini", DIR = QC_DIR, sample = SAMPLES)
    params: dir=QC_DIR
    shell: "iu-gen-configs {input} -o {params.dir}"

rule qc:
    version: 1.0
    input: QC_DIR + "/{sample}.ini"
    output: 
        r1= QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq.gz", 
        r2= QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq.gz"
    threads: 4
    shell: "iiu-filter-quality-minoche {input}"

rule megahit:
    version: 1.0
    input:
        r1= QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq.gz", 
        r2= QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq.gz"
    params:
        dir = ASSEMBLY_DIR + "/{sample}_TEMP",
        MIN_CONTIG_LENGTH_FOR_ASSEMBLY = config["MIN_CONTIG_LENGTH_FOR_ASSEMBLY"],
        memory_portion_usage_for_assembly = config["memory_portion_usage_for_assembly"]
    output: temp(ASSEMBLY_DIR + "/{sample}_TEMP")
    threads: 11
    shell: "megahit -1 {input.r1} -2 {input.r2} --min-contig-len {params.MIN_CONTIG_LENGTH_FOR_ASSEMBLY} -m {params.memory_portion_usage_for_assembly} -o {params.dir} -t {threads} "

rule reformat_fasta:
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
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: CONTIGS_DIR + "/{sample}-contigs.db"
    threads: 5
    shell: "anvi-gen-contigs-database -f {input} -o {output}"

rule anvi_run_hmms:
    """ Run anvi-run-hmms"""
    version: 1.0
    input: rules.gen_contigs_db.output
    output: touch("anvi_run_hmms-{sample}.done")
    threads: 20
    shell: "anvi-run-hmms -c {input} -T {threads}"

rule bowtie_build:
    """ Run bowtie-build on the contigs fasta"""
    # consider runnig this as a shadow rule
    version: 1.0
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: touch("%s/{sample}/{sample}-contigs" % MAPPING_DIR) #I touch this file because the files created have different suffix
    threads: 4
    shell: "bowtie2-build {input} {output}" 

rule bowtie:
    """ Run mapping with bowtie2,  sort and convert to bam with samtools"""
    version: 1.0
    input:
        build_output = rules.bowtie_build.output,
        r1 = rules.megahit.input.r1,
        r2 = rules.megahit.input.r2,
    output: temp("%s/{sample}/{sample}.sam" % MAPPING_DIR)
    params: dir = MAPPING_DIR + "/{sample}"
    threads: 10
    shell: "bowtie2 --threads {threads} -x {input.build_output} -1 {input.r1} -2 {input.r2} --no-unal -S {output}"

rule samtools_view:
    """ sort sam file with samtools and create a RAW.bam file"""
    version: 1.0
    input: rules.bowtie.output
    output: temp("%s/{sample}/{sample}-RAW.bam" % MAPPING_DIR)
    threads: 4
    shell: "samtools view -F 4 -bS {input} > {output}"

rule anvi_init_bam:
    """ run anvi-init-bam on RAW bam file to create a bam file ready for anvi-profile"""
    version: 1.0 # later we can decide if we want the version to use the version of anvi'o
    input: rules.samtools_view.output
    output: "%s/{sample}/{sample}.bam" % MAPPING_DIR
    threads: 4
    shell: "anvi-init-bam {input} -o {output}"

rule anvi_profile:
    """ run anvi-profile on the bam file"""
    version: 1.0
    input:
        bam = rules.anvi_init_bam.output,
        contigs = rules.gen_contigs_db.output,
        hmms = rules.anvi_run_hmms.output # this is here just so snakemake would run the hmms
    output: "%s/{sample}" % PROFILE_DIR
    params:
        MIN_CONTIG_SIZE_FOR_PROFILE_DB = config["MIN_CONTIG_SIZE_FOR_PROFILE_DB"],
        CLUSTER_CONTIGS = config["CLUSTER_CONTIGS"],
    threads: 5
    shell: "anvi-profile -i {input.bam} -c {input.contigs} -o {output} -M {params.MIN_CONTIG_SIZE_FOR_PROFILE_DB} -T {threads} --cluster-contigs"

