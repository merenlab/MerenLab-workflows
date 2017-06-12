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


os.makedirs('00_LOGS', exist_ok=True)
configfile: "config.json"

# loading the samples.txt file 
samples_txt_file = config["samples_txt"]
fastq_files = u.get_TAB_delimited_file_as_dictionary(samples_txt_file)
SAMPLES = set(fastq_files.keys())


rule all:
    input: expand("03_CONTIGS/{sample}-contigs.db", sample=SAMPLES)

rule megahit:
    version: 1.0
    input:
        r1="01_QC/{sample}-QUALITY_PASSED_R1.fastq.gz", 
        r2="01_QC/{sample}-QUALITY_PASSED_R2.fastq.gz"
    threads: 11
    params: dir="02_ASSEMBLY/{sample}_TEMP"
    output: temp("02_ASSEMBLY/{sample}_TEMP")
    shell: "megahit -1 {input.r1} -2 {input.r2} --min-contig-len 1000 -m 0.4 -o {params.dir} -t {threads} "

rule reformat_fasta:
    version: 1.0
    input:
        "02_ASSEMBLY/{sample}_TEMP"
    output:
        contig=protected("02_ASSEMBLY/{sample}/{sample}-contigs.fa"),
        report="02_ASSEMBLY/{sample}/{sample}-reformat-report.txt"
    shell: "anvi-script-reformat-fasta {input}/final.contigs.fa -o {output.contig} -r {output.report} --simplify-names --prefix {wildcards.sample}"

if config["remove_human_contamination"] == "yes":
    # These rules will only run if the user asked for removal of Human contamination
    rule remove_human_dna_using_centrifuge:
        """ this is just a placeholder for now """
        version: 1.0
        input: "02_ASSEMBLY/{sample}/{sample}-contigs.fa"
        output: "02_ASSEMBLY/{sample}/{sample}-contigs-filtered.fa"
        shell: "touch {output}"

rule gen_contigs_db:
    """ Generates a contigs database using anvi-gen-contigs-database """
    version: 1.0
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: "03_CONTIGS/{sample}-contigs.db"
    shell: "anvi-gen-contigs-database -f {input} -o {output}"

rule bowtie_build:
    """ Run bowtie-build on the contigs fasta"""
    # consider runnig this as a shadow rule
    version: 1.0
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: "04_MAPPING/{sample}-contigs"
    shell: "bowtie-build {input} {output}" 

rule bowtie:
    """ Run mapping with bowtie2,  sort and convert to bam with samtools"""
    version: 1.0
    input: rules.bowtie_build.output
    output: {sample}.bam
    params: 
        threads = {cluster.n},
        r1 = rules.megahit.input.r1
        r2 = rules.megahit.input.r2
        sam = "{sample}.sam"
        raw_bam = "{sample}-RAW.bam"
        dir = "04_MAPPING/{sample}"
    shadow: "shallow" # By making this rule a shadow rule, we don't need to cleanup the sam file and the raw.bam file.
    shell:
    """
    bowtie2 --threads {params.threads} -x {input} -1 {params.r1} -2 {params.r2} --no-unal -S {params.dir}/{params.sam}
    samtools view -F 4 -bS {params.dir}/{params.sam} > {params.dir}/{params.raw_bam}
    """

rule samtools_view:
    """ sort sam file with samtools"""

