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
                -n {threads} -log {log}' --jobs 4 --latency-wait 100 -p 

    Note on rule order: whenever the order of rule execution was ambiguous
        mypreferred approach was to use the rule dependencies. See:
        http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-dependencies

    Note on cluster configuration: because multiple rules require the 
    number of threads as input (for example anvi-profil, megahit), and I
    couldn't find a way to make the number of threads from the
    cluster.config file available within rules, then instead I define the 
    number of threads within each rule. I'm aware it's less elegant than
    having all cluster configuration in the cluster.json file, and would
    love to learn about an alternative solution if you have one.

    Note on log files: in order for the stdout and stderr to be written
    into log files, I have added `&>> {log}` to each shell command. if
    running on a cluster, I suggested including something like this in
    your `--cluster` command: `--log {log}`.
'''
import os
import anvio
import pandas as pd
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
MERGE_DIR = "06_MERGED"

#If it doesn't already exist then create a 00_LOGS folder
os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

# The config file contains many essential configurations for the workflow
configfile: "config.json"

# loading the samples.txt file
samples_txt_file = config["samples_txt"]

# getting the samples information (names, [group], path to r1, path to r2) from samples.txt
samples_information = pd.read_csv(samples_txt_file, sep='\t', index_col=False)
# get a list of the sample names
sample_names = list(samples_information['sample'])

# Collecting information regarding groups.
if "group" in samples_information.columns:
    # if groups were specified then members of a groups will be co-assembled.
    group_names = list(samples_information['group'].unique())
    # creating a dictionary with groups as keys and number of samples in
    # the groups as values
    group_sizes = samples_information['group'].value_counts().to_dict()
else:
    # If not groups were specified then each sample would be assembled 
    # separately
    samples_information['group'] = samples_information['sample']
    group_names = list(sample_names)
    group_sizes = dict.fromkeys(group_names,1)
    

if not os.path.isfile(QC_DIR + "/path-to-raw-fastq-files.txt"):
    # only create the file if it doesn't always exist.
    # if we create the file every time, then it snakemake would run from the begininnig every time
    samples_information.to_csv(QC_DIR + "/path-to-raw-fastq-files.txt", sep='\t', columns=['sample','r1','r2'],index=False)


rule all:
    '''
        The final product of the workflow is an anvi'o merged profile directory
        for each group
    '''
    input: expand("{DIR}/{group}/PROFILE.db", DIR=MERGE_DIR, group=group_names)


rule gen_configs:
    '''
        Generating a config file for each sample. Notice that this step
        is ran only once and generates the config files for all samples
    '''
    version: 1.0
    log: LOGS_DIR + "/gen_configs.log"
    input: ancient(QC_DIR + "/path-to-raw-fastq-files.txt")
    output: expand("{DIR}/{sample}.ini", DIR=QC_DIR, sample=sample_names)
    params: dir=QC_DIR
    shell: "iu-gen-configs {input} -o {params.dir} &>> {log}"


rule qc:
    ''' Run QC using iu-filter-quality-minoche '''
    version: 1.0
    log: LOGS_DIR + "/{sample}-qc.log"
    # making the config file as "ancient" so QC wouldn't run just because
    # a new config file was produced.
    input: ancient(QC_DIR + "/{sample}.ini")
    output: 
        r1 = QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq",
        r2 = QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq"
    threads: 4
    shell: "iu-filter-quality-minoche {input} --ignore-deflines &>> {log}"


rule gzip_fastas:
    ''' Compressing the quality controlled fastq files'''
    version: 1.0
    log: LOGS_DIR + "/{sample}-{R}-gzip.log"
    input: QC_DIR + "/{sample}-QUALITY_PASSED_{R}.fastq"
    output: QC_DIR + "/{sample}-QUALITY_PASSED_{R}.fastq.gz"
    shell: "gzip {input} &>> {log}"


def input_for_megahit(wildcards):
    '''
        Creating a dictionary containing the input files for megahit.
        This could have also been done with a lambda expression, but for
        easier readability it was done in a function.
    '''
    r1 = expand("{DIR}/{sample}-QUALITY_PASSED_R1.fastq.gz", DIR=QC_DIR, sample=list(samples_information[samples_information["group"] == wildcards.group]["sample"]))
    r2 = expand("{DIR}/{sample}-QUALITY_PASSED_R2.fastq.gz", DIR=QC_DIR, sample=list(samples_information[samples_information["group"] == wildcards.group]["sample"]))
    return {'r1': r1, 'r2': r2}


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
    log: LOGS_DIR + "/{group}-megahit.log"
    input: unpack(input_for_megahit)
    params:
        # the minimum length for contig (smaller contigs will be discarded)
        MIN_CONTIG_LENGTH_FOR_ASSEMBLY = config["MIN_CONTIG_LENGTH_FOR_ASSEMBLY"],
        # portion of total memory to use by megahit
        memory_portion_usage_for_assembly = config["memory_portion_usage_for_assembly"]
    # output folder for megahit is temporary
    # TODO: maybe change to shaddow, because with current configuration, if a job is canceled then all
    # the files that were created stay there.
    output: temp(ASSEMBLY_DIR + "/{group}_TEMP")
    threads: 11
    run:
        r1 = ','.join(input.r1)
        r2 = ','.join(input.r2)
        
        cmd = "megahit -1 %s -2 %s" % (r1, r2) + \
            " --min-contig-len {params.MIN_CONTIG_LENGTH_FOR_ASSEMBLY}" + \
            " -m {params.memory_portion_usage_for_assembly}" + \
            " -o {output}" + \
            " -t {threads}" + \
            " &>> {log}"
        print(cmd)
        shell(cmd)


rule touch_megahit_output:
    '''
        Since the output of the megahit rule is a folder (see the comments
        for the rule above), this rule is here to move the final assembly
        fasta to the final assembly folder. This allows later rules to be
        ignorant of the fact that the megahit rule output is a folder.
        This way if in the future we will want to use a different assembler
        that would work well with the snakemake way then we wouldnt have
        to change downstream rules.
    '''
    log: LOGS_DIR + "/{group}-touch_megahit_output.log"
    input:
        dir = ASSEMBLY_DIR + "/{group}_TEMP"
    output:
        contigs = temp(ASSEMBLY_DIR + "/{group}/final.contigs.fa")
    shell:
        "mv {input.dir}/final.contigs.fa {output.contigs}"



rule reformat_fasta:
    '''
        Reformating the headers of the contigs fasta files in order to
        give contigs meaningful names; so that if the group name is
        'MYSAMPLE01', the contigs would look like this:
        > MYSAMPLE01_000000000001
        > MYSAMPLE01_000000000002
    '''
    version: 1.0
    log: LOGS_DIR + "/{group}-reformat_fasta.log"
    input:
        contigs = ASSEMBLY_DIR + "/{group}/final.contigs.fa"
    output:
        contig = protected(ASSEMBLY_DIR + "/{group}/{group}-contigs.fa"),
        report = ASSEMBLY_DIR + "/{group}/{group}-reformat-report.txt"
    params: prefix = "{group}"
    wrapper:
        "wrappers/reformat-fasta"


if config["remove_human_contamination"] == "yes":
    # These rules will only run if the user asked for removal of Human contamination
    rule remove_human_dna_using_centrifuge:
        """ this is just a placeholder for now """
        version: 1.0
        log: LOGS_DIR + "/{group}-remove-human-dna-using-centrifuge.log"
        input: ASSEMBLY_DIR + "/{group}/{group}-contigs.fa"
        output: ASSEMBLY_DIR + "/{group}/{group}-contigs-filtered.fa"
        shell: "touch {output} &>> {log}"


rule gen_contigs_db:
    """ Generates a contigs database using anvi-gen-contigs-database """
    # Setting the version to the same as that of the contigs__version in anvi'o
    version: anvio.__contigs__version__
    log: LOGS_DIR + "/{group}-gen_contigs_db.log"
    # depending on whether human contamination using centrifuge was done
    # or not, the input to this rule will be the raw assembly or the 
    # filtered.
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    output: CONTIGS_DIR + "/{group}-contigs.db"
    threads: 5
    shell: "anvi-gen-contigs-database -f {input} -o {output} &>> {log}"


if config["assign_taxonomy_with_centrifuge"] == "yes":
    # If the user wants taxonomy to be assigned with centrifuge
    # then these following rules would run.
    rule export_gene_calls:
        ''' Export gene calls and use for centrifuge'''
        version: 1.0
        log: LOGS_DIR + "/{group}-export_gene_calls.log"
        # marking the input as ancient in order to ignore timestamps.
        input: ancient(rules.gen_contigs_db.output)
        # output is temporary. No need to keep this file.
        output: temp(CONTIGS_DIR + "/{group}-gene-calls.fa")
        shell: "anvi-get-dna-sequences-for-gene-calls -c {input} -o {output} &>> {log}"


    rule run_centrifuge:
        ''' Run centrifuge on the exported gene calls of the contigs.db'''
        version: 1.0
        log: LOGS_DIR + "/{group}-run_centrifuge.log"
        input: rules.export_gene_calls.output
        output:
            hits = CONTIGS_DIR + "/{group}-centrifuge_hits.tsv",
            report = CONTIGS_DIR + "/{group}-centrifuge_report.tsv"
        params: db=config['centrifuge']['db']
        threads: 5
        shell: "centrifuge -f -x {params.db} {input} -S {output.hits} --report-file {output.report} --threads {threads} &>> {log}"


    rule import_taxonomy:
        ''' Run anvi-import-taxonomy '''
        version: 1.0
        log: LOGS_DIR + "/{group}-import_taxonomy.log"
        input:
            hits = rules.run_centrifuge.output.hits,
            report = rules.run_centrifuge.output.report,
            # marking the contigs.db as ancient in order to ignore timestamps.
            contigs = ancient(rules.gen_contigs_db.output)
        # using a flag file because no file is created by this rule.
        # for more information see:
        # http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#flag-files
        output: touch(CONTIGS_DIR + "/{group}-anvi_import_taxonomy.done")
        params: parser = "centrifuge"
        shell: "anvi-import-taxonomy -c {input.contigs} -i {input.report} {input.hits} -p {params.parser} &>> {log}"


rule anvi_run_hmms:
    """ Run anvi-run-hmms"""
    # TODO: add rule for running hmms for ribosomal genes and import
    # their new gene calls. 
    version: 1.0
    log: LOGS_DIR + "/{group}-anvi_run_hmms.log"
    # if the user requested to run taxonomy using centrifuge, then this
    # will be ran only after centrifuge finished. Otherwise, this rule
    # will run after anvi-gen-contigs-database
    # marking the input as ancient in order to ignore timestamps.
    input: ancient(rules.gen_contigs_db.output)
    # using a snakemake flag file as an output since no file is generated
    # by the rule.
    output: touch(CONTIGS_DIR + "/anvi_run_hmms-{group}.done")
    threads: 20
    shell: "anvi-run-hmms -c {input} -T {threads} &>> {log}"


rule bowtie_build:
    """ Run bowtie-build on the contigs fasta"""
    # TODO: consider runnig this as a shadow rule
    version: 1.0
    log: LOGS_DIR + "/{group}-bowtie_build.log"
    input: rules.remove_human_dna_using_centrifuge.output if config["remove_human_contamination"] == "yes" else rules.reformat_fasta.output.contig
    # I touch this file because the files created have different suffix
    output: touch("%s/{group}/{group}-contigs" % MAPPING_DIR) 
    threads: 4
    shell: "bowtie2-build {input} {output} &>> {log}"


rule bowtie:
    """ Run mapping with bowtie2,  sort and convert to bam with samtools"""
    version: 1.0
    log: LOGS_DIR + "/{group}-{sample}-bowtie.log"
    input:
        build_output = lambda wildcards: expand(MAPPING_DIR + "/{group}/{group}-contigs", group=list(samples_information[samples_information["sample"] == wildcards.sample]["group"])),
        r1 = QC_DIR + "/{sample}-QUALITY_PASSED_R1.fastq.gz",
        r2 = QC_DIR + "/{sample}-QUALITY_PASSED_R2.fastq.gz"
    # setting the output as temp, since we only want to keep the bam file.
    output: temp("%s/{group}/{sample}.sam" % MAPPING_DIR)
    params: dir = MAPPING_DIR + "/{sample}"
    threads: 10
    shell: "bowtie2 --threads {threads} -x {input.build_output} -1 {input.r1} -2 {input.r2} --no-unal -S {output} &>> {log}"


rule samtools_view:
    """ sort sam file with samtools and create a RAW.bam file"""
    version: 1.0
    log: LOGS_DIR + "/{group}-{sample}-samtools_view.log"
    input: rules.bowtie.output
    # output as temp. we only keep the final bam file
    output: temp("%s/{group}/{sample}-RAW.bam" % MAPPING_DIR)
    threads: 4
    shell: "samtools view -F 4 -bS {input} > {output} &>> {log}"


rule anvi_init_bam:
    """
        run anvi-init-bam on RAW bam file to create a bam file ready for
        anvi-profile.
    """
    version: 1.0 # later we can decide if we want the version to use the version of anvi'o
    log: LOGS_DIR + "/{group}-{sample}-anvi_init_bam.log"
    input: rules.samtools_view.output
    output:
        bam = "%s/{group}/{sample}.bam" % MAPPING_DIR,
        bai = "%s/{group}/{sample}.bam.bai" % MAPPING_DIR
    threads: 4
    shell: "anvi-init-bam {input} -o {output.bam} &>> {log}"


rule anvi_profile:
    """ run anvi-profile on the bam file"""
    # setting the rule version to be as the version of the profile database of anvi'o
    version: anvio.__profile__version__
    log: LOGS_DIR + "/{group}-{sample}-anvi_profile.log"
    input:
        bam = "%s/{group}/{sample}.bam" % MAPPING_DIR,
        # TODO: add option to profile all to all (all samples to all contigs)
        # marking the contigs.db as ancient in order to ignore timestamps.
        contigs = ancient(lambda wildcards: CONTIGS_DIR + "/%s-contigs.db" % samples_information[samples_information["sample"] == wildcards.sample]["group"].values[0]),
        # this is here just so snakemake would run the taxonomy before running this rule
        taxonomy = rules.import_taxonomy.output if config["assign_taxonomy_with_centrifuge"] == "yes" else rules.anvi_init_bam.output,
        # this is here just so snakemake would run the hmms before running this rule
        hmms = rules.anvi_run_hmms.output 
    output:
        profile = "%s/{group}/{sample}/PROFILE.db" % PROFILE_DIR,
        aux = PROFILE_DIR + "/{group}/{sample}/AUXILIARY-DATA.h5"
    params:
        # minimal length of contig to include in the profiling
        MIN_CONTIG_SIZE_FOR_PROFILE_DB = config["MIN_CONTIG_SIZE_FOR_PROFILE_DB"],
        # if profiling to individual assembly then clustering contigs
        # see --cluster-contigs in the help manu of anvi-profile
        cluster_contigs = lambda wildcards: '--cluster-contigs' if group_sizes[wildcards.group] == 1 else '',
        name = "{sample}",
        profile_AA = "--profile-AA-frequencies" if config["profile_AA"] == "yes" else "",
        output_dir = PROFILE_DIR + "/{group}/{sample}"
    threads: 5
    shell: "anvi-profile -i {input.bam} -c {input.contigs} -o {params.output_dir} -M {params.MIN_CONTIG_SIZE_FOR_PROFILE_DB} -S {params.name} -T {threads} --overwrite-output-destinations {params.cluster_contigs} {params.profile_AA} &>> {log}"

rule anvi_merge:
    '''
        If there are multiple profiles mapped to the same contigs database,
        then merge these profiles. For individual profile, create a symlink
        to the profile database. The purpose is to have one folder in
        which for every contigs database there is a profile database (or
        a symlink to a profile database) that could be used together for
        anvi-interactive.
    '''
    version: anvio.__profile__version__
    log: LOGS_DIR + "/{group}-anvi_merge.log"
    # The input are all profile databases that belong to the same group
    input:
        # marking the contigs.db as ancient in order to ignore timestamps.
        contigs = ancient(rules.gen_contigs_db.output),
        profiles = lambda wildcards: expand(PROFILE_DIR + "/{group}/{sample}/PROFILE.db", sample=list(samples_information[samples_information['group'] == wildcards.group]['sample']), group=wildcards.group)
    output:
        profile = MERGE_DIR + "/{group}/PROFILE.db",
        aux = MERGE_DIR + "/{group}/AUXILIARY-DATA.h5"
    threads: 1
    params:
        output_dir = MERGE_DIR + "/{group}",
        name = "{group}",
        profile_dir = PROFILE_DIR + "/{group}"
    run:
        # using run instead of shell so we can choose the appropriate shell command.
        # In accordance with: https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency#comment-29348196
        if group_sizes[wildcards.group] == 1:
            # for individual assemblies, create a symlink to the profile database
            shell("ln -s {params.profile_dir}/*/* -t {params.output_dir} &>> {log}")
            shell("touch -h {params.profile_dir}/*/*")
        else:
            shell("anvi-merge {input.profiles} -o {params.output_dir} -c {input.contigs} -S {params.name} --overwrite-output-destinations &>> {log}")

