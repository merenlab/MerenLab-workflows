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
        "memory": "0.4",
        "min_contig_len": "1000",
        "min_contig_length": "2500",
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
import anvio.dbops as dbops
import anvio.terminal as terminal

from anvio.errors import ConfigError

__author__ = "Alon Shaiber"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Alon Shaiber"
__email__ = "alon.shaiber@gmail.com"

run = terminal.Run()
progress = terminal.Progress()

# The config file contains many essential configurations for the workflow
configfile: "config.json"

# Setting the names of all directories
dir_list = ["LOGS_DIR", "QC_DIR", "ASSEMBLY_DIR", "CONTIGS_DIR", "MAPPING_DIR", "PROFILE_DIR", "MERGE_DIR"]
dir_names = ["00_LOGS", "01_QC", "02_ASSEMBLY", "03_CONTIGS", "04_MAPPING", "05_ANVIO_PROFILE", "06_MERGED"]
dirs_dict = dict(zip(dir_list, dir_names))






def A(_list, d, default_value = ""):
    '''
        A helper function to make sense of config details.
        string_list is a list of strings (or a single string)
        d is a dictionary

        this function checks if the strings in x are nested values in y.
        For example if x = ['a','b','c'] then this function checkes if the
        value y['a']['b']['c'] exists, if it does then it is returned
    '''
    # converting to list for the cases of only one item
    if type(_list) is not list:
        _list = [_list]
    while _list:
        a = _list.pop(0)
        if a in d:
            d = d[a]
        else:
            return default_value
    return d

# a helper function to get the user defined number of threads for a rule
def T(rule_name, N=1): return A([rule_name,'threads'], config, default_value=N)

for d in A("output_dirs", config):
    # renaming folders according to the config file, if the user specified.
    if d not in dir_list:
        # making sure the user is asking to rename an existing folder.
        raise ConfigError("You define a name for the directory '%s' in your "\
                          "config file, but the only available folders are: "\
                          "%s" % (d, dir_list))

    dirs_dict[d] = A(d,config["output_dirs"])

# setting configuration for optional steps
run_remove_human_dna_using_centrifuge = A(["remove_human_dna_using_centrifuge", "run"], config)
# default is NOT running taxonomy with centrifuge
run_taxonomy_with_centrifuge = A(["run_centrifuge", "run"], config)
# default is running anvi_run_hmms
run_anvi_run_hmms = A(["anvi_run_hmms", "run"], config, default_value=True)
# Sanity check for centrifuge db
if run_taxonomy_with_centrifuge:
    if not A(["run_centrifuge", "db"], config):
        raise ConfigError("If you plan to run centrifuge, then you must "\
                          "provide a path for the centrifuge db in the "\
                          "config file. See documentation for more details.")


# loading the samples.txt file
# The default samples file is samples.txt
samples_txt_file = A("samples_txt", config, default_value="samples.txt")
# getting the samples information (names, [group], path to r1, path to r2) from samples.txt
samples_information = pd.read_csv(samples_txt_file, sep='\t', index_col=False)
# get a list of the sample names
sample_names = list(samples_information['sample'])

# if no groups are supplied then group names are sample names
group_names = sample_names

if 'references_txt' in config:
    # if the user supplied a reference.txt file, then there is no need to
    # create an assembly (see documentation for 'reference-mode')
    references_txt_file = config["references_txt"]
    references_information = pd.read_csv(references_txt_file, sep='\t', index_col=0).to_dict(orient='index')
    group_names = list(references_information.keys())

# Collecting information regarding groups.
if "group" in samples_information.columns:
    # if groups were specified then members of a groups will be co-assembled.
    group_names = list(samples_information['group'].unique())
    # creating a dictionary with groups as keys and number of samples in
    # the groups as values
    group_sizes = samples_information['group'].value_counts().to_dict()

    if 'references_txt' in config:
        # sanity check to see that groups specified in samples.txt match
        # the names of references.
        mismatch = set(group_names) - set(references_information.keys())
        if mismatch:
            raise ConfigError("Group names specified in the samples.txt \
                               file must match the names of references \
                               in the reference.txt file. These are the \
                               mismatches: %s" % mismatch)

else:
    if 'references_txt' in config:
        # if the user didn't provide a group column in the samples.txt,
        # in reference mode the default is 'all_against_all'.
        config['all_against_all'] = True
    else:
        # if not groups were specified then each sample would be assembled
        # separately
        samples_information['group'] = samples_information['sample']
        group_names = list(sample_names)
        group_sizes = dict.fromkeys(group_names,1)

if A('all_against_all', config) :
    # in all_against_all, the size of each group is as big as the number
    # of samples.
    group_sizes = dict.fromkeys(group_names,len(sample_names))


rule all:
    '''
        The final product of the workflow is an anvi'o merged profile directory
        for each group
    '''
    input: expand("{DIR}/{group}/PROFILE.db", DIR=dirs_dict["MERGE_DIR"], group=group_names)


rule gen_input_for_gen_configs:
    log: dirs_dict["LOGS_DIR"] + "/gen_input_for_gen_configs"
    output: dirs_dict["QC_DIR"] + "/path-to-raw-fastq-files.txt"
    threads: T('gen_input_for_gen_configs')
    resources: nodes = T('gen_input_for_gen_configs'),
    run:
        samples_information.to_csv(output, sep='\t', columns=['sample','r1','r2'],index=False)


rule gen_configs:
    '''
        Generating a config file for each sample. Notice that this step
        is ran only once and generates the config files for all samples
    '''
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/gen_configs.log"
    # the input file is marked as 'ancient' so snakemake wouldn't run it
    # just because a new path-to-raw-fastq-files.txt file was created.
    input: ancient(dirs_dict["QC_DIR"] + "/path-to-raw-fastq-files.txt")
    output: temp(expand("{DIR}/{sample}.ini", DIR=dirs_dict["QC_DIR"], sample=sample_names))
    params: dir=dirs_dict["QC_DIR"]
    threads: T('gen_configs')
    resources: nodes = T('gen_configs'),
    shell: "iu-gen-configs {input} -o {params.dir} &>> {log}"


rule qc:
    ''' Run QC using iu-filter-quality-minoche '''
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{sample}-qc.log"
    # making the config file as "ancient" so QC wouldn't run just because
    # a new config file was produced.
    input: ancient(dirs_dict["QC_DIR"] + "/{sample}.ini")
    output:
        r1 = dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_R1.fastq",
        r2 = dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_R2.fastq"
    threads: T('qc', 4)
    resources: nodes = T('qc', 4),
    shell: "iu-filter-quality-minoche {input} --ignore-deflines &>> {log}"


rule gzip_fastas:
    ''' Compressing the quality controlled fastq files'''
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{sample}-{R}-gzip.log"
    input: dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_{R}.fastq"
    output: dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_{R}.fastq.gz"
    threads: T('gzip_fastas')
    resources: nodes = T('gzip_fastas'),
    shell: "gzip {input} &>> {log}"


def input_for_megahit(wildcards):
    '''
        Creating a dictionary containing the input files for megahit.
        This could have also been done with a lambda expression, but for
        easier readability it was done in a function.
    '''
    r1 = expand("{DIR}/{sample}-QUALITY_PASSED_R1.fastq.gz", DIR=dirs_dict["QC_DIR"], sample=list(samples_information[samples_information["group"] == wildcards.group]["sample"]))
    r2 = expand("{DIR}/{sample}-QUALITY_PASSED_R2.fastq.gz", DIR=dirs_dict["QC_DIR"], sample=list(samples_information[samples_information["group"] == wildcards.group]["sample"]))
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
    log: dirs_dict["LOGS_DIR"] + "/{group}-megahit.log"
    input: unpack(input_for_megahit)
    params:
        # the minimum length for contig (smaller contigs will be discarded)
        min_contig_len = int(A(["megahit", "min_contig_len"], config, default_value="1000")),
        # portion of total memory to use by megahit
        memory = float(A(["megahit", "memory"], config, default_value=0.4))
    # output folder for megahit is temporary (using the snakemake temp())
    # TODO: maybe change to shaddow, because with current configuration, if a job is canceled then all
    # the files that were created stay there.
    output: temp(dirs_dict["ASSEMBLY_DIR"] + "/{group}_TEMP")
    threads: T('megahit', 11)
    resources: nodes = T('megahit', 11),
    run:
        r1 = ','.join(input.r1)
        r2 = ','.join(input.r2)

        cmd = "megahit -1 %s -2 %s" % (r1, r2) + \
            " --min-contig-len {params.min_contig_len}" + \
            " -m {params.memory}" + \
            " -o {output}" + \
            " -t {threads}" + \
            " &>> {log}"
        print("Running: %s" % cmd)
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
    log: dirs_dict["LOGS_DIR"] + "/{group}-touch_megahit_output.log"
    input:
        dir = dirs_dict["ASSEMBLY_DIR"] + "/{group}_TEMP"
    output:
        contigs = temp(dirs_dict["ASSEMBLY_DIR"] + "/{group}/final.contigs.fa")
    threads: T('touch_megahit_output')
    resources: nodes = T('touch_megahit_output'),
    shell:
        "mv {input.dir}/final.contigs.fa {output.contigs}"


def input_for_reformant_fasta(wildcards):
    '''define the input for the rule reformat fasta.'''

    if 'references_txt' in config:
        # in 'reference mode' the input is the reference fasta
        contig = references_information[wildcards.group]['path']
    else:
        contig = dirs_dict["ASSEMBLY_DIR"] + "/%s/final.contigs.fa" % wildcards.group

    return contig


rule reformat_fasta:
    '''
        Reformating the headers of the contigs fasta files in order to
        give contigs meaningful names; so that if the group name is
        'MYSAMPLE01', the contigs would look like this:
        > MYSAMPLE01_000000000001
        > MYSAMPLE01_000000000002
    '''
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-reformat_fasta.log"
    input:
        contigs = input_for_reformant_fasta
    output:
        # write protecting the contig fasta file using protected() because
        # runnig the assembly is probably the most time consuming step and
        # we don't want anyone accidentaly deleting or changing this file.
        contig = protected(dirs_dict["ASSEMBLY_DIR"] + "/{group}/{group}-contigs.fa"),
        report = dirs_dict["ASSEMBLY_DIR"] + "/{group}/{group}-reformat-report.txt"
    params: prefix = "{group}"
    threads: T('reformat_fasta')
    resources: nodes = T('reformat_fasta'),
    wrapper:
        # Notice that path to wrapper is relative to the workdir (if you
        # want an absolute path, use 'file://' instead of 'file:')
        "file:wrappers/reformat-fasta"


if run_remove_human_dna_using_centrifuge:
    # These rules will only run if the user asked for removal of Human contamination
    rule remove_human_dna_using_centrifuge:
        """ this is just a placeholder for now """
        version: 1.0
        log: dirs_dict["LOGS_DIR"] + "/{group}-remove-human-dna-using-centrifuge.log"
        input: dirs_dict["ASSEMBLY_DIR"] + "/{group}/{group}-contigs.fa"
        output: dirs_dict["ASSEMBLY_DIR"] + "/{group}/{group}-contigs-filtered.fa"
        threads: T('remove_human_dna_using_centrifuge')
        resources: nodes = T('remove_human_dna_using_centrifuge'),
        shell: "touch {output} &>> {log}"


rule gen_contigs_db:
    """ Generates a contigs database using anvi-gen-contigs-database """
    # Setting the version to the same as that of the contigs__version in anvi'o
    version: anvio.__contigs__version__
    log: dirs_dict["LOGS_DIR"] + "/{group}-gen_contigs_db.log"
    # depending on whether human contamination using centrifuge was done
    # or not, the input to this rule will be the raw assembly or the
    # filtered.
    input: rules.remove_human_dna_using_centrifuge.output if run_remove_human_dna_using_centrifuge else rules.reformat_fasta.output.contig
    output:
        db = dirs_dict["CONTIGS_DIR"] + "/{group}-contigs.db",
        aux = dirs_dict["CONTIGS_DIR"] + "/{group}-contigs.h5"
    threads: T('gen_contigs_db', 5)
    resources: nodes = T('gen_contigs_db', 5),
    shell: "anvi-gen-contigs-database -f {input} -o {output.db} &>> {log}"


if run_taxonomy_with_centrifuge:
    # If the user wants taxonomy to be assigned with centrifuge
    # then these following rules would run.
    rule export_gene_calls:
        ''' Export gene calls and use for centrifuge'''
        version: 1.0
        log: dirs_dict["LOGS_DIR"] + "/{group}-export_gene_calls.log"
        # marking the input as ancient in order to ignore timestamps.
        input: ancient(rules.gen_contigs_db.output.db)
        # output is temporary. No need to keep this file.
        output: temp(dirs_dict["CONTIGS_DIR"] + "/{group}-gene-calls.fa")
        threads: T('run_taxonomy_with_centrifuge')
        resources: nodes = T('run_taxonomy_with_centrifuge'),
        shell: "anvi-get-dna-sequences-for-gene-calls -c {input} -o {output} &>> {log}"


    rule run_centrifuge:
        ''' Run centrifuge on the exported gene calls of the contigs.db'''
        version: 1.0
        log: dirs_dict["LOGS_DIR"] + "/{group}-run_centrifuge.log"
        input: rules.export_gene_calls.output
        output:
            hits = dirs_dict["CONTIGS_DIR"] + "/{group}-centrifuge_hits.tsv",
            report = dirs_dict["CONTIGS_DIR"] + "/{group}-centrifuge_report.tsv"
        params: db=config['centrifuge']['db']
        threads: T('run_centrifuge', 5)
        resources: nodes = T('run_centrifuge', 5),
        shell: "centrifuge -f -x {params.db} {input} -S {output.hits} --report-file {output.report} --threads {threads} &>> {log}"


    rule import_taxonomy:
        ''' Run anvi-import-taxonomy '''
        version: 1.0
        log: dirs_dict["LOGS_DIR"] + "/{group}-import_taxonomy.log"
        input:
            hits = rules.run_centrifuge.output.hits,
            report = rules.run_centrifuge.output.report,
            # marking the contigs.db as ancient in order to ignore timestamps.
            contigs = ancient(rules.gen_contigs_db.output.db)
        # using a flag file because no file is created by this rule.
        # for more information see:
        # http://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#flag-files
        output: touch(dirs_dict["CONTIGS_DIR"] + "/{group}-anvi_import_taxonomy.done")
        params: parser = "centrifuge"
        threads: T('import_taxonomy')
        resources: nodes = T('import_taxonomy'),
        shell: "anvi-import-taxonomy -c {input.contigs} -i {input.report} {input.hits} -p {params.parser} &>> {log}"


if run_anvi_run_hmms:
    rule anvi_run_hmms:
        """ Run anvi-run-hmms"""
        # TODO: add rule for running hmms for ribosomal genes and import
        # their new gene calls.
        version: 1.0
        log: dirs_dict["LOGS_DIR"] + "/{group}-anvi_run_hmms.log"
        # if the user requested to run taxonomy using centrifuge, then this
        # will be ran only after centrifuge finished. Otherwise, this rule
        # will run after anvi-gen-contigs-database
        # marking the input as ancient in order to ignore timestamps.
        input: ancient(rules.gen_contigs_db.output.db)
        # using a snakemake flag file as an output since no file is generated
        # by the rule.
        output: touch(dirs_dict["CONTIGS_DIR"] + "/anvi_run_hmms-{group}.done")
        threads: T('run_anvi_run_hmms', 20)
        resources: nodes = T('run_anvi_run_hmms', 20),
        shell: "anvi-run-hmms -c {input} -T {threads} &>> {log}"


rule bowtie_build:
    """ Run bowtie-build on the contigs fasta"""
    # TODO: consider runnig this as a shadow rule
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-bowtie_build.log"
    input: rules.remove_human_dna_using_centrifuge.output if run_remove_human_dna_using_centrifuge else rules.reformat_fasta.output.contig
    # I touch this file because the files created have different suffix
    output:
        o1 = expand(dirs_dict["MAPPING_DIR"] + "/{group}/{group}-contigs" + '.{i}.bt2', i=[1,2,3,4], group="{group}"),
        o2 = expand(dirs_dict["MAPPING_DIR"] + "/{group}/{group}-contigs" + '.rev.{i}.bt2', i=[1,2], group="{group}")
    params:
        prefix = dirs_dict["MAPPING_DIR"] + "/{group}/{group}-contigs"
    threads: T('bowtie_build', 4)
    resources: nodes = T('bowtie_build', 4),
    shell: "bowtie2-build {input} {params.prefix} &>> {log}"


rule bowtie:
    """ Run mapping with bowtie2,  sort and convert to bam with samtools"""
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{sample}-bowtie.log"
    input:
        build_output = rules.bowtie_build.output,
        r1 = dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_R1.fastq.gz",
        r2 = dirs_dict["QC_DIR"] + "/{sample}-QUALITY_PASSED_R2.fastq.gz"
    # setting the output as temp, since we only want to keep the bam file.
    output: temp(dirs_dict["MAPPING_DIR"] + "/{group}/{sample}.sam")
    params:
        dir = dirs_dict["MAPPING_DIR"] + "/{sample}",
        bowtie_build_prefix = rules.bowtie_build.params.prefix
    threads: T('bowtie', 10)
    resources: nodes = T('bowtie', 10),
    shell: "bowtie2 --threads {threads} -x {params.bowtie_build_prefix} -1 {input.r1} -2 {input.r2} --no-unal -S {output} &>> {log}"


rule samtools_view:
    """ sort sam file with samtools and create a RAW.bam file"""
    version: 1.0
    log: dirs_dict["LOGS_DIR"] + "/{group}-{sample}-samtools_view.log"
    input: rules.bowtie.output
    # output as temp. we only keep the final bam file
    output: temp(dirs_dict["MAPPING_DIR"] + "/{group}/{sample}-RAW.bam")
    threads: T('samtools_view', 4)
    resources: nodes = T('samtools_view', 4),
    shell: "samtools view -F 4 -bS {input} -o {output} &>> {log}"


rule anvi_init_bam:
    """
        run anvi-init-bam on RAW bam file to create a bam file ready for
        anvi-profile.
    """
    version: 1.0 # later we can decide if we want the version to use the version of anvi'o
    log: dirs_dict["LOGS_DIR"] + "/{group}-{sample}-anvi_init_bam.log"
    input: rules.samtools_view.output
    output:
        bam = dirs_dict["MAPPING_DIR"] + "/{group}/{sample}.bam",
        bai = dirs_dict["MAPPING_DIR"] + "/{group}/{sample}.bam.bai"
    threads: T('anvi_init_bam', 4)
    resources: nodes = T('anvi_init_bam', 4),
    shell: "anvi-init-bam {input} -o {output.bam} &>> {log}"


rule anvi_profile:
    """ run anvi-profile on the bam file"""
    # setting the rule version to be as the version of the profile database of anvi'o
    version: anvio.__profile__version__
    log: dirs_dict["LOGS_DIR"] + "/{group}-{sample}-anvi_profile.log"
    input:
        bam = dirs_dict["MAPPING_DIR"] + "/{group}/{sample}.bam",
        # TODO: add option to profile all to all (all samples to all contigs)
        # marking the contigs.db as ancient in order to ignore timestamps.
        contigs = ancient(dirs_dict["CONTIGS_DIR"] + "/{group}-contigs.db")
    output:
        profile = dirs_dict["PROFILE_DIR"] + "/{group}/{sample}/PROFILE.db",
        aux = dirs_dict["PROFILE_DIR"] + "/{group}/{sample}/AUXILIARY-DATA.h5",
        runlog = dirs_dict["PROFILE_DIR"] + "/{group}/{sample}/RUNLOG.txt"
    params:
        # minimal length of contig to include in the profiling
        # if not specified in the config file then default is 2,500.
        min_contig_length = A(["anvi_profile", "min_contig_length"], config, default_value=2500),
        # if profiling to individual assembly then clustering contigs
        # see --cluster-contigs in the help manu of anvi-profile
        cluster_contigs = lambda wildcards: '--cluster-contigs' if group_sizes[wildcards.group] == 1 else '',
        name = "{sample}",
        profile_AA = "--profile-AA-frequencies" if A(["anvi_profile", "profile_AA"], config)  else "",
        output_dir = dirs_dict["PROFILE_DIR"] + "/{group}/{sample}"
    threads: T('anvi_profile', 5)
    resources: nodes = T('anvi_profile', 5),
    shell: "anvi-profile -i {input.bam} -c {input.contigs} -o {params.output_dir} -M {params.min_contig_length} -S {params.name} -T {threads} --overwrite-output-destinations {params.cluster_contigs} {params.profile_AA} &>> {log}"


def input_for_anvi_merge(wildcards):
    '''
        Create dictionary as input for rule anvi_merge.
        The reason we need a function as an input is to allow the user
        to choose between an option of an "all against all" vs. "normal"
        modes. See the documentation to learn more about the difference
        between these modes.
    '''

    if A('all_against_all', config):
        # If the user specified 'all against all' in the configs file
        # the end product would be a merge of all samples per group
        profiles = expand(dirs_dict["PROFILE_DIR"] + "/{group}/{sample}/PROFILE.db", sample=list(samples_information['sample']), group=wildcards.group)

    else:
        # The default behaviour is to only merge (and hence map and profile)
        # together samples that belong to the same group.
        profiles = expand(dirs_dict["PROFILE_DIR"] + "/{group}/{sample}/PROFILE.db", sample=list(samples_information[samples_information['group'] == wildcards.group]['sample']), group=wildcards.group)

    return profiles


def create_fake_output_files(_message, output):
    # creating "fake" output files with an informative message for
    # user.
    for o in output:
        with open(o, 'w') as f:
            f.write(_message + '\n')


def remove_empty_profile_databases(profiles, group):
    '''remove profiles that recruited zero reads from the metagenome.'''

    empty_profiles = []
    progress.new("Checking for empty profile databases")
    for p in profiles:
        db = dbops.ProfileDatabase(p)
        n = next(iter(db.meta['total_reads_mapped'].values()))
        if n == 0:
            # this profile is empty so we can't include it in the merged profile.
            profiles.remove(p)
            empty_profiles.append(p)
    progress.end()

    return profiles

    if not profiles:
        # if there are no profiles to merge then notify the user
        run.warning('It seems that all your profiles are empty for the \
                     contigs database: %s.db. And so cannot be merged.' \
                     % group)

    run.info('Number of non-empty profile databases', len(profiles))
    run.info('Number of empty profile databases', len(empty_profiles))
    if len(empty_profiles) > 0:
        run.info('The following databases are empty: %s' % empty_profiles,)


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
    log: dirs_dict["LOGS_DIR"] + "/{group}-anvi_merge.log"
    # The input are all profile databases that belong to the same group
    input:
        # marking the contigs.db as ancient in order to ignore timestamps.
        contigs = ancient(rules.gen_contigs_db.output.db),
        profiles = input_for_anvi_merge,
        # this is here just so snakemake would run the taxonomy before running this rule
        taxonomy = rules.import_taxonomy.output if run_taxonomy_with_centrifuge else ancient(rules.gen_contigs_db.output.db),
        # this is here just so snakemake would run the hmms before running this rule
        hmms = rules.anvi_run_hmms.output if run_anvi_run_hmms else ancient(rules.gen_contigs_db.output.db)
    output:
        profile = dirs_dict["MERGE_DIR"] + "/{group}/PROFILE.db",
        aux = dirs_dict["MERGE_DIR"] + "/{group}/AUXILIARY-DATA.h5",
        runlog = dirs_dict["MERGE_DIR"] + "/{group}/RUNLOG.txt"
    threads: T('anvi_merge')
    resources: nodes = T('anvi_merge'),
    params:
        output_dir = dirs_dict["MERGE_DIR"] + "/{group}",
        name = "{group}",
        profile_dir = dirs_dict["PROFILE_DIR"] + "/{group}"
    run:
        # using run instead of shell so we can choose the appropriate shell command.
        # In accordance with: https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency#comment-29348196

        # remove empty profile databases
        input.profiles = remove_empty_profile_databases(input.profiles, wildcards.group)

        if not input.profiles:
            # there are no profiles to merge.
            # this should only happen if all profiles were empty.
            _message = "Nothing to merge for %s. This should " \
                       "only happen if all profiles were empty " \
                       "(you can check the log file: {log} to see " \
                       "if that is indeed the case). " \
                       "This file was created just so that your workflow " \
                       "would continue with no error (snakemake expects " \
                       "to find these output files and if we don't create " \
                       "them, then it will be upset). As we see it, " \
                       "there is no reason to throw an error here, since " \
                       "you mapped your metagenome to some fasta files " \
                       "and you got your answer: whatever you have in " \
                       "your fasta file is not represented in your  " \
                       "metagenomes. Feel free to contact us if you think " \
                       "That this is our fault. sincerely, Meren Lab" \
                       % wildcards.group
            # creating the expected output files for the rule
            create_fake_output_files(_message, output)

        elif group_sizes[wildcards.group] == 1:
            # for individual assemblies, create a symlink to the profile database
            #shell("ln -s {params.profile_dir}/*/* -t {params.output_dir} &>> {log}")
            #shell("touch -h {params.profile_dir}/*/*")

            # Still waiting to get an answer on this issue:
            # https://groups.google.com/d/msg/snakemake/zU_wkfZ7YCs/GZP0Z_RoAgAJ
            # Until then, I will just create fake file so snakemake is happy
            _message = "Only one file was profiles with %s so there " \
                       "is nothing to merge. But don't worry, you can " \
                       "still use anvi-interacite with the single profile " \
                       "database that is here: %s" \
                       % (wildcards.group, input.profiles)
            create_fake_output_files(_message, output)

        elif len(input.profiles) == 1:
            # if only one sample is not empty, but the group size was
            # bigger than 1 then it means that --cluster-contigs was
            # not performed during anvi-profile.
            _message = "Only one samlpe had reads recruited to %s " \
                       "and hence merging could not occur." \
                       % wildcards.group
            create_fake_output_files(_message, output)

        else:
            shell("anvi-merge {input.profiles} -o {params.output_dir} -c {input.contigs} -S {params.name} --overwrite-output-destinations &>> {log}")

