# Snakemake workflow for assembly based metagenomics

The majority of the steps used in this pipeline are extensively described in the [anvi'o user tutorial for metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). But this pipeline includes also the first steps that are not described in the anvi'o user tutorial for metagenomic workflow. The entering point to this pipeline are the unprocessed raw reads of a collection (or a single) metagenomes, and the output of the pipline is an anvi'o merged profile database ready for refinement of bins (or whatever it is that you want to do with it).
The pipline includes the following steps:

1. QC of the metagenomes using [illumina-utils](https://github.com/merenlab/illumina-utils/).
2. (Co-)Assembly using [megahit](https://github.com/voutcn/megahit).
3. Generating an anvi'o CONTIGS database.
4. Store HMM hits in the CONTIGS database using [anvi-run-hmms](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms).
5. Run [centrifuge](https://ccb.jhu.edu/software/centrifuge/) and import taxonomy to the CONTIGS database using [anvi-import-taxonomy](http://merenlab.org/2016/06/18/importing-taxonomy/).
6. Mapping short reads using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
7. Profiling the individual bam files using [anvi-profile](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile).
8. Merging the individual profile databases using [anvi-merge](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-merge).

A directed acyclic graph (DAG) describing the workflow for a mock dataset could be seen below:

![alt text][mock_files_for_merenlab_metagenomics_pipeline/mock-dag.png]


If you want to create a DAG for your dataset, simply run:

```
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --dag | dot -Tsvg > dag.svg
```

### Using `dot` on mac

If you are on mac, you can use `dot` by installing `graphviz`, simply run `brew install graphviz`.

## Standard usage

All you need is a bunch of fastq files and a `samples.txt` file to specify the names of your samples, which group they belong to (optional),
and the path to the pair-end fastq files. For now, only pair end fastq files are supported.

An example for a properly formatted `.txt` file is available [here](mock_files_for_merenlab_metagenomics_pipeline/samples.txt).
If you prefer a different name for your samples text file, then you can provide a different name in the config file (see below).
If nothing was provided in the config file then the default is `samples.txt`.

## The config file

To make changes easy and accessible for the user, we tried our best to make all relevant configurations available
to the user using a `json` formatted config file, and thus avoiding the need to change the Snakefile. For an example 
config file go [here](mock_files_for_merenlab_metagenomics_pipeline/config.json).

### Output directories

The following directories are expected to exist at the end of the workflow:
"00_LOGS", "01_QC", "02_ASSEMBLY", "03_CONTIGS", "04_MAPPING", "05_ANVIO_PROFILE", "06_MERGED"

Don't like these names? You can specify what is the name of the folder, by providing the following information in the config file:

```
    "output_dirs":{
        "LOGS_DIR" : "00_MY_beAuTiFul_LOGS",
        "QC_DIR" : "BEST_QC_DIR_EVER",
        "CONTIGS_DIR" : "/absolute/path/to/my/contigs/dir",
        "MAPPING_DIR" : "make/mapping/great/again"
    }
```

You can change all or just some of the names of these folders.
And you can provide an absolute path or a relative path.

### Optional steps

The following steps are only optional (the default behaviour is to run all optional steps):
1. Assigning taxonomy with centrifuge.
2. Running hmm profiles on the contigs database.

To choose not to run these steps add these lines in your config file:
```
    "run_centrifuge":{
        "run": "False"
    },
    "anvi_run_hmms":{
        "run": "False"
    }
```

### The "all against all" option

The default behaviour for this workflow is to create a contigs database for each _group_ and map (and profile, and merge) the samples that belong to that _group_. If you wish to map all samples to all contigs, use the `all_against_all` option in the config file:

```
    "all_against_all: "True"
```

For those of you who are learning `snakemake`, you might be surprised of how easy the switch between the modes is. All we need to do is tell the `anvi_merge` rule that we want all samples merged for each _group_, and snakemake immediatly infers that it needs to also run the extra mapping, and profiling steps. *Thank you snakemake!* (says everyone).

An updated DAG for the workflow for our mock data is available below:

![alt text][mock_files_for_merenlab_metagenomics_pipeline/mock-dag-all-against-all.png]

A little more of a mess! But also has a beauty to it :-).