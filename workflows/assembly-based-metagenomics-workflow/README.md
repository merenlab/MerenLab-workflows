# Snakemake workflow for assembly based metagenomics

**Important note**: this pipeline was evaluated using snakemake version 3.12.0. If you are using an older version, then we suggest upgrading to the newest version.

The majority of the steps used in this pipeline are extensively described in the [anvi'o user tutorial for metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/). But this pipeline includes also the first steps that are not described in the anvi'o user tutorial for metagenomic workflow. The entering point to this pipeline are the unprocessed raw reads of a collection (or a single) metagenomes, and the output of the pipline is an anvi'o merged profile database ready for refinement of bins (or whatever it is that you want to do with it).
The pipline includes the following steps:

1. QC of the metagenomes using [illumina-utils](https://github.com/merenlab/illumina-utils/).
2. (Co-)Assembly using [megahit](https://github.com/voutcn/megahit).
3. Generating an anvi'o CONTIGS database.
4. Store HMM hits in the CONTIGS database using [anvi-run-hmms](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms) (optional step).
5. Run [centrifuge](https://ccb.jhu.edu/software/centrifuge/) and import taxonomy to the CONTIGS database using [anvi-import-taxonomy](http://merenlab.org/2016/06/18/importing-taxonomy/) (optional step).
6. Mapping short reads using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
7. Profiling the individual bam files using [anvi-profile](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile).
8. Merging the individual profile databases using [anvi-merge](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-merge).

A directed acyclic graph (DAG) describing the workflow for a mock dataset could be seen below:

![Alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag.png?raw=true "mock-dag")


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

## Running on a cluster

When submitting to a cluster, you can utilize the [snakemake cluster execution](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution). Notice that the number of threads per rule could be changed using the `config.json` file (and not by using the [cluster configuration](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) file). For more details, refer to the documentation of the configuration file below.

When submitting a workflow to a cluster, snakemake requires you to limit the number of jobs using `--jobs`. If you prefer to limit the number of threads that would be used by your workflow (for example, if you share your cluster with others and you don't want to take control of all threads), then we made use of the snakemake built-in `resources`. You can set the number of jobs to your limit (or to a very big number like this one: 99834 if you dont care), and use `--resources nodes=30`, if you wish to only use 30 threads. We used the word `nodes` so that to not confuse with the reserved word `threads` in snakemake.

### A note on cluster-config

This note is here mainly for documentation of the code, and for those of you who are interested in snakemake. The reason we decided not to use the [cluster configuration](http://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) file to control the number of threads per rule, is becuase certain software require the number of threads as an input (for example `megahit` and `anvi-profile`), but the cluster config file is not available for shell commands within snakemake rules. To bypass this issue we simply put the threads configuration in the `config.json`, thus available for the user to modify.

## The config file

To make changes easy and accessible for the user, we tried our best to make all relevant configurations available
to the user using a `json` formatted config file, and thus avoiding the need to change the Snakefile. For an example 
config file go [here](mock_files_for_merenlab_metagenomics_pipeline/config.json). There are some general configurations, and there are step specific configurations.

## General configurations

### Output directories

The following directories are expected to exist at the end of the workflow:
"00\_LOGS", "01\_QC", "02\_ASSEMBLY", "03\_CONTIGS", "04\_MAPPING", "05\_ANVIO_PROFILE", "06\_MERGED"

Don't like these names? You can specify what is the name of the folder, by providing the following information in the config file:

```
    "output_dirs":{
        "LOGS_DIR" : "00_MY_beAuTiFul_LOGS",
        "QC_DIR" : "BEST_QC_DIR_EVER",
        "ASSEMBLY_DIR" : "assemblies",
        "CONTIGS_DIR" : "/absolute/path/to/my/contigs/dir",
        "MAPPING_DIR" : "relative/path/to/my/mapping/dir",
        "PROFILE_DIR": "/I/already/did/my/profiling/and/this/is/where/you/can/find/it/",
        "MERGE_DIR": "06_Keep_Calm_and_Merge_On"
    }
```

When using "references mode" (see below) the default name for the `ASSEMBLY_DIR` is `02_REFERENCE_FASTA` (since in references mode no assembly is performed in the workflow, and the fasta files could be not assemblies, but rather full genomes). In order to change it, you may use `"REFERENCES_DIR"`, i.e.:

```
    "output_dirs":{
    	"REFERENCES_DIR" : "02_REF"
    }
```

You can change all or just some of the names of these folders.
And you can provide an absolute path or a relative path.

### The "all against all" option

The default behaviour for this workflow is to create a contigs database for each _group_ and map (and profile, and merge) the samples that belong to that _group_. If you wish to map all samples to all contigs, use the `all_against_all` option in the config file:

```
    "all_against_all: "True"
```

For those of you who are learning `snakemake`, you might be surprised of how easy the switch between the modes is. All we need to do is tell the `anvi_merge` rule that we want all samples merged for each _group_, and snakemake immediatly infers that it needs to also run the extra mapping, and profiling steps. *Thank you snakemake!* (says everyone).

An updated DAG for the workflow for our mock data is available below:

![alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag-all-against-all.png?raw=true "mock-dag-all-against-all")

A little more of a mess! But also has a beauty to it :-).

### Define the number of threads per rule in a cluster

In order to change the number of threads per rule when running on a cluster, the following structure should be used: 

```
	"rule_name":
		"threads": number_of_threads
```

The following defaults have been set:

**rule**|**threads**
:-----:|:-----:
qc|2
megahit|11
gen\_contigs\_db|5
run\_centrifuge|5
run\_anvi\_run\_hmms|20
bowtie\_build|4
bowtie|10
samtools\_view|4
anvi\_init\_bam|4
anvi\_profile|5

All other rules use 1 thread by default.

### Optional steps

The following steps are only optional:

1. Assigning taxonomy with centrifuge (default is **not** running).
2. Running hmm profiles on the contigs database (default is **running**).

For more details refer to the specific documentation for these steps below.

## Step-specific configurations 

Some of the steps in the workflow have parameters with defaults that could be changed. We tried to keep things flexible and accessible for the user, but we know we didn't do everything possible. If there is something that you want to have access to and is not possible, please create an issue on [github](https://github.com/merenlab/MerenLab-workflows/issues). Or, better yet, make those changes and send us a pull request. We plan to do a better job to let you access in a flexible form all the parameters of each step, and if this is of special interest to you, you can refer to the note below regarding wrappers.

The step-specific configurations in the `config.json` file always have the following structure:
```
	"step_name":{
		"configurable_parameter": "value"
	}
```

Notice that everything has to have quotation marks (to be compatible with the JSON format).

### `qc`

If you already performed QC, and so wish to skip qc, then simply add this to your config file:

```
	"qc": {
		"run": false
	}
```

A nice trick worth knowing: if you only want to qc your files and then compress them (and not do anything else), simply invoke the workflow with the following command:

```
snakemake --snakefile merenlab-metagenomics-pipeline.snakefile --until gzip_fastqs
```

To understand this better, refer to the snakemake documentation.

### `reformat_fasta`

In "references mode", you may choose to skip this step, and keep your contigs names. In order to do so, add this to your config file:

```json
	"reformat_fasta": {
		"run": false
	}
```

In assembly mode, this rule is always excecuted.

### `megahit`

The following parameters are available:

`memory` (see `-m/--memory` in the megahit documentation) - The default is 0.4.

`min_contig_len` (`--min-contig-len`) - default is 1,000.

### `run_centrifuge`

`run` - could get values of `true` or `false` (all lower case!) - to configure whether to run centrifuge or not. The default is `false`.

`db` - if you choose run centrifuge, you **must** provide the path to the database (for example `$CENTRIFUGE_BASE/p+h+v/p+h+v`).

### `run_anvi_hmms`

`run` - could get values of `true` or `false` (all lower case!) - to configure whether to run hmms or not. The default is `true`.

### `anvi_profile`

`min_contig_length` - see anvi-profile documentation for `--min-contig-length`. The default is going with the default of `anvi-profile` (which is 2,500).

### example

So let's say I want to run centrifuge, I don't want to run hmms, and I want my minimum contig length for megahit and anvi-profile to be 500 and and 3,000 respectively. Then my config file would like like this:

```
{
	"run_centrifuge":{
		"run": true,
		"db": "$CENTRIFUGE_BASE/p+h+v/p+h+v"
	},
	"run_anvi_hmms":{
		"run": false
	},
	"anvi_profile:{
		"min_contig_length": 3000,
		"threads": 10
	},
	"megahit":{
		"min_contig_len": 500
	}
}
```
## Estimating occurence of population genomes in metagenomes

Along with assembly-based metagenomics, we often use anvi'o to explore the occurence of population genomes accross metagenomes. You can see a nice example of that here: [Please insert a nice example here. Probably the blog about DWH thingy](link-to-nice-example).
In that case, what you have is a bunch of fastq files (metagenomes) and fasta files (reference genomes), and all you need to do is to let the workflow know where to find these files, using to `.txt` files: `samples.txt`, and `references.txt`. The `samples.txt` stays as before, but this time the `group` column will specify for each sample, which reference should be used. If the `samples.txt` files doesn't have a `group` column, then an "all against all" mode would be provoked. Below you can see how the DAG looks like for this mode:

![alt text](mock_files_for_merenlab_metagenomics_pipeline/mock-dag-references-mode.png?raw=true "mock-dag-references-mode")

## wrappers

add a note about transferring to wrappers
