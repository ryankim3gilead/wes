# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

# Introduction

Sarek is a workflow designed to detect germline and somatic variants on whole genome, whole exome, or targeted sequencing data.

Initially designed for human and mouse, it can work on any species if a reference genome is available.
Sarek is designed to handle single samples, such as single-normal or single-tumor samples, and tumor-normal pairs including additional relapses.

# Running the pipeline

## Quickstart

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/sarek --input samplesheet.csv --outdir <OUTDIR> --genome GATK.GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
> The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/sarek -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
input: 'data'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Input: Sample sheet configurations

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 3 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

If necessary, a tumor sample can be associated to a normal sample as a pair, if specified with the same `patient` ID, a different `sample`, and the respective `status`.
An additional tumor sample (such as a relapse for example), can be added if specified with the same `patient` ID, a different `sample`, and the `status` value `1`.

Sarek will output results in a different directory for _each sample_.
If multiple samples IDs are specified in the CSV file, Sarek will consider all files to be from different samples.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample and tool configuration (or normal/tumor pair if applicable).

### Overview: Samplesheet Columns

| Column    | Description                                                                                                                                                                                                                                                                                                                       |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `patient` | **Custom patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. normal and tumor). <br /> _Required_                                                                                                                                                      |
| `sex`     | **Sex chromosomes of the patient**; i.e. XX, XY..., only used for Copy-Number Variation analysis in a tumor/pair<br /> _Optional, Default: `NA`_                                                                                                                                                                                  |
| `status`  | **Normal/tumor status of sample**; can be `0` (normal) or `1` (tumor).<br /> _Optional, Default: `0`_                                                                                                                                                                                                                             |
| `sample`  | **Custom sample ID** for each tumor and normal sample; more than one tumor sample for each subject is possible, i.e. a tumor and a relapse; samples can have multiple lanes for which the _same_ ID must be used to merge them later (see also `lane`). Sample IDs must be unique for unique biological samples <br /> _Required_ |
| `lane`    | Lane ID, used when the `sample` is multiplexed on several lanes. Must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character <br /> _Required for `--step mapping`_                                                                                 |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `bam`     | Full path to (u)BAM file                                                                                                                                                                                                                                                                                                          |
| `bai`     | Full path to BAM index file                                                                                                                                                                                                                                                                                                       |
| `cram`    | Full path to CRAM file                                                                                                                                                                                                                                                                                                            |
| `crai`    | Full path to CRAM index file                                                                                                                                                                                                                                                                                                      |
| `table`   | Full path to recalibration table file                                                                                                                                                                                                                                                                                             |
| `vcf`     | Full path to vcf file                                                                                                                                                                                                                                                                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Start with mapping (`--step mapping` [default])

This step can be started either from FastQ files or (u)BAMs. The CSV must contain at least the columns `patient`, `sample`, `lane`, and either `fastq_1/fastq_2` or `bam`.

#### Examples

Minimal config file:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

```bash
patient,sample,lane,bam
patient1,test_sample,lane_1,test.bam
```

In this example, the sample is multiplexed over three lanes:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,test_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,test_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
```

```bash
patient,sample,lane,bam
patient1,test_sample,1,test_L001.bam
patient1,test_sample,2,test_L002.bam
patient1,test_sample,3,test_L003.bam
```

#### Full samplesheet

In this example, all possible columns are used. There are three lanes for the normal sample, two for the tumor sample, and one for the relapse sample, including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,lane,fastq_1,fastq_2
patient1,XX,0,normal_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,XX,0,normal_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,XX,0,normal_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
patient1,XX,1,tumor_sample,lane_1,test2_L001_1.fastq.gz,test2_L001_2.fastq.gz
patient1,XX,1,tumor_sample,lane_2,test2_L002_1.fastq.gz,test2_L002_2.fastq.gz
patient1,XX,1,relapse_sample,lane_1,test3_L001_1.fastq.gz,test3_L001_2.fastq.gz
```

```bash
patient,sex,status,sample,lane,bam
patient1,XX,0,normal_sample,lane_1,test_L001.bam
patient1,XX,0,normal_sample,lane_2,test_L002.bam
patient1,XX,0,normal_sample,lane_3,test_L003.bam
patient1,XX,1,tumor_sample,lane_1,test2_L001.bam
patient1,XX,1,tumor_sample,lane_2,test2_L002.bam
patient1,XX,1,relapse_sample,lane_1,test3_L001.bam
```

### Start with duplicate marking (`--step markduplicates`)

#### Duplicate Marking

For starting from duplicate marking, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`

> **NB:** When using [GATK4 MarkduplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) reads should be name-sorted for efficient execution

Example:

```bash
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/mapped.csv` if in a previous run `--save_mapped` was set and will automatically be used as an input when specifying the parameter `--step markduplicates`. Otherwise this file will need to be manually generated.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_mapped.bam,test_mapped.bam.bai
patient1,XX,1,tumor_sample,test2_mapped.bam,test2_mapped.bam.bai
patient1,XX,1,relapse_sample,test3_mapped.bam,test3_mapped.bam.bai
```

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

### Start with preparing the recalibration tables (`--step prepare_recalibration`)

For starting directly from preparing the recalibration tables, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```bash
patient,sample,bam,bai
patient1,test_sample,test_md.bam,test_md.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_md.cram,test_md.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates_no_table.csv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_md.bam,test_md.bam.bai
patient1,XX,1,tumor_sample,test2_md.bam,test2_md.bam.bai
patient1,XX,1,relapse_sample,test3_md.bam,test3_md.bam.bai
```

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_md.cram,test_md.cram.crai
patient1,XX,1,tumor_sample,test2_md.cram,test2_md.cram.crai
patient1,XX,1,relapse_sample,test3_md.cram,test3_md.cram.crai
```

### Start with base quality score recalibration (`--step recalibrate`)

For starting from base quality score recalibration the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai`, `table` or `patient`, `sample`, `cram`, `crai`, `table` containing the paths to _non-recalibrated CRAM/BAM_ files and the associated recalibration table.

Example:

```bash
patient,sample,bam,bai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

```bash
patient,sample,cram,crai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates.csv` and will automatically be used as an input when specifying the parameter `--step recalibrate`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,cram,crai,table
patient1,XX,0,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai,test2.table
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai,test3.table
```

### Start with variant calling (`--step variant_calling`)

For starting from the variant calling step, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```bash
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/recalibrated.csv` and will automatically be used as an input when specifying the parameter `--step variant_calling`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

### Start with annotation (`--step annotate`)

For starting from the annotation step, the CSV file must contain at least the columns `patient`, `sample`, `vcf`.

As Sarek will use [bgzip](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html) to compress and index the annotated VCF files, it expects the input VCF files to be sorted and compressed.

Example:

```bash
patient,sample,vcf
patient1,test_sample,test.vcf.gz
```

The Sarek-generated CSV file is stored under `results/csv/variantcalled.csv` and will automatically be used as an input when specifying the parameter `--step annotation`.

#### Full samplesheet

In this example, all possible columns are used including the `variantcaller` information per sample:

```bash
patient,sample,variantcaller,vcf
test,sample3,strelka,sample3.variants.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.diploid_sv.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.somatic_sv.vcf.gz
```

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/sarek
```

## Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `3.1.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 3.1.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

# Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

## `-profile`

Use this parameter to choose a configuration profile.
Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of `Docker` or `Singularity` containers for full pipeline reproducibility, however when this is not possible, `Conda` is also supported.

The pipeline also dynamically loads configurations from [github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

## `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

# Custom configuration

## Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

# Troubleshooting & FAQ

## How to test the pipeline

When using default parameters only, sarek runs preprocessing and `Strelka2`.
This is reflected in the default test profile:

```bash
nextflow run nf-core/sarek -r 3.2.1 -profile test,<container/institute> --outdir results
```

Expected run output:

```bash
[85/6b7739] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM1_INDEX (genome.fasta)                                                [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                                               -
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                                                           -
[22/cf54a8] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (genome.fasta)                               [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                                                           -
[28/dad25a] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX (genome.fasta)                                               [100%] 1 of 1 ✔
[23/3fe964] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_DBSNP (dbsnp_146.hg38.vcf)                                            [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                                                     -
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_KNOWN_SNPS                                                            -
[14/26e286] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS (mills_and_1000G.indels.vcf)                             [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_PON                                                                   -
[76/04d107] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:CREATE_INTERVALS_BED (genome.interval_list)                              [100%] 1 of 1 ✔
[d4/f97174] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:GATK4_INTERVALLISTTOBED (genome)                                         [100%] 1 of 1 ✔
[70/82ba3c] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr22_1-40001)                          [100%] 1 of 1 ✔
[d4/c2d0c4] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (genome)                              [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                                                  -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_UNMAP                                              -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_MAP                                                -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_UNMAP                                                -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_MERGE_UNMAP                                                   -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_UNMAP                                                    -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_MAP                                                      -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:CAT_FASTQ                                                              -
[c4/f59e5a] process > NFCORE_SAREK:SAREK:FASTQC (test-test_L1)                                                                      [100%] 1 of 1 ✔
[0b/c5a999] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:BWAMEM1_MEM (test)                                         [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:BWAMEM2_MEM                                                -
[-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:DRAGMAP_ALIGN                                              -
[c7/664cd1] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test)                                             [100%] 1 of 1 ✔
[13/bc73b6] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:INDEX_MARKDUPLICATES (test)                                             [100%] 1 of 1 ✔
[2a/99608e] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (test)                         [100%] 1 of 1 ✔
[f2/0420ca] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (test)                               [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM                                                                                -
[eb/46945a] process > NFCORE_SAREK:SAREK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (test)                                         [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS                                               -
[ec/2377d4] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:GATK4_APPLYBQSR (test)                                                       [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM                                         -
[88/3af664] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (test)                                  [100%] 1 of 1 ✔
[f4/828fde] process > NFCORE_SAREK:SAREK:CRAM_QC_RECAL:SAMTOOLS_STATS (test)                                                        [100%] 1 of 1 ✔
[fb/a9d66f] process > NFCORE_SAREK:SAREK:CRAM_QC_RECAL:MOSDEPTH (test)                                                              [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM_RECAL                                                                          -
[ef/026185] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:STRELKA_SINGLE (test)  [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA          -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA_GENOME   -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:STRELKA_SINGLE       -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA        -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA_GENOME -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:STRELKA_SOMATIC        -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:MERGE_STRELKA_INDELS   -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:MERGE_STRELKA_SNVS     -
[bc/f3f5cf] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (test)                                             [100%] 1 of 1 ✔
[21/8d4f02] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (test)                                        [100%] 1 of 1 ✔
[36/957fba] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (test)                                         [100%] 1 of 1 ✔
[70/a8e064] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_SUMMARY (test)                                           [100%] 1 of 1 ✔
[36/e35b1b] process > NFCORE_SAREK:SAREK:CUSTOM_DUMPSOFTWAREVERSIONS (1)                                                            [100%] 1 of 1 ✔
[3f/3c3356] process > NFCORE_SAREK:SAREK:MULTIQC                                                                                    [100%] 1 of 1 ✔
-[nf-core/sarek] Pipeline completed successfully-
Completed at: 09-Jun-2023 13:46:31
Duration    : 1m 50s
CPU hours   : (a few seconds)
Succeeded   : 27
```

The pipeline comes with a number of possible paths and tools that can be used.

Due to the small test data size, unfortunately not everything can be tested from top-to-bottom, but often is done by utilizing the pipeline's `--step` parameter.

For more extensive testing purpose, we have the `test_cache` profile that contain the same data, but on which the path to the reference and input files can be changed using the `--test_data_base` params.

Annotation is generally tested separately from the remaining workflow, since we use references for `C.elegans`, while the remaining tests are run on downsampled human data.

```bash
nextflow run nf-core/sarek -r 3.2.1 -profile test_cache,<container/institute> --outdir results --tools snpeff --step annotation
```

If you are interested in any of the other tests that are run on every code change or would like to run them yourself, you can take a look at `tests/<filename>.yml`.
For each entry the respective nextflow command run and the expected output is specified.

Some of the currently, available test profiles:

| Test profile    | Run command                                                                           |
| :-------------- | :------------------------------------------------------------------------------------ |
| annotation      | `nextflow run main.nf -profile test_cache,annotation,docker --tools snpeff,vep,merge` |
| no_intervals    | `nextflow run main.nf -profile test_cache,no_intervals,docker`                        |
| targeted        | `nextflow run main.nf -profile test_cache,targeted,docker`                            |
| tools_germline  | `nextflow run main.nf -profile test_cache,tools_germline,docker --tools strelka`      |
| tools_tumoronly | `nextflow run main.nf -profile test_cache,tools_tumoronly,docker --tools strelka`     |
| tools_somatic   | `nextflow run main.nf -profile test_cache,tools_somatic,docker --tools strelka`       |
| trimming        | `nextflow run main.nf -profile test_cache,trim_fastq,docker`                          |
| umi             | `nextflow run main.nf -profile test_cache,umi,docker`                                 |
| use_gatk_spark  | `nextflow run main.nf -profile test_cache,use_gatk_spark,docker`                      |

If you are interested in any of the other profiles that are used, you can take a look at `conf/test/<filename>.config`.

## How can the different steps be used

Sarek can be started at different points in the analysis by setting the parameter `--step`. Once started at a certain point, the pipeline runs through all the following steps without additional intervention. For example when starting from `--step mapping` (set by default) and `--tools strelka,vep`, the input reads will be aligned, duplicate marked, recalibrated, variant called with Strelka, and finally VEP will annotate the called variants.

## Which variant calling tool is implemented for which data type?

This list is by no means exhaustive and it will depend on the specific analysis you would like to run. This is a suggestion based on the individual docs of the tools specifically for human genomes and a garden-variety sequencing run as well as what has been added to the pipeline.

| Tool                                                                                                    | WGS | WES |  Panel |  Normal | Tumor | Somatic |
| :------------------------------------------------------------------------------------------------------ | :-: | :-: | :----: | :-----: | :---: | :-----: |
| [DeepVariant](https://github.com/google/deepvariant)                                                    |  x  |  x  |   x    |    x    |   -   |    -    |
| [FreeBayes](https://github.com/ekg/freebayes)                                                           |  x  |  x  |   x    |    x    |   x   |    x    |
| [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) |  x  |  x  |   x    |    x    |   -   |    -    |
| [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2)                 |  x  |  x  |   x    |    -    |   x   |    x    |
| [mpileup](https://www.htslib.org/doc/samtools-mpileup.html)                                             |  x  |  x  |   x    |    x    |   x   |    -    |
| [Strelka2](https://github.com/Illumina/strelka)                                                         |  x  |  x  |   x    |    x    |   x   |    x    |
| [Manta](https://github.com/Illumina/manta)                                                              |  x  |  x  |   x    |    x    |   x   |    x    |
| [TIDDIT](https://github.com/SciLifeLab/TIDDIT)                                                          |  x  |  x  |   x    |    x    |   x   |    x    |
| [ASCAT](https://github.com/VanLoo-lab/ascat)                                                            |  x  |  x  |   -    |    -    |   -   |    x    |
| [CNVKit](https://cnvkit.readthedocs.io/en/stable/)                                                      |  x  |  x  |   -    |    x    |   x   |    x    |
| [Control-FREEC](https://github.com/BoevaLab/FREEC)                                                      |  x  |  x  |   x    |    -    |   x   |    x    |
| [MSIsensorPro](https://github.com/xjtu-omics/msisensor-pro)                                             |  x  |  x  |   x    |    -    |   -   |    x    |

## How to run ASCAT with whole-exome sequencing data?

ASCAT runs out of the box on whole genome sequencing data using iGenomes resources. While the ASCAT implementation in sarek is capable of running with whole-exome sequencing data, the needed references are currently not provided with the igenomes.config. According to the [developers](https://github.com/VanLoo-lab/ascat/issues/97) of ASCAT, loci and allele files (one file per chromosome) can be downloaded directly from the [Battenberg repository](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52).

Please note that:

- Row names (for GC and RT correction files) should be `${chr}_${position}` (there is no SNP/probe ID for HTS data).
- All row names in GC and RT correction files should also appear in the loci files
- Loci and allele files must contain the same set of SNPs
- ASCAT developers strongly recommend using a BED file for WES/TS data. This prevents considering SNPs covered by off-target reads that would add noise to log/BAF tracks.
- The total number of GC correction loci in a sample must be at least 10% of the number of loci with logR values. If the number of GC correction loci is too small compared to the total number of loci, ASCAT will throw an error.

From 'Reference files' https://github.com/VanLoo-lab/ascat:

> For WES and targeted sequencing, we recommend using the reference files (loci, allele and logR correction files) as part of the Battenberg package. Because they require a high-resolution input, our reference files for WGS are not suitable for WES and targeted sequencing. For WES, loci and allele files from the Battenberg package can be fed into ascat.prepareHTS. For targeted sequencing, allele files from the Battenberg package can be fed into ascat.prepareTargetedSeq, which will generate cleaned loci and allele files that can be fed into ascat.prepareHTS.

### How to generate ASCAT resources for exome or targeted sequencing

1. Fetch the GC content correction and replication timing (RT) correction files from the [Dropbox links provided by the ASCAT developers](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS) and intersect the SNP coordinates with the exome target coordinates. If the target file has 'chr' prefixes, make a copy with these removed first. Extract the GC and RT information for only the on target SNPs and zip the results.

```bash
sed -e 's/chr//' targets_with_chr.bed > targets.bed

for t in GC RT
do
  unzip ${t}_G1000_hg38.zip

  cut -f 1-3 ${t}_G1000_hg38.txt > ascat_${t}_snps_hg38.txt
  tail -n +2 ascat_${t}_snps_hg38.txt | awk '{ print $2 "\t" $3-1 "\t" $3 "\t" $1 }' > ascat_${t}_snps_hg38.bed
  bedtools intersect -a ascat_${t}_snps_hg38.bed -b targets.bed | awk '{ print $1 "_" $3 }' > ascat_${t}_snps_on_target_hg38.txt

  head -n 1 ${t}_G1000_hg38.txt > ${t}_G1000_on_target_hg38.txt
  grep -f ascat_${t}_snps_on_target_hg38.txt ${t}_G1000_hg38.txt >> ${t}_G1000_on_target_hg38.txt
  zip ${t}_G1000_on_target_hg38.zip ${t}_G1000_on_target_hg38.txt

  rm ${t}_G1000_hg38.zip
done
```

2. Download the Battenberg 1000G loci and alleles files. The steps below follow downloading from the [Battenberg repository at the Oxford University Research Archive](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52). The files are also available via Dropbox links from the same page as the GC and RT correction files above.

```bash
wget https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/files/rt435gd52w
mv rt345gd52w battenberg.zip
jar xf battenberg.zip

unzip 1000G_loci_hg38_chr.zip
cd 1000G_loci_hg38
mkdir battenberg_alleles_on_target_hg38
mv *allele* battenberg_alleles_on_target_hg38/
mkdir battenberg_loci_on_target_hg38
mv *loci* battenberg_loci_on_target_hg38/
```

3. Copy the `targets_with_chr.bed` and `GC_G1000_on_target_hg38.txt` files into the newly created `battenberg_loci_on_target_hg38` folder before running the next set of steps. ASCAT generates a list of GC correction loci with sufficient coverage in a sample, then intersects that with the list of all loci with tumour logR values in that sample. If the intersection is <10% the size of the latter, it will fail with an error. Because the Battenberg loci/allele sets are very dense, subsetting to on-target regions is still too many loci. This script ensures that all SNPs with GC correction information are included in the loci list, plus a random sample of another 30% of all on target loci. You may need to vary this proportion depending on your set of targets. A good rule of thumb is that the size of your GC correction loci list should be about 15% the size of your total loci list. This allows for a margin of error.

```bash
cd battenberg_loci_on_target_hg38/
rm *chrstring*
rm 1kg.phase3.v5a_GRCh38nounref_loci_chr23.txt
for i in {1..22} X
do
   awk '{ print $1 "\t" $2-1 "\t" $2 }' 1kg.phase3.v5a_GRCh38nounref_loci_chr${i}.txt > chr${i}.bed
   grep "^${i}_" GC_G1000_on_target_hg38.txt | awk '{ print "chr" $1 }' > chr${i}.txt
   bedtools intersect -a chr${i}.bed -b targets_with_chr.bed | awk '{ print $1 "_" $3 }' > chr${i}_on_target.txt
   n=`wc -l chr${i}_on_target.txt | awk '{ print $1 }'`
   count=$((n * 3 / 10))
   grep -xf chr${i}.txt chr${i}_on_target.txt > chr${i}.temp
   shuf -n $count chr${i}_on_target.txt >> chr${i}.temp
   sort -n -k2 -t '_' chr${i}.temp | uniq | awk 'BEGIN { FS="_" } ; { print $1 "\t" $2 }' > battenberg_loci_on_target_hg38_chr${i}.txt
done
zip battenberg_loci_on_target_hg38.zip battenberg_loci_on_target_hg38_chr*.txt
```

4. Extract the alleles for the same set of SNPs. Uses a short R script defined below.

```bash
cd ../battenberg_alleles_on_target_hg38/
rm 1kg.phase3.v5a_GRCh38nounref_allele_index_chr23.txt
for i in {1..22} X
do
  Rscript intersect_ascat_alleles.R ../battenberg_loci_on_target_hg38/battenberg_loci_on_target_hg38_chr${i}.txt \
    1kg.phase3.v5a_GRCh38nounref_allele_index_chr${i}.txt battenberg_alleles_on_target_hg38_chr${i}.txt
done
zip battenberg_alleles_on_target_hg38.zip battenberg_alleles_on_target_hg38_chr*.txt
```

Rscript `intersect_ascat_alleles.R`

```bash
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

loci = read.table(args[1], header=F, sep="\t", stringsAsFactors=F)
alleles = read.table(args[2], header=T, sep="\t", stringsAsFactors=F)

i = intersect(loci$V2, alleles$position)

out = subset(alleles, alleles$position %in% i)
write.table(out, args[3], col.names=T, row.names=F, quote=F, sep="\t")
```

5. Move or copy all of the zip files you've created to a suitable location. Specify these in your parameters, e.g.

```json
{
  "ascat_alleles": "/path/to/battenberg_alleles_on_target_hg38.zip",
  "ascat_loci": "/path/to/battenberg_loci_on_target_hg38.zip",
  "ascat_loci_gc": "/path/to/GC_G1000_on_target_hg38.zip",
  "ascat_loci_rt": "/path/to/RT_G1000_on_target_hg38.zip"
}
```

## What are the bwa/bwa-mem2 parameters?

For mapping, sarek follows the parameter suggestions provided in this [paper](https://www.nature.com/articles/s41467-018-06159-4):

`-K 100000000` : for deterministic pipeline results, for more info see [here](https://github.com/CCDG/Pipeline-Standardization/issues/2)

`-Y`: force soft-clipping rather than default hard-clipping of supplementary alignments

In addition, currently the mismatch penalty for reads with tumor status in the sample sheet are mapped with a mismatch penalty of `-B 3`.

## How to manage scatter/gathering (parallelization with-in each sample)

While Nextflow ensures all samples are run in parallel, the pipeline can split input files for each sample into smaller chunks which are processes in parallel.
This speeds up analysis for individual chunks, but might occupy more storage space.

Therefore, the different scatter/gather options can be set by the user:

### Split Fastq files

By default, the input fastq files are split into smaller chunks with FASTP, mapped in parallel, and then merged and duplicate marked. This can be customized by setting the parameter `--split_fastq`.
This parameter determines how many reads are within each split. Setting it to `0` will turn of any splitting and only one mapping process is run per input fastq file.

> FastP creates as many chunks as CPUs are specified (by default 12) and subdivides them further, if the number of reads in a chunk is larger then the value specified in `--split_fastq`. Thus, the parameter `--split_fastq` is an upper bound, e.g. if 1/12th of the Fastq file exceeds the provided value another fastq file will be generated.

### Intervals for Base Quality Score Recalibration and Variantcalling

The pipeline can parallelize base quality score recalibration and variant calling across genomic chunks of roughly similar sizes.
For this, a bed file containing genomic regions of interest is used, it's the intervals file.
By default, the intervals file for WGS used is the one provided by GATK (details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals-)).
When running targeted analysis, it is recommended to use the bed file containing the targeted regions.

The amount of scatter/gathering can be customized by adjusting the parameter `--nucleotides_per_second`.

> **NB:** The _same_ intervals are processed regardless of the number of groups. The number of groups however determines over how many compute nodes the analysis is scattered on.

The default value is `200000`, increasing this value will _reduce_ the number of groups that are processed in parallel.
Generally, smaller numbers of groups (each group has more regions), the slower the processing, and less storage space is consumed.
In particular, in cloud computing setting it is often advisable to reduce the number of groups to be run in parallel to reduce data staging steps.

## How to create a panel-of-normals for Mutect2

For a detailed tutorial on how to create a panel-of-normals, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132).

## Spark related issues

If you have problems running processes that make use of Spark such as `MarkDuplicates`.
You are probably experiencing issues with the limit of open files in your system.
You can check your current limit by typing the following:

```bash
ulimit -n
```

The default limit size is usually 1024 which is quite low to run Spark jobs.
In order to increase the size limit permanently you can:

Edit the file `/etc/security/limits.conf` and add the lines:

```bash
*     soft   nofile  65535
*     hard   nofile  65535
```

Edit the file `/etc/sysctl.conf` and add the line:

```bash
fs.file-max = 65535
```

Edit the file `/etc/sysconfig/docker` and add the new limits to OPTIONS like this:

```bash
OPTIONS=”—default-ulimit nofile=65535:65535"
```

Re-start your session.

Note that the way to increase the open file limit in your system may be slightly different or require additional steps.

### Cannot delete work folder when using docker + Spark

Currently, when running spark-based tools in combination with docker, it is required to set `docker.userEmulation = false`. This can unfortunately causes permission issues when `work/` is being written with root permissions. In case this happens, you might need to configure docker to run without `userEmulation` (see [here](https://github.com/Midnighter/nf-core-adr/blob/main/docs/adr/0008-refrain-from-using-docker-useremulation-in-nextflow.md)).

## How to handle UMIs

Sarek can process UMI-reads, using [fgbio](http://fulcrumgenomics.github.io/fgbio/tools/latest/) tools.

In order to use reads containing UMI tags as your initial input, you need to include `--umi_read_structure [structure]` in your parameters.

This will enable pre-processing of the reads and UMI consensus reads calling, which will then be used to continue the workflow from the mapping steps. For post-UMI processing depending on the experimental setup, duplicate marking and base quality recalibration can be skipped with [`--skip_tools`].

### UMI Read Structure

This parameter is a string, which follows a [convention](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) to describe the structure of the umi.

As an example: if your reads contain a UMI only on the forward read, the string can only represent one structure (i.e. "2M11S+T"); should your reads contain a UMI on both reas, the string will contain two structures separated by a blank space (i.e. "2M11S+T 2M11S+T"); should your reads contain a UMI only on the reverse read, your structure must represent the template only for the forward read and template plus UMI for the reverse read (i.e. +T 12M11S+T). Please do refer to FGBIO documentation for more details, as providing the correct structure is essential and specific to the UMI kit used.

### Limitations and future updates

Recent updates to Samtools have been introduced, which can speed-up performance of fgbio tools used in this workflow.
The current workflow does not handle duplex UMIs (i.e. where opposite strands of a duplex molecule have been tagged with a different UMI), and best practices have been proposed to process this type of data.
Both changes will be implemented in a future release.

## How to run sarek when no(t all) reference files are in igenomes

For common genomes, such as GRCh38 and GRCh37, the pipeline is shipped with (almost) all necessary reference files. However, sometimes it is necessary to use custom references for some or all files:

### No igenomes reference files are used

If none of your required genome files are in igenomes, `--igenomes_ignore` must be set to ignore any igenomes input and `--genome null`. The `fasta` file is the only required input file and must be provided to run the pipeline. All other possible reference file can be provided in addition. For details, see the paramter documentation.

Minimal example for custom genomes:

```bash
nextflow run nf-core/sarek --genome null --igenomes_ignore --fasta <custom.fasta>
```

### Overwrite specific reference files

If you don't want to use some of the provided reference genomes, they can be overwritten by either providing a new file or setting the respective file parameter to `false`, if it should be ignored:

Example for using a custom known indels file:

```bash
nextflow run nf-core/sarek --known_indels <my_known_indels.vcf.gz> --genome GRCh38.GATK
```

Example for not using known indels, but all other provided reference file:

```bash
nextflow run nf-core/sarek --known_indels false --genome GRCh38.GATK
```

### Where do the used reference genomes originate from

For GATK.GRCh38 the links for each reference file and the corresponding processes that use them is listed below. For GATK.GRCh37 the files originate from the same sources:

| File                  | Tools                                                                                                                                                                                                                                                                                                                                                                                                                                                | Origin                                                                                                                | Docs                                                                                 |
| :-------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------- |
| ascat_alleles         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/uouszfktzgoqfy7/G1000_alleles_hg38.zip                                                      | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci            | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_gc         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_rt         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/xlp99uneqh6nh6p/RT_G1000_hg38.zip                                                           | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| bwa                   | bwa-mem                                                                                                                                                                                                                                                                                                                                                                                                                                              | bwa index -p bwa/${fasta.baseName} $fasta                                                                             |                                                                                      |
| bwamem2               | bwa-mem2                                                                                                                                                                                                                                                                                                                                                                                                                                             | bwa-mem2 index -p bwamem2/${fasta} $fasta                                                                             |                                                                                      |
| dragmap               | DragMap                                                                                                                                                                                                                                                                                                                                                                                                                                              | dragen-os --build-hash-table true --ht-reference $fasta --output-directory dragmap                                    |                                                                                      |
| dbsnp                 | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| dbsnp_tbi             | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| dict                  | Baserecalibrator(Spark), CNNScoreVariant, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, MarkDulpicates(Spark), MergeVCFs, Mutect2, Variantrecalibrator                                                                                                                                                                                               | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta                 | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, snpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta_fai             | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, snpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| germline_resource     | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                           | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| germline_resource_tbi | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                           | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| intervals             | ApplyBQSR(Spark), ASCAT, Baserecalibrator(Spark), BCFTools, CNNScoreVariants, ControlFREEC, Deepvariant, FilterVariantTranches, FreeBayes, GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, Strelka, mpileup, MSISensorPro, Mutect2, VCFTools                                                                                                                                                                                                      | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels          | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                       | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels_tbi      | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                       | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps            | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps_tbi        | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |
| mappability           | ControlFREEC                                                                                                                                                                                                                                                                                                                                                                                                                                         | http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip                                                                | http://boevalab.inf.ethz.ch/FREEC/tutorial.html                                      |
| pon                   | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                              | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |
| pon_tbi               | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                              | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |

## How to customise snpeff and vep annotation

### Using the nf-core containers with pre-downloaded cache

For common genomes, it is already configured within the [igenomes.config](https://github.com/nf-core/sarek/blob/master/conf/igenomes.config) file, so nothing to be done there.

Note: These containers are only created for some species and some cache/tools versions combinations (cf DockerHub tags for these containers [`nfcore/snpeff`](https://hub.docker.com/r/nfcore/snpeff/tags) and [`nfcore/vep`](https://hub.docker.com/r/nfcore/vep/tags).

These containers can be quite huge especially for human, it is recommended to use annotation cache on a path if possible

### Create containers with pre-downloaded cache

For each tool, an helper script `build.sh` can be found at the root of the tool folder in the nf-core module repo ([snpeff](https://github.com/nf-core/modules/tree/master/modules/nf-core/snpeff) and [ensemblvep](https://github.com/nf-core/modules/tree/master/modules/nf-core/ensemblvep)), and can be adapted for your usage.

### Use Sarek to download cache and annotate in one go

Use the params `--download_cache`, and specify with `--tools` for which annotation tool you need to download the cache (`snpeff` and or `vep`)

Sarek will automatically download the cache, use the biocontainers container for said tools, and use it to annotate any vcfs produced.

### Only download cache

Using the params `--build_only_index` allow for only downloading the cache for the specified tools.

### Location for the cache

Cache can be downloaded in the specified `--outdir_cache` location.
Else, it will be downloaded in `cache/` in the specified `--outdir` location.

To download cache on a cloud infrastructure, an absolute path is needed.

Params `--snpeff_cache` and `--vep_cache` are to used to specify the locations to the root of the annotation cache folder.

For example this is what can be seen when cache has been downloaded for `GATK.GRCh38` and `WBcell235` for both tools using the default values in the [igenomes.config](https://github.com/nf-core/sarek/blob/master/conf/igenomes.config) file:

```bash
ls /data/snpeff_cache /data/vep_cache/*
/data/snpeff_cache:
GRCh38.105
WBcel235.105

/data/vep_cache/caenorhabditis_elegans:
106_WBcel235
/data/vep_cache/homo_sapiens:
106_GRCh38
```

### Change cache version and species

By default all is specified in the [igenomes.config](https://github.com/nf-core/sarek/blob/master/conf/igenomes.config) file.
Explanation can be found for all params in the documentation:

- [snpeff_db](https://nf-co.re/sarek/latest/parameters#snpeff_db)
- [snpeff_genome](https://nf-co.re/sarek/latest/parameters#snpeff_genome)
- [vep_genome](https://nf-co.re/sarek/latest/parameters#vep_genome)
- [vep_species](https://nf-co.re/sarek/latest/parameters#vep_species)
- [vep_cache_version](https://nf-co.re/sarek/latest/parameters#vep_cache_version)

With the previous example of `GRCh38`, these are the values that were used for these params:

```bash
snpeff_db         = '105'
snpeff_genome     = 'GRCh38'
vep_genome        = 'GRCh38'
vep_species       = 'homo_sapiens'
vep_cache_version = '106'
```

### Usage recommendation with AWS iGenomes

Annotation cache is a resource separated from AWS iGenomes, which as its own structure and a frequent update cycle.
So it is not recommended to put any annotation cache in your local AWS iGenomes folder.

A classical organisation could be:

```bash
/data/igenomes/
/data/cache/ensemblvep
/data/cache/snpeff
```

which can then be used this way in sarek:

```bash
nextflow run nf-core/sarek \\
    --igenomes_base /data/igenomes/ \\
    --snpeff_cache /data/cache/snpeff/ \\
    --vep_cache /data/cache/ensemblvep/ \\
    ...
```

Or similarly on the cloud:

```bash
s3://data/igenomes/
s3://data/cache/ensemblvep
s3://data/cache/snpeff
```

which can then be used this way in sarek:

```bash
nextflow run nf-core/sarek \\
    --igenomes_base s3://data/igenomes/ \\
    --snpeff_cache s3://data/cache/snpeff/ \\
    --vep_cache s3://data/cache/ensemblvep/ \\
    ...
```

These params can be specified in a config file or in a profile using the params scope, or even in a json or a yaml file using the `-params-file` nextflow option.

### Using VEP plugins

#### dbnsfp

Enable with `--vep_dbnsfp`. The following parameters are mandatory:

- `--dbnsfp`, to specify the path to the dbNSFP processed file.
- `--dbnsfp_tbi`, to specify the path to the dbNSFP tabix indexed file.

The following parameters are optionnal:

- `--dbnsfp_consequence`, to filter/limit outputs to a specific effect of the variant.
  - The set of consequence terms is defined by the Sequence Ontology and an overview of those used in VEP can be found [here](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).
  - If one wants to filter using several consequences, then separate those by using '&' (i.e. `--dbnsfp_consequence '3_prime_UTR_variant&intron_variant'`.",
- `--dbnsfp_fields`, to retrieve individual values from the dbNSFP file.
  - The values correspond to the name of the columns in the dbNSFP file and are separated by comma.
  - The column names might differ between the different dbNSFP versions. Please check the Readme.txt file, which is provided with the dbNSFP file, to obtain the correct column names. The Readme file contains also a short description of the provided values and the version of the tools used to generate them.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp).

#### LOFTEE

Enable with `--vep_loftee`.

For more details, see [here](https://github.com/konradjk/loftee).

#### SpliceAi

Enable with `--vep_spliceai`. The following parameters are mandatory:

- `--spliceai_snv`, to specify the path to SpliceAI raw scores snv file.
- `--spliceai_snv_tbi`, to specify the path to SpliceAI raw scores snv tabix indexed file.
- `--spliceai_indel`, to specify the path to SpliceAI raw scores indel file.
- `--spliceai_indel_tbi`, to specify the path to SpliceAI raw scores indel tabix indexed file.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceai).

#### SpliceRegions

Enable with `--vep_spliceregion`.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceregion) and [here](https://www.ensembl.info/2018/10/26/cool-stuff-the-vep-can-do-splice-site-variant-annotation/)."

## Requested resources for the tools

Resource requests are difficult to generalize and are often dependent on input data size. Currently, the number of cpus and memory requested by default were adapted from tests on 5 ICGC paired whole-genome sequencing samples with approximately 40X and 80X depth.
For targeted data analysis, this is overshooting by a lot. In this case resources for each process can be limited by either setting `--max_memory` and `-max_cpus` or tailoring the request by process name as described [here](#resource-requests). If you are using sarek for a certain data type regulary, and would like to make these requests available to others on your system, an institution-specific, pipeline-specific config file can be added [here](https://github.com/nf-core/configs/tree/master/conf/pipeline/sarek).

## MultiQC related issues

### Plots for SnpEff are missing

When plots are missing, it is possible that the fasta and the custom SnpEff database are not matching https://pcingola.github.io/SnpEff/se_faq/#error_chromosome_not_found-details.
The SnpEff completes without throwing an error causing nextflow to complete successfully. An indication for the error are these lines in the `.command` files:

```text
ERRORS: Some errors were detected
Error type      Number of errors
ERROR_CHROMOSOME_NOT_FOUND      17522411
```

## How to set up sarek to use sentieon

Sarek is currently not supporting sentieon. It is planned for the upcoming release 3.3. In the meantime, please revert to the last release 2.7.2.
