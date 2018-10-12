# METAGenOmic Analysis Workflow CLI

A CLI that aims to automate common metagenome pipeline steps, with a subsequent focus on viral downstream analysis.  The base steps include read quality control, assembly, binning, and gene calling, and should work for most any metagenomic analysis with a little tweaking of the parameters.

# Workflow Overview

## General Steps and Opinionated Software List

1. Quality Control: [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2. Assembly: [MegaHit](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)
3. Binning: [Maxbin 2.0](http://sourceforge.net/projects/maxbin/)
4. Gene Calling: [Prodigal](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) | [EMIRGE (ribosomal reconstruction)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44)
5. Viral Analysis: [VirSorter](https://peerj.com/articles/985/?utm_source=TrendMD&utm_campaign=PeerJ_TrendMD_0&utm_medium=TrendMD)
6. Annotation: [Prokka](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)

## Resources

1. Micro-Phage Interaction Database: [MVP](http://mvp.medgenius.info/home)

# Setup
The following steps assume you have a python and pip install.  If you don't and potentially are using a server with unprivileged access, you may want to consider a user install of Anaconda.  Some rationale and steps to do this [can be found here](https://github.com/DuhaimeLab/Virus_Knowledgebase/wiki/Flux-Helper#personal-environment) in the Personal Environment section.

## Github Install
```
$ git clone https://github.com/jordangumm/metago.git
$ cd metago && ./build.sh
$ pip install -e .
```
## Singularity Install

> In Development!

# Usage

## Base Command
The `metago` command is your interface to a myriad of workflow commands.  It requires fastq files to be organized in Illumina fashion, that is in the form of `Run_[RUNID]/Sample_[SAMPLEID]/[SAMPLEID].fastq`.  Run-based commands target a run directory and will process every sample automatically.  Sample-based commands target a single sample or fastq file.  The quality control step interleaves fastqs, so ensure you run your fastq sample files through that step first if you want to be able to leverage downstream commands.

`$ metago --help`
```
Usage: metago [OPTIONS] COMMAND [ARGS]...

  Metago Command Line Interface

  Note: Use absolute paths to all files and directories

Options:
  -o, --output TEXT
  --flux / --no-flux
  -a, --account TEXT
  -p, --ppn INTEGER
  -m, --mem TEXT
  -w, --walltime TEXT
  --help               Show this message and exit.

Commands:
  run_assembly        Assemble run sample reads
  run_mapping         Read mapping of run to reference
  run_minhash         Minhash compare fastqs in run
  run_pseudoalign     Read assignment of run to reference
  run_qc              Quality control of Illumina run
  sample_assembly     Assembly sample reads
  sample_mapping      Read mapping of reference to sample
  sample_pseudoalign  Read assignment of sample to reference
  sample_qc           Quality control of Illumina sample
```

## Secondary Commands
The workflow commands require their own arguments and options.

`$ metago sample_qc --help`
```
Usage: metago sample_qc [OPTIONS] SAMPLE_DP

  Quality control of Illumina sample

Options:
  --help  Show this message and exit.
```

Make sure to provide arguments and options at the appropriate command level.  `metago` expects resource and output information, while secondary commands require more data-specific information, like what fastq to process or analyze.  Below is an example command that specifies an output path and 4 hour walltime limit for quality controlling a sample.

`$ metago -o /scratch/analysis/Sample_1234 -w 4:00:00 sample_qc /nfs/longterm_storage/Sample_1234`
