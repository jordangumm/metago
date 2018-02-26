# Virus Omics Analysis Workflow

A workflow that aims to automate common viral metagenome pipeline steps.  This includes read quality control, assembly, gene calling, and other basic analysis routines.

# Workflow Overview

## General Steps

1. Quality Control: [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
2. Assembly: [MegaHit](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)
3. Gene Calling: [Prodigal](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119) | [EMIRGE (ribosomal reconstruction)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44)
4. Viral Analysis: [VirSorter](https://peerj.com/articles/985/?utm_source=TrendMD&utm_campaign=PeerJ_TrendMD_0&utm_medium=TrendMD) | [VirHostMatcher](https://academic.oup.com/nar/article/45/1/39/2605663)
5. Annotation: [Prokka](https://academic.oup.com/bioinformatics/article/30/14/2068/2390517)

## Resources

1. Micro-Phage Interaction Database: [MVP](http://mvp.medgenius.info/home)
