# Simulating populations and short read data

The purpose of this pipeline is to simulate data for testing imputation methods on low-coverage sequencing (LCS) data, such as STITCH, GLIMPSE, or QUILT. The pipeline performs SLIM4 simulation, and uses the simulated genotypes to simulate Illumina reads at high and low coverage that would be used as reference panel/ground truth and target, respectively, for imputation.

This Snakemake pipeline runs entirely on docker image and conda environment.

## Table of content

- [SLIM4 population simulation](#slim4-population-simulation)
- [Simulation of low-coverage sequencing data](#simulation-of-low-coverage-sequencing-data)
- [Simulation of high-coverage sequencing data](#simulation-of-high-coverage-sequencing-data)
- [Software requirements](#software-requirements)
- [Software included in the pipeline](#software-included-in-the-pipeline)
- [Run the pipeline](#run-the-pipeline)

## SLiM4 population simulation 

The only input file of this pipeline is a SLiM recipe. A template for SLiM recipe is [here](./test/test.slim). Model settings and events can be modified at dispersal, however output file of *outputVCF()* or *outputVCFSample()* must be blank as they will be printed in Snakemake output file system. Output of population simulation step includes a VCF file and the reference sequence (ancestral sequence) in a fasta file. VCF data will be converted into single fasta files for each individual, by [vcf2fasta](https://github.com/vcflib/vcflib/blob/master/doc/vcf2fasta.md), for read simulation.

## Simulation of low-coverage sequencing data

There are three main steps to simulate LCS data.

1. Read simulation performed by [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art). Sequencing platform, depth, read lengths, etc. can be modified in the config file(a template for config file is [here](config.yaml)).

2. Fastq files as output of ART will be mapped against the reference sequence (output of SLiM4), by bwa mem. SAM files are converted into BAM files by samtools view, with filtering parameters to be defined by the user in the config file, then sorted by samtools sort.

3. Duplicated reads are marked by picard tools.

## Simulation of high-coverage sequencing data

HCS data will be used as reference panel and/or ground truth for imputation. Simulation steps include the three steps as described for simulation of LCS data, followed by the following steps:

4. Variant calling, performed by bcftools mpileup. Only biallelic SNPs will be retained.

5. SNP filtering by bcftools filter. Filtering parameters are specified on the config file.

## Software requirements
- Snakemake (>= 8.10.0, recommended)
- Miniconda3/conda
- Singularity >= 3.0.0

## Software included in the pipeline
- SLiM 4.3
- samtools and bcftools 1.21
- ART-MountRainier (lastest)
- bwa (0.7.18)
- picard (3.3.0)
- vcflib (lastest)

## Run the pipeline

To run the pipeline:
1. Prepare a SLiM recipe file, which simulates a population and output a VCF file (see this [example](./test/test.slim))
2. Prepare a config file with this [template](config.yaml). Path to the SLiM recipe, and parameters for [ART read simulation](https://www.niehs.nih.gov/research/resources/software/biostatistics/art), and [samtools view](https://www.htslib.org/doc/samtools-view.html) need to be specify.
3. Run snakemake: `snakemake --sdm apptainer conda --configfile [config file]`

After downloading the pipeline, a test run can be executed with the existing SLiM recipe and config file with the following commands.
```
git clone https://github.com/vibaotram/simulation.git
cd simulation
snakemake --sdm apptainer conda --configfile config/config.yaml
```