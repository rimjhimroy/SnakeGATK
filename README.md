# Snakemake workflow: dna-seq-gatk-variant-calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.7.1-brightgreen.svg)](https://snakemake.readthedocs.io)

This Snakemake pipeline implements the [GATK best-practices workflow](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) for calling small genomic variants from already created bam files.

## Author

* Rimjhim Roy Choudhury (https://rimjhimroy.github.io)

This workflow is a **modified fork** of the [dna-seq-gatk-variant-calling](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling) workflow by [Johannes Köster](https://koesterlab.github.io)  

## Usage

#### Step 1: Obtain a copy of this workflow

1. Fork a the repository.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 3: Install conda and Snakemake

```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh 
conda config --set auto_activate_base False
conda update conda
conda create -n snakemake -c bioconda snakemake=5.7.4
conda activate snakemake
```


#### Step 3: Configure workflow

1. Edit the samples.tsv and units.tsv to accomodate the provided bam files.

2. Download ***recaliberation files***:  
    • [dbsnp](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz) as set of known variants.  
    • Also got: [hapmap](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz); [1000G_omni](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz) and [1000G_phase1](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz) and [Mills_and_1000G_gold_standard.indels](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz) for VQSR.
    • Put all the above files in `ref/`folder.

3. Please run initialize rule first

    snakemake --use-conda -s rules/index.smk  


**NOTE**: I was not able to perforrm vqsr on indels. It gave me error "java.lang.IllegalArgumentException: No data found." So I went for hardfiltering the indels. But the vqsr rule for indell is provided in `rules/vqsr_indels.smk` and can be plugged in after a few tweaks. 

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

Execute in a cluster with SLURM

    sbatch snakemake.sh

## Map of workflow

<img src="dag/dag.svg" alt="Directed acyclic graph (DAG) of the workflow"  width="4200" height="3000">
