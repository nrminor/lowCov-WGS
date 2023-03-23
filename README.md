# lowCov-WGS
A NextFlow Pipeline inspired by [Lou et al. 2021 (Molecular Ecology)](https://onlinelibrary.wiley.com/doi/10.1111/mec.16077). Its purpose is to collate and streamline tools used to analyze WGS data with low sequencing effort per sample.

_DISCLAIMER: This pipeline is in the early stages of development and should not be used by other researchers._

### What this pipeline will run:
The pipeline starts with routine bioinformatic processing and QC using open-source tools, which include:
- Merging paired reads with `bbmerge.sh`
- Orienting reads with `vsearch`
- Filtering low quality reads and trimming adapters with `fastp` and `Trimmomatic`
- Generating read quality reports with `FASTQC` and `MultiQC`
- Mapping reads to a reference sequence of your choice with `bwa mem`.
- Assessing depth with `mosdepth` and plotting it per-sample and per-reference-scaffold with `R ggplot2`
- Inferring genotype likelihoods for each species with `ANGSD`
- Calling variants for each sample with `bcftools`, filtering them to reliable SNPs with `vcftools`, and merging them per-species and library-prep-type with `bcftools`

From there, the pipeline can perform a number of what Lou et al. 2021 refer to as individual-level analyses, most of which skew toward using genotype-likelihoods instead of hard variant calls. These include:
- Population structure with `STRUCTURE`
- Nonparametric population structure with `PCANGSD` and `FASTPCA`
- Selection scan with `Ohana`
- GWAS with `ANGSD_GWAS`
- Linkage disequilibrium scan with `NGSLD`
- Runs of homozygosity scan with `bcftools roh`

The pipeline can also perform so-called population level analyses, again mostly with genotype likelihoods, including:
- Allele frequency assessment with `ANGSD -af`
- Site frequency spectrum (SFS) inference with `ANGSD -sfs` and visualization with `R`
- Demographic history reconstruction with `stairwayplot`
- Fst estimation with `ANGSD FST` and `NGSTOOLS FST`
- Fst probability assessment with `vcflib`

All of these tools will be prepackaged and versioned in a Docker image associated with this workflow, which means the only tools you will need to install yourself are Git, NextFlow, Java, and either Docker or Singularity.

### Example invocation
```
nextflow run lowCov-WGS.nf \
--samplesheet resources/nrm_test.csv \
--low_disk_mode true \
--reference resources/GCA_014549065.1_CCar_1.0_genomic.fna
```