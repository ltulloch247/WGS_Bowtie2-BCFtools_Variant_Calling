# **Whole-genome sequence data processing and variant calling pipeline**



## **Overview**



This repository contains a bioinformatics pipeline for processing whole-genome sequencing (WGS) data, performing copy number variation (CNV) analysis and identifying candidate genetic variants.

The pipeline is designed to be applicable to any organism with a suitable reference genome and corresponding annotation files. Although developed and initially applied to *Leishmania* drug-resistance studies, it is not restricted to this organism.



The pipeline performs the following analyses:

1\. **Read alignment:** Aligns paired-end WGS reads to a reference genome using Bowtie2 and generates BAM files.

2\. **Copy number variation analysis:** Estimates gene and chromosome copy numbers from the aligned sequencing reads.

3\. **Variant calling:** Uses BCFtools mpileup to identify genetic variants and generate VCF files.

4\. **Variant filtering and annotation:** Applies quality filters, including parental-line filtering based on QUAL, DP and GQ, and annotates variants using snpEff.

5\. **Results output:** Produces Excel files containing variant and CNV data for further filtering, analysis and presentation.



The pipeline includes detailed instructions within the script covering input file organisation, configuration, execution and interpretation of the output files.



## **Workflow**



1\. Bowtie2 alignment

2\. CNV analysis

3\. BCFtools variant calling

4\. snpEff annotation and data filtering

5\. Excel output and results presentation



## **Requirements (input files)**



1\. WGS FASTQ.gz data

2\. Reference genome (FASTA and GFF annotations)

3\. Variant calling script (run\_variants\_BCFtools\_v1.sh - included in this repository)

4\. Excel templates:

&#x20;   - Genomes.xlsx - reference genome annotation template, included in this repository

&#x20;   - Variants.Template\_v10.2.1.xlsx - template file that results are written into, included in this repository
5\. A conda environment (named bcftools_env) is required, along with several packages - the script contains additional information about these dependencies 



## **Output**



When run correctly, the script will write results into the variants template file, creating an Excel output named according to the experiment, reference genome and data contained.



Data output to the Excel results file:

1. **FASTQ\_metadata** - WGS filename and MD5 checksums to unambiguously show which data were processed. Simplifies deposition of WGS to nucleotide repository.
2. **Genome** - reference genome GFF annotation derived from the corresponding GFF file provided in Excel format.
3. **Bam\_QC** - information on the number of sequence reads per sample, % mapped to the reference genome, and mean coverage.
4. **Gene\_CNV** - gene and chromosome copy number variation analysis.
5. **BCFtools** - variant calls, quality scores and annotations.
6. **Run\_log** - a copy of the master log file, retaining information on script version, configuration and execution.



## **Acknowledgements**



This pipeline was created at the University of Dundee thanks to funding from the Wellcome Trust

