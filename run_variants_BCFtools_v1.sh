#!/bin/bash

######################################################################################################################################################################
SCRIPT=run_variants_BCFtools_v1.sh       # If modifying this script name, also change the name in the command to run the script and log below!			     #
#                                                                                                                                                                    #
# Overview                                                                                                                                                           #
#                                                                                                                                                                    #
# This script aligns WGS data to a reference genome, performs copy number variation analysis and variant calling.						     #
# If using this script, please cite xxxxxxx															     #
#                                                                                                                                                                    #
# User configuration begins on line 134 - scroll down to update data if everything else is already setup.							     #
#                                                                                                                                                                    #
# Detailed overview 																		     #
#                                                                                                                                                                    #
#   Part 1) Aligns WGS reads to a reference genome with Bowti2, outputting binary alignment (BAM) files;                                                             #
#   Part 2) Calculates copy number variation (CNV) using the Bowtie2-generated BAM files;                                                                            #
#   Part 3) Performs variant calling with BCFtools mpileup using the Bowtie2-generated BAM files, outputting a variant call format (VCF) file;                       #
#   Part 4) Filters poor quality parental calls from the VCF, based on read depth (DP), genotype qualtiy (GQ) and overall quality (QUAL).                            #
#     - results are output to Excel sheet ${PROJECT}.${GENOME}.Bowtie2.CNV.BCF.xlsx for additional filtering and analysis.                                           #
#                                                                                                                                                                    #
#  Note that filters in part 4 of this script only filter out variants when the QUAL, DP and GQ in the parental line are below the cutoff:			     #
#     - this prevents variants with poor scores in some lines from 'hiding' the same real mutations in other lines;						     #
#     - this prevents gene deletions in some lines from 'hiding' the real variants in other lines;								     #
#     - the Excel results sheet allows variants to be inspected and filters can be made more stringent and include mutant filtering there.                           #
#     - results are output to Excel sheet ${PROJECT}.${GENOME}.OVarFlow.CNV.OV.BCF.DEEP.xlsx for additional filtering and analysis.                                  #
#                                                                                                                                                                    #
# Files required to run this script (and file naming rules):													     #
#                                                                                                                                                                    #
#   1) Paired-end WGS read files (FASTQ files)															     #
#	- each sample must have only ONE paired set named *_R1.fastq.gz and *_R2.fastq.gz, where *=sample name)							     #
#            - rename FASTQ prefix to match the BGI sample name (eg, LdRES-6690-WT / LdRES-6690-R1 etc)   				                             #
#            - replace _1 and _2 (for paired read 1 and paired erad 2) to _R1 and _R2                   						             #
#            - rename extension from .fq.gz to .fastq.gz    													     #
#            - examples of suitable names (best practise naming convention: SpeciesExperiment-Drug/Gene-Line_R1/2.fastq.gz):					     #					     
#		   - LdRES-6690-WT_R1.fastq.gz    /    LdRES-6690-WT_R2.fastq.gz	(L. donovani Resistance to 6690	- WT)					     #
#		   - LdRES-6690-R1_R1.fastq.gz    /    LdRES-6690-R1_R2.fastq.gz	(L. donovani Resistance to 6690	- Resistant line 1)			     #
#		   - LdRES-6690-R3_R1.fastq.gz    /    LdRES-6690-R3_R2.fastq.gz	(L. donovani Resistance to 6690	- Resistant line 3)			     #
#		   - LdXPER-DLP-WT_R1.fastq.gz    /    LdXPER-DLP-WT_R2.fastq.gz	(L. donovani CRISPR edit of DLP	- WT)					     #
#		   - LdXPER-DLP-DKO_R1.fastq.gz   /    LdXPER-DLP-DKO_R2.fastq.gz	(L. donovani CRISPR edit of DLP	- DLP DKO)				     #
#		   - LdXPER-DLP-G606V_R1.fastq.gz /    LdXPER-DLP-G606V__R2.fastq.gz	(L. donovani CRISPR edit of DLP	- DLP G606V)				     #
#            - if more than one pair exist per sample, concentate _R1 reads to single file and concentate _R2 reads to single file    	                             #
#            - if a directory only contains FASTQ reads from one sample, all of the above can be achieved with the following:    	                             #
#                   cat *_1.fq.gz > LdRES-6690-WT_R1.fastq.gz && cat *_2.fq.gz > LdRES-6690-WT_R2.fastq.gz      # but obviously change the prefix!                   #
#            - All FASTQ files must be in the FASTQ_INPUT_DIR directroy with no subfolders					    	                             #
#                                                                                                                                                                    #
#   2) Reference genome:																	     #
#	- The reference genome (fasta format) must be acompanied by the corresponding annotations file (gff format)						     #
#	- Can be readily downloaded,for example from TriTrypdb													     #
#       - Reference fasta:																	     #
#            - rename .fasta extension to .fa				(rename with mv *.fasta *.fa)								     #
#            - a compressed fasta (*.fa.gz) is also required		(compress with "bgzip -c *.fa > *.fa.gz")						     #
#       - Reference gff:																	     #
#            - a sorted gff is required with sequential annotation	(sort with "(grep '^#' *.gff; grep -v '^#' *.gff | sort -k1,1 -k4,4n) > *.sorted.gff")	     #
#       - All four reference genome files must be in REFERENCE_INPUT_DIR 											     #
#            - this script copies the files from my reference genome folder to REFERENCE_INPUT_DIR (you can change this)					     #
#                                                                                                                                                                    #
#   3) Excel file containing sorted annotated GFFs 			(Genomes.xls and located in EXCEL_TEMPLATES)						     #                                                                                                                                                                    
#   4) Excel template into which all results will be inserted 		(Variants.Template_v10.2.xls and located in EXCEL_TEMPLATES)				     #
#                                                                                                                                                                    #
# Files and folders map																		     #
#										     										     #
#   INPUT_PATH					# Folder where the this script is run from									     #
#	|										     									     #
#	|	Input files/Folders																     #
#	|										    									     #
#	|--- run_variants_BCFtools_v1.sh	# This script to run alignment, CNV analysis and variant calling						     #
# 	|--- FASTQ_INPUT_DIR			# Folder containing *_R1.fastq.gz and *_R2.fastq.gz files for every sample (without subfolders)			     #
# 	|--- REFERENCE_INPUT_DIR		# Folder to which reference genome will be copied								     #
# 	|--- OLD_GVCF_FILES			# Leave empty													     #
#	|--- EXCEL_TEMPLATES			# Folder containing Excel templates (Genomes.xlsx and Variants.Template_v10.2.xlsx)				     #
#       |													 						     #
#       |																			     #
# 	|	Output files/folders																     #
# 	|																			     #
# 	|---${PROJECT}.${GENOME}																     #
# 	 	|													 					     #
# 		|--- Bowtie2_BAM					# Folder containing Bowtie2-generated BAM files.					     #
# 		|--- Bowtie2_CNV					# Folder for CNV calculation using Bowtie2-aligned data.				     #
# 		|--- Bowtie2_VCF					# Folder containing raw and filtered VCFs from Bowtie2-aligned data.    		     #
# 		|--- ${PROJECT}.${GENOME}.Bowtie.CNV_BCF.xlsx		# Bowtie2/BCFtools final results Excel file.						     #                                                                                                                                                        
# 		|																		     #
# 	 	|--- Excel_out 						# Folder with all Excel outputs. Not needed once run is complete			     #
# 	 	|--- Genome						# Folder with reference genome, to be copied over with final results			     #
# 	 	|--- LOGS 						# Folder for log files. Not needed once run is complete					     #
#																				     #
#																				     #
#      WGS					# ${WGS_PATH} - Folder where I have reference genome, singularity images and snpeff installed			     #
#	|																			     #
# 	|--- genomes/${GENOME%_v68}/${REF}	# Folder with REF.fa.gz and REF.sorted.gff 	# Note,"%_v68" means without the version extension		     #
#	|--- snpEff				# Folder for SNP Effects (SnpEff installed here)								     #
#		|---- snpEff.jar		#														     #
#		|---- snpEff.config		# Your reference genme will be added to this file by this script						     #
#		|---- data/${GENOME}		# genes.gff, sequences.fa and snpEffectPredictor.bin will be added here automatically				     #
#                                                                                                                                                                    #
# Packages required to run this script:																     #
#                                                                                                                                                                    #
#   To run this script you need a conda environment called "bcftools_env" containing the following packages: 						             #
#        channels:                                                                                                                                                   #
#         - bioconda                                                                                                                                                 #
#          - conda-forge                                                                                                                                             #
#          - defaults                                                                                                                                                #
#        dependencies:                                                                                                                                               #
#          - bcftools v1.23                                                                                                                                          #
#          - samtools v1.23                                                                                                                                          #
#          - subread v2.1.1                                                                                                                                          #
#          - python v3.11.15                                                                                                                                         #
#          - pandas v3.0.1                                                                                                                                           #
#          - numpy v2.4.2                                                                                                                                            #
#          - openpyxl v3.1.5                                                                                                                                         #
#          - snpeff v5.1d                                                                                                                                            #
#          - bowtie v2.5.5                                                                                                                                           # 
#                                                                                                                                                                    #
# To use this script:                                                                                                                                                #
#   1) Copy FASTQ data, Reference genome and Excel templates to the folders shown in the map above.								     #
#   2) Update info in the User configuration section below (from line 131):											     #
#          - If setting up fo the first time, enter all info													     #
#          - If running with new data, you only need to update:													     #
#                 1) Project name (eg, SpeciesRES-Drug). All outputs will be prefixed by this.									     #
#                 2) Samples																	     #
#                 3) Reference genome																     #
#                 4) Parts to run																     #									
#																				     #
# Submit to cluster with:																	     #
#     "sbatch [cluster-specific parameters] --output=run_variants_BCFtools_v1.log run_variants_BCFtools_v1.sh" (slurm)						     #
# Check the script is running with "squeue --me"     														     #
# Alternatively, run in your shell with ./run_variants_BCFtools_v1.sh (but not recommended as it takes a long time and will likely timeout).			     #                                                                                                  
# Cancel a slurm job with "scancel "[job_number]                                                                                                                     #
#                                                                                                                                                                    #
# Lindsay Tulloch 09-September-2026                                                                                                                                  #
######################################################################################################################################################################


########################################################
##               USER Configuration 1                 ##
##          Define the following EVERY run            ##
########################################################

# ---- Parent folder where the script will run from ----
WGS_PATH="[enter the path to your WGS folder here]"
INPUT_PATH="${WGS_PATH}/processing/variant_calling"		# This is where your script will run, raplace [path] for the actusl path


# ---- Project Name ---- (Usage: SpRES-Drug). Results will be output to a folder with this name ----
PROJECT="LdRES-DNDI6690"	# usage: speciesEXPERIMENT-drug # examples: LdRES-DNDI6690 (L. donovani RESISTANCE GENERATION - DNDI-6690) or LdCRI-DLP (L. donovani RESISTANCE GENERATION - Dynamin-Like Protein)


# ---- Samples ----
PARENT="LdRES-6690-WT"	# Put the parental line here (WT, or the parent that has been mutated)
MUT1="LdRES-6690-R1"	#	
MUT2="LdRES-6690-R2"	#	
MUT3="LdRES-6690-R3"	#	
MUT4="LdRES-6690-R4"	#	
MUT5="None"	#	
MUT6="None"	# This will not be processed in the final Excel sheet (V10.1 and earlier)
#
# Notes: 1) Enter the sample names to match the prefixes of your _R1.fastq.gz file names (ideally use the sample names submitted to BGI, unless theya re not informative)
#        2) Leave unused samples as "None"
#        3) The maximum number of cell lines that can be processed is 7; Excel sheet V10.1 and earlier only really processes 6
#        4) This sample info (as well as reference genome info below) will be used to create the "samples_and_read_groups.csv" so you don't need to (reduces errors).


# ---- FASTQ directory ---- # Keep these in the folder that OVarFlow requires them to be in so as not to duplicate (they are very large files) 
FASTQ_DIR="$INPUT_PATH/FASTQ_INPUT_DIR"



# ---- Reference genome (uncomment ONE) ----
REF="${WGS_PATH}/genomes/LdBPK282A1/TriTrypDB-68_LdonovaniBPK282A1_GK" GENOME="LdBPK282A1_v68"
# Any reference genome can be used here, remember to update the GENOME parameter to match the tab name in the Geonomes.xls Excel template


# ---- Ploidy  (uncomment ONE) ----
# PLOIDY=1    # Use for haploid cells (eg, Plasmodium, Cryptococcus) # But actually use 2 for now while testing
PLOIDY=2    # Use for diploid cells (eg, Leishmania, Trypanosoma)


# ---- Options to skip parts (for troubleshooting) ----
SKIP_PART_0="No"	# "Yes" will skip part 0 (Preparation of files for the pipeline)
SKIP_PART_1=""	# "Yes" will skip part 1 (WGS alignment with Bowtie2 and Bowtie2_BAM generation)
SKIP_PART_2=""	# "Yes" will skip part 2 (CNV calculation from Bowtie2-generated BAMs and output to Excel)
SKIP_PART_3=""	# "Yes" will skip part 3 (Variant calling of Bowtie2-generated BAMs with BCFtools)
SKIP_PART_4=""	# "Yes" will skip part 4 (Filtering Bowtie2/BCFtools VCF and output to Excel)
# Note Parts 1-4 must be performed sequentially to generate a final Bowtie2/BCFtools output. Once all parts have been completed, earlier parts can be skipped, for example if changing filtering parameters.


# ---- VCF filtering options (select filtered VCF to be written to Excel output) ----
# FINVCF=5	# VCF after QUAL, DP and GQ filtering (minimum filtering to remove poor calls), removal of 'variants' that are identical in all lines including the WT/parental line, and removal of WT≠0/0 variants
FINVCF=6	# above VCF containing only variants inside mRNA	# Reccommended for polycisronic organisms only


########################################################
##               USER Configuration 2                 ##
##    Define the following the first time you run     ##
##  But you may wish to alter QUAL/DP/GQ and re-run   ##
########################################################


# ---- Excel files (templates used for the final output) ----
EXCEL_TEMP="$INPUT_PATH/EXCEL_TEMPLATES"								# Location of all Excel input templates required for this script
EXCEL_TEMPLATE="${EXCEL_TEMP}/Variants.Template_v10.2.1.xlsx"						# Template file to which FASTQ info, reference genome, Bam_QC, CNV and VCF data wil be added
EXCEL_GENOMES="${EXCEL_TEMP}/Genomes.xlsx"								# Genome annotation Excel file from which to copy the genome to the above

# ---- snpEff paths (used for mutation annotation) ----
SNPEFF_HOME="${WGS_PATH}/snpEff/"									# Location of the snpEff installation
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"								# location of the snpEff.config (with your genome - you may have to add your genom eto the config file)
SNPEFF_DATA="${SNPEFF_HOME}/data"									# snpEff data (gff, fa and bin) must be located here, ie, ${SNPEFF_HOME}/data/PfDd2_v68/xxx.fa (etc)
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"									# snpEff.jar must be located here

# ---- Bowtie2 settings ----
BOWTIE2_THREADS=40

# ---- BCFtools settings ----
CALL_MODE=cv		# use "cv" for consensus caller or "mv" for multiallelic caller 		# Recommend using cv
MIN_QUAL_BCF=30	    	# Min QUAL for subsequent filtering of BCFtools VCF				# Recommend using 30 here, but set to 100 in the Excel sheet
MIN_DEPTH_BCF=20	# Min Parent Read Depth (DP) for subsequent filtering of BCFtools VCF		# Recommend using 20
MIN_GQ_BCF=70		# Min Parent Genotype Quality (GQ) for filtering of BCFtools VCF		# Recommend using 99 here to ensure the parental call is correct
# Note that filters are only applied to the parental line in this script. Filters can be increased and applied to mutants in the Excel output.

# ---- Number of threads (for parallel processing on cluster)
THREADS=40


############################################################
##   ** No need to modify anything beyond this point **   ##
##         as these paths float with INPUT_PATH           ##
############################################################



############################################################
##               Paths, files and destinations            ##
############################################################

# ---- Reference genome and samples ----
REFERENCE_FASTA="${REF}.fa"
REFERENCE_GFF="${REF}.sorted.gff"


# ---- Folders ----
PROJECT_OUT="$INPUT_PATH/${PROJECT}.${GENOME}"							# Folder to contain all results of data analysis for this project & genome
#
LOGS_OUT="$PROJECT_OUT/LOGS"									# Directory for log files from 
SNPeff_DIR="$PROJECT_OUT/SNPeff"								# Directory where the effects of SNPs will be processed
#
BOW_BAM_OUT="$PROJECT_OUT/Bowtie2_BAM"								# Directory for Bowtie2-aligned BAM files
BOW_VCF_OUT="$PROJECT_OUT/Bowtie2_VCF"								# Directory for  VCF files from Bowtie2 alignment
BOW_CNV_OUT="$PROJECT_OUT/Bowtie2_CNV"								# Directory for  CNV processing from Bowtie2 alignment
#
EXCEL_OUT="$PROJECT_OUT/Excel_output"								# Directory where Excel results outputs will be put
GENOME_DIR="$PROJECT_OUT/Genome"								# Directory where reference genome fasta and gff will be copied to for reference


# ---- Log files ----
BOWTIE_BUILD_LOG="${LOGS_OUT}/bowtie2_build.log"						# Log file for Bowtie2 build
MPILEUP_BOW_LOG="${LOGS_OUT}/mpileup_Bowtie2.log"						# Log file for Mpileup with Bowtie2-aligned BAMs (Part 3)


# ---- Bowtie2 output ----
BOWTIE2_INDEX="${BOW_BAM_OUT}/${GENOME}_bowtie2"


# ---- List of CSV files generated ----
BAM_LIST_BOW="${BOW_BAM_OUT}/Bowtie2_Bam_list.csv"						# List of BAM files (with path) output by Bowtie12
BAM_QC_REPORT_BOW="$BOW_BAM_OUT/bam_qc_summary.csv"						# BAM QC report from Bowtie2-aligned reads


# ---- Bowtie2/BCFtools VCF files (created in this order) ----
BOW_BCF_RAW_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.0.Raw.vcf.gz"		# Bowtie2-aligned/BCFtools-called raw VCF output, copied to VCF output directory
BOW_BCF_ANNOTATED_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.1.Annotated.vcf.gz"	# above VCF after mutation annotation in step 5.1
BOW_BCF_REORDERED_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.2.Reordered.vcf.gz"	# above VCF with samples reordered so WT comes first, created in step 5.2
BOW_BCF_FILTERED_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.3.Filtered.vcf.gz"	# above VCF after QUAL, DP and GQ filtering in step 5.3
BOW_BCF_DIFF_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.4.Different.vcf.gz"	# above VCF where >=1 cell line must be different to the others (including WT), created in step 5.4
BOW_BCF_DIFF_WT0_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.5.WT00.vcf.gz"	# above VCF where WT=0/0, created in step 5.5
BOW_BCF_MRNA_VCF="${BOW_VCF_OUT}/${PROJECT}.${GENOME}.Bowtie_BCFtools.6.mRNA.vcf.gz"		# above VCF containing only variants inside mRNA, created in step 2.6


# ---- Excel files output (in this order) ----
EXCEL_FASTQ="${EXCEL_OUT}/${PROJECT}.FASTQ.xlsx"						# Excel template with FASTQ data added to 'FASTQ_metadata' tab
EXCEL_GENOME="${EXCEL_OUT}/${PROJECT}.${GENOME}.xlsx"						# Above with genome added to 'Genome' tab
#
EXCEL_BOW_BAMqc="${EXCEL_OUT}/${PROJECT}.${GENOME}.Bowtie.Bamqc.xlsx"				# 'FASTQ/Genome' Excel sheet (above) with Bowtie2 BAM-QC added to 'Bam_QC' tab
EXCEL_BOW_CNV="${EXCEL_OUT}/${PROJECT}.${GENOME}.Bowtie.CNV.xlsx"				# Above with RPKM data added to 'CNV' tab
EXCEL_BOW_BCF="${EXCEL_OUT}/${PROJECT}.${GENOME}.Bowtie.CNV_BCF.xlsx"				# Above with BCFtools calls data added to 'BCFtools' tab


# ---- snpEff files (input, intermediate and output) ----
SNPEFF_GENOME="$GENOME"
SNPEFF_HTML="${SNPeff_DIR}/snpEff_summary.html"
SNPEFF_LOG="${LOGS_OUT}/SnpEff_build.log"


############################################################
##                     CLUSTER SETTINGS                   ##
############################################################

#$ -adds l_hard local_free 200G
#$ -mods l_hard m_mem_free 30G
#$ -adds l_hard avx 1
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 40
#$ -N Variants
#$ -o Variants_errors_$JOB_ID

set -euo pipefail


############################################################
##                  Log of user settings                  ##
############################################################

echo " 
***********************************************************************************************************
Master log file for $PROJECT
***********************************************************************************************************

----- $(date) : Running $SCRIPT
                                     Job ID : $SLURM_JOB_ID
                                     Running in : $INPUT_PATH
                                     Results folder : $PROJECT_OUT
                                     From user configuration:
                                        Project name : $PROJECT
                                        Samples : $PARENT $MUT1 $MUT2 $MUT3 $MUT4 $MUT5 $MUT6
                                        Reference genome : $REF
                                     If details are incorrect, kill the job with *scancel $SLURM_JOB_ID*
                                     Parts to perform:"
if [[ "$SKIP_PART_0" != "Yes" ]]; then
     echo "                                        0. Prepare files for the pipeline"
fi
if [[ "$SKIP_PART_1" != "Yes" ]]; then
     echo "                                        1. WGS alignment and variant calling with Bowtie2"
fi
if [[ "$SKIP_PART_2" != "Yes" ]]; then
     echo "                                        2. Gene and chromosome CNV calculation from Bowtie2-aligned reads and output to Excel"
fi
if [[ "$SKIP_PART_3" != "Yes" ]]; then
     echo "                                        3. Variant calling of Bowtie2-aligned reads with BCFtools"
fi
if [[ "$SKIP_PART_4" != "Yes" ]]; then
     echo "                                        4. Filtering of the Bowtie2-aligned, BCFtools-called VCF"
fi
echo " "



############################################################
##                        PART 0                          ##
##             Preparing files for the pipeline           ##
############################################################

echo "----- $(date) : Starting part 0 - Preparing files for the pipeline (8 jobs to complete, estimated time ~2 min):"

############################################################
## 0.1 Activate bcftools conda environment                ##
############################################################

echo "                                        $(date) : 0.1: Activating conda environment bcftools_env"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate bcftools_env


############################################################
## 0.2 Copy current script to Project_out folder          ##
############################################################

echo "                                        $(date) : 0.2: Copying this script to PROJECT_OUT folder"
mkdir -p "${PROJECT_OUT}"			# Creating PROJECT_OUT folder where the results will be saved
cp $INPUT_PATH/$SCRIPT "${PROJECT_OUT}/."	# Copies this script to the output folder for re-running or so it is clear what variables were used


############################################################
## 0.3 Copy reference genome to REFERENCE folder          ##
############################################################

echo "                                        $(date) : 0.3: Copying reference genome files to ./REFERENCE_INPUT_DIR"
cp ${REF}.fa REFERENCE_INPUT_DIR/.
cp ${REF}.fa.fai REFERENCE_INPUT_DIR/.
cp ${REF}.fa.gz REFERENCE_INPUT_DIR/.
cp ${REF}.sorted.gff REFERENCE_INPUT_DIR/.

# The following copies the reference genomefasta and gff to a genomes folder so it is with the other data
mkdir -p "$GENOME_DIR"
cp ${REF}.fa $GENOME_DIR/.
cp ${REF}.sorted.gff $GENOME_DIR/.

############################################################
## 0.4 Create "samples_and_read_groups.csv"               ##
############################################################

# This creates the samples_and_read_groups.csv that OVarFlow reads as input

echo "                                        $(date) : 0.4: Creating samples_and_read_groups.csv"
SAMPLES=("$PARENT" "$MUT1" "$MUT2" "$MUT3" "$MUT4" "$MUT5" "$MUT6")
cat > "$INPUT_PATH/samples_and_read_groups.csv" <<EOF
Reference Sequence:,${REF##*/}.fa.gz,,,,,
Reference Annotation:,${REF##*/}.sorted.gff,,,,,
,,,,,,
Min sequence length:,2000,,,,,
,,,,,,
old gvcf to include:,,,,,,
,,,,,,
forward reads,reverse reads,ID,PL - platform technology,CN - sequencing center,LB - library name,SM - uniq sample name
EOF

# Add only defined samples
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ "$SAMPLE" != "None" ]]; then
        echo "${SAMPLE}_R1.fastq.gz,${SAMPLE}_R2.fastq.gz,id_${SAMPLE},illumina,ENA,lib_${SAMPLE},${SAMPLE}" >> "$INPUT_PATH/samples_and_read_groups.csv"
    fi
done


############################################################
## 0.5 Generate "Bam_list.csv" for Bowtie2 and OVarFlow   ##
############################################################

# Creates two Bam_lists (one for Bowtie and one for OVarFlow) for CNV analysis and variant calling

echo "                                        $(date) : 0.5: Creating Bowtie2_BAM/Bam_list.csv"

mkdir -p "$BOW_BAM_OUT"

# CSV output
echo "sample_name,Bam_file" > "$BAM_LIST_BOW"	# generate Bam list for Bowtie2-generated Bam files (includes path)

for SAMPLE in "${SAMPLES[@]}"; do
    if [[ "$SAMPLE" != "None" ]]; then
        R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
        if [[ ! -f "$R1" || ! -f "$R2" ]]; then
            echo "ERROR: FASTQ files missing for $SAMPLE"
            exit 1
        fi
        # Just use the actual sample name, no prefix
        echo "${SAMPLE},"${BOW_BAM_OUT}"/${SAMPLE}.bam" >> "$BAM_LIST_BOW"
    fi
done


################################################################
##  0.6 Generate "FASTQ_metadata.csv"                         ##
################################################################

# This scetion creates a CSV file with metadata extracted from the fastq files:
#   - outputs md5sum so that one can later be certain of the fastq files used
#   - outputs other info that is useful for submission of the data to the NCBI SRA

echo "                                        $(date) : 0.6: Creating FASTQ_metadata.csv"

OUTPUT="${PROJECT_OUT}/FASTQ_metadata.csv"

echo "Sample,FASTQ_file,md5sum" > "$OUTPUT"

SAMPLES=("$PARENT" "$MUT1" "$MUT2" "$MUT3" "$MUT4" "$MUT5" "$MUT6")

for SAMPLE in "${SAMPLES[@]}"
do
    [[ "$SAMPLE" == "None" ]] && continue

    for READ in R1 R2
    do
        FASTQ="${FASTQ_DIR}/${SAMPLE}_${READ}.fastq.gz"

        [[ ! -f "$FASTQ" ]] && continue

        FILE=$(basename "$FASTQ")
        MD5=$(md5sum "$FASTQ" | awk '{print $1}')

        echo "${SAMPLE},${FILE},${MD5}" >> "$OUTPUT"

    done
done

# This next section copies the FASQ_metadata.csv into the Excel template


############################################################
## 0.7 Add FASTQ_metedata, Genome and Run_log to Excel out##
############################################################

echo "                                        $(date) : 0.7: Adding FASTQ metadata, Genome and Run_log to Excel output (~2 min)"

mkdir -p "$EXCEL_OUT"

python3 << EOF

import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

# File locations
source_csv = "$PROJECT_OUT/FASTQ_metadata.csv"
source_genome = "$EXCEL_GENOMES"
genome_sheet = "$GENOME"

input_excel = "$EXCEL_TEMPLATE"
output_excel = "$EXCEL_GENOME"

log_file = "$INPUT_PATH/${SCRIPT%.sh}.log"

MAX_GENOME_ROWS = 17000
MAX_GENOME_COLS = 20


# Load Excel template
wb = load_workbook(input_excel)


#=================================
# Add FASTQ_metadata
#=================================

df = pd.read_csv(source_csv)

if "FASTQ_metadata" in wb.sheetnames:
    ws = wb["FASTQ_metadata"]
else:
    ws = wb.create_sheet("FASTQ_metadata")

# Clear existing contents but retain formatting
for row in ws.iter_rows():
    for cell in row:
        cell.value = None

# Write header
for col_idx, col_name in enumerate(df.columns, start=1):
    ws.cell(
        row=1,
        column=col_idx,
        value=col_name
    )

# Write data
for row_idx, row in enumerate(df.values, start=2):
    for col_idx, value in enumerate(row, start=1):
        ws.cell(
            row=row_idx,
            column=col_idx,
            value=value
        )


#=================================
# Add Genome
#=================================

src = load_workbook(
    source_genome,
    read_only=True,
    data_only=True
)

if genome_sheet not in src.sheetnames:
    raise Exception(
        f"Sheet '{genome_sheet}' not found in {source_genome}"
    )

src_ws = src[genome_sheet]

if "Genome" in wb.sheetnames:
    dst_ws = wb["Genome"]
else:
    dst_ws = wb.create_sheet("Genome")

# Clear existing contents but retain formatting
for row in dst_ws.iter_rows():
    for cell in row:
        cell.value = None

# Copy genome data
row_count = 0

for r_idx, row in enumerate(
    src_ws.iter_rows(values_only=True),
    start=1
):

    if r_idx > MAX_GENOME_ROWS:
        break

    for c_idx, value in enumerate(
        row[:MAX_GENOME_COLS],
        start=1
    ):
        dst_ws.cell(
            row=r_idx,
            column=c_idx,
            value=value
        )

    row_count += 1

src.close()


#=================================
# Add Run_log
#=================================

log_ws = wb["Run_log"]

white_fill = PatternFill(
    fill_type="solid",
    fgColor="FFFFFF"
)

# Clear existing contents but retain formatting
for row in log_ws.iter_rows():
    for cell in row:
        cell.value = None

# Read current log
if os.path.exists(log_file):

    with open(log_file, "r") as f:
        log_lines = f.readlines()

    # Put each log line into a separate Excel row
    for row_idx, line in enumerate(log_lines, start=1):

        cell = log_ws.cell(
            row=row_idx,
            column=1
        )

        cell.value = line.rstrip("\r\n")
        cell.data_type = "s"
        cell.fill = white_fill

else:

    log_ws["A1"] = "Log file not found: " + log_file
    log_ws["A1"].fill = white_fill

# Set column width
log_ws.column_dimensions["A"].width = 150


#=================================
# Saving Excel workbook
#=================================

wb.save(output_excel)

EOF


############################################################
## 0.8: Build SnpEff database                             ##
############################################################

SNPEFF_GENOME_DIR="${SNPEFF_DATA}/${GENOME}"

if [[ -f "${SNPEFF_GENOME_DIR}/snpEffectPredictor.bin" ]]; then
echo "                                        $(date) : 0.8: SnpEff database already built!  Skipping build"
else

echo "                                        $(date) : 0.8: Building SnpEff database:"

mkdir -p "$SNPEFF_GENOME_DIR"

# ---- Step 2.0.1: Copy GFF ----
echo "                                           $(date) : 0.8.1: Copying GFF..."
cp "$REFERENCE_GFF" "${SNPEFF_GENOME_DIR}/genes.gff"

# ---- Step 2.0.2: Prepare FASTA ----
echo "                                           $(date) : 0.8.2: Preparing FASTA..."
cp "${REFERENCE_FASTA}" "${SNPEFF_GENOME_DIR}/sequences.fa"

# ---- Step 2.0.3: Add genome to snpEff.config if missing ----
if ! grep -q "^${GENOME}\.genome" "$SNPEFF_CONFIG"; then
    echo "                                             $(date) : 0.8.3: Adding genome to snpEff.config..."
    echo "${GENOME}.genome : ${GENOME}" >> "$SNPEFF_CONFIG"
fi

# ---- Step 2.0.4: Build SnpEff database ----
echo "                                           $(date) : 0.8.4: Building SnpEff database..."

{
    echo "----- $(date) : SnpEff build log for ${GENOME} ====="
    echo "Reference GFF : $REFERENCE_GFF"
    echo "Reference FASTA : $REFERENCE_FASTA"
    echo "SnpEff genome : $GENOME"
    echo ""
    
    java -Xmx8g -jar "$SNPEFF_JAR" build \
        -gff3 \
        -noCheckCds \
        -noCheckProtein \
        -v "$GENOME"
} > "$SNPEFF_LOG" 2>&1

if [[ $? -ne 0 ]]; then
    echo "ERROR: SnpEff build failed"
    exit 1
fi

# ---- Step 2.0.5: Verify output ----
if [[ ! -f "${SNPEFF_GENOME_DIR}/snpEffectPredictor.bin" ]]; then
    echo "ERROR: snpEffectPredictor.bin not created"
    exit 1
fi

echo "                                           $(date) :        SnpEff database built successfully"

# ---- Optional: compress FASTA (not required, but tidy) ----
bgzip -f "${SNPEFF_GENOME_DIR}/sequences.fa"

# ---- Optional: remove per-chromosome .bin files (not required, but tidy) ----
rm -f "${SNPEFF_GENOME_DIR}"/sequence.*.bin

fi

echo "----- $(date) : Part 0. Pipeline file preparation complete"
echo "                                        Part 0 can now be skipped if the pipeline is rerun or crashes from here onwards"
echo " "



############################################################
##                        PART 1                          ##
##                Bowtie2 read alignment                  ##
############################################################


############################################################
##                   PART 1 skip option                   ##
############################################################

if [[ "$SKIP_PART_1" == "Yes" ]]; then
     echo "----- $(date) : You chose to skip Part 1 manually - skipping Bowtie-alignment of WGS reads."
     echo " "
else

echo "----- $(date) : Starting Part 1. WGS Alignment with Bowtie2 and BAM production. 3 jobs, estimated time ~10 min per sample:"

############################################################
## 1.1 Run Bowtie2 read alignment and output Bams         ##
############################################################

echo "                                        $(date) : 1.1: Starting Bowtie2:"
# This part runs Bowtie2 to align WGS data to the selected reference genome and outputs results to $PROJECT_OUT

mkdir -p "$LOGS_OUT"

# Build Bowtie2 index if it doesn't already exist
if [[ ! -f "${BOWTIE2_INDEX}.1.bt2" && ! -f "${BOWTIE2_INDEX}.1.bt2l" ]]; then

    echo "                                           $(date) : 1.1.1: Building Bowtie2 index"

    bowtie2-build \
        "${REF}.fa" \
        "$BOWTIE2_INDEX" > "$BOWTIE_BUILD_LOG" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Bowtie2 index construction failed"
        exit 1
    fi
else
    echo "                                           $(date) : 1.1.1: Bowtie2 index already exists - skipping build"
fi


# Align each sample
for SAMPLE in "${SAMPLES[@]}"; do

    if [[ "$SAMPLE" != "None" ]]; then

        R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
        BAM="${BOW_BAM_OUT}/${SAMPLE}.bam"
        BOWTIE2_LOG="${LOGS_OUT}/${SAMPLE}.Bowtie2.log"
        RG_ID="id_${SAMPLE}"

        echo "                                           $(date) : 1.1.2: Aligning ${SAMPLE}"

        bowtie2 \
            --very-sensitive \
            -p "$BOWTIE2_THREADS" \
            -x "$BOWTIE2_INDEX" \
            -1 "$R1" \
            -2 "$R2" \
            --rg-id "$RG_ID" \
            --rg "SM:${SAMPLE}" \
            --rg "PL:illumina" \
            --rg "CN:ENA" \
            --rg "LB:lib_${SAMPLE}" \
            2> "$BOWTIE2_LOG" \
        | samtools view -b - 2>> "$BOWTIE2_LOG" \
        | samtools sort -@ "$BOWTIE2_THREADS" -o "$BAM" - 2>> "$BOWTIE2_LOG"

        if [[ ${PIPESTATUS[0]} -ne 0 || ${PIPESTATUS[1]} -ne 0 || ${PIPESTATUS[2]} -ne 0 ]]; then
            echo "ERROR: Bowtie2/samtools alignment failed for ${SAMPLE}"
            exit 1
        fi

        samtools index -@ "$BOWTIE2_THREADS" "$BAM"

    fi

done


############################################################
## 1.2  Generating read and alignment QC report           ##
############################################################

echo "                                        $(date) : 1.2: Generating read and alignment QC report for all samples"

# Write header
echo "sample,total_reads,mapped_reads,mapped_percent,mean_coverage" > "$BAM_QC_REPORT_BOW"

# Loop through BAM files
for bam in "$BOW_BAM_OUT"/*.bam; do

    sample=$(basename "$bam" .bam)

    # Read and alignment statistics
    stats=$(samtools flagstat "$bam")

    total=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped=$(echo "$stats" | grep "mapped (" | head -n1 | awk '{print $1}')
    percent=$(echo "$stats" | grep "mapped (" | head -n1 | grep -oP '\(\K[0-9.]+')

    # Calculate genome-wide mean coverage
    mean_cov=$(samtools coverage "$bam" | \
        awk 'NR > 1 {
            bases += $3
            depth += $7 * $3
        }
        END {
            if (bases > 0)
                printf "%.2f", depth/bases
            else
                print "NA"
        }')

    echo "${sample},${total},${mapped},${percent},${mean_cov}" >> "$BAM_QC_REPORT_BOW"

done

############################################################
## 1.3  Adding Alignment QC to Excel template             ##
############################################################

echo "                                        $(date) : 1.3: Adding alignment QC report to Excel output (BAMqc tab)"

python3 << EOF

import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

source_csv = "$BOW_BAM_OUT/bam_qc_summary.csv"
input_excel = "$EXCEL_GENOME"
output_excel = "$EXCEL_BOW_BAMqc"

log_file = "$INPUT_PATH/${SCRIPT%.sh}.log"


#=================================
# Add Bam_QC data
#=================================

df = pd.read_csv(source_csv)

wb = load_workbook(input_excel)

# Remove existing Bam_QC sheet if present
if "Bam_QC" in wb.sheetnames:
    del wb["Bam_QC"]

ws = wb.create_sheet("Bam_QC")

# Write header
for col_idx, col_name in enumerate(df.columns, start=1):
    ws.cell(
        row=1,
        column=col_idx,
        value=col_name
    )

# Write data
for row_idx, row in enumerate(df.values, start=2):
    for col_idx, value in enumerate(row, start=1):
        ws.cell(
            row=row_idx,
            column=col_idx,
            value=value
        )


#=================================
# Update Run_log
#=================================

log_ws = wb["Run_log"]

white_fill = PatternFill(
    fill_type="solid",
    fgColor="FFFFFF"
)

# Clear existing contents but retain formatting
for row in log_ws.iter_rows():
    for cell in row:
        cell.value = None

# Read current log
if os.path.exists(log_file):

    with open(log_file, "r") as f:
        log_lines = f.readlines()

    # Put each log line into a separate Excel row
    for row_idx, line in enumerate(log_lines, start=1):

        cell = log_ws.cell(
            row=row_idx,
            column=1
        )

        cell.value = line.rstrip("\r\n")
        cell.data_type = "s"
        cell.fill = white_fill

else:

    log_ws["A1"] = "Log file not found: " + log_file
    log_ws["A1"].fill = white_fill

# Set column width
log_ws.column_dimensions["A"].width = 150


#=================================
# Saving Excel workbook
#=================================

wb.save(output_excel)

EOF


echo "----- $(date) : Part 1. WGS Alignment and BAM production with Bowtie2 compelete"
echo "                                        Part 1 can now be skipped if the pipeline is rerun or crashes from here onwards"
echo " "

fi ####### ... End of user-defined skip option if SKIP_PART_1="Yes"


###################################################################
##                            PART 2                             ##
##            Gene and chromosome CNV determination              ##
###################################################################


###################################################################
##                      PART 2 skip option                       ##
###################################################################

if [[ "$SKIP_PART_2" == "Yes" ]]; then
     echo "----- $(date) : You chose to skip Part 2 manually - skipping gene and chromosome CNV determination."
     echo " "
else

echo "----- $(date) : Starting Part 2: Gene and chromosome CNV determination from Bowtie aignment. 4 jobs, estimated time ~5 min "

###################################################################
## 2.1: count reads per gene                                     ##
###################################################################

echo "                                        $(date) : 2.1: Calculating gene read counts with featureCounts - for errors, check the featureounts log:"
echo "                                                                             - $LOGS_OUT/Feature_counts.log"

mkdir -p "$BOW_CNV_OUT"

BAM_LIST=()

# Read BAM paths from CSV (skip header)
while IFS=',' read -r SAMPLE BAM; do
    if [[ "$SAMPLE" != "sample_name" ]]; then
        BAM_LIST+=("$BAM")
    fi
done < "$BAM_LIST_BOW"

#The following two lines are for debugging only
#echo "Ordered BAM list for RPKM:"
#printf "%s\n" "${BAM_LIST[@]}"

# Run featureCounts
featureCounts \
    -T ${THREADS:-1} \
    -F GFF \
    -t CDS \
    -g gene_id \
    -p \
    -a "$REFERENCE_GFF" \
    -o "$BOW_CNV_OUT/gene_counts.txt" \
    "${BAM_LIST[@]}" \
    > "$LOGS_OUT/Bowtie2_Feature_counts.log" 2>&1


###################################################################
## 2.2: calculate total mapped reads                             ##
###################################################################

echo "                                        $(date) : 2.2: Counting mapped reads"

TOTAL_READS="$BOW_CNV_OUT/total_reads.txt"
> "$TOTAL_READS"

# Use the same BAM_LIST to guarantee PARENT first
for BAM in "${BAM_LIST[@]}"; do
    SAMPLE=$(basename "$BAM")
    SAMPLE=${SAMPLE%.bam}   # strip extension cleanly

    COUNT=$(samtools view -c -F 260 "$BAM")
    echo -e "$SAMPLE\t$COUNT" >> "$TOTAL_READS"
done

###################################################################
## 2.3: calculate RPKM                                           ##
###################################################################

echo "                                        $(date) : 2.3: Calculating RPKM "

python3 << EOF
import pandas as pd

counts = pd.read_csv("$BOW_CNV_OUT/gene_counts.txt", sep="\t", comment="#")

# Read total reads
totals = {}
with open("$BOW_CNV_OUT/total_reads.txt") as f:
    for line in f:
        s, c = line.strip().split("\t")
        totals[s] = int(c)

# Build samples list in order (Parent first)
samples = ["${PARENT}", "${MUT1}", "${MUT2}", "${MUT3}", "${MUT4}", "${MUT5}", "${MUT6}"]
samples = [s for s in samples if s != "None"]

# Initialize RPKM dataframe
rpkm = pd.DataFrame()
rpkm["GeneID"] = counts["Geneid"]
lengths = counts["Length"]

for sample in samples:
    # Match the correct column from featureCounts output
    bam_col = [c for c in counts.columns if sample in c][0]
    total = totals[sample]
    rpkm[sample] = (1e9 * counts[bam_col]) / (total * lengths)

rpkm.to_csv("$BOW_CNV_OUT/gene_RPKM.csv", index=False)
EOF


###################################################################
## 2.4: copy RPKM data to final Excel file                       ##
###################################################################

echo "                                        $(date) : 2.4: Copying RPKM data to Excel output (CNV tab) "

python3 << EOF

import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill
from openpyxl.utils import column_index_from_string
from openpyxl.utils.dataframe import dataframe_to_rows


# --- Paths ---
csv_file = "$BOW_CNV_OUT/gene_RPKM.csv"
excel_file_in = "${EXCEL_BOW_BAMqc}"
excel_file_out = "${EXCEL_BOW_CNV}"
sheet_name = "Gene_CNV"

log_file = "$INPUT_PATH/${SCRIPT%.sh}.log"


# --- Destination start cell ---
start_col = column_index_from_string("AF")
start_row = 30


# --- Load CSV ---
df = pd.read_csv(csv_file)


# --- Load workbook ---
wb = load_workbook(excel_file_in)


# --- Get or create worksheet ---
if sheet_name in wb.sheetnames:
    ws = wb[sheet_name]
else:
    ws = wb.create_sheet(sheet_name)


# --- Clear previous data region (AB:AI from row 30 downward) ---
end_col = column_index_from_string("AI")

for r in range(start_row, ws.max_row + 1):
    for c in range(start_col, end_col + 1):
        ws.cell(row=r, column=c).value = None


# --- Write dataframe ---
for r_idx, row in enumerate(
    dataframe_to_rows(df, index=False, header=True)
):
    for c_idx, value in enumerate(row):

        ws.cell(
            row=start_row + r_idx,
            column=start_col + c_idx,
            value=value
        )


#=================================
# Update Run_log
#=================================

log_ws = wb["Run_log"]

white_fill = PatternFill(
    fill_type="solid",
    fgColor="FFFFFF"
)

# Clear existing contents but retain formatting
for row in log_ws.iter_rows():
    for cell in row:
        cell.value = None

# Read current log
if os.path.exists(log_file):

    with open(log_file, "r") as f:
        log_lines = f.readlines()

    # Put each log line into a separate Excel row
    for row_idx, line in enumerate(log_lines, start=1):

        cell = log_ws.cell(
            row=row_idx,
            column=1
        )

        cell.value = line.rstrip("\r\n")
        cell.data_type = "s"
        cell.fill = white_fill

else:

    log_ws["A1"] = "Log file not found: " + log_file
    log_ws["A1"].fill = white_fill

# Set column width
log_ws.column_dimensions["A"].width = 150


#=================================
# Save workbook once
#=================================

wb.save(excel_file_out)

EOF
echo "                                        $(date) : Part 2: Gene and chromosome CNV determination from Bowtie2 alignmnt complete"
echo "                                        $(date) : If pipeline crashes in next step, you can skip part 2 and rerun"
echo " "

fi ####### ... End of user-defined skip option if SKIP_PART_2="Yes"


############################################################
##                         PART 3                         ##
##      Variant calling with mpileup and BCFtools         ##
############################################################

############################################################
##                    PART 3 skip option                  ##
############################################################

if [[ "$SKIP_PART_3" == "Yes" ]]; then
     echo "----- $(date) : You chose to skip Part 3 manually - skipping variant calling with BCFtools"
    echo " "
else

echo "----- $(date) : Starting Part 3. Variant calling with BCFtools (using Bowtie2 alignments). 1 job, estimated time = ~20 min per sample:"

mkdir -p "$BOW_VCF_OUT"

############################################################
##  3.1: mpileup & variant calling                        ##
############################################################

echo "                                        $(date) : 3.1: mpileup with variant calling"

BAM_LIST=()

# Read BAM paths from CSV (skip header) - this t=is to ensure the correct order of the VCF
while IFS=',' read -r SAMPLE BAM; do
    if [[ "$SAMPLE" != "sample_name" ]]; then
        BAM_LIST+=("$BAM")
    fi
done < "$BAM_LIST_BOW"

#The following two lines are for debugging only
#echo "Ordered BAM list for BCFtools:"
#printf "%s\n" "${BAM_LIST[@]}"

bcftools mpileup \
    -f "${REF}.fa" \
    -a FORMAT/DP,FORMAT/AD \
    --max-depth 10000 \
    --threads 40 \
    -Ou \
    "${BAM_LIST[@]}" 2> "$MPILEUP_BOW_LOG" | \
bcftools call \
    -"$CALL_MODE" \
    -f GQ \
    --ploidy "$PLOIDY" \
    -Oz \
    -o "$BOW_BCF_RAW_VCF" 2>> "$MPILEUP_BOW_LOG"

bcftools index -f "$BOW_BCF_RAW_VCF"

echo "                                        $(date) : Part 3. Variant calling of Bowtie2 alignment with BCFtools complete"
echo "                                        Part 3 can now be skipped if the pipeline is rerun or crashes from here onwards"
echo " "

fi ####### ... End of user-defined skip option if SKIP_PART_3="Yes"



############################################################
##                         PART 4                         ##
##              VCF annotation and filtering              ##
############################################################

############################################################
##                    PART 4 skip option                  ##
############################################################

if [[ "$SKIP_PART_4" == "Yes" ]]; then
     echo "----- $(date) : You chose to skip Part 4 manually. Skipping VCF annotation and filtering"
     echo " "
else

echo "----- $(date) : Starting Part 4. Bowtie2/BCFtools VCF annotation and filtering. 7 jobs, estimated time = 2 min:"

############################################################
## 4.1: Annotate variants with SnpEff                     ##
############################################################

echo "                                        $(date) : 4.1: Annotating BCFtools variants with SnpEff"

mkdir -p "${SNPeff_DIR}"

BOW_BCF_ANNOTATED_VCF_TMP="${BOW_VCF_OUT}/annotated.tmp.vcf"

# Check that the input VCF exists
if [[ ! -f "$BOW_BCF_RAW_VCF" ]]; then
    echo "ERROR: Input VCF not found: $BOW_BCF_RAW_VCF"
    exit 1
fi

# Run SnpEff and capture stdout/stderr
#echo "----- $(date) : Step 5.1.1: Running SnpEff annotation..."
java -Xmx8g -jar "${SNPEFF_JAR}" \
    -c "$SNPEFF_CONFIG" \
    -dataDir "$SNPEFF_DATA" \
    -v \
    -stats "$SNPEFF_HTML" \
    "$SNPEFF_GENOME" \
    "$BOW_BCF_RAW_VCF" \
    > "$BOW_BCF_ANNOTATED_VCF_TMP" 2> "${LOGS_OUT}/BOW_BCF_snpEff.log"

# Check that SnpEff completed successfully and output exists
if [[ $? -ne 0 ]]; then
    echo "ERROR: SnpEff failed. See log: ${LOGS_OUT}/BOW_BCF_snpEff.log"
    exit 1
fi
if [[ ! -s "$BOW_BCF_ANNOTATED_VCF_TMP" ]]; then
    echo "ERROR: SnpEff produced an empty VCF: $BOW_BCF_ANNOTATED_VCF_TMP"
    exit 1
fi

# Compress with bgzip AFTER SnpEff completes
#echo "----- $(date) : Step 5.1.2: Compressing annotated VCF with bgzip..."
bgzip -f "$BOW_BCF_ANNOTATED_VCF_TMP"
if [[ $? -ne 0 ]]; then
    echo "ERROR: bgzip failed"
    exit 1
fi
# Rename to final annotated VCF
mv "${BOW_BCF_ANNOTATED_VCF_TMP}.gz" "$BOW_BCF_ANNOTATED_VCF"

# Index with tabix
#echo "----- $(date) : Step 5.1.3: Indexing annotated VCF..."
tabix -f "$BOW_BCF_ANNOTATED_VCF"
if [[ $? -ne 0 ]]; then
    echo "ERROR: tabix indexing failed"
    exit 1
fi

# Verify final VCF is readable
zcat "$BOW_BCF_ANNOTATED_VCF" | tail -5 > /dev/null
if [[ $? -ne 0 ]]; then
    echo "ERROR: Final annotated VCF is corrupted"
    exit 1
fi

############################################################
## 4.2: Reorder variants for correct format               ##
############################################################

echo "                                        $(date) : 4.2: Reordering variants for correct format"

# Build sample order string (comma-separated)
ORDERED_SAMPLES=""

for s in "$PARENT" "$MUT1" "$MUT2" "$MUT3" "$MUT4" "$MUT5" "$MUT6"; do
    if [[ "$s" != "None" ]]; then
        if [[ -z "$ORDERED_SAMPLES" ]]; then
            ORDERED_SAMPLES="$s"
        else
            ORDERED_SAMPLES="${ORDERED_SAMPLES},${s}"
        fi
    fi
done


# Reorder samples properly
bcftools view \
    -s "$ORDERED_SAMPLES" \
    -Oz -o "$BOW_BCF_REORDERED_VCF" \
    "$BOW_BCF_ANNOTATED_VCF"

bcftools index -f "$BOW_BCF_REORDERED_VCF"


############################################################
## 4.3: Filter variants by depth & quality                ##
############################################################

echo "                                        $(date) : 4.3: Filtering BCFtools variants by QUAL ($MIN_QUAL_BCF), DP ($MIN_DEPTH_BCF) and GQ ($MIN_GQ_BCF)"

bcftools filter -i "FORMAT/DP[0]>=$MIN_DEPTH_BCF && QUAL>=$MIN_QUAL_BCF && FORMAT/GQ[0]>=$MIN_GQ_BCF" -Oz -o "$BOW_BCF_FILTERED_VCF" "$BOW_BCF_REORDERED_VCF"		# This filter applies only to the first sample (PARENT)

bcftools index -f "$BOW_BCF_FILTERED_VCF"


############################################################
## 4.4: Keep variants different from PARENT               ##
############################################################

echo "                                        $(date) : 4.4: Removing "variants" detected in parental line (PARENT)"

TMP_GT=$(mktemp)
TMP_KEEP=$(mktemp)

# Extract genotypes (from FILTERED VCF!)
bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$BOW_BCF_FILTERED_VCF" > "$TMP_GT"

# Identify positions where any sample differs from WT
awk '
{
    parent=$3; keep=0;
    for(i=4;i<=NF;i++){
        if($i != parent){keep=1; break}
    }
    if(keep) print $1"\t"$2
}' "$TMP_GT" > "$TMP_KEEP"

bcftools view -T "$TMP_KEEP" -Oz -o "$BOW_BCF_DIFF_VCF" "$BOW_BCF_FILTERED_VCF"
bcftools index -f "$BOW_BCF_DIFF_VCF"

rm -f "$TMP_GT" "$TMP_KEEP"

############################################################
## 4.5: Keep variants where PARENT is 0/0                 ##
############################################################

echo "                                        $(date) : 4.5: Removing variants where PARENT ≠ 0/0"

bcftools view -i 'GT[0]="0/0"' -Oz -o "$BOW_BCF_DIFF_WT0_VCF" "$BOW_BCF_DIFF_VCF"
bcftools index -f "$BOW_BCF_DIFF_WT0_VCF"


############################################################
## 4.6: Keep variants within mRNA regions only            ##
############################################################

echo "                                        $(date) : 4.6: Removing variants outwith mRNA regions"

# ---- Define files ----
MRNA_BED="${BOW_VCF_OUT}/mRNA_regions.bed"

# ---- Step 4.6.1: Extract mRNA features from GFF ----

awk '$3 == "mRNA"' "$REFERENCE_GFF" | \
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' > "$MRNA_BED"

# Check file is not empty
if [[ ! -s "$MRNA_BED" ]]; then
    echo "ERROR: No mRNA features found in GFF"
    exit 1
fi

# ---- Step 4.6.2: Sort BED file (important for bcftools) ----
sort -k1,1 -k2,2n "$MRNA_BED" -o "$MRNA_BED"


# ---- Step 4.6.3: Filter VCF ----

bcftools view \
    -R "$MRNA_BED" \
    -Oz -o "$BOW_BCF_MRNA_VCF" \
    "$BOW_BCF_DIFF_WT0_VCF"

# ---- Step 4.6.4: Index ----
bcftools index -f "$BOW_BCF_MRNA_VCF"


###################################################################
##  4.7: copying VCF data to Excel template                      ##
###################################################################

if [[ "$FINVCF" -eq 5 ]]; then
    echo "                                        $(date) : 4.7: Exporting filtered BCFtools VCF 5 data to Excel output (BCFtools tab)"
    gunzip -f "$BOW_BCF_DIFF_WT0_VCF"
    tsv_file="${BOW_BCF_DIFF_WT0_VCF%.gz}"
else
    echo "                                        $(date) : 4.7: Exporting filtered BCFtools VCF 6 data to Excel output (BCFtools tab)"
    gunzip -f "$BOW_BCF_MRNA_VCF"
    tsv_file="${BOW_BCF_MRNA_VCF%.gz}"
fi

python3 << EOF

import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill


# --- Excel files ---
template_file = "${EXCEL_BOW_CNV}"
excel_file = "${EXCEL_BOW_BCF}"

sheet_name = "BCFtools"
start_row, start_col = 1, 1

log_file = "$INPUT_PATH/${SCRIPT%.sh}.log"


#=================================
# Add BCFtools VCF data
#=================================

with open("${tsv_file}", 'r') as f:
    lines = f.readlines()

header_lines = [
    line.strip()
    for line in lines
    if line.startswith('#')
]

data_lines = [
    line.strip()
    for line in lines
    if not line.startswith('#')
]

df_data = pd.DataFrame(
    [line.split('\t')[:16] for line in data_lines]
)

df_header = pd.DataFrame([[line] for line in header_lines])

df_combined = pd.concat(
    [df_header, df_data],
    ignore_index=True
)


# --- Load workbook ---
wb = load_workbook(template_file)

ws = wb[sheet_name]


# --- Write VCF data ---
for r_idx, row in enumerate(
    df_combined.values,
    start=start_row
):

    for c_idx, value in enumerate(
        row,
        start=start_col
    ):

        ws.cell(
            row=r_idx,
            column=c_idx,
            value=value
        )


#=================================
# Update Run_log
#=================================

log_ws = wb["Run_log"]

white_fill = PatternFill(
    fill_type="solid",
    fgColor="FFFFFF"
)

# Clear existing contents but retain formatting
for row in log_ws.iter_rows():
    for cell in row:
        cell.value = None

# Read current log
if os.path.exists(log_file):

    with open(log_file, "r") as f:
        log_lines = f.readlines()

    # Put each log line into a separate Excel row
    for row_idx, line in enumerate(log_lines, start=1):

        cell = log_ws.cell(
            row=row_idx,
            column=1
        )

        cell.value = line.rstrip("\r\n")
        cell.data_type = "s"
        cell.fill = white_fill

else:

    log_ws["A1"] = "Log file not found: " + log_file
    log_ws["A1"].fill = white_fill

# Set column width
log_ws.column_dimensions["A"].width = 150


#=================================
# Remove unused variant-calling tabs
#=================================

for sheet in ["OVarFlow", "DeepVariant"]:
    if sheet in wb.sheetnames:
        del wb[sheet]

#=================================
# Save workbook once
#=================================

wb.save(excel_file)

EOF

cp $EXCEL_BOW_BCF $PROJECT_OUT/.

echo "                                        $(date) : Part 4. Bowtie2/BCFtools VCF annotation and filtering complete "
echo "                                        Part 4 can now be skipped if the pipeline is rerun or crashes from here onwards"
echo " "
echo "                                        Parts 1-4: Bowtie2-alignment and BCFtools variant calling complete - inspect results in:"
echo "                                          - $EXCEL_BOW_BCF"
echo " "

fi ####### ... End of user-defined skip option if SKIP_PART_4="Yes"


###################################################################
##  Pipeline complete - report following at end of master log    ##
###################################################################

echo " 
***********************************************************************************************************
----- $(date) : $SCRIPT pipeline complete 
***********************************************************************************************************

Check the following results files to make sure that the data has been processed correctly
   - Bowtie2-aligned output:  $EXCEL_BOW_BCF

You can rerun sections without needing to rerun the entire script, for example:
   - you can rerun VCF filtering without needing to re-aling or recall;
But you will need to rerun all subsequent steps so that the Excel output is updated accordingly.

Once you're satisfied that processing is complete, do the following:
   - Copy the following *results folder* to permanent storage (JAF drive, RITseq and OE analysis\Raw data):
        - $PROJECT_OUT
   - Copy the FASTQ_INPUT_DIR folder (containing raw WGS data) into the *results folder* you just made in permanent storage;
   - Delete the *results folder* and FASTQ.gz files from the cluster (when verified that your copy is sound) to free up disk space."

# ---- copy log master log file to output directory
cp $INPUT_PATH/${SCRIPT%.sh}.log $PROJECT_OUT/.
