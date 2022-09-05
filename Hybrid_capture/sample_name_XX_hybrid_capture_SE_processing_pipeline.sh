#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --job-name='sample_name_XX_HybridCapture'
#SBATCH --output='sample_name_XX_HybridCapture.log'
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task='2'
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH -m block
/bin/hostname


genome=hg38

# path to the directory with raw (unmapped) bam files
local_raw_data_location=SPECIFY_local_raw_data_location

# path to the annotation file with sample names
sample_names=SPECIFY_sample_names

# path to the mapping results directory
mapping_results_dir=SPECIFY_mapping_results_dir

sample_name=sample_name_XX

echo ${sample_name}



############################################
# Generate Fastq files, fastqc, trimming   #
############################################
date

#Target to produce: fastq file *fq.gz

### convert raw reads to fastq format
module load BamTools/2.5.1-foss-2018b

bamtools convert -format fastq -in ${local_raw_data_location}/${sample_name}.bam | gzip > ${mapping_results_dir}/fastq/${sample_name}.fq.gz



#Target to produce trimmed file *trimmed.fq.gz
module load Trimmomatic/0.39-Java-11

`java -Xmx60000m -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -phred33 -threads 2 ${mapping_results_dir}/fastq/${sample_name}.fq.gz ${mapping_results_dir}/fastq/${sample_name}_trimmed.fq.gz ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`



#Target to produce fastQC report
module load FastQC/0.11.9-Java-11

`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_trimmed.fq.gz`



###########################
### Alignment using bwa ###
###########################

module load bio/BWA/0.7.17-GCC-9.3.0

bwa aln -t 4 /research/lab_winter/resources/reference_files/genomes/hg38/indices_for_BWA/hg38.fa ${mapping_results_dir}/fastq/${sample_name}_trimmed.fq.gz > ${mapping_results_dir}/bam/${sample_name}_aln.bwa

bwa samse /research/lab_winter/resources/reference_files/genomes/hg38/indices_for_BWA/hg38.fa ${mapping_results_dir}/bam/${sample_name}_aln.bwa ${mapping_results_dir}/fastq/${sample_name}_trimmed.fq.gz > ${mapping_results_dir}/bam/${sample_name}_aln.sam


######################################################
### Mark duplicate reads, convert to .bam and sort ###
######################################################

module load bio/picard/2.25.1-Java-11

### clean sam (unmapped reads removed --> otherwise give errors downstream)
java -jar -Xmx2g $EBROOTPICARD/picard.jar CleanSam I=${mapping_results_dir}/bam/${sample_name}_aln.sam O=${mapping_results_dir}/bam/${sample_name}_clean.sam

### sort and convert to bam
java -jar -Xmx2g $EBROOTPICARD/picard.jar SortSam I=${mapping_results_dir}/bam/${sample_name}_clean.sam O=${mapping_results_dir}/bam/${sample_name}_sort.bam SORT_ORDER=coordinate

### Mark duplicates
java -jar -Xmx2g $EBROOTPICARD/picard.jar MarkDuplicates I=${mapping_results_dir}/bam/${sample_name}_sort.bam O=${mapping_results_dir}/bam/${sample_name}_dedup.bam  METRICS_FILE=${mapping_results_dir}/bam/${sample_name}_metrics.txt



### add read groups (RG) for GATK to accept the bam files
#RGID(String)   Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
#RGLB (String)  Read Group library Required.
#RGPL (String)  Read Group platform (e.g. illumina, solid) Required.
#RGPU (String)  Read Group platform unit (eg. run barcode) Required.
#RGSM (String)  Read Group sample name Required


java -jar -Xmx2g $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${mapping_results_dir}/bam/${sample_name}_dedup.bam O=${mapping_results_dir}/bam/${sample_name}_dedup_RG.bam RGID=4 RGLB=HC RGPL=illumina RGPU=${sample_name} RGSM=${sample_name} VALIDATION_STRINGENCY=LENIENT

### Build bam index
java -jar -Xmx2g $EBROOTPICARD/picard.jar BuildBamIndex I=${mapping_results_dir}/bam/${sample_name}_dedup_RG.bam



date

####################