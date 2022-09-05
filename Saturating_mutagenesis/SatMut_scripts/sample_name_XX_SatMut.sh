#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --job-name='sample_name_XX_SatMut'
#SBATCH --output='sample_name_XX_SatMut.log'
#SBATCH --mem-per-cpu=100G
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH -m block
/bin/hostname

module load FastQC/0.11.9-Java-11
module load multigrep/20220218

module load bio/BWA/0.7.17-GCC-9.3.0
module load bio/SAMtools/1.10-GCC-9.3.0
module load Trimmomatic/0.39-Java-11
module load bio/picard/2.25.1-Java-11


genome=hg38

# path to the directory with SatMut input files
directory_with_SatMut_input_files=SPECIFY_directory_with_SatMut_input_files

# path to the directory with SatMut scripts
directory_with_SatMut_scripts=SPECIFY_directory_with_SatMut_scripts

# path to the directory with raw (unmapped) bam files
local_raw_data_location=SPECIFY_local_raw_data_location

# path to the annotation file with sample names
sample_names=SPECIFY_sample_names

# path to the mapping results directory
mapping_results_dir=SPECIFY_mapping_results_dir

sample_name=sample_name_XX

echo ${sample_name}

date

#Target to produce: fastq file
fastq_out=$(echo ${mapping_results_dir}/fastq/${sample_name}_R1.fastq)

`samtools view ${local_raw_data_location}/${sample_name}.bam | awk -v fastq_out=${fastq_out} '{ print "@"$1"\n"$10"\n+\n"$11 > fastq_out; }'`

`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1.fastq`


# Trimming
# Target to produce: *trimmed.fastq

# Apply HEADCROP:5 and MINLEN:31
`java -Xmx60000m -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -phred33 -threads 2 ${mapping_results_dir}/fastq/${sample_name}_R1.fastq ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq HEADCROP:5 ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:31`



`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq`
date


############################################
# Mapping reads to the VHL/CRBN, sorting
############################################

gene=$(cat ${sample_names} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f3|sort -u)
orf=$(cat ${sample_names} | multigrep.sh -f 1 -w -p ${sample_name}|cut -f10|sort -u)

bwa aln -n 5 ${directory_with_SatMut_input_files}/${gene}.fa ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq > ${mapping_results_dir}/bam/${sample_name}_R1_trimmed.sai

bwa samse ${directory_with_SatMut_input_files}/${gene}.fa ${mapping_results_dir}/bam/${sample_name}_R1_trimmed.sai ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq > ${mapping_results_dir}/bam/${sample_name}_R1_trimmed.bam


java -jar -Xmx2g $EBROOTPICARD/picard.jar SortSam I=${mapping_results_dir}/bam/${sample_name}_R1_trimmed.bam O=${mapping_results_dir}/bam/${sample_name}_R1_trimmed.sorted.bam SORT_ORDER=coordinate



############################################
# Call mutations
############################################

echo call mutations for sample: ${sample_name}, gene: ${gene}, orf: ${orf}

module load bio/GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8


java -jar -Xmx2g /cm/shared/specific/apps/GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8/gatk-package-4.1.8.1-local.jar AnalyzeSaturationMutagenesis -I ${mapping_results_dir}/bam/${sample_name}_R1_trimmed.sorted.bam -R ${directory_with_SatMut_input_files}/${gene}.fa --orf ${orf} -O ${mapping_results_dir}/mut_calling/${sample_name}.gatk --paired-mode false



# Within the reference, the ORFs are as follows: VHL = 1278-1919, CRBN = 1278-2606


############################################
# Generate a freqency matrix
############################################

echo generate a frequency matrix for sample: ${sample_name}, gene: ${gene}, orf: ${orf}

module load R/4.1.0-foss-2021a

R --vanilla --args ${mapping_results_dir}/mut_calling/${sample_name}.gatk.aaCounts ${mapping_results_dir}/wt_residues_${gene}.t.txt ${mapping_results_dir}/mut_calling/Freq_aaCounts/${sample_name}.freq.txt < ${directory_with_SatMut_scripts}/SatMut_freq_matrix.R



