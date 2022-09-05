#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name='sample_name_XX_MuTect2'
#SBATCH --output='sample_name_XX_MuTect2.log'
#SBATCH --mem-per-cpu=100G
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH -m block
/bin/hostname



genome=hg38


# path to the mapping results directory
mapping_results_dir=SPECIFY_mapping_results_dir

sample_name=sample_name_XX

sample_name_control=SPECIFY_sample_name_control

echo ${sample_name}

echo ${sample_name_control}


input_tumor=${sample_name}
input_normal=${sample_name_control}

date



#############################
# Call variants
#############################

module load GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8

java -jar -Xmx2g /cm/shared/specific/apps/GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8/gatk-package-4.1.8.1-local.jar Mutect2 -R /research/lab_winter/resources/reference_files/genomes/hg38/indices_for_BWA/hg38.fa -I ${mapping_results_dir}/bam/${input_tumor}_dedup_RG.bam -I ${mapping_results_dir}/bam/${input_normal}_dedup_RG.bam -normal ${input_normal} -O ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.vcf.gz 

java -jar -Xmx2g /cm/shared/specific/apps/GATK/4.1.8.1-GCCcore-9.3.0-Java-1.8/gatk-package-4.1.8.1-local.jar FilterMutectCalls -R /research/lab_winter/resources/reference_files/genomes/hg38/indices_for_BWA/hg38.fa -V ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.vcf.gz -O ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vcf.gz 



#############################
# Annotate variants
#############################

module load VEP/103.1-GCC-10.2.0

vep \
    --allele_number --allow_non_variant \
    --assembly GRCh38 \
    --cache --dir_cache /nobackup/lab_bsf/resources/VEP --offline \
    --dir_plugins /nobackup/lab_bsf/resources/VEP/Plugins_103 \
    --dont_skip --everything --failed 1 \
    --fasta /nobackup/lab_bsf/resources/VEP/homo_sapiens_merged/103_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --flag_pick_allele_gene --force_overwrite --format vcf \
    --gencode_basic --hgvsg --exclude_predicted \
    --plugin CADD,/nobackup/lab_bsf/resources/CADD/hg38/1.6/whole_genome_SNVs.tsv.gz,/nobackup/lab_bsf/resources/CADD/hg38/1.6/gnomad.genomes.r3.0.indel.tsv.gz \
    --custom /nobackup/lab_bsf/resources/ClinVar/hg38/clinvar_20210213.vcf.gz,ClnV,vcf,exact,0,ALLELEID,ORIGIN,CLNSIG,CLNREVSTAT,CLNDN \
    --species homo_sapiens --merged \
    --vcf --fork 16 \
    --input_file ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vcf.gz \
    --output_file ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vep.vcf \
    --stats_file ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vep.html \
    --warning_file ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vep.warning.txt;



# split gene annotation in a vcf file
/research/lab_winter/users/himrichova/resources/scripts/split_gene_annotation_in_vcf.sh ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vep.vcf ${mapping_results_dir}/MuTect2_variants/${input_tumor}.vs.${input_normal}.variants.filter.vep_split.vcf


