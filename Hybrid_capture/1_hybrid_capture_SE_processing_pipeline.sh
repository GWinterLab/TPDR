###################################################################################
# Process sequencing reads of AML cell lines cultured for 6 months in ARV-771
# - from raw bam files to mutations
###################################################################################

# path to a directory with raw (unmapped) bam files
local_raw_data_location=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/11_SM_and_Hybrid_capture_FINAL/2_Hybrid_capture/1_Hybrid_capture_KBM7/raw_bam_files

# path to a sample annotation file
sample_names=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/11_SM_and_Hybrid_capture_FINAL/2_Hybrid_capture/1_Hybrid_capture_KBM7/sample_annotation.txt 

# path to a directory where results will be generated
mapping_results_dir=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/11_SM_and_Hybrid_capture_FINAL/2_Hybrid_capture/1_Hybrid_capture_KBM7/mapping_results



# create directories for data processing results:
mkdir -p ${mapping_results_dir}
mkdir -p ${mapping_results_dir}/fastq
mkdir -p ${mapping_results_dir}/fastqc
mkdir -p ${mapping_results_dir}/bam
mkdir -p ${mapping_results_dir}/MuTect2_variants



############################################
# Hybrid capture data preprocessing
############################################

cd ${mapping_results_dir}

# copy a hybrid capture preprocessing pipeline to the working directory
cp /research/lab_winter/users/himrichova/resources/scripts/12_SatMut/Hybrid_capture_data_analysis/sample_name_XX_hybrid_capture_SE_processing_pipeline.sh .



# adapt the script per each sample 
for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' $sample_names)
do
echo $sample_name
sed 's/sample_name_XX/'${sample_name}'/g' sample_name_XX_hybrid_capture_SE_processing_pipeline.sh | sed 's+SPECIFY_local_raw_data_location+'${local_raw_data_location}'+g' | sed 's+SPECIFY_sample_names+'${sample_names}'+g' | sed 's+SPECIFY_mapping_results_dir+'${mapping_results_dir}'+g' > ${sample_name}_hybrid_capture_SE_processing_pipeline.sh
done



# submit a job per each sample
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names} )
do
echo $sample_name
sbatch ${sample_name}_hybrid_capture_SE_processing_pipeline.sh
done






############################################
# Call mutations with MuTect2
############################################



# copy a mutations-calling pipeline to the working directory
cp /research/lab_winter/users/himrichova/resources/scripts/12_SatMut/Hybrid_capture_data_analysis/sample_name_XX_MuTect2.SE_KBM7.sh .



# adapt the script per each sample 
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names} |grep -v wt)
do

sample_name_control=wt

echo $sample_name
sed 's/sample_name_XX/'${sample_name}'/g' sample_name_XX_MuTect2.SE_KBM7.sh | sed 's+SPECIFY_sample_name_control+'${sample_name_control}'+g'  | sed 's+SPECIFY_mapping_results_dir+'${mapping_results_dir}'+g' > ${sample_name}_MuTect2.SE_KBM7.sh
done



# submit a job per each sample
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names} |grep -v wt)
do
echo $sample_name
sbatch ${sample_name}_MuTect2.SE_KBM7.sh
done



############################################