###################################################################################################
# Deep mutational scanning pipeline - from raw data (bam files) preprocessing to heatmaps
###################################################################################################

# specify a path to the directory with SatMut input files:
directory_with_SatMut_input_files=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/SatMut_input_files

# specify a path to the directory with SatMut scripts:
directory_with_SatMut_scripts=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/SatMut_scripts

# specify a path to the directory with raw (unmapped) bam files:
local_raw_data_location=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/11_SM_and_Hybrid_capture_FINAL/1_SM_VHL_and_CRBN/raw_bam_files

# specify a sample annotation file:
sample_names=${directory_with_SatMut_input_files}/sample_annotation.txt

# specify your output directory:
mapping_results_dir=/research/lab_winter/users/himrichova/1_Projects/12_TPDR/12_SM_and_Hybrid_capture_rerun_test/1_SM_VHL_and_CRBN/mapping_results





# create directories for data processing results

mkdir -p ${mapping_results_dir}
mkdir -p ${mapping_results_dir}/fastq
mkdir -p ${mapping_results_dir}/fastqc
mkdir -p ${mapping_results_dir}/bam
mkdir -p ${mapping_results_dir}/mut_calling
mkdir -p ${mapping_results_dir}/mut_calling/Freq_aaCounts



cd ${mapping_results_dir}

# create wt residues files, including the number per each residue: 
cat ${directory_with_SatMut_input_files}/wt_residues_VHL.txt |tr "\t" "\n" |awk -F "\t" '{print NR "\t" $1}'|sed 's/*/X/g' > ${mapping_results_dir}/wt_residues_VHL.t.txt 

cat ${directory_with_SatMut_input_files}/wt_residues_CRBN.txt |tr "\t" "\n" |awk -F "\t" '{print NR "\t" $1}'|sed 's/*/X/g' > ${mapping_results_dir}/wt_residues_CRBN.t.txt 




# copy a sample_name_XX_SatMut.sh script to your directory:

cp ${directory_with_SatMut_scripts}/sample_name_XX_SatMut.sh .

 
# adapt the script per each sample - add information from the sample annotation file:
for sample_name in $( awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
sed 's/sample_name_XX/'${sample_name}'/g' sample_name_XX_SatMut.sh | sed 's+SPECIFY_local_raw_data_location+'${local_raw_data_location}'+g' | sed 's+SPECIFY_sample_names+'${sample_names}'+g' | sed 's+SPECIFY_mapping_results_dir+'${mapping_results_dir}'+g' | sed 's+SPECIFY_directory_with_SatMut_input_files+'${directory_with_SatMut_input_files}'+g' | sed 's+SPECIFY_directory_with_SatMut_scripts+'${directory_with_SatMut_scripts}'+g' > ${sample_name}_SatMut.sh
done



# submit a job per each sample:
for sample_name in $(awk -F "\t" '{if(NR>1){print $1}}' ${sample_names})
do
echo $sample_name
sbatch ${sample_name}_SatMut.sh
done



module load R/4.2.1-foss-2022a

# generate Fold-change matrices Drug/DMSO:

R --vanilla --args ${mapping_results_dir}/mut_calling/Freq_aaCounts ${directory_with_SatMut_input_files} < ${directory_with_SatMut_scripts}/SM_all_samples_replicates_FC_matrix.R


# generate heatmaps

R --vanilla --args ${mapping_results_dir}/mut_calling/Freq_aaCounts ${directory_with_SatMut_input_files} < ${directory_with_SatMut_scripts}/SM_heatmaps_VHL_CRBN.R


# generate correlation plots
R --vanilla --args ${mapping_results_dir}/mut_calling/Freq_aaCounts ${directory_with_SatMut_input_files} < ${directory_with_SatMut_scripts}/SM_correlations_VHL_CRBN.R


