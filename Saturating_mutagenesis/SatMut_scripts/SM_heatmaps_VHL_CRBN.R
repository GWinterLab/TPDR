rm(list=ls())

list.files()
require("pheatmap")
library(gridExtra)

info_text="minmax_norm_heatmaps_"

working_directory <- as.character(commandArgs()[4])
directory_with_SatMut_input_files <- as.character(commandArgs()[5])



setwd(paste(working_directory,sep=""))

annotation <- read.table(paste(directory_with_SatMut_input_files,"/aa_characterization.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE)
aa_Char <- read.table(paste(directory_with_SatMut_input_files,"/chars_ordered.withX.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)

wt_residues_VHL <- read.table(paste(directory_with_SatMut_input_files,"/VHL_wt_with_position.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)
wt_residues_CRBN <- read.table(paste(directory_with_SatMut_input_files,"/CRBN_wt_with_position.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)

selected_residues_VHL=c(64:122,136:151) #VHL
selected_residues_CRBN=c(54:63,86,98:104,145:156,347:361,365,368:403,414:418) #CRBN





#Colors

library(RColorBrewer)

selected_colormap='PuOr'



#######################################
# Load minmax-normalized FC-matrices:
#######################################

FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix <- read.table("FC_matrix__average_replicates.VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV50_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_ARV50_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV5_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_ARV5_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_2ndBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_2ndBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_MZ500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_MZ500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")


FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_CC885_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive <- read.table("FC_matrix__average_replicates.CRBN_SM_CC885_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_CC90009_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive <- read.table("FC_matrix__average_replicates.CRBN_SM_CC90009_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET57_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET57_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")




FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_vs_DMSO_FC_matrix.10000.selectRes.normalized_x_pos.txt")



#FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_1stBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos_neg.txt")

#FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_1stBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.normalized_x_pos.txt")






plot_list=list()

# Specify which conditions should be visualized:
for(sample_comparison_number in c(1,2,6,17:32,36,38,40,42,44,46,48,50,52,54,56,62))  
{
  print(sample_comparison_number)
  
  
  
  if(sample_comparison_number==1)
  { 
    
BET_degraders_mean_pos<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive + FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive+FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive + FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive)/4   
SMARCA_degraders_mean_pos<-(FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive)/1
FC_FC_matrix<-(BET_degraders_mean_pos - SMARCA_degraders_mean_pos)

color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)

comparison_name="VHL: mean(positive log2FC of 4 BET degraders) - mean(positive log2FC of a SMARCA degrader) (500uM:ARV-771,MZ1,DAT548; 2uM:DAT551,ACBI1)"
comparison_file_name="VHL_mean_pos_log2FC_4_BET_degraders_vs_mean_pos_log2FC_1_SMARCA_degraders"

selected_color_scale_values="NEG_POS"

selected_gene="VHL"

#neutral_white_up_to=0.5
neutral_white_up_to=0.2
}
  
  
  
  
  if(sample_comparison_number==2)
  { 
    
    all_degraders_mean_pos<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive + FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive+FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive+FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive+FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive)/5   
    FC_FC_matrix<-all_degraders_mean_pos
    
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    comparison_name="VHL: mean of positive log2FC values (500uM:ARV-771,MZ1,DAT548; 2uM:DAT551,ACBI1)"
    comparison_file_name="VHL_mean_of_positive_log2FC_values_5_degraders"
    
    selected_color_scale_values="POSITIVE"
    
    
    selected_gene="VHL"
    
    #neutral_white_up_to=1
    neutral_white_up_to=0.4
    }
  
  
  if(sample_comparison_number==3)
  { 
    
    all_degraders_mean<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix + FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix+FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix+FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix+FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix)/5   
    FC_FC_matrix<-all_degraders_mean
    
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    comparison_name="VHL: mean of all log2FC values (500uM:ARV-771,MZ1,DAT548; 2uM:DAT551,ACBI1)"
    comparison_file_name="VHL_mean_of_all_log2FC_values_5_degraders"
    
    
    selected_color_scale_values="NEG_POS"

    
    selected_gene="VHL"
    
    neutral_white_up_to=1
    }

  
  
  
  
 
  if(sample_comparison_number==5)
  { 
    all_degraders_mean<-(FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix+FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix+FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix+FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix)/4
    
    FC_FC_matrix<-all_degraders_mean
    
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    comparison_name="CRBN: mean of all log2FC values (4 degraders: CC885,CC90009,dBET57,dBET6)"
    comparison_file_name="CRBN_mean_of_all_log2FC_values_4_degraders"
    
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"

    neutral_white_up_to=0.5
  }
  
  
  if(sample_comparison_number==6)
  { 
    all_degraders_mean_pos<-(FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive+FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive+FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive+FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive)/4
    
    FC_FC_matrix<-all_degraders_mean_pos
    
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    comparison_name="CRBN: mean of positive log2FC values (4 degraders: CC885,CC90009,dBET57,dBET6)"
    comparison_file_name="CRBN_mean_of_pos_log2FC_values_4_degraders"
    
    
    selected_color_scale_values="POSITIVE"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.2
  }
  
  
 
  
  
  ###
  
 
  
  
  
  
  
  #####
  # diff CRBN
  
 
  if(sample_comparison_number==17)
  { 
    FC_FC_matrix<-(FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive - FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: CC885 - CC90009"
    comparison_file_name="CRBN_diff_CC885_vs_CC90009"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  
  if(sample_comparison_number==18)
  { 
    FC_FC_matrix<-(FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive - FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: CC885 - dBET57"
    comparison_file_name="CRBN_diff_CC885_vs_dBET57"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  
  if(sample_comparison_number==19)
  { 
    FC_FC_matrix<-(FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive - FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: CC885 - dBET6"
    comparison_file_name="CRBN_diff_CC885_vs_dBET6"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==20)
  { 
    FC_FC_matrix<-(FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive - FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: CC90009 - dBET57"
    comparison_file_name="CRBN_diff_CC90009_vs_dBET57"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==21)
  { 
    FC_FC_matrix<-(FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive - FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: CC90009 - dBET6"
    comparison_file_name="CRBN_diff_CC90009_vs_dBET6"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    #neutral_white_up_to=0.5 
    neutral_white_up_to=0.2
  }
  
  
  if(sample_comparison_number==22)
  { 
    FC_FC_matrix<-(FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive - FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive)
    comparison_name="CRBN: dBET57 - dBET6"
    comparison_file_name="CRBN_diff_dBET57_vs_dBET6"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
    
    
    
    
    ### VHL_SM_rep1rep2:
    
    if(sample_comparison_number==23)
    { 
      FC_FC_matrix<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive - FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ARV_500nM - MZ1_500nM"
      
      comparison_file_name="VHL_diff_log2FC_ARV_500nM_vs_log2FC_MZ1_500nM"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      #neutral_white_up_to=0.5
      neutral_white_up_to=0.2
    }
    
    
    
    if(sample_comparison_number==24)
    {
      FC_FC_matrix<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive - FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ARV_500nM - ACBI1_2uM"
      comparison_file_name="VHL_diff_log2FC_ARV_500nM_vs_log2FC_ACBI1_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
      
    }
    

    
    
    if(sample_comparison_number==25)
    {
      FC_FC_matrix<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive - FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ARV_500nM - DAT548_500nM"
      comparison_file_name="VHL_diff_log2FC_ARV_500nM_vs_log2FC_DAT548_500nM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    
    
    
    
    
    if(sample_comparison_number==26)
    {
      FC_FC_matrix<-(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive - FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ARV_500nM - DAT551_2uM"
      comparison_file_name="VHL_diff_log2FC_ARV_500nM_vs_log2FC_DAT551_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
      neutral_white_up_to=0.2
      }
    

    
    
    
    if(sample_comparison_number==27)
    {
      FC_FC_matrix<-(FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive - FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: MZ1_500nM - ACBI1_2uM"
      comparison_file_name="VHL_diff_log2FC_MZ1_500nM_vs_log2FC_ACBI1_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    

    
    
    if(sample_comparison_number==28)
    {
      FC_FC_matrix<-(FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive - FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: MZ1_500nM - DAT548_500nM"
      comparison_file_name="VHL_diff_log2FC_MZ1_500nM_vs_log2FC_DAT548_500nM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    
    
    
    if(sample_comparison_number==29)
    {
      FC_FC_matrix<-(FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive - FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: MZ1_500nM - DAT551_2uM"
      comparison_file_name="VHL_diff_log2FC_MZ1_500nM_vs_log2FC_DAT551_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    


    
    
    if(sample_comparison_number==30)
    {
      FC_FC_matrix<-(FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive - FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ACBI1_2uM - DAT548_500nM"
      comparison_file_name="VHL_diff_log2FC_ACBI1_2uM_7d_vs_log2FC_DAT548_500nM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    
    
    
    if(sample_comparison_number==31)
    {
      FC_FC_matrix<-(FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive - FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: ACBI1_2uM - DAT551_2uM"
      comparison_file_name="VHL_diff_log2FC_ACBI1_2uM_7d_vs_log2FC_DAT551_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    
    if(sample_comparison_number==32)
    {
      FC_FC_matrix<-(FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive - FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive)
      comparison_name="VHL: DAT548_500nM - DAT551_2uM"
      comparison_file_name="VHL_diff_log2FC_DAT548_500nM_7d_vs_log2FC_DAT551_2uM_7d"
        color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    
    ###
    
   
    
    
    
#### CRBN drug vs DMSO
  
  
  
  
  if(sample_comparison_number==35)
  {
    FC_FC_matrix<-FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix
    comparison_name="CRBN: CC885 vs DMSO"
    comparison_file_name="CRBN_NegPosValues_CC885_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==36)
  {
    FC_FC_matrix<-FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive
    comparison_name="CRBN: CC885 vs DMSO"
    comparison_file_name="CRBN_PosValues_CC885_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.2
  }
  
  
  if(sample_comparison_number==37)
  {
    FC_FC_matrix<-FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix
    comparison_name="CRBN: CC90009 vs DMSO"
    comparison_file_name="CRBN_NegPosValues_CC90009_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==38)
  {
    FC_FC_matrix<-FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive
    comparison_name="CRBN: CC90009 vs DMSO"
    comparison_file_name="CRBN_PosValues_CC90009_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==39)
  {
    FC_FC_matrix<-FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix
    comparison_name="CRBN: dBET57 vs DMSO"
    comparison_file_name="CRBN_NegPosValues_dBET57_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==40)
  {
    FC_FC_matrix<-FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive
    comparison_name="CRBN: dBET57 vs DMSO"
    comparison_file_name="CRBN_PosValues_dBET57_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="CRBN"
    
    #neutral_white_up_to=0.5
    neutral_white_up_to=0.2
  }
    
  
  if(sample_comparison_number==41)
  {
    FC_FC_matrix<-FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix
    comparison_name="CRBN: dBET6 vs DMSO"
    comparison_file_name="CRBN_NegPosValues_dBET6_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==42)
  {
    FC_FC_matrix<-FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive
    comparison_name="CRBN: dBET6 vs DMSO"
    comparison_file_name="CRBN_PosValues_dBET6_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.2
  }
  
#### VHL individual drug vs DMSO
  
  if(sample_comparison_number==43)
  {
    FC_FC_matrix<-FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix
    comparison_name="VHL: ARV771_500 vs DMSO"
    comparison_file_name="VHL_NegPosValues_ARV500_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    #selected_color_scale_values="POS"
    selected_gene="VHL"
    
    #neutral_white_up_to=0.5
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==44)
  {
    FC_FC_matrix<-FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive
    comparison_name="VHL: ARV771_500 vs DMSO"
    comparison_file_name="VHL_PosValues_ARV500_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="VHL"
    
    #neutral_white_up_to=0.5
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==45)
  {
    FC_FC_matrix<-FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix
    comparison_name="VHL: MZ1 vs DMSO"
    comparison_file_name="VHL_NegPosValues_MZ1_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==46)
  {
    FC_FC_matrix<-FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive
    comparison_name="VHL: MZ1 vs DMSO"
    comparison_file_name="VHL_PosValues_MZ1_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==47)
  {
    FC_FC_matrix<-FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix
    comparison_name="VHL: DAT548 vs DMSO"
    comparison_file_name="VHL_NegPosValues_DAT548_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==48)
  {
    FC_FC_matrix<-FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive
    comparison_name="VHL: DAT548 vs DMSO"
    comparison_file_name="VHL_PosValues_DAT548_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==49)
  {
    FC_FC_matrix<-FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix
    comparison_name="VHL: DAT551 vs DMSO"
    comparison_file_name="VHL_NegPosValues_DAT551_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==50)
  {
    FC_FC_matrix<-FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive
    comparison_name="VHL: DAT551 vs DMSO"
    comparison_file_name="VHL_PosValues_DAT551_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.2
  }
  
  if(sample_comparison_number==51)
  {
    FC_FC_matrix<-FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix
    comparison_name="VHL: ACBI1 vs DMSO"
    comparison_file_name="VHL_NegPosValues_ACBI1_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="NEG_POS"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.5
  }
  
  if(sample_comparison_number==52)
  {
    FC_FC_matrix<-FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive
    comparison_name="VHL: ACBI1 vs DMSO"
    comparison_file_name="VHL_PosValues_ACBI1_vs_DMSO"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="VHL"
    
    neutral_white_up_to=0.2
  }
  
    

    if(sample_comparison_number==53)
    {
      FC_FC_matrix<-FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix
      comparison_name="VHL: ARV771_50 vs DMSO"
      comparison_file_name="VHL_NegPosValues_ARV50_vs_DMSO"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    if(sample_comparison_number==54)
    {
      FC_FC_matrix<-FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix_positive
      comparison_name="VHL: ARV771_50 vs DMSO"
      comparison_file_name="VHL_PosValues_ARV50_vs_DMSO"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="POSITIVE"
      
      selected_gene="VHL"
      
      #neutral_white_up_to=0.5
      neutral_white_up_to=0.2
    }
    
    if(sample_comparison_number==55)
    {
      FC_FC_matrix<-FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix
      comparison_name="VHL: ARV771_5 vs DMSO"
      comparison_file_name="VHL_NegPosValues_ARV5_vs_DMSO"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="NEG_POS"
      
      selected_gene="VHL"
      
      neutral_white_up_to=0.5
    }
    
    if(sample_comparison_number==56)
    {
      FC_FC_matrix<-FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix_positive
      comparison_name="VHL: ARV771_5 vs DMSO"
      comparison_file_name="VHL_PosValues_ARV5_vs_DMSO"
      color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
      
      selected_color_scale_values="POSITIVE"
      
      selected_gene="VHL"
      
      #neutral_white_up_to=0.5
      neutral_white_up_to=0.2
    }
      
     
        
        ####
        # VHL - 1st vs 2nd batch DMSO
        
        if(sample_comparison_number==61)
        {
          FC_FC_matrix<-FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix
          comparison_name="VHL: ARV771_500 vs DMSO2nd"
          comparison_file_name="VHL_NegPosValues_ARV500_vs_DMSO2nd"
          color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
          
          selected_color_scale_values="NEG_POS"
          
          selected_gene="VHL"
          
          neutral_white_up_to=0.5
        }
        
        if(sample_comparison_number==62)
        {
          FC_FC_matrix<-FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive
          comparison_name="VHL: ARV771_500 vs DMSO2nd"
          comparison_file_name="VHL_PosValues_ARV500_vs_DMSO2nd"
          color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
          
          selected_color_scale_values="POSITIVE"
          
          selected_gene="VHL"
          
          neutral_white_up_to=0.2
        }
        
        
       
        
        ####
        # VHL - DMSO 1st vs 2nd batch DMSO
        
        if(sample_comparison_number==63)
        {
          FC_FC_matrix<-FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix
          comparison_name="VHL: DMSO1st vs DMSO2nd"
          comparison_file_name="VHL_NegPosValues_DMSO1st_vs_DMSO2nd"
          color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
          
          selected_color_scale_values="NEG_POS"
          
          selected_gene="VHL"
          
          neutral_white_up_to=0.5
        }
        
        if(sample_comparison_number==64)
        {
          FC_FC_matrix<-FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix_positive
          comparison_name="VHL: DMSO1st vs DMSO2nd"
          comparison_file_name="VHL_PosValues_DMSO1st_vs_DMSO2nd"
          color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
          
          selected_color_scale_values="POSITIVE"
          
          selected_gene="VHL"
          
          neutral_white_up_to=0.2
          
          
        } 
  
  
  
  
  if(sample_comparison_number==65)
  {
    FC_FC_matrix<-FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive
    comparison_name="CRBN: dBET6 vs DMSO (1stBatch)"
    comparison_file_name="CRBN_PosValues_dBET6_vs_DMSO_1stBatch"
    color_scale_max_value=max(abs(FC_FC_matrix),na.rm=T)
    
    selected_color_scale_values="POSITIVE"
    
    selected_gene="CRBN"
    
    neutral_white_up_to=0.5
    
    
  } 
############################################################ 
# Values, colors
############################################################
  #neutral_white_up_to=0.2

  color_scale_max_value=color_scale_max_value+0.01      
   
   # set_breaks_color_scale and colors for heatmaps with POSITIVE values only:
   set_breaks_color_scale_POSITIVE=seq(0, color_scale_max_value, 0.01)
   
   # set_breaks_color_scale and colors for heatmaps with ALL values:
   set_breaks_color_scale_NEG_POS=seq(-color_scale_max_value, color_scale_max_value, 0.01)
   
   if(selected_gene=="VHL")
   { 
   colors_NEG<-brewer.pal(n=11,name = selected_colormap)[1:6]
   colors_POS<-brewer.pal(n=11,name = selected_colormap)[6:11]
  
   wt_residues<-wt_residues_VHL
   selected_residues<-selected_residues_VHL
   
   
    }
   
   
   if(selected_gene=="CRBN")
   { 
     colors_NEG<-brewer.pal(n=11,name = selected_colormap)[1:6]
     colors_POS<-brewer.pal(n=11,name = selected_colormap)[6:11]

   
   wt_residues<-wt_residues_CRBN
   selected_residues<-selected_residues_CRBN
   }
   
   
   set_breaks_color_scale_NEG=seq(-color_scale_max_value, -neutral_white_up_to, 0.01)
   set_breaks_color_scale_POS=seq(neutral_white_up_to, color_scale_max_value, 0.01)
   set_breaks_color_scale_NEUTRAL=seq(-neutral_white_up_to, neutral_white_up_to, 0.01)
   set_breaks_color_scale_POS_NEUTRAL=seq(0, neutral_white_up_to, 0.01)
   
   
   my_palette_NEG <- colorRampPalette(c(colors_NEG))(n = length(set_breaks_color_scale_NEG)-1)
   my_palette_POS <- colorRampPalette(c(colors_POS))(n = length(set_breaks_color_scale_POS)-1)
   my_palette_NEUTRAL <- colorRampPalette("#FFFFFF")(n = length(set_breaks_color_scale_NEUTRAL)-1)
   my_palette_POS_NEUTRAL <- colorRampPalette("#FFFFFF")(n = length(set_breaks_color_scale_POS_NEUTRAL)-1)
   
   my_palette_NEG_POS <- c(my_palette_NEG,my_palette_NEUTRAL,my_palette_POS)
   
   my_palette_POSITIVE <- c(my_palette_POS_NEUTRAL,my_palette_POS)
   
   
   
   if(selected_color_scale_values=="NEG_POS")
   { 
   set_breaks_color_scale=set_breaks_color_scale_NEG_POS
   my_palette=my_palette_NEG_POS
   }
   
   
   if(selected_color_scale_values=="POSITIVE")
   { 
     set_breaks_color_scale=set_breaks_color_scale_POSITIVE
     my_palette=my_palette_POSITIVE
   }
  
############################################################ 
# Heatmap and stackBar
############################################################

y=data.matrix(FC_FC_matrix)[,]
y=data.matrix(FC_FC_matrix)
y2<-y[,c(as.character(aa_Char$V1))]

write.table(y2,paste(info_text,"_",comparison_file_name,"_FC_matrix.",".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)

mean_value_per_residue<-apply(y2,1,mean,na.rm=T)
median_value_per_residue<-apply(y2,1,median,na.rm=T)

write.table(mean_value_per_residue,paste(info_text,"_",comparison_file_name,"_MEAN_value_per_residue",".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=F)
write.table(median_value_per_residue,paste(info_text,"_",comparison_file_name,"_MEDIAN_value_per_residue",".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=F)

H_positive_lfc<-pheatmap(t(y2)
            ,na_col="gray98"
            ,border_color = "white"
            ,col = my_palette
            ,Rowv=NA,Colv=NA
            ,legend=T
            ,cellwidth=9,cellheight=10
            ,cluster_rows=FALSE
            ,cluster_cols=FALSE
            ,ylab="Amino acid substitution"
            ,margins =c(3,3),cexCol=1
            ,scale = "none"
            ,na.rm=F
            #,main=comparison_name
            ,labels_row=paste(aa_Char$V1, " (",aa_Char$V2,")",sep="")
            ,labels_col = paste(t(wt_residues)[selected_residues],sep="")
            ,breaks=set_breaks_color_scale
)

plot_H_positive_lfc = H_positive_lfc[[4]]

dev.off()

sorted_heatmap<-c()
for(i in c(1:dim(y2)[1]))
{
  sorted_y2_i <- sort(y2[i,], decreasing = FALSE, na.last = FALSE)
  sorted_heatmap<-cbind(sorted_heatmap,sorted_y2_i)
}


H_positive_lfc_stackBar<-pheatmap((sorted_heatmap)
                                  ,na_col="gray98"
                                  ,border_color = NA
                                  ,col = my_palette
                                  ,Rowv=NA,Colv=NA
                                  ,cellwidth=9
                                  ,cellheight=1.5
                                  ,legend=T
                                  ,cluster_rows=FALSE
                                  ,cluster_cols=FALSE
                                  ,ylab="Amino acid substitution"
                                  ,margins =c(3,3),cexCol=1
                                  ,scale = "none"
                                  ,na.rm=F,main=paste(comparison_name,sep="")
                                  ,breaks=set_breaks_color_scale
                                  ,labels_row=paste(aa_Char$V1, " (",aa_Char$V2,")",sep="")
                                  ,labels_col = ""
)

plot_H_positive_lfc_stackBar = H_positive_lfc_stackBar[[4]]

dev.off()


plot_list[["a"]]=plot_H_positive_lfc_stackBar
plot_list[["b"]]=plot_H_positive_lfc
pdf(paste(info_text,"_",comparison_file_name,"_FC_matrix.",".pdf",sep="")
    ,width=12,height=7)
g<-do.call(grid.arrange,plot_list)
dev.off()

  }
 

