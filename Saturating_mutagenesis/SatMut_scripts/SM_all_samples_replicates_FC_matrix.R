rm(list=ls())
require("pheatmap")
require("ggpubr")
library(gridExtra)

working_directory <- as.character(commandArgs()[4])
directory_with_SatMut_input_files <- as.character(commandArgs()[5])

setwd(paste(working_directory,sep=""))

annotation <- read.table(paste(directory_with_SatMut_input_files,"/aa_characterization.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE)
aa_Char <- read.table(paste(directory_with_SatMut_input_files,"/chars_ordered.withX.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)

wt_residues_VHL <- read.table(paste(directory_with_SatMut_input_files,"/VHL_wt_with_position.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)
wt_residues_CRBN <- read.table(paste(directory_with_SatMut_input_files,"/CRBN_wt_with_position.txt",sep=""),sep="\t", header=F, dec=".", fill = TRUE)




selected_residues_VHL=c(64:122,136:151) #VHL
selected_residues_CRBN=c(54:63,86,98:104,145:156,347:361,365,368:403,414:418) #CRBN



info_text<-"FC_matrix_"

#cut-off for filtered matrix:
test_alternative="two.sided"
signif_cutoff_to_show_residues_in_heatmaps=0.05



total_reads_denominator=10000 #final threshold to be used
fdr_n=1




for(sample_comparison_number in c(1:17))
{
  print(sample_comparison_number)

  #1
  if(sample_comparison_number==1)
  {
    comparison_name="VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO"
    
    dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
    treat_rep1_name<-"S_20_SM_r1_ACBI1_2_7d_S66511"
    
    dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
    treat_rep2_name<-"S_9_SM_r2_ACBI1_2_7d_S66512"
    
    
    selected_gene="VHL"

  }

  
  #2
  if(sample_comparison_number==2)
  {
    comparison_name="VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO"
    
    
    dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
    treat_rep1_name<-"S_13_SM_r1_DAT548_500_7d_S66503"
    
    dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
    treat_rep2_name<-"S_2_SM_r2_DAT548_500_7d_S66504"
    
    selected_gene="VHL"
  }
  
  #3
  if(sample_comparison_number==3)
  {
    comparison_name="VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO"
    
    
    dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
    treat_rep1_name<-"S_17_SM_r1_DAT551_2_7d_S66507"
    
    dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
    treat_rep2_name<-"S_6_SM_r2_DAT551_2_7d_S66508"
    
    selected_gene="VHL"
  }
  
  
      
      if(sample_comparison_number==4)
      {
        comparison_name="VHL_ARV500_vs_DMSO_1stBatch"
        
        
        dmso_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"
        treat_rep1_name<-"S_47_VHL_SM_r1_ARV500_S60246"
        
        
        dmso_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"
        treat_rep2_name<-"S_54_VHL_SM_r2_ARV500_S60257"
        
        selected_gene="VHL"
      }
      
      
      if(sample_comparison_number==5)
      {
        comparison_name="VHL_ARV500_vs_DMSO_2ndBatch"
        
        
        dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
        treat_rep1_name<-"S_47_VHL_SM_r1_ARV500_S60246"
        
        
        dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
        treat_rep2_name<-"S_54_VHL_SM_r2_ARV500_S60257"
        
        selected_gene="VHL"
      }
      
      
      
      if(sample_comparison_number==6)
      {
        comparison_name="VHL_MZ500_vs_DMSO_1stBatch"
        
        
        dmso_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"
        treat_rep1_name<-"S_50_VHL_SM_r1_MZ500_S60241"        
        
        dmso_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"
        treat_rep2_name<-"S_57_VHL_SM_r2_MZ500_S60259"
        
        selected_gene="VHL"
      }
      
      
      if(sample_comparison_number==7)
      {
        comparison_name="VHL_MZ500_vs_DMSO_2ndBatch"
        
        
        dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
        treat_rep1_name<-"S_50_VHL_SM_r1_MZ500_S60241"
        
        dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
        treat_rep2_name<-"S_57_VHL_SM_r2_MZ500_S60259"
        
        selected_gene="VHL"
      }
      
      
      
      if(sample_comparison_number==8)
      {
        comparison_name="VHL_ARV50_vs_DMSO_1stBatch"
        
        
        dmso_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"
        treat_rep1_name<-"S_46_VHL_SM_r1_ARV50_S60245"
        
        
        dmso_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"
        treat_rep2_name<-"S_53_VHL_SM_r2_ARV50_S60244"
        
        selected_gene="VHL"
      }
      
      if(sample_comparison_number==9)
      {
        comparison_name="VHL_ARV5_vs_DMSO_1stBatch"
        
        
        dmso_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"
        treat_rep1_name<-"S_45_VHL_SM_r1_ARV5_S60236"
        
        
        dmso_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"
        treat_rep2_name<-"S_52_VHL_SM_r2_ARV5_S60243"
        
        selected_gene="VHL"
      }
      
      
  
  
  if(sample_comparison_number==10)
  {
    comparison_name="VHL_DMSO_1stBatch_vs_DMSO_2ndBatch"
    
    
    dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
    treat_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"

    
    dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
    treat_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"

    selected_gene="VHL"
  }
  
  
  
  if(sample_comparison_number==11)
  {
    comparison_name="VHL_DMSO_2ndBatch_vs_DMSO_1stBatch"
    
    
    dmso_rep1_name<-"S_44_VHL_SM_r1_DMSO_S60235"
    treat_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"

    
    dmso_rep2_name<-"S_51_VHL_SM_r2_DMSO_S60242"
    treat_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"

    selected_gene="VHL"
  }
  
  
  if(sample_comparison_number==12)
  {
    comparison_name="VHL_DMSO_2ndBatch_vs_Library"
    
    
    dmso_rep1_name<-"P1_SM_r1_Lib_S66497"
    treat_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"

    
    dmso_rep2_name<-"P2_SM_r2_Lib_S66498"
    treat_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"

    selected_gene="VHL"
  }
  
  
  
  if(sample_comparison_number==13)
  {
    comparison_name="VHL_Library_vs_DMSO_2ndBatch"
    
    
    dmso_rep1_name<-"S_11_SM_r1_DMSO_7d_S66509"
    treat_rep1_name<-"P1_SM_r1_Lib_S66497"

    
    dmso_rep2_name<-"S_10_SM_r2_DMSO_7d_S66510"
    treat_rep2_name<-"P2_SM_r2_Lib_S66498"

    selected_gene="VHL"
  }
  
  

      
      
      
      if(sample_comparison_number==14)
      {
        comparison_name="CRBN_SM_CC885_vs_DMSO"
        
        dmso_rep1_name<-"DMSO_r1_S75013"
        treat_rep1_name<-"CC885_r1_S75011"
        
        dmso_rep2_name<-"DMSO_r2_S75015"
        treat_rep2_name<-"CC885_r2_S75008"
        
        dmso_rep3_name<-"DMSO_r3_S75017"
        treat_rep3_name<-"CC885_r3_S75009"
        
        
        selected_gene="CRBN"
        
      }
      
      
      
      
      if(sample_comparison_number==15)
      {
        comparison_name="CRBN_SM_CC90009_vs_DMSO"
        
        dmso_rep1_name<-"DMSO_r1_S75013"
        treat_rep1_name<-"CC90009_r1_S75024"
        
        dmso_rep2_name<-"DMSO_r2_S75015"
        treat_rep2_name<-"CC90009_r2_S75025"
        
        dmso_rep3_name<-"DMSO_r3_S75017"
        treat_rep3_name<-"CC90009_r3_S75022"
        
        
        selected_gene="CRBN"
        
      }  
      
      
      

      if(sample_comparison_number==16)
      {
        comparison_name="CRBN_SM_dBET57_vs_DMSO"
        
        dmso_rep1_name<-"DMSO_r1_S75013"
        treat_rep1_name<-"dBET57_r1_S75012"
        
        dmso_rep2_name<-"DMSO_r2_S75015"
        treat_rep2_name<-"dBET57_r2_S75014"
        
        dmso_rep3_name<-"DMSO_r3_S75017"
        treat_rep3_name<-"dBET57_r3_S75010"
        
        
        selected_gene="CRBN"
        
      }
      

      if(sample_comparison_number==17)
      {
        comparison_name="CRBN_SM_dBET6_vs_DMSO"
        
        dmso_rep1_name<-"DMSO_r1_S75013"
        treat_rep1_name<-"dBET6_r1_S75019"
        
        dmso_rep2_name<-"DMSO_r2_S75015"
        treat_rep2_name<-"dBET6_r2_S75016"
        
        dmso_rep3_name<-"DMSO_r3_S75017"
        treat_rep3_name<-"dBET6_r3_S75018"
        
        
        selected_gene="CRBN"
        
        
      }
      
      

  if(sample_comparison_number==18)
  {
    comparison_name="CRBN_SM_dBET6_1stBatch_vs_DMSO_1stBatch"
    
    dmso_rep1_name<-"S_58_CRBN_SM_r1_DMSO"
    treat_rep1_name<-"S_61_CRBN_SM_r1_dBET500_S60226"
    
    dmso_rep2_name<-"S_58_CRBN_SM_r1_DMSO"
    treat_rep2_name<-"S_61_CRBN_SM_r1_dBET500_S60226"
    
    dmso_rep3_name<-"S_58_CRBN_SM_r1_DMSO"
    treat_rep3_name<-"S_61_CRBN_SM_r1_dBET500_S60226"
    
    
    selected_gene="CRBN"
    

    
  }
  
      if(selected_gene=="VHL")
      { 
        
        wt_residues<-wt_residues_VHL
        selected_residues<-selected_residues_VHL
        
        
      }
      
      
      if(selected_gene=="CRBN")
      { 
        
        wt_residues<-wt_residues_CRBN
        selected_residues<-selected_residues_CRBN
      }
      

      
      
      
      if(selected_gene=="VHL")
      {      
      
  dmso_rep1<-c()
  dmso_rep1_raw<-c()
  treat_rep1<-c()
  dmso_rep2<-c()
  dmso_rep2_raw<-c()
  treat_rep2<-c()
  
  wt_counts_dmso_rep1<-c()
  
dmso_rep1=as.data.frame(read.table(paste(dmso_rep1_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
dmso_rep1_raw=as.data.frame(read.table(paste("../",dmso_rep1_name,".gatk.aaCounts",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))

treat_rep1=as.data.frame(read.table(paste(treat_rep1_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))



dmso_rep2=as.data.frame(read.table(paste(dmso_rep2_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
dmso_rep2_raw=as.data.frame(read.table(paste("../",dmso_rep2_name,".gatk.aaCounts",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))


treat_rep2=as.data.frame(read.table(paste(treat_rep2_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))



wt_counts_dmso_rep1<-dmso_rep1$counts_wt_res


all_values<-c()
FC_values<-c()
FDR_values<-c()

for(i in c(1:nrow(dmso_rep1)))

  { 
  
  #print(i)
  
  for(j in c(1:(ncol(dmso_rep1[1:21]))))
  {
    #print(j)
    cutoff_rep1<-dmso_rep1[i,22]/total_reads_denominator
    cutoff_rep2<-dmso_rep2[i,22]/total_reads_denominator
    
    wt_counts_dmso_rep1[i]
    
    

      
    if(dmso_rep1_raw[i,j] < cutoff_rep1 | dmso_rep2_raw[i,j] < cutoff_rep2 )
    {fdr_value <- NA 
    log10_p <- NA
    log10_fdr <- NA
    mean_of_x<-NA
    mean_of_y<-NA
    log2_fold_change <- NA
    mean_x_and_y <- NA
    }
      
    if(dmso_rep1_raw[i,j] >= cutoff_rep1 & dmso_rep2_raw[i,j] >= cutoff_rep2 )
    {

      ttest<-t.test(c((dmso_rep1[i,j]+0.00001),(dmso_rep2[i,j]+0.00001)),c((treat_rep1[i,j]+0.00001),(treat_rep2[i,j]+0.00001)),alternative=test_alternative)


      
      log2_fold_change<-log2(mean(c((treat_rep1[i,j]+0.00001),(treat_rep2[i,j]+0.00001))) / mean(c((dmso_rep1[i,j]+0.00001),(dmso_rep2[i,j]+0.00001))) )
      
      p_value<-ttest$p.value
      
      log10_p <- -log10(p_value)
    
      
      fdr_value<-p.adjust(p_value,method="fdr",n=fdr_n)
      log10_fdr <- -log10(fdr_value)
      
      mean_of_x<-ttest$estimate[[1]]
      mean_of_y<-ttest$estimate[[2]]
      mean_x_and_y <- ttest$estimate[[1]]+ttest$estimate[[2]]
    }
    
    
    FC_values<-rbind(FC_values,cbind(i,j,mean_of_y, mean_of_x,log2_fold_change))
    FDR_values<-rbind(FDR_values,cbind(i,j,mean_of_y, mean_of_x,fdr_value))
    all_values<-rbind(all_values,cbind(i,j,mean_of_y, mean_of_x,mean_x_and_y,log2_fold_change,log10_p,log10_fdr,fdr_value))
  
    }
}








head(all_values)
dim(all_values)

signif_cutoff_to_show_residues_in_heatmaps
all_filtered_for_signif<-rbind( cbind(all_values[which(all_values[,9]<=signif_cutoff_to_show_residues_in_heatmaps),],all_values[which(all_values[,9]<=signif_cutoff_to_show_residues_in_heatmaps),6]),  cbind(all_values[which(all_values[,9]>signif_cutoff_to_show_residues_in_heatmaps | is.na(all_values[,9]) ),],NA))

sort_all_filtered_for_signif<-all_filtered_for_signif[order(as.numeric(all_filtered_for_signif[,1]),as.numeric(all_filtered_for_signif[,2])) , ]
head(sort_all_filtered_for_signif)


filtered_FC_matrix<-matrix(sort_all_filtered_for_signif[,10], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)

colnames(filtered_FC_matrix)<-colnames(dmso_rep1[1:21])
rownames(filtered_FC_matrix)<-rownames(dmso_rep1)

filtered_FC_matrix<-filtered_FC_matrix[selected_residues,]
head(filtered_FC_matrix)



write.table(filtered_FC_matrix,paste(info_text,"_average_replicates.",comparison_name,"_filtered_FC_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)

write.table(all_values,paste(info_text,"_all_values.",comparison_name,"_",total_reads_denominator,".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)







FC_matrix<-matrix(FC_values[,5], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)

colnames(FC_matrix)<-colnames(dmso_rep1[1:21])
rownames(FC_matrix)<-rownames(dmso_rep1)

FC_matrix<-FC_matrix[selected_residues,]

#to be used:
write.table(FC_matrix,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)


FC_matrix_positive <-FC_matrix
FC_matrix_positive[FC_matrix_positive<0] <- 0

#Normalized Data pos_neg
x_pos_neg=FC_matrix
x_pos_only<-x_pos_neg
x_pos_only[x_pos_only<0] <- 0

x_neg_only<-x_pos_neg
x_neg_only[x_neg_only>0] <- 0
x_neg_only <- abs(x_neg_only)

normalized_x_pos_only = (x_pos_only-min(x_pos_only,na.rm=T))/(max(x_pos_only,na.rm=T)-min(x_pos_only,na.rm=T))
normalized_x_neg_only = (x_neg_only-min(x_neg_only,na.rm=T))/(max(x_neg_only,na.rm=T)-min(x_neg_only,na.rm=T))

normalized_x_pos_neg<- normalized_x_pos_only-normalized_x_neg_only
write.table(normalized_x_pos_neg,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.normalized_x_pos_neg.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)


#Normalized Data pos
x_pos=FC_matrix_positive
normalized_x_pos = (x_pos-min(x_pos,na.rm=T))/(max(x_pos,na.rm=T)-min(x_pos,na.rm=T))

write.table(normalized_x_pos,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.normalized_x_pos.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)







Mean_of_x_matrix <- matrix(FC_values[,4], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)

colnames(Mean_of_x_matrix)<-colnames(dmso_rep1[1:21])
rownames(Mean_of_x_matrix)<-rownames(dmso_rep1)

Mean_of_x_matrix<-Mean_of_x_matrix[selected_residues,]

write.table(Mean_of_x_matrix,paste(info_text,"_average_replicates.",dmso_rep1_name,dmso_rep2_name,"_log2_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)


Mean_of_y_matrix <- matrix(FC_values[,3], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)

colnames(Mean_of_y_matrix)<-colnames(dmso_rep1[1:21])
rownames(Mean_of_y_matrix)<-rownames(dmso_rep1)

Mean_of_y_matrix<-Mean_of_y_matrix[selected_residues,]


write.table(Mean_of_y_matrix,paste(info_text,"_average_replicates.",treat_rep1_name, treat_rep2_name,"_log2_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)


  }


  
  if(selected_gene=="CRBN")
  {      
    
    dmso_rep1<-c()
    dmso_rep1_raw<-c()
    treat_rep1<-c()
    dmso_rep2<-c()
    dmso_rep2_raw<-c()
    treat_rep2<-c()
    dmso_rep3<-c()
    dmso_rep3_raw<-c()
    treat_rep3<-c()
    
    wt_counts_dmso_rep1<-c()
    
    dmso_rep1=as.data.frame(read.table(paste(dmso_rep1_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    dmso_rep1_raw=as.data.frame(read.table(paste("../",dmso_rep1_name,".gatk.aaCounts",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    treat_rep1=as.data.frame(read.table(paste(treat_rep1_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    
    
    dmso_rep2=as.data.frame(read.table(paste(dmso_rep2_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    dmso_rep2_raw=as.data.frame(read.table(paste("../",dmso_rep2_name,".gatk.aaCounts",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    
    treat_rep2=as.data.frame(read.table(paste(treat_rep2_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    dmso_rep3=as.data.frame(read.table(paste(dmso_rep3_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    dmso_rep3_raw=as.data.frame(read.table(paste("../",dmso_rep3_name,".gatk.aaCounts",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    
    treat_rep3=as.data.frame(read.table(paste(treat_rep3_name,".freq.txt",sep=""),sep="\t", header=TRUE, dec=".", fill = TRUE))
    
    
    
    
    wt_counts_dmso_rep1<-dmso_rep1$counts_wt_res
    
    
    all_values<-c()
    FC_values<-c()
    FDR_values<-c()
    
    for(i in c(1:nrow(dmso_rep1)))
      
    { 
      
      #print(i)
      
      for(j in c(1:(ncol(dmso_rep1[1:21]))))
      {
        #print(j)
        cutoff_rep1<-dmso_rep1[i,22]/total_reads_denominator
        cutoff_rep2<-dmso_rep2[i,22]/total_reads_denominator
        cutoff_rep3<-dmso_rep3[i,22]/total_reads_denominator
        
        wt_counts_dmso_rep1[i]
        
        
        
        
        if(dmso_rep1_raw[i,j] < cutoff_rep1 | dmso_rep2_raw[i,j] < cutoff_rep2 | dmso_rep3_raw[i,j] < cutoff_rep3)
        {fdr_value <- NA 
        log10_p <- NA
        log10_fdr <- NA
        mean_of_x<-NA
        mean_of_y<-NA
        log2_fold_change <- NA
        mean_x_and_y <- NA
        }
        
        if(dmso_rep1_raw[i,j] >= cutoff_rep1 & dmso_rep2_raw[i,j] >= cutoff_rep2 & dmso_rep3_raw[i,j] >= cutoff_rep3)
        {

          # t.test to be used
          ttest<-t.test(c(log2(dmso_rep1[i,j]+0.00001),log2(dmso_rep2[i,j]+0.00001),log2(dmso_rep3[i,j]+0.00001)),c(log2(treat_rep1[i,j]+0.00001),log2(treat_rep2[i,j]+0.00001),log2(treat_rep3[i,j]+0.00001)),alternative=test_alternative)
          
          
          
          log2_fold_change<-log2(mean(c((treat_rep1[i,j]+0.00001),(treat_rep2[i,j]+0.00001),(treat_rep3[i,j]+0.00001))) / mean(c((dmso_rep1[i,j]+0.00001),(dmso_rep2[i,j]+0.00001),(dmso_rep3[i,j]+0.00001))) )
          
          p_value<-ttest$p.value
          
          log10_p <- -log10(p_value)
          
          
          fdr_value<-p.adjust(p_value,method="fdr",n=fdr_n)
          log10_fdr <- -log10(fdr_value)
          
          mean_of_x<-ttest$estimate[[1]]
          mean_of_y<-ttest$estimate[[2]]
          mean_x_and_y <- ttest$estimate[[1]]+ttest$estimate[[2]]
        }
        
        
        FC_values<-rbind(FC_values,cbind(i,j,mean_of_y, mean_of_x,log2_fold_change))
        FDR_values<-rbind(FDR_values,cbind(i,j,mean_of_y, mean_of_x,fdr_value))
        all_values<-rbind(all_values,cbind(i,j,mean_of_y, mean_of_x,mean_x_and_y,log2_fold_change,log10_p,log10_fdr,fdr_value))
        
      }
    }
    
    
    
    
    
    
    
    
    

    head(all_values)
    dim(all_values)

        signif_cutoff_to_show_residues_in_heatmaps
    all_filtered_for_signif<-rbind( cbind(all_values[which(all_values[,9]<=signif_cutoff_to_show_residues_in_heatmaps),],all_values[which(all_values[,9]<=signif_cutoff_to_show_residues_in_heatmaps),6]),  cbind(all_values[which(all_values[,9]>signif_cutoff_to_show_residues_in_heatmaps | is.na(all_values[,9]) ),],NA))

    sort_all_filtered_for_signif<-all_filtered_for_signif[order(as.numeric(all_filtered_for_signif[,1]),as.numeric(all_filtered_for_signif[,2])) , ]
    head(sort_all_filtered_for_signif)
    

    filtered_FC_matrix<-matrix(sort_all_filtered_for_signif[,10], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)
    colnames(filtered_FC_matrix)<-colnames(dmso_rep1[1:21])
    rownames(filtered_FC_matrix)<-rownames(dmso_rep1)
    
    filtered_FC_matrix<-filtered_FC_matrix[selected_residues,]
    head(filtered_FC_matrix)
    
    
    
    write.table(filtered_FC_matrix,paste(info_text,"_average_replicates.",comparison_name,"_filtered_FC_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    write.table(all_values,paste(info_text,"_all_values.",comparison_name,"_",total_reads_denominator,".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    #filtered_FC_matrix
    
    
    
    
    
    
    FC_matrix<-matrix(FC_values[,5], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)
    
    colnames(FC_matrix)<-colnames(dmso_rep1[1:21])
    rownames(FC_matrix)<-rownames(dmso_rep1)
    
    FC_matrix<-FC_matrix[selected_residues,]
    
    #to be used:
    write.table(FC_matrix,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    FC_matrix_positive <-FC_matrix
    FC_matrix_positive[FC_matrix_positive<0] <- 0
    
    #Normalized Data pos_neg
    x_pos_neg=FC_matrix
    x_pos_only<-x_pos_neg
    x_pos_only[x_pos_only<0] <- 0
    
    x_neg_only<-x_pos_neg
    x_neg_only[x_neg_only>0] <- 0
    x_neg_only <- abs(x_neg_only)
    
    normalized_x_pos_only = (x_pos_only-min(x_pos_only,na.rm=T))/(max(x_pos_only,na.rm=T)-min(x_pos_only,na.rm=T))
    normalized_x_neg_only = (x_neg_only-min(x_neg_only,na.rm=T))/(max(x_neg_only,na.rm=T)-min(x_neg_only,na.rm=T))
    
    normalized_x_pos_neg<- normalized_x_pos_only-normalized_x_neg_only
    write.table(normalized_x_pos_neg,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.normalized_x_pos_neg.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    
    #Normalized Data pos
    x_pos=FC_matrix_positive
    normalized_x_pos = (x_pos-min(x_pos,na.rm=T))/(max(x_pos,na.rm=T)-min(x_pos,na.rm=T))
    
    write.table(normalized_x_pos,paste(info_text,"_average_replicates.",comparison_name,"_FC_matrix.",total_reads_denominator,".selectRes.normalized_x_pos.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    
    
    
    Mean_of_x_matrix <- matrix(FC_values[,4], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)
    
    colnames(Mean_of_x_matrix)<-colnames(dmso_rep1[1:21])
    rownames(Mean_of_x_matrix)<-rownames(dmso_rep1)
    
    Mean_of_x_matrix<-Mean_of_x_matrix[selected_residues,]
    
    write.table(Mean_of_x_matrix,paste(info_text,"_average_replicates.",dmso_rep1_name,dmso_rep2_name,dmso_rep3_name,"_log2_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    
    Mean_of_y_matrix <- matrix(FC_values[,3], nrow=nrow(dmso_rep1), ncol=ncol(dmso_rep1[1:21]),byrow = TRUE)
    
    colnames(Mean_of_y_matrix)<-colnames(dmso_rep1[1:21])
    rownames(Mean_of_y_matrix)<-rownames(dmso_rep1)
    
    Mean_of_y_matrix<-Mean_of_y_matrix[selected_residues,]
    
    
    write.table(Mean_of_y_matrix,paste(info_text,"_average_replicates.",treat_rep1_name, treat_rep2_name,treat_rep3_name,"_log2_matrix.",total_reads_denominator,".selectRes.txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)
    
    
  }
}



