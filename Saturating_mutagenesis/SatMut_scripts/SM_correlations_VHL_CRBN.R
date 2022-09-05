rm(list=ls())

list.files()
require("pheatmap")
library(gridExtra)

working_directory <- as.character(commandArgs()[4])
directory_with_SatMut_input_files <- as.character(commandArgs()[5])

setwd(paste(working_directory,sep=""))








FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix <- read.table("FC_matrix__average_replicates.VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix_positive <-FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix
FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix_positive[FC_VHL_DMSO_2ndBatch_vs_DMSO_1stBatch_FC_matrix_positive<0] <- 0

FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive <-FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix
FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive[FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive<0] <- 0

FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV50_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix_positive <-FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix
FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix_positive[FC_3_VHL_SM_r1r2_ARV50_vs_DMSO_FDR_matrix_positive<0] <- 0

FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV5_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix_positive <-FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix
FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix_positive[FC_5_VHL_SM_r1r2_ARV5_vs_DMSO_FDR_matrix_positive<0] <- 0



FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_ARV500_vs_DMSO_2ndBatch_FC_matrix.10000.selectRes.txt")
FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive <-FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix
FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive[FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive<0] <- 0




FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_MZ500_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive <- FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix
FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive[FC_2_VHL_SM_r1r2_MZ500_vs_DMSO_FDR_matrix_positive<0] <- 0



FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive <- FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix
FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive[FC_3a_VHL_SM_r1r2_ACBI1_2uM_7d_vs_DMSO_FDR_matrix_positive<0] <- 0


FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive <- FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix
FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive[FC_4_VHL_SM_r1r2_DAT548_500nM_7d_vs_DMSO_FDR_matrix_positive<0] <- 0


FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix <- read.table("FC_matrix__average_replicates.VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive <- FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix
FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive[FC_5_VHL_SM_r1r2_DAT551_2uM_7d_vs_DMSO_FDR_matrix_positive<0] <- 0



FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_CC885_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive<-FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix
FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive[FC_2_CRBN_SM_CC885_vs_DMSO_FC_matrix_positive<0] <- 0

FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_CC90009_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive<-FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix
FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive[FC_3_CRBN_SM_CC90009_vs_DMSO_FC_matrix_positive<0] <- 0

FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET57_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive<-FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix
FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive[FC_4_CRBN_SM_dBET57_vs_DMSO_FC_matrix_positive<0] <- 0

FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_vs_DMSO_FC_matrix.10000.selectRes.txt")
FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive<-FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix
FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive[FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive<0] <- 0


FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix <- read.table("FC_matrix__average_replicates.CRBN_SM_dBET6_1stBatch_vs_DMSO_1stBatch_FC_matrix.10000.selectRes.txt")
FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive<-FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix
FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive[FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive<0] <- 0








####################################################################################
# Calculate correlations between selected conditions


###############################################
# ARV500/DMSO_1st vs ARV500/DMSO_2nd
###############################################


sample1<-as.vector(as.matrix(FC_1_VHL_SM_r1r2_ARV500_vs_DMSO_FDR_matrix_positive[,1:21]))

sample2<-as.vector(as.matrix(FC_1b_VHL_SM_r1r2_ARV500_vs_DMSO2nd_FDR_matrix_positive[,1:21]))



my_data<-as.data.frame(cbind(sample1,sample2))

cor.test(sample1,sample2,use="complete.obs",method="spearman",alternative = "two.sided", exact=FALSE)



library("ggpubr")

pdf("correlation_ARV500_vs_DMSO1st_or_DMSO2nd.pdf")
ggscatter(my_data, x = "sample1", y = "sample2", 
          add = "reg.line", conf.int = T, 
          color="black",size=1,fill = "lightgray",
          cor.coef = T, cor.method = "spearman",
          xlab = "Log2FC (ARV500 / DMSO_1st)"
          , ylab = "Log2FC (ARV500 / DMSO_2nd)"
          ,cor.coef.size = 7
          ,cor.coef.coord = c(0, 2.5)
          )


dev.off()



shuffled_data3<-as.data.frame(cbind(sample(sample1),sample(sample2) ))

pdf("correlation_ARV500_vs_DMSO1st_or_DMSO2nd.shuffled.pdf")

ggscatter(shuffled_data3, x = "V1", y = "V2", 
          add = "reg.line", conf.int = T, 
          color="black",size=1,fill = "lightgray",
          cor.coef = T, cor.method = "spearman",
          xlab = "Log2FC (ARV500 / DMSO_1st) (shuffled)"
          , ylab = "Log2FC (ARV500 / DMSO_2nd) (shuffled)"
          ,cor.coef.size = 7
          ,cor.coef.coord = c(0, 2.5)
          )

dev.off()





###############################################################
# dBET6_1st_vs_DMSO_1st vs dBET6_2nd_vs_DMSO_2nd
###############################################################

sample1<-as.vector(as.matrix(FC_6_CRBN_SM_dBET6_vs_DMSO_1stBatch_FC_matrix_positive[,1:21]))

sample2<-as.vector(as.matrix(FC_5_CRBN_SM_dBET6_vs_DMSO_FC_matrix_positive[,1:21]))



my_data<-as.data.frame(cbind(sample1,sample2))

cor.test(sample1,sample2,use="complete.obs",method="spearman",alternative = "two.sided")



library("ggpubr")

pdf("correlation_dBET6vsDMSO1stBatch_vs_dBET6vsDMSO2ndBatch.pdf")
ggscatter(my_data, x = "sample1", y = "sample2", 
          add = "reg.line", conf.int = T, 
          color="black",size=1,fill = "lightgray",
          cor.coef = T, cor.method = "spearman",
          xlab = "Log2FC (dBET6_1st / DMSO_1st)"
          , ylab = "Log2FC (dBET6_3rd / DMSO_3rd)"
          ,cor.coef.size = 7
          ,cor.coef.coord = c(0, 2.5)
)


dev.off()



shuffled_data3<-as.data.frame(cbind(sample(sample1),sample(sample2) ))

pdf("correlation_dBET6vsDMSO1stBatch_vs_dBET6vsDMSO2ndBatch.shuffled.pdf")

ggscatter(shuffled_data3, x = "V1", y = "V2", 
          add = "reg.line", conf.int = T, 
          color="black",size=1,fill = "lightgray",
          cor.coef = T, cor.method = "spearman",
          xlab = "Log2FC (dBET6_1st / DMSO_1st) (shuffled)"
          , ylab = "Log2FC (dBET6_3rd / DMSO_3rd) (shuffled)"
          ,cor.coef.size = 7
          ,cor.coef.coord = c(0, 2.5)
)

dev.off()
