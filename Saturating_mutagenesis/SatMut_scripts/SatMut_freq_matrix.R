rm(list=ls())

################################################################################################################################################################################################
# Saturation mutations - generate a frequency mut. matrix from count matrix
#
# R --vanilla --args ${mapping_results_dir}/mut_calling/${sample_name}.gatk.aaCounts ${mapping_results_dir}/wt_residues_${gene}.t.txt < /data/groups/lab_winter/himrichova/resources/scripts/SatMut_freq_matrix.R
################################################################################################################################################################################################



input_aaCounts_file_name <- as.character(commandArgs()[4])
wt_residue_file_name <- as.character(commandArgs()[5])
output_name <- as.character(commandArgs()[6])

infile_aaCounts=as.data.frame(read.table(input_aaCounts_file_name,sep="\t", header=TRUE, dec=".", fill = TRUE, row.names=NULL))

wt_residues_file=as.data.frame(read.table(wt_residue_file_name,sep="\t", header=F, dec=".", fill = TRUE, row.names=NULL))



freq_aaCounts<-c()

COUNTS <- apply(as.matrix(infile_aaCounts),1,sum)

for(i in c(1:nrow(infile_aaCounts)))

{
print(i)

total_counts<-COUNTS[i]

wt_res<-as.character(wt_residues_file[which(wt_residues_file$V1==i),2])

counts_wt_res<-infile_aaCounts[i,which(colnames(infile_aaCounts)==wt_res)]

total_counts_minus_wt_res<-COUNTS[i]-counts_wt_res

# divie by total number of reads per 10k
total_counts_per_10k<-total_counts/10000

freq_row<-cbind(infile_aaCounts[i,]/total_counts,total_counts,counts_wt_res,total_counts_minus_wt_res,total_counts_per_10k,i)

freq_aaCounts<-rbind(freq_aaCounts,freq_row)

}

write.table(freq_aaCounts, file=output_name,quote=FALSE,sep="\t",row.names=T,col.names=T)



q()


