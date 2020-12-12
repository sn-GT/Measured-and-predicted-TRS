#######################################################################################################
# USER INPUTS REQUIRED
# ADD name of tissue, enter one of the following [Rectum/Colon/Cortex/Muscle/Blood]
reftranscriptome="Rectum"
# Set path to where you saved this folder on your system
home_path="/gpfs/scratch1/3/snagpal3/PPTRS_main"

#########################################################################################################

# pptrs_weights
pval <- read.table(paste0(home_path,"/Data/PPTRS_gene_weights/",reftranscriptome,"_gene_weights.txt"), header = T)

filenames <- as.character(pval$Filename)

pptrs <- read.table(paste0(home_path,"/Results/Predicted_GE/",reftranscriptome,"/", filenames[1], "_GE.txt"), header =T)
pptrs <- data.frame(pptrs[,"FID"])
names(pptrs)[1] <- "FID"

for(file in filenames){
	print(file)

	predictedge <- read.table(paste0(home_path,"/Results/Predicted_GE/",reftranscriptome,"/", file, "_GE.txt"), header =T)
	predictedge <- predictedge[,c("FID", "Predicted_Exp")]
	
	beta <- pval[pval$Filename == file, "Beta"]

	predictedge[,2] <- beta*predictedge[,2]
	names(predictedge)[2] <- file

	pptrs <- cbind(pptrs,predictedge[,2])
}


pptrs$WeightedPPTRS <- rowSums(pptrs[,-1])
names(pptrs)[1] <- "FID"
final_pptrs <- pptrs[,c("FID", "WeightedPPTRS")]

write.table(final_pptrs, paste0(home_path,"/Results/PPTRS/",reftranscriptome,"_PPTRS.txt"), row.names=F, col.names=T, quote = F)
