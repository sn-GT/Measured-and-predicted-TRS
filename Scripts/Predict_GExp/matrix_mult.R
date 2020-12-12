library("optparse")
library(data.table)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--home_path"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--reference_transcriptome"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

print("entered R")
print(opt$reference_transcriptome)
print(opt$file)

effectsize <- read.table(paste0(opt$home_path,"/Data/DPR_ciseQTL_weights/",opt$reference_transcriptome, "/weights/", opt$file,".param.txt"), header = T)
effectsize$rs <- gsub(':.*(.*)','\\1', effectsize$rs)
effectsize$BETA <- effectsize$b + effectsize$beta
effectsize_sub <- effectsize[,c("rs", "BETA")]

recodeA <- fread(paste0(opt$home_path, "/Data/Other/raw_genotypes/", opt$reference_transcriptome, "/",opt$file, ".raw"))

names(recodeA) <- gsub('_.*(.*)','\\1', names(recodeA))

exp <- recodeA[,c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]

rsids <- colnames(recodeA[,-(1:6)])
geno <- data.frame(t(recodeA[,-(1:6)]))
geno$rs <- rsids

beta_dos <- merge(effectsize_sub, geno, by = "rs", sort = F)

k <- which(is.na(beta_dos), arr.ind=TRUE)
beta_dos[k] <- rowMeans(beta_dos[,-(1:2)], na.rm=TRUE)[k[,1]]

beta_matrix <- as.matrix(beta_dos$BETA)
geno_matrix <- as.matrix(t(beta_dos[,-(1:2)]))

# Matrix multiplication
gene_expression <- data.frame(geno_matrix %*% beta_matrix)
names(gene_expression)[1] <- "Predicted_Exp"
exp_df <- cbind(exp, gene_expression)
exp_df <- exp_df[,-(3:6)]

write.table(exp_df, paste0(opt$home_path, "/Results/Predicted_GE/",opt$reference_transcriptome,"/",opt$file,"_GE.txt"), row.names = F, col.names = T, quote = F)

print("done")