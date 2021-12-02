
# This script written by Angela Mo, Georgia Tech,  please cite while using it
# Paper: https://pubmed.ncbi.nlm.nih.gov/34450030/

#PROTECT- differential expression analysis & calculation of PC1
#Fig 1

#Load raw counts, phenotype
expr_orig<-read.table("0_raw_counts_448.txt",header=T,sep="\t",row.names=1)
phen_orig<-read.table("0_phen_351.txt",header=T,sep="\t")
genes_w_expr<-read.table("0_genes_w_expr.txt",header=F,sep="\t")

#restrict to genotyped individuals
phen<-phen_orig
expr<-expr_orig[rownames(expr_orig) %in% genes_w_expr$V1,colnames(expr_orig) %in% phen$X]
phen<-phen[order(phen$X),]
expr<-expr[,order(colnames(expr))]

library(edgeR)
#edgeR plain logcpm
y <- DGEList(counts=expr)
keep <- rowSums(cpm(y)>2) >= nrow(phen)/2
y <- calcNormFactors(y)
logcpm_plain <- cpm(y, prior.count=2, log=TRUE)
write.table(logcpm_plain,"logcpm_plain.txt",sep="\t",quote=F)

#edgeR normalized counts by bl/wk52
group <- factor(phen$Visit.Number)
y <- DGEList(counts=expr, group=group)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("BL","FU")
head(design)

keep <- rowSums(cpm(y)>2) >= nrow(design)/2
y <- y[keep,]

y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=c(1,-1))
detags <- rownames(topTags(lrt, p=0.05, adjust="BH"))
summary(decideTestsDGE(lrt, p=0.05, adjust="BH"))
toptags_table<-topTags(lrt,n=15000)
write.table(toptags_table,"de_bl_fu.txt",sep="\t",quote=F)
logcpm_time <- cpm(y, prior.count=2, log=TRUE)
write.table(logcpm_time,"logcpm_time.txt",sep="\t",quote=F)

#edgeR normalized counts by colectomy
group <- factor(phen$COLECTOMY_WK52_KM)
y <- DGEList(counts=expr, group=group)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("NoCol","Col")
head(design)

keep <- rowSums(cpm(y)>2) >= nrow(design)/2
y <- y[keep,]

y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=c(1,-1))
detags <- rownames(topTags(lrt, p=0.05, adjust="BH"))
summary(decideTestsDGE(lrt, p=0.05, adjust="BH"))
toptags_table<-topTags(lrt,n=15000)
write.table(toptags_table,"de_colectomy_allgenotyped.txt",sep="\t",quote=F)
logcpm_col <- cpm(y, prior.count=2, log=TRUE)
write.table(logcpm_col,"logcpm_col.txt",sep="\t",quote=F)


############################
##### 2F.B) FIT SVA            

library(sva)

### A) go for leek's advice: protect it in SVA, and then fit it as covariate
### https://support.bioconductor.org/p/54719/
mod = model.matrix(~as.factor(INITIAL_TRT_C4)+as.factor(Visit.Number)+as.factor(COLECTOMY_WK52_KM)+as.factor(FEMALE)+as.factor(RACE_NONWHITE), data=phen)
mod0 = model.matrix(~as.factor(FEMALE)+as.factor(RACE_NONWHITE),data=phen)
edata <- logcpm+1;
## This has 23 SV but protects disease!!
#n.sv = num.sv(edata,mod0);   ## because we dont want to protect subtype and gender
#n.sv = num.sv(edata,mod);   ## because we dont want to protect subtype and gender
#svobj = sva(as.matrix(edata),mod,mod0,n.sv=n.sv)      
##### I try this and in the end it fits 16 if I don't tell it how many to remove n.sv
n.sv = num.sv(edata,mod,method="leek")
svobj = sva(as.matrix(edata),mod,mod0)
sv.est <- svobj$sv;
name.print <- paste('02f_b_EstimatedSV_bigmodel2_update.txt',sep='')
write.table(round(sv.est,3),name.print,quote=F,append=F,col.names=F,row.names=F,sep='\t')

### ANOVA of SV with disease
name.print <- paste('02f_b_EstimatedSV_ANOVA_Dis_bigmodel2.pdf',sep='');
pdf(name.print);
par(mfrow=c(3,5),mar=c(5,1.5,5,1.5))
for(i in 1:ncol(sv.est)){
  pvalue <- summary(aov(sv.est[,i]~phen$INITIAL_TRT_C4))[[1]][["Pr(>F)"]][1];
  if(pvalue > 0.05){  
    boxplot(sv.est[,i]~phen$INITIAL_TRT_C4,cex.axis=1,las=2,main=paste('SV',i,'; P=',round(pvalue,2),sep=''));
  }  
  if(pvalue <= 0.05){  
    pvalue <- format(pvalue,scientific=T,digits=1); 
    #legend('top',paste('SV',i,'; P=',round(pvalue,2),sep=''),bty='n',col='red');
    boxplot(sv.est[,i]~phen$INITIAL_TRT_C4,cex.axis=1,las=2,main=paste('SV',i,'; P=',pvalue,sep=''),sep='');
  }  
}
dev.off();

### ANOVA of SV with timepoint
name.print <- paste('02f_b_EstimatedSV_ANOVA_Timepoint_bigmodel2.pdf',sep='');
pdf(name.print);
par(mfrow=c(3,5),mar=c(5,1.5,5,1.5))
for(i in 1:ncol(sv.est)){
  pvalue <- summary(aov(sv.est[,i]~phen$Visit.Number))[[1]][["Pr(>F)"]][1];
  if(pvalue > 0.05){  
    boxplot(sv.est[,i]~phen$Visit.Number,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',round(pvalue,2),sep=''));
  }  
  if(pvalue <= 0.05){  
    pvalue <- format(pvalue,scientific=T,digits=1); 
    #legend('top',paste('SV',i,'; P=',round(pvalue,2),sep=''),bty='n',col='red');
    boxplot(sv.est[,i]~phen$Visit.Number,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',pvalue,sep=''),sep='');
  }  
}
dev.off();

### ANOVA of SV with colectomy
name.print <- paste('02f_b_EstimatedSV_ANOVA_Colectomy_bigmodel2.pdf',sep='');
pdf(name.print);
par(mfrow=c(3,5),mar=c(5,1.5,5,1.5))
for(i in 1:ncol(sv.est)){
  pvalue <- summary(aov(sv.est[,i]~phen$COLECTOMY_WK52_KM))[[1]][["Pr(>F)"]][1];
  if(pvalue > 0.05){  
    boxplot(sv.est[,i]~phen$COLECTOMY_WK52_KM,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',round(pvalue,2),sep=''));
  }  
  if(pvalue <= 0.05){  
    pvalue <- format(pvalue,scientific=T,digits=1); 
    #legend('top',paste('SV',i,'; P=',round(pvalue,2),sep=''),bty='n',col='red');
    boxplot(sv.est[,i]~phen$COLECTOMY_WK52_KM,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',pvalue,sep=''),sep='');
  }  
}
dev.off();

### ANOVA of SV with gender
name.print <- paste('02f_b_EstimatedSV_ANOVA_Gender_bigmodel2.pdf',sep='');
pdf(name.print);
par(mfrow=c(3,5),mar=c(5,1.5,5,1.5))
for(i in 1:ncol(sv.est)){
  pvalue <- summary(aov(sv.est[,i]~phen$FEMALE))[[1]][["Pr(>F)"]][1];
  if(pvalue > 0.05){  
    boxplot(sv.est[,i]~phen$FEMALE,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',round(pvalue,2),sep=''));
  }  
  if(pvalue <= 0.05){  
    pvalue <- format(pvalue,scientific=T,digits=1); 
    #legend('top',paste('SV',i,'; P=',round(pvalue,2),sep=''),bty='n',col='red');
    boxplot(sv.est[,i]~phen$FEMALE,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',pvalue,sep=''),sep='');
  }  
}
dev.off();

### ANOVA of SV with race
name.print <- paste('02f_b_EstimatedSV_ANOVA_Race_bigmodel2.pdf',sep='');
pdf(name.print);
par(mfrow=c(3,5),mar=c(5,1.5,5,1.5))
for(i in 1:ncol(sv.est)){
  pvalue <- summary(aov(sv.est[,i]~phen$RACE_NONWHITE))[[1]][["Pr(>F)"]][1];
  if(pvalue > 0.05){  
    boxplot(sv.est[,i]~phen$RACE_NONWHITE,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',round(pvalue,2),sep=''));
  }  
  if(pvalue <= 0.05){  
    pvalue <- format(pvalue,scientific=T,digits=1); 
    #legend('top',paste('SV',i,'; P=',round(pvalue,2),sep=''),bty='n',col='red');
    boxplot(sv.est[,i]~phen$RACE_NONWHITE,cex.axis=0.8,las=2,main=paste('SV',i,'; P=',pvalue,sep=''),sep='');
  }  
}
dev.off();


##SNM
library(snm);

## PREPARING THE INPUT INTEGER FILE (1)
int.var <- as.data.frame(1:dim(expr)[2],ncol=1); 
colnames(int.var) <- 'Array';    
int.var$Array <- as.factor(int.var$Array);

## PREPARING THE ADJUSTMENT FILE, by PROCESSING DATE (2)
sv <- read.table('02f_b_EstimatedSV_bigmodel2.txt',header=F);
sv_keep<-sv[,-c(5,6,8,9,19)]
snm.adj <- as.data.frame(sv_keep);          
colnames(snm.adj) <- paste('SV_keep',c(1:ncol(sv_keep)),sep=''); 
row.names(snm.adj) <- colnames(expr);
snm.adj <- model.matrix(~.,snm.adj);

## PREPARING THE BIOLOGICAL FILE: Disease Yes/No + Types of Disease (3)
snm.bio <- data.frame(init_trt=phen$INITIAL_TRT_C4,time=phen$Visit.Number,colectomy=phen$COLECTOMY_WK52_KM);  
colnames(snm.bio) <- c('Init_trt','Timepoint','Colectomy');
snm.bio <- model.matrix(~., snm.bio);

## PREPARING THE EXPRESSION FILE (4)
snm.raw <- as.matrix(edata);    

## RUN SNM
snmR.exposerum1 <- snm(snm.raw,bio.var=snm.bio,adj.var=snm.adj,int.var=int.var,rm.adj=T,num.iter=1);
snmR.exposerum2 <- snm(snm.raw,bio.var=snm.bio,adj.var=snm.adj,int.var=int.var,rm.adj=T,num.iter=1)
snmR.exposerum3 <- snm(snm.raw,bio.var=snm.bio,adj.var=snm.adj,int.var=int.var,rm.adj=T,num.iter=1)
expr_snm1 <- snmR.exposerum1$norm.dat;
expr_snm2 <- snmR.exposerum2$norm.dat;
expr_snm3 <- snmR.exposerum3$norm.dat;
expr_snm1 <- round(expr_snm1,3);
expr_snm2 <- round(expr_snm2,3);
expr_snm3 <- round(expr_snm3,3);
colnames(expr_snm1) <- phen$X;
colnames(expr_snm2) <- phen$X;
colnames(expr_snm3) <- phen$X;
name.print1 <- paste('0_sva_snm_ver_nonrem.txt',sep='');
name.print2 <- paste('0_sva_snm_ver_rem.txt',sep='');
name.print3 <- paste('0_sva_snm_ver_3.txt',sep='');
write.table(expr_snm1,name.print1,quote=F,sep="\t",row.names=T,col.names=T)
write.table(expr_snm2,name.print2,quote=F,sep="\t",row.names=T,col.names=T)
write.table(expr_snm3,name.print3,quote=F,sep="\t",row.names=T,col.names=T)

library(NMF);
logcpm_gender<-logcpm[rownames(logcpm) %in% c('RPS4Y1','EIF1AY','DDX3Y','KDM5D','XIST')]
annCol <- list(phen$FEMALE);
names(annCol) <- 'Gender';
annColors <- list(c('orange','lightblue','grey'));
pdf('02c_Heatmap_CheckGender_Log2AndRaw.pdf');
raw_heatmap<-aheatmap(as.data.frame(logcpm_gender),annCol=annCol,annColors=annColors);
aheatmap(as.data.frame(snm_expr_gender),annCol=annCol,annColors=annColors);
dev.off();

#Remove 50081, 10355 gender /10587, 50471 PCA outliers
snm_expr_final<-expr_snm3[,-which(colnames(expr_snm3) %in% c('BPY50081','BPY10355','BPY10587','BPY50471'))]
write.table(snm_expr_final,"0_sva_snm_final.txt",quote=F,sep="\t",row.names=T,col.names=T)
phen_update4<-phen[-which(phen$X %in% c('BPY50081','BPY10355','BPY10587','BPY50471')),]

### Inverse quantile normalization per gene
inv.norm <- function(expresion){
  z.val <- rep(0,length(expresion));
  z.val[order(expresion)] <- qnorm(seq(from=1,to=length(expresion))/(1+length(expresion)));
  return(z.val);
}
expr.quant.z <- t(apply(expr_snm1,1,inv.norm));  
colnames(expr.quant.z) <- colnames(expr_snm1);

name.print <- paste('02f_b_ConcatenatedExpr_RPKM_15K_FitSV_Subtype_QuantAndInvNormAgain.txt',sep='');
write.table(round(expr.quant.z,3),name.print,quote=F,append=F,row.names=F,col.names=T,sep='\t');

phen_baseline<-phen[phen$Visit.Number == "CO:01 Screening/Baseline",]
phen_wk52<-phen[phen$Visit.Number == "CO:08 52wk FU",]
expr_baseline<-expr.quant.z[,colnames(expr.quant.z) %in% phen_baseline$X]
expr_wk52<-expr.quant.z[,colnames(expr.quant.z) %in% phen_wk52$X]


expr_baseline<-as.data.frame(t(expr_baseline))
expr_baseline<-expr_baseline[order(rownames(expr_baseline)),]
phen_baseline<-phen_baseline[order(phen_baseline$X),]
expr_baseline<-as.data.frame(t(expr_baseline))
write.table(round(expr_baseline,3),"0_sva_snm_baseline_expr.txt",quote=F,append=F,row.names=T,col.names=T,sep='\t')


expr_wk52<-as.data.frame(t(expr_wk52))
expr_wk52<-expr_wk52[order(rownames(expr_wk52)),]
phen_wk52<-phen_wk52[order(phen_wk52$X),]
expr_wk52<-as.data.frame(t(expr_wk52))
write.table(round(expr_wk52,3),"0_sva_snm_wk52_expr.txt",quote=F,append=F,row.names=T,col.names=T,sep='\t')

write.table(bl_wk52like_expr,"0_sva_snm_baseline_173_expr.txt",sep="\t",quote=F)
write.table(wk52_wk52like_expr,"0_sva_snm_wk52_77_expr.txt",sep="\t",quote=F)



#voom DEG with adjusted logcpm
library(edgeR)

phen<-phen_update4
expr_snm1<-read.table("0_sva_snm_final.txt",sep="\t",row.names=1,header=T)

colectomy<-factor(phen$COLECTOMY_WK52_KM)

group<-interaction(race,dx)
group<-factor(colectomy)

d0 <- DGEList(counts=expr_snm1+3)
d <- d0 
dim(d)

#plotMDS(d,col=as.numeric(group),pch=1)
#black, red, green, blue, teal, pink 
#AA.CD EU.CD AA.CTRL EU.CTRL AA.UC EU.UC
#AA.FALSE EU.FALSE AA.TRUE EU.TRUE

mm<-model.matrix(~0+group)
y <- voom(d,mm)
fit <- lmFit(y)
head(coef(fit))

contr <- makeContrasts(group0 - group1, levels = colnames(coef(fit))) #[1] 6920

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "DEG_NoCol_vs_Col_SVA_fin.txt", row.names = F, sep = "\t", quote = F)

################

setwd("C:/Users/amo3/Desktop/PROTECT/Paper/Figures")
load("C:/Users/amo3/Desktop/PROTECT/eQTL_redo/env1.RData")
load("env1.RData")

library(showtext)
font_add("Arial","arial.ttf")

top.table$direction<-c(rep(NA,12577))
for (i in 1:12577){
  top.table[i,8]<-ifelse(top.table[i,2] > 0 & top.table[i,6] < 4*10^-6,"b_positive",ifelse(top.table[i,2] < 0 & top.table[i,6] < 4*10^-6,"c_negative","a_nonsignificant"))
}
top.table$direction<-factor(top.table$direction)
top.table<-top.table[order(top.table$direction),]

pdf(file="Updated_voom_volcano.pdf",height=8,width=8)
par(family="Arial",mar=c(8,8,4,4))

p_vol<-ggplot(data=top.table, aes(x=logFC,y=-log10(adj.P.Val),color=direction)) +
  geom_vline(xintercept=0,color="grey20",linetype="dotted") +
  geom_hline(yintercept=-log10(4*10^-6),color="grey20",linetype="dotted") +
  geom_point(shape=1,size=2.5) +
  scale_color_manual(values=c("grey80","orangered1","royalblue3")) +
  theme_bw() +
  theme(plot.title=element_text(size=20),text=element_text(color="black",size=20),legend.position = "none",plot.margin = unit(c(1,1,1,1),"lines"),axis.text.x=element_text(size=20,color="black"),axis.text.y=element_text(size=20,color="black"),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  labs(x="Log Difference",y="-Log10 Adj. P Value")
p_vol

dev.off();

pdf(file="Updated_voom_volcano_2021_2.pdf",height=8,width=8)
par(family="Arial",mar=c(8,8,4,4))

p_vol<-ggplot(data=top.table, aes(x=logFC,y=-log10(P.Value),color=direction)) +
  geom_vline(xintercept=0,color="grey20",linetype="dotted") +
  geom_hline(yintercept=-log10(0.05),color="grey20",linetype="dotted") +
  geom_point(shape=1,size=2.5) +
  scale_color_manual(values=c("grey80","orangered1","royalblue3")) +
  theme_bw() +
  theme(plot.title=element_text(size=20),text=element_text(color="black",size=20),legend.position = "none",plot.margin = unit(c(1,1,1,1),"lines"),axis.text.x=element_text(size=20,color="black"),axis.text.y=element_text(size=20,color="black"),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  labs(x="Log Difference",y="-Log10 P Value")
p_vol

dev.off();

with(top.table, plot(logFC, -log10(adj.P.Val), pch=1, col="grey85", main="Colectomy Gene Expression", ylim=c(0,115),xlim=c(-1,1),cex=1.2,lwd=1.5,xlab="Log Difference", ylab = "-Log10 Adj. P Value"))
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)

points(top.table[top.table$adj.P.Val < 0.05 & top.table$logFC > 0,2],-log10(top.table[top.table$adj.P.Val < 0.05 & top.table$logFC > 0,6]),pch=1,col="orangered1",cex=1.2,lwd=1.5)
points(top.table[top.table$adj.P.Val < 0.05 & top.table$logFC < 0,2],-log10(top.table[top.table$adj.P.Val < 0.05 & top.table$logFC < 0,6]),pch=1,col="royalblue3",cex=1.2,lwd=1.5)

abline(h=-log10(0.05), lwd = 2, lty=2)
abline(v=0, lwd = 2, lty=2)

text(-0.45,115,labels="Upregulated in Col")
text(0.45,115,labels="Downregulated in Col")

dev.off()

#Extract DEG & calculate PC1
load("C:/Users/amo3/Desktop/PROTECT/eQTL_redo/env1.RData")

top_150_genes<-top.table[1:150,]
top_500_genes<-top.table[1:500,]

top_150_genes_expr<-expr_snm1[rownames(expr_snm1) %in% top_150_genes$Gene,]
top_500_genes_expr<-expr_snm1[rownames(expr_snm1) %in% top_150_genes$Gene,]
smillie_epi_genes_expr<-expr_snm1[rownames(expr_snm1) %in% smillie_epi_unique$Gene,]
smillie_imm_genes_expr<-expr_snm1[rownames(expr_snm1) %in% smillie_imm_unique$Gene,]

smillie_epi_genes_expr200<-expr_snm1[rownames(expr_snm1) %in% smillie_epi_unique_200$Gene,]
smillie_imm_genes_expr200<-expr_snm1[rownames(expr_snm1) %in% smillie_imm_unique_200$Gene,]

print_150<-t(top_150_genes_expr)
print_500<-t(top_500_genes_expr)
print_epi<-t(smillie_epi_genes_expr)
print_imm<-t(smillie_imm_genes_expr)

print_epi200<-t(smillie_epi_genes_expr200)
print_imm200<-t(smillie_imm_genes_expr200)

write.table(print_150,"top_150_genes_expr_for_pc.txt",sep="\t",quote=F,row.names=T)
write.table(print_500,"top_500_genes_expr_for_pc.txt",sep="\t",quote=F,row.names=T)
write.table(phen,"top_for_pc_phen.txt",sep="\t",quote=F,row.names=T)

write.table(top.table,"gsea_genes_list.txt",sep="\t",quote=F)

library(ggplot2)
library(ggpubr)
library(showtext)
font_add("Arial","arial.ttf")

pca_150<-prcomp(print_150, rank. = 1,scale. = T)
pc1_scores_150<-as.vector(pca_150$x,)
phen$PC1_colectomy_score<-pc1_scores_150

pca_epi<-prcomp(print_epi, rank. = 1,scale. = T)
pc1_scores_epi<-as.vector(pca_epi$x,)
phen$PC1_epi<-pc1_scores_epi

pca_imm<-prcomp(print_imm, rank. = 1,scale. = T)
pc1_scores_imm<-as.vector(pca_imm$x,)
phen$PC1_imm<-pc1_scores_imm

pca_epi200<-prcomp(print_epi200, rank. = 1,scale. = T)
pc1_scores_epi200<-as.vector(pca_epi200$x,)
phen$PC1_epi200<-pc1_scores_epi200

pca_imm200<-prcomp(print_imm200, rank. = 1,scale. = T)
pc1_scores_imm200<-as.vector(pca_imm200$x,)
phen$PC1_imm200<-pc1_scores_imm200

phen$colectomy_status<-as.factor(phen$COLECTOMY_WK52_KM)
phen$time<-as.factor(phen$Visit.Number)
phen$init_trt<-as.factor(phen$INITIAL_TRT_C4)
phen$pc1<-as.numeric(pc1_scores_150$PC1)
phen$COLECTOMY_YR3<-as.factor(conv_col_3yr)
phen$COLECTOMY_stepwise<-as.factor(as.numeric(phen$COLECTOMY_YR3)+as.numeric(phen$COLECTOMY_WK52_KM))

phen_bl<-phen[phen$Visit.Number == "CO:01 Screening/Baseline",]

pdf(file="Boxplot_time.pdf",width=16,height=4)
par(family="Arial",mar=c(8,8,4,4))

p_col<-ggplot(data=phen_bl, aes(x=as.factor(COLECTOMY_WK52_KM),y=as.numeric(PC1),color=as.factor(COLECTOMY_WK52_KM))) +
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#F5793A","#A95AA1")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Colectomy Status by Week 52",x=element_blank(),y="GE-Based Risk Score") + 
  scale_x_discrete(labels=c("0" = "No Colectomy", "1" = "Colectomy"))
p_col

#phen_bl$init_trt<-factor(phen_bl$init_trt,levels(phen_bl$init_trt)[c(1,3,2)])
p_init<-ggplot(data=phen_bl, aes(x=init_trt,y=pc1,color=init_trt)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#63ACBE","#601A4A","#EE442F")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Initial Treatment",x=element_blank(),y="GE-Based Risk Score") +
  scale_x_discrete(labels=c("5ASA" = "5-ASA", "CS-Oral" = "Oral CS", "CS-IV" = "IV CS"))
#p_init

p_time<-ggplot(data=phen, aes(x=time,y=pc1,color=Both_timepoints)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#85C0F9","#85C0F9","#0F2080")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Timepoint",x=element_blank(),y="GE-Based Risk Score") + 
  scale_x_discrete(labels=c("CO:01 Screening/Baseline" = "Baseline", "CO:08 52wk FU" = "Week 52"))
p_time

ggarrange(p_time,p_time,p_time,ncol=3,nrow=1,align="h")

dev.off();


p_col3<-ggplot(data=phen_bl, aes(x=COLECTOMY_YR3,y=pc1,color=COLECTOMY_YR3)) +
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#F5793A","#A95AA1")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Colectomy Status by Week 52",x=element_blank(),y="GE-Based Risk Score") + 
  scale_x_discrete(labels=c("0" = "No Colectomy", "1" = "Colectomy"))
p_col3

pdf(file="Supp_colyr3.pdf",width=6,height=4)
par(family="Arial",mar=c(8,8,4,4))

p_col4<-ggplot(data=phen_bl, aes(x=COLECTOMY_stepwise,y=pc1,color=COLECTOMY_stepwise)) +
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#63ACBE","#601A4A","#EE442F")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Colectomy Status",x=element_blank(),y="GE-Based Risk Score") + 
  scale_x_discrete(labels=c("1" = "No Colectomy", "2" = "Colectomy\nYear 3", "3" = "Colectomy\nYear 1"))
print(p_col4)

dev.off();

pdf(file="Supp_PC1_epi_imm.pdf",width=10,height=4)

p_time_epi<-ggplot(data=phen, aes(x=time,y=PC1_epi,color=time)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#0F2080","#85C0F9")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Timepoint",x=element_blank(),y="PC1 Epithelial GE") + 
  scale_x_discrete(labels=c("CO:01 Screening/Baseline" = "Baseline", "CO:08 52wk FU" = "Week 52"))
p_time_epi

p_time_imm<-ggplot(data=phen, aes(x=time,y=PC1_imm,color=time)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#0F2080","#85C0F9")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Timepoint",x=element_blank(),y="PC1 Immune GE") + 
  scale_x_discrete(labels=c("CO:01 Screening/Baseline" = "Baseline", "CO:08 52wk FU" = "Week 52"))
p_time_imm

p_time_epi<-ggplot(data=phen, aes(x=time,y=-PC1_epi200,color=time)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#0F2080","#85C0F9")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Timepoint",x=element_blank(),y="PC1 Epithelial GE") + 
  scale_x_discrete(labels=c("CO:01 Screening/Baseline" = "Baseline", "CO:08 52wk FU" = "Week 52"))
p_time_epi

p_time_imm<-ggplot(data=phen, aes(x=time,y=-PC1_imm200,color=time)) + 
  geom_jitter(alpha=1,width=0.15) +
  scale_color_manual(values=c("#0F2080","#85C0F9")) +
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6,linetype="dotted", color="grey20",alpha=0) + 
  geom_boxplot(width=0.4,outlier.shape=NA,lwd=0.6, color="grey20",coef=0,alpha=0) +
  theme(plot.title=element_text(size=16),text=element_text(color="black",size=16),legend.position = "none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.line = element_line(color="black"),panel.background=element_blank(),plot.margin = unit(c(0.5,0.5,0.5,0.5),"lines"),axis.text.x=element_text(size=16,color="black"),axis.text.y=element_text(size=16,color="black")) + 
  labs(title="Timepoint",x=element_blank(),y="PC1 Immune GE") + 
  scale_x_discrete(labels=c("CO:01 Screening/Baseline" = "Baseline", "CO:08 52wk FU" = "Week 52"))
p_time_imm

ggarrange(p_time_epi,p_time_imm,ncol=2,nrow=1,align="h")

dev.off()