TRS from measured gene expression performed by Angela Mo, Georgia Tech.

Please cite paper if using the idea/code: https://pubmed.ncbi.nlm.nih.gov/34450030/ 
Mo, Nagpal, et al. Stratification of risk of progression to colectomy in ulcerative colitis via measured and predicted gene expression (2021).

Data: Bulk RNA-seq data for PROTECT https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150961 

PROTECT_deg.R Script for differential expression analysis, normalization, SVA, adjusting for covariates.

### TRS Computation
Transcriptional risk score (TRS) is the summation of polarized gene expression (observed) for genes found at the intersection of GWAS and eQTL. In the original paper (Marigorta, et al. PMID: 28805827), the list of genes associated with inflammatory bowel disease used for computing TRS were obtained with SMR and coloc analysis. You can find published studies to obtain the list of genes for your phenotype and use them to compute TRS. To compute TRS, you can use R to: 
	- First, obtain dataset with measured gene expression for your list of genes
    - Convert them to z-scores in R using dataset<-scale(gene_expr)
    - Compute the magnitude of differential expression,  
        case-controlstatus ~ zscore-gene-expression, in R using test_model<-glm(case_status ~ ., data=dataset)
    - For each gene, you will get a weight (beta) from glm in the above step using test_model$coefficients
    - Using these weights, finally, TRS is the weighted sum of gene expression, 
        TRS = Sum [beta x [z-score-gene-expression] ] for i=1 to n;
        dataset$predictor<- I(-0.04*dataset$APEH) + I(0.03*dataset$CDC42SE2) + ...)

