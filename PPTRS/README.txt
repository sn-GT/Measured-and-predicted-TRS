TWAS-based PPTRS performed by Sini Nagpal, Georgia Tech (sini.nagpal@gatech.edu)

Please cite paper if using the idea/code: https://pubmed.ncbi.nlm.nih.gov/34450030/ 
Mo, Nagpal, et al. Stratification of risk of progression to colectomy in ulcerative colitis via measured and predicted gene expression (2021)

Softwares needed:
- Plink2
- R 
- Linux to run shell script

R packages: 
- Install R packages prior to running the scripts, using install.packages("package_name") in R. 
  Only 2 R packages are required 1) data.table 2) opt.parse
		
Genotype data on which to compute PPTRS: 
- Convert your genotype data into standard plink format: bed/bim/fam 
- ALl analysis was done using hg19 positions and list of snps given in the above file. Convert your data (rsIDs and positions) to hg19 using liftOver prior to running the scripts.
- All the SNPs used in this analyses with their chr, position in hg19 are given in file: https://www.dropbox.com/s/zvqnbbtm0h2jw8s/All_SNPs_hg19_used_in_this_analysis.txt?dl=0
  You can check the overlap of snps with your data.
- Missing SNPs in your data: If ony few SNPs are missing, it should not matter, since each gene consist of 4-8k SNPs, so missing out a few of them with small effect should not matter. In case of many missing SNPs, you could pick SNP proxies and use the DPR cis-eQTL weight for the proxy SNP. However, that would require more coding work of screening snps within each gene and replacing the SNP with its proxy.

- Ancestry: Since cis-eqtl weights are computed from European subset, you can subset your data based on European ancestry.

#########################################################################################################################

Weights provided to compute PPTRS:

cis-eQTL SNP weights computed from DPR have been provided for:
- rectal gene expression (PROTECT)
- colon gene expression (GTEx)
- cortex gene expression (GTEx)
- muscle gene expression (GTEx)


TWAS weights of UC vs Control from UK Biobank provided for: 
- 820 genes of rectum
- 1097 genes of colon
- 1075 genes of cortex
- 777 genes of muscle

####################### CODING PIPELINE ##################################################################################

There are 2 steps of computing PPTRS:

Step1: Add paths and other user inputs, install R packages in PredcitGE.sh and run ./PredictGE.sh on your server. 
Step2: Add paths in the pptrs.r and run the R script pptrs.r on your server
Please note that for both of these steps, run the script as nohup or pbs job, to avoid exiting/logging out of shell/terminal while the script is running.
Each script will take some time to run depending upon the sample size of your data. 


Coding pipeline explained in detail: 
1) Predict Gene Expression from genotype data (Scripts/Predict_GExp/PredictGE.sh)
	- USER INPUTS  REQUIRED in PredictGE.sh:
	 	please add paths to softwares, PPTRS folder before running the script. 
	 	keep_ids contains the list of IDs (in fam format) for your data, you can use this to subset European IDs. 
	 	Run script for one tissue at a time by choosing reference_transcriptome as one of ["Rectum"/ "Colon"/ "Cortex"/ "Muscle"]. 
	 	Reference transcriptome/tissue determines that cis-eQTL weights derived from which tissue is to be used to predict gene expression.
	- Step1: Convert plink files into raw genotype matrix having allele count.
	- Step2: Once the genotype matrix is obtained using plink, R script matrix_multi.R is run. This script works as follows:
		For each gene, SNPs within +/- 1MB from gene start and end were taken. So each gene consists of about 4-8k SNPs. 
		For each gene, the cis-eqtl weights are taken from Data/DPR_ciseQTL_weights, corresponding genotype matrix generated in previous step is taken from Data/Other/raw_genotype. 
		Then matrix multiplication of cis-eqtl weights with genotype matrix is performed to get the predicted gene expression. 
		PredictedGeneExpression = cis-eqtl-weights x genotype-matrix    [...for each gene]
	-	To run the script, change user inputs and then type command: ./PredictGE.sh
	- Predicted gene expression is written into folder: Results/Predicted_GE/.
	- To ensure the script has successfully run, you can check that files with predicted gene expression for your data must be generated at Results/Predicted_GE/. There will be total 820 files for rectum, 1097 for colon, 1075 for cortex and 777 for muscle.  

2) Compute PPTRS (Scripts/PPTRS/pptrs.R):
  - USER INPUTS REQUIRED: set home path and reference trascriptome before running.
  - Now we have the predicted gene expression from first step. 
  - We have the PPTRS gene weights in Data/PPTRS_gene_weights. 
  - To compute PPTRS, we simply take the sum of the predicted gene expression, weighted by their TWAS weights; where twas weights are computed from logisitc regression of UCstatus ~ PredictedGE in UKbiobank, for genes with TWAS p < 0.05
  - PPTRS results are written in folder: Results/PPTRS/

###########################################################################################################################################################################

Contact: Please reach out to sini.nagpal@gatech.edu for any questions.

###########################################################################################################################################################################
