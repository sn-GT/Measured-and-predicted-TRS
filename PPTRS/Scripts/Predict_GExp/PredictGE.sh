#!/bin/bash

echo "entered"

while getopts "i:" arg; do
  case $arg in
    i)
      i=$OPTARG
      ;;
  esac
done

################### USER INPUTS REQUIRED ##########################
# ADD name of tissue, enter one of the following [Rectum/Colon/Cortex/Muscle/Blood]
reference_transcriptome="Rectum"

# ADD full path where you have saved the folder named PPTRS
home_path="/storage/scratch1/3/snagpal3/PPTRS"

# Path to PLINK2
plink2_path='/storage/coda1/p-ggibson3/0/snagpal3/rich_project_bio-gibson/Tools/plink2_Mar19'

# Genotype data in plink format (bed/bim/fam file path)
plink_input_bfile='/storage/coda1/p-ggibson3/0/snagpal3/rich_project_bio-gibson/Datasets/PROTECT/Genotypes/0_chr_imp_merged'

# Add path to IDs of your genotype data which you want to keep (subset European IDs here), the file should be in fam format
ids='/storage/coda1/p-ggibson3/0/snagpal3/rich_project_bio-gibson/Projects/IBD/GEXP_Imputation_PROTECT/keep_ids.txt'

# Path to Rscript on your system
R='/usr/local/pace-apps/spack/packages/0.12/linux-rhel7-x86_64/gcc-4.8.5/r-3.6.0-xrs45zl5vgbulifddk4swnfyesguxsqt/bin/Rscript'

#####################################################################

# Number of genes for each tissue
if [ $reference_transcriptome == "Rectum" ]
then
   n=820
elif [ $reference_transcriptome == "Colon" ]
then
   n=1097
elif [ $reference_transcriptome == "Cortex" ]
then
   n=1075
elif [ $reference_transcriptome == "Muscle" ]
then
   n=777   
else
   n=122
fi


# allele code 
if [ $reference_transcriptome == "Rectum" ]
then
   allele_code=${home_path}/Data/Other/allele_code/PROTECT_rectum_allele_code.txt
elif [ $reference_transcriptome == "Blood" ]
then
   allele_code=${home_path}/Data/Other/allele_code/CAGE_Blood_allele_code.txt   
else
   allele_code=${home_path}/Data/Other/allele_code/Gtex_allele_code.txt
fi

####################################################################
mkdir -p ${home_path}/Data/Other/raw_genotypes/${reference_transcriptome}

# Path for file list for Rectum
gene_annot=${home_path}/Data/Other/file_list/${reference_transcriptome}_files.txt

start=1
for (( i=$start; i<=$n; i++ )) do

file=$(awk -v i=$i 'NR==i{print $1}' "$gene_annot")

# snp list for each gene to be extracted
snp_path=${home_path}/Data/DPR_ciseQTL_weights/${reference_transcriptome}/snp_list/${file}_snp.txt

if [ $reference_transcriptome == "Rectum" ]
then
  echo "$(cut -d ':' -f 1 $snp_path)" > $snp_path
fi

# Get genotype matrix from plink files (recodeA)
${plink2_path} --bfile ${plink_input_bfile} --extract ${snp_path} --keep ${ids} --export A --export-allele ${allele_code}  --out ${home_path}/Data/Other/raw_genotypes/${reference_transcriptome}/${file}  

# Predict gene expression with Rscript matrix_multi.R R script
${R} --vanilla ${home_path}/Scripts/Predict_GExp/matrix_mult.R -f $file -p $home_path -r $reference_transcriptome &

done

cd ${home_path}/Data/Other/raw_genotypes/${reference_transcriptome}/
rm *