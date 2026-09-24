#-----------------------------------------------------------
# Set working directory
# This is the folder where input files are stored
# and output files will be saved
#-----------------------------------------------------------
setwd("E:/Final_Data_GBC_Paper/MR_11_Sept_2026/Phewas_Results")


#-----------------------------------------------------------
# Load required libraries
#-----------------------------------------------------------

# data.table -> used for fast reading of large files
library(data.table)

# ieugwasr -> contains functions for PheWAS and OpenGWAS access
library(ieugwasr)

# dplyr -> used for data manipulation and filtering
library(dplyr)


#-----------------------------------------------------------
# Read the GWAS input file
#-----------------------------------------------------------
data <- fread("GBC_GWAS_summary_data_8_lead_snps.csv")
dim(data)
head(data)
cat(colnames(data),sep = "\n")


#-----------------------------------------------------------
# Check smallest and largest p-values
# after filtering
#-----------------------------------------------------------
min(data$`P - value`)
max(data$`P - value`)


#-----------------------------------------------------------
# Extract only SNP column from filtered data
#-----------------------------------------------------------
snps <- data$rsID


#-----------------------------------------------------------
# Display first few SNPs
#-----------------------------------------------------------
head(snps)


#-----------------------------------------------------------
# Count total number of SNPs
#-----------------------------------------------------------
length(snps)


#-----------------------------------------------------------
# Set OpenGWAS authentication token
# This token allows access to OpenGWAS database
#-----------------------------------------------------------
Sys.setenv(OPENGWAS_JWT = "Open_GWAS_Token")


#-----------------------------------------------------------
# Run PheWAS analysis
#
# variants = SNP list
# pval = significance threshold
# batch = database batch selection
#
# batch = c() means search across all available datasets
#-----------------------------------------------------------
result2 <- phewas(variants = snps, pval = 1e-05, batch = c("ukb-a", "ukb-b", "ukb-d", "ukb-e"))


#-----------------------------------------------------------
# View PheWAS result
#-----------------------------------------------------------
result2


#-----------------------------------------------------------
# Check dimensions of result dataset
#-----------------------------------------------------------
dim(result2)


#-----------------------------------------------------------
# Save final PheWAS result into CSV file
#-----------------------------------------------------------

write.csv(result2,"Phewas_results_of_8_lead_snps.csv",row.names = FALSE)
