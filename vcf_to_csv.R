setwd("E:/GWAS_catalog_Age_at_menarche_EFO_0004703/ukb-b-3768_Ben_Elsworth/")

library(data.table)
#install.packages("R.utils")
library(R.utils)

# Read VCF
read_data <- fread("ukb-b-3768.vcf.gz",skip = "#CHROM",header = TRUE,sep = "\t")
head(read_data)
dim(read_data)
# Rename first column
setnames(read_data, "#CHROM", "CHROM")

# Split FORMAT names
format_names <- unlist(strsplit(read_data$FORMAT[1], ":"))

# Split sample column values
split_values <- tstrsplit(read_data[[10]], ":")

# Convert to dataframe
split_df <- as.data.table(split_values)

# Assign column names
setnames(split_df, format_names)

# Combine with original VCF columns
final_data <- cbind(read_data, split_df)

# Optional rename
setnames(final_data,
         old = c("ES", "SE", "LP", "AF","ID"),
         new = c("BETA", "SE", "LOGP", "Allele_Frequency","rsid"))

# Save CSV
data_converted <- fwrite(final_data, "ukb-b-3768_full_by_R_code.csv")



######################################################################################################################################################################

data <- fread("ieu-b-5136_full_by_R_code.csv")
head(data)
dim(data)
