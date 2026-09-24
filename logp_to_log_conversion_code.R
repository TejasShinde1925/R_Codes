
setwd("E:/GWAS_catalog_Age_at_menarche_EFO_0004703/ukb-b-3768_Ben_Elsworth/")

library(data.table)

#-----------------------------------------------------
# Read GWAS file (exposure dataset)
#-----------------------------------------------------
result_g <- fread("ukb-b-3768_full_by_R_code.csv")
dim(result_g)
cat(colnames(result_g),sep = "\n")
head(result_g)
#----------------------------------------------------------
#convert in P
#----------------------------------------------------------
result_g$P <- 10^(-result_g$LOGP)
write.csv(result_g, "ukb-b-3768_Pvalue.csv", row.names = FALSE)
head(result_g)
dim(result_g)
