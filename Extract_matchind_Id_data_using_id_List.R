#############################################################################################
# Find matching snps from GBC Meta Gwas Data
#############################################################################################

id_list <- read.csv("id_list.csv", stringsAsFactors = FALSE)
reference <- fread("G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv")
dim(reference)
min(reference$P)
max(reference$P)

result <- reference[reference$SNP %in% id_list$SNP, ]

dim(result)
write.csv(result, "matched_rsid_data_GBC.csv", row.names = FALSE)
names(result)
