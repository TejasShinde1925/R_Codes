#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")
library(devtools)
library(TwoSampleMR)

setwd("E:/Final_Data_GBC_Paper/MR_11_Sept_2026/With_SAS/With_clumping")


#############################################################################################
#Extract instruments from Opengwas
#############################################################################################

gs_exp_dat <- fread("exp15_sm.csv")
dim(gs_exp_dat)
min(gs_exp_dat$pval)
max(gs_exp_dat$pval)
names(gs_exp_dat)


#############################################################################################
# Exposrure data
#############################################################################################
gs_exp_dat <- read_exposure_data(
  filename = "exp15_sm.csv",
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
)

#write.csv(gs_exp_dat, "EXPOSURE_DATA.csv", row.names = FALSE)
#dim(gs_exp_dat)

#############################################################################################
# Clumping
############################################################################################# 
#gs_exp_dat <- clump_data(gs_exp_dat)

gs_exp_dat <- clump_data(
  gs_exp_dat,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  #pop = "SAS",
  bfile = "E:/Final_Data_GBC_Paper/MR_11_Sept_2026/With_SAS/With_clumping/1000_genome_SAS/1000G_SAS",
  plink_bin = "C:/Users/admin/Desktop/plink_win64_20250806/plink.exe"
)

write.csv(gs_exp_dat,"clumped_data.csv", row.names=FALSE)
dim(gs_exp_dat)
min(gs_exp_dat$pval.exposure)
max(gs_exp_dat$pval.exposure)

gs_exp_dat <- fread("clumped_data.csv")
dim(gs_exp_dat)
min(gs_exp_dat$pval.exposure)
max(gs_exp_dat$pval.exposure)

write.csv(gs_exp_dat,"EXPOSURE_DATA.csv", row.names=FALSE)

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

############################################################################################# 
out_data_GBC_meta <- fread("matched_rsid_data_GBC.csv")
dim(out_data_GBC_meta)
names(out_data_GBC_meta)
############################################################################################# 
gbc_out_dat <- read_outcome_data(
  filename = "matched_rsid_data_GBC.csv",
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1_MinorAllele",
  other_allele_col = "A2_MajorAllele",
  eaf_col = "A1_MinorAlleleFreq",
  pval_col = "P",
)

write.csv(gbc_out_dat, "OUTCOME_DATA.CSV", row.names = FALSE)
dim(gbc_out_dat)
min(gbc_out_dat$pval.outcome)
max(gbc_out_dat$pval.outcome)

dat <- harmonise_data(
  exposure_dat = gs_exp_dat, 
  outcome_dat = gbc_out_dat
)

dat <- harmonise_data(gs_exp_dat, gbc_out_dat, action=2)
write.csv(dat, "HARMONISED_ACTION_2.CSV", row.names = FALSE)



harmonise_data(gs_exp_dat, gbc_out_dat, action=3)

res <- mr(dat)
res
write.csv(res, file = "mr_results.csv", row.names=T)
write.csv(generate_odds_ratios(res), file = "generate_odds_ratios.csv", row.names=T)
mr_heterogeneity(dat)
write.csv(mr_heterogeneity(dat), file = "mr_heterogeneity.csv", row.names=T)
res_single <- mr_singlesnp(dat)
write.csv(mr_singlesnp(dat), file = "res_single.csv", row.names=T)
res_loo <- mr_leaveoneout(dat)
write.csv(mr_leaveoneout(dat), file = "leaveoneout.csv", row.names=T)
mr_pleiotropy_test(dat)
write.csv(mr_pleiotropy_test(dat), file = "mr_pleiotropy_test.csv", row.names=T)


res <- mr(dat)
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_two_sample_ml"))
p2 <- mr_forest_plot(res_single)
p2[[1]]
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

mr_report(dat)


###############################################################################################################
#
#
###############################################################################################################
#Assess outliers
#Requires following output from the MR base (2SMR) package: 
# dat = harmonized dataset for exposure and outcome: 
#dat <- harmonise_data(exposure_dat, outcome_dat = chd_out_dat

# res_single = single SNP associations
# res_single <- mr_singlesnp(dat)

#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/25, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#2 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)
write.csv(generate_odds_ratios(mr_results2), file = "mr_odds_outliers2.csv", row.names=T)


results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nooutliers.txt")

#MR presso 


devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
write.csv(mr_presso, file = "mr_presso.csv", row.names=T)
