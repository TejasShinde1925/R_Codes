setwd("C:/Users/admin/Desktop/CYTOBAND_CODE")

##############################
# STEP 1: Install packages (run only once)
##############################
install.packages("BiocManager")

BiocManager::install(c("biomaRt", "GenomicRanges", "AnnotationHub"))


##############################
# STEP 2: Load libraries
##############################
library(biomaRt)
library(GenomicRanges)
library(AnnotationHub)


##############################
# STEP 3: Connect to Ensembl GRCh37 (IMPORTANT)
##############################
mart <- useEnsembl(
  biomart = "snp",
  dataset = "hsapiens_snp",
  GRCh = 37
)


##############################
# STEP 4: Input your SNPs
##############################
snps <- c("rs1059542","rs111280611","rs2312959")


##############################
# STEP 5: Get chromosome + position
##############################
pos <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start"),
  filters = "snp_filter",
  values = snps,
  mart = mart
)

print(pos)


##############################
# STEP 6: Load hg19 cytoband
##############################
ah <- AnnotationHub()

# Search hg19 cytoband
query(ah, "cytoBand hg19")

# Use correct ID (commonly AH5086)
cyto <- ah[["AH5012"]]

pos$chr_name <- paste0("chr", pos$chr_name)




##############################
# STEP 7: Convert SNPs to GRanges
##############################
#gr <- GRanges(
#  seqnames = paste0("chr", pos$chr_name),
#  ranges = IRanges(
#    start = pos$chrom_start,
#    end = pos$chrom_start
#  )
#)
gr <- GRanges(
  seqnames = pos$chr_name,
  ranges = IRanges(start = pos$chrom_start, end = pos$chrom_start)
)


##############################
# STEP 8: Find overlaps (cytoband)
##############################
hits <- findOverlaps(gr, cyto)


##############################
# STEP 9: Add cytoband to result
##############################
pos$cytoband <- NA
pos$cytoband[queryHits(hits)] <- mcols(cyto)$name[subjectHits(hits)]


##############################
# STEP 10: Final output
##############################
print(pos)


##############################
# STEP 11: Save output (optional)
##############################
write.csv(pos, "cytoband_output_hg19.csv", row.names = FALSE)
