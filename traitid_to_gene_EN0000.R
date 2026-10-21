setwd("E:/ebiomedicine_phewas_MR_25_May_2026/Phewas")



#install.packages("BiocManager")
#BiocManager::install("biomaRt")

# Install packages (run once)
install.packages("BiocManager")
BiocManager::install("biomaRt")

# Load library
library(biomaRt)

# Read your CSV file
data <- read.csv("trait_id.csv")

# Connect to Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene names
gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ),
  
  filters = "ensembl_gene_id",
  
  values = data$trait,
  
  mart = mart
)

# View result
head(gene_info)

# Save output
write.csv(gene_info, "gene_name_output.csv", row.names = FALSE)


