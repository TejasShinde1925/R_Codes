# ============================================================
# GBC PheWAS VOLCANO PLOT
# 8 Lead SNPs
# ============================================================
#
# BIOLOGICAL CATEGORIES SHOWN:
#
# 1. Hepatobiliary
# 2. Metabolic/Lipids
# 3. Hematological
#
# IMPORTANT:
# Inflammatory/Autoimmune traits are REMOVED completely.
# Traits classified as "Other" are also REMOVED completely.
#
# Labels:
# ONLY distinct traits with P < 1 × 10^-8
#
# Special clinical endpoints:
# K80  = Cholelithiasis
# J18.3 = Cholecystectomy
#
# ============================================================


# ============================================================
# 1. SET WORKING DIRECTORY
# ============================================================

setwd(
  "E:/Final_Data_GBC_Paper/Volcano_plot"
)


# ============================================================
# 2. LOAD LIBRARIES
# ============================================================

library(readxl)
library(data.table)
library(ggplot2)
library(ggrepel)
library(grid)


# ============================================================
# 3. READ DATA
# ============================================================

data <- fread(
  "phewas_result_of_GBC_8_lead_snps_data.csv"
)


# ============================================================
# 4. CHECK COLUMN NAMES
# ============================================================

cat(
  "\n==============================================\n"
)

cat(
  "COLUMN NAMES IN DATASET\n"
)

cat(
  "==============================================\n\n"
)

print(
  colnames(data)
)


# ============================================================
# 5. CONVERT REQUIRED COLUMNS
# ============================================================

data$OR <- as.numeric(
  data$OR
)

data$pval <- as.numeric(
  data$p
)


# ============================================================
# 6. REMOVE INVALID VALUES
# ============================================================

data <- data[
  !is.na(OR) &
    !is.na(pval) &
    pval > 0 &
    OR > 0
]


# ============================================================
# 7. CALCULATE -log10(P)
# ============================================================

data$logP <- -log10(
  data$pval
)


# ============================================================
# 8. DEFINE MAIN SIGNIFICANCE THRESHOLD
# ============================================================

# Main PheWAS threshold:
#
# P = 1 × 10^-5

p_threshold <- 1e-5

threshold <- -log10(
  p_threshold
)


# ============================================================
# 9. DEFINE LABEL THRESHOLD
# ============================================================

# Only traits with:
#
# P < 1 × 10^-8
#
# will be labelled.

label_p_threshold <- 1e-8

label_threshold <- -log10(
  label_p_threshold
)


# ============================================================
# 10. CREATE RISK / PROTECTIVE GROUP
# ============================================================

data$Group <- "Normal"


# ------------------------------------------------------------
# Risk
# ------------------------------------------------------------

data$Group[
  data$pval < p_threshold &
    data$OR > 1
] <- "Risk"


# ------------------------------------------------------------
# Protective
# ------------------------------------------------------------

data$Group[
  data$pval < p_threshold &
    data$OR < 1
] <- "Protective"


data$Group <- factor(
  data$Group,
  
  levels = c(
    "Normal",
    "Protective",
    "Risk"
  )
)


# ============================================================
# 11. CREATE BIOLOGICAL CATEGORY
# ============================================================

# Initially classify everything as Other.
#
# Later, Other and Inflammatory/Autoimmune traits
# will be removed.

data$Category <- "Other"


# ============================================================
# 12. HEPATOBILIARY
# ============================================================

data$Category[
  grepl(
    
    paste(
      c(
        "gallbladder",
        "gall bladder",
        "biliary",
        "cholelithiasis",
        "gall stone",
        "gallstone",
        "cholecystectomy",
        "cholecyst",
        "K80",
        "J18.3",
        "alanine aminotransferase",
        "alanine transaminase",
        "ALT",
        "aspartate aminotransferase",
        "aspartate transaminase",
        "AST",
        "alkaline phosphatase",
        "ALP",
        "bilirubin",
        "liver",
        "hepatic"
      ),
      
      collapse = "|"
    ),
    
    data$trait,
    
    ignore.case = TRUE
    
  )
] <- "Hepatobiliary"


# ============================================================
# 13. METABOLIC / LIPIDS
# ============================================================

data$Category[
  data$Category == "Other" &
    
    grepl(
      
      paste(
        c(
          "LDL",
          "cholesterol",
          "apolipoprotein",
          "apoliprotein",
          "ApoB",
          "triglyceride",
          "HDL",
          "lipid",
          "body mass index",
          "BMI",
          "body weight",
          "weight",
          "fat mass",
          "trunk fat",
          "waist circumference",
          "hip circumference",
          "adiposity",
          "urate",
          "uric acid",
          "diabetes",
          "glyc"
        ),
        
        collapse = "|"
      ),
      
      data$trait,
      
      ignore.case = TRUE
      
    )
] <- "Metabolic/Lipids"


# ============================================================
# 14. HEMATOLOGICAL
# ============================================================

data$Category[
  data$Category == "Other" &
    
    grepl(
      
      paste(
        c(
          "platelet",
          "haemoglobin",
          "hemoglobin",
          "red blood",
          "white blood",
          "erythrocyte",
          "leukocyte",
          "lymphocyte",
          "neutrophil",
          "monocyte",
          "eosinophil",
          "basophil",
          "mean corpuscular",
          "haematocrit",
          "hematocrit",
          "blood cell"
        ),
        
        collapse = "|"
      ),
      
      data$trait,
      
      ignore.case = TRUE
      
    )
] <- "Hematological"


# ============================================================
# 15. INFLAMMATORY / AUTOIMMUNE
# ============================================================

# This classification is retained here so that the
# original classification logic remains unchanged.
#
# These traits will be removed in Step 16.

data$Category[
  data$Category == "Other" &
    
    grepl(
      
      paste(
        c(
          "rheumatoid arthritis",
          "arthritis",
          "autoimmune",
          "inflammatory",
          "inflammation",
          "hyperthyroidism",
          "thyrotoxicosis",
          "hypothyroidism",
          "myxoedema",
          "thyroid",
          "lupus",
          "crohn",
          "ulcerative colitis",
          "psoriasis",
          "immune"
        ),
        
        collapse = "|"
      ),
      
      data$trait,
      
      ignore.case = TRUE
      
    )
] <- "Inflammatory/Autoimmune"


# ============================================================
# 16. REMOVE INFLAMMATORY / AUTOIMMUNE
# ============================================================
#
# THIS IS THE ONLY MAIN CHANGE.
#
# All Inflammatory/Autoimmune traits are removed completely.
#
# ============================================================

data <- data[
  Category != "Inflammatory/Autoimmune"
]


# ============================================================
# 17. REMOVE OTHER
# ============================================================
#
# Remove all traits that do not belong to the remaining
# three biological categories.
#
# ============================================================

data <- data[
  Category != "Other"
]


# ============================================================
# 18. CONVERT CATEGORY TO FACTOR
# ============================================================

data$Category <- factor(
  
  data$Category,
  
  levels = c(
    "Hepatobiliary",
    "Metabolic/Lipids",
    "Hematological"
  )
)


# ============================================================
# 19. SHOW CATEGORY COUNTS
# ============================================================

cat(
  "\n==============================================\n"
)

cat(
  "FINAL BIOLOGICAL CATEGORY COUNTS\n"
)

cat(
  "==============================================\n\n"
)

print(
  table(
    data$Category
  )
)


# ============================================================
# 20. CREATE SPECIAL CLINICAL ENDPOINT TYPES
# ============================================================
#
# IMPORTANT:
# "Other" here is NOT a biological category.
#
# It simply means:
# not K80 and not J18.3.
#
# ============================================================

data$Trait_Type <- "Other"


# ============================================================
# 21. CHOLELITHIASIS / K80
# ============================================================

data$Trait_Type[
  grepl(
    
    "K80|cholelithiasis|gall stone|gallstone",
    
    data$trait,
    
    ignore.case = TRUE
    
  )
] <- "Cholelithiasis (K80)"


# ============================================================
# 22. CHOLECYSTECTOMY / J18.3
# ============================================================

data$Trait_Type[
  grepl(
    
    "J18\\.3|cholecystectomy|gall bladder removal|gallbladder removal",
    
    data$trait,
    
    ignore.case = TRUE
    
  )
] <- "Cholecystectomy (J18.3)"


# ============================================================
# 23. AUTOMATIC LABEL SELECTION
# ============================================================

# ONLY:
#
# P < 1 × 10^-8
#
# will be labelled.

label_data <- data[
  pval < label_p_threshold
]


# ============================================================
# 24. SORT BY P-VALUE
# ============================================================

label_data <- label_data[
  order(pval)
]


# ============================================================
# 25. KEEP ONE LABEL PER DISTINCT TRAIT
# ============================================================

label_data <- label_data[
  !duplicated(trait)
]


# ============================================================
# 26. PRINT LABELLED TRAITS
# ============================================================

cat(
  "\n==============================================\n"
)

cat(
  "AUTOMATIC LABEL SELECTION\n"
)

cat(
  "==============================================\n"
)

cat(
  "Label criterion: P < 1 × 10^-8\n"
)

cat(
  "Number of distinct labelled traits: ",
  
  nrow(label_data),
  
  "\n",
  
  sep = ""
)

cat(
  "==============================================\n\n"
)


print(
  label_data[
    ,
    .(
      Trait = trait,
      OR = OR,
      P_value = pval,
      Minus_log10_P = logP,
      Category = Category,
      Trait_Type = Trait_Type
    )
  ]
)


# ============================================================
# 27. MAXIMUM Y VALUE
# ============================================================

max_logP <- max(
  data$logP,
  na.rm = TRUE
)


# ============================================================
# 28. CREATE VOLCANO PLOT
# ============================================================

p <- ggplot(
  
  data,
  
  aes(
    x = OR,
    y = logP
  )
  
) +
  
  
  # ==========================================================
# ALL POINTS
# ==========================================================

geom_point(
  
  aes(
    color = Category,
    shape = Trait_Type
  ),
  
  size = 2.7,
  
  alpha = 0.80
  
) +
  
  
  # ==========================================================
# P = 1 × 10^-5
# ==========================================================

geom_hline(
  
  yintercept = threshold,
  
  linetype = "dashed",
  
  linewidth = 0.65,
  
  color = "grey35"
  
) +
  
  
  # ==========================================================
# OR = 1
# ==========================================================

geom_vline(
  
  xintercept = 1,
  
  linetype = "dashed",
  
  linewidth = 0.65,
  
  color = "grey35"
  
) +
  
  
  # ==========================================================
# P = 1 × 10^-8 LABEL THRESHOLD
# ==========================================================

geom_hline(
  
  yintercept = label_threshold,
  
  linetype = "dotted",
  
  linewidth = 0.55,
  
  color = "grey55"
  
) +
  
  
  # ==========================================================
# AUTOMATIC LABELS
#
# ONLY P < 1 × 10^-8
#
# ONE LABEL PER DISTINCT TRAIT
# ==========================================================

geom_label_repel(
  
  data = label_data,
  
  aes(
    
    label = trait,
    
    fill = Category
    
  ),
  
  # --------------------------------------------------------
  # LABEL SIZE
  # --------------------------------------------------------
  
  size = 5.2,
  
  color = "white",
  
  fontface = "plain",
  
  
  # --------------------------------------------------------
  # LABEL BOX
  # --------------------------------------------------------
  
  label.padding = unit(
    
    0.30,
    
    "lines"
    
  ),
  
  label.r = unit(
    
    0.12,
    
    "lines"
    
  ),
  
  
  # --------------------------------------------------------
  # LABEL SPACING
  # --------------------------------------------------------
  
  box.padding = 0.75,
  
  point.padding = 0.30,
  
  
  # --------------------------------------------------------
  # CONNECTOR
  # --------------------------------------------------------
  
  segment.color = "grey35",
  
  segment.size = 0.40,
  
  min.segment.length = 0,
  
  
  # --------------------------------------------------------
  # REPULSION
  # --------------------------------------------------------
  
  force = 5,
  
  force_pull = 0.15,
  
  direction = "both",
  
  
  # --------------------------------------------------------
  # SHOW ALL QUALIFYING LABELS
  # --------------------------------------------------------
  
  max.overlaps = Inf,
  
  max.time = 120,
  
  max.iter = 50000,
  
  seed = 123,
  
  show.legend = FALSE
  
) +
  
  
  # ==========================================================
# PUBLICATION-FRIENDLY COLORS
# ==========================================================

scale_color_manual(
  
  name = "Biological Category",
  
  values = c(
    
    "Hepatobiliary" =
      "#D55E00",
    
    "Metabolic/Lipids" =
      "#009E73",
    
    "Hematological" =
      "#0072B2"
    
  )
  
) +
  
  
  # ==========================================================
# LABEL BOX COLORS
# ==========================================================

scale_fill_manual(
  
  values = c(
    
    "Hepatobiliary" =
      "#D55E00",
    
    "Metabolic/Lipids" =
      "#009E73",
    
    "Hematological" =
      "#0072B2"
    
  ),
  
  guide = "none"
  
) +
  
  
  # ==========================================================
# CLINICAL ENDPOINT SHAPES
# ==========================================================

scale_shape_manual(
  
  name = "Clinical Endpoint",
  
  values = c(
    
    "Other" = 16,
    
    "Cholelithiasis (K80)" = 17,
    
    "Cholecystectomy (J18.3)" = 15
    
  )
  
) +
  
  
  # ==========================================================
# X AXIS
# ==========================================================

scale_x_continuous(
  
  limits = c(
    
    0.2,
    
    1.8
    
  ),
  
  breaks = seq(
    
    0.2,
    
    1.8,
    
    0.2
    
  ),
  
  expand = c(
    
    0,
    
    0
    
  )
  
) +
  
  
  # ==========================================================
# Y AXIS
# ==========================================================

scale_y_continuous(
  
  limits = c(
    
    0,
    
    max_logP + 8
    
  ),
  
  breaks = seq(
    
    0,
    
    ceiling(max_logP / 20) * 20,
    
    20
    
  ),
  
  expand = expansion(
    
    mult = c(
      
      0.02,
      
      0.04
      
    )
    
  )
  
) +
  
  
  # ==========================================================
# AXIS LABELS
# ==========================================================

labs(
  
  x =
    "Odds Ratio (OR)",
  
  y =
    expression(
      -log[10](P~value)
    )
  
) +
  
  
  # ==========================================================
# THEME
# ==========================================================

theme_classic(
  
  base_size = 14
  
) +
  
  theme(
    
    # --------------------------------------------------------
    # AXIS TITLE
    # --------------------------------------------------------
    
    axis.title =
      
      element_text(
        
        face = "bold",
        
        size = 15
        
      ),
    
    
    # --------------------------------------------------------
    # AXIS TEXT
    # --------------------------------------------------------
    
    axis.text =
      
      element_text(
        
        size = 12,
        
        color = "black"
        
      ),
    
    
    # --------------------------------------------------------
    # LEGEND TITLE
    # --------------------------------------------------------
    
    legend.title =
      
      element_text(
        
        face = "bold",
        
        size = 13
        
      ),
    
    
    # --------------------------------------------------------
    # LEGEND TEXT
    # --------------------------------------------------------
    
    legend.text =
      
      element_text(
        
        size = 11
        
      ),
    
    
    # --------------------------------------------------------
    # LEGEND POSITION
    # --------------------------------------------------------
    
    legend.position =
      
      "right",
    
    
    # --------------------------------------------------------
    # BOX AROUND LEGENDS
    # --------------------------------------------------------
    
    legend.background =
      
      element_rect(
        
        fill = "white",
        
        color = "black",
        
        linewidth = 0.7
        
      ),
    
    
    # --------------------------------------------------------
    # LEGEND KEY
    # --------------------------------------------------------
    
    legend.key =
      
      element_blank(),
    
    
    # --------------------------------------------------------
    # LEGEND MARGIN
    # --------------------------------------------------------
    
    legend.margin =
      
      margin(
        
        8,
        
        10,
        
        8,
        
        10
        
      ),
    
    
    # --------------------------------------------------------
    # PLOT BORDER
    # --------------------------------------------------------
    
    panel.border =
      
      element_rect(
        
        color = "black",
        
        fill = NA,
        
        linewidth = 0.6
        
      ),
    
    
    # --------------------------------------------------------
    # REMOVE GRID
    # --------------------------------------------------------
    
    panel.grid.major =
      
      element_blank(),
    
    panel.grid.minor =
      
      element_blank(),
    
    
    # --------------------------------------------------------
    # PLOT MARGIN
    # --------------------------------------------------------
    
    plot.margin =
      
      margin(
        
        20,
        
        30,
        
        20,
        
        20
        
      )
    
  ) +
  
  
  # ==========================================================
# DO NOT CLIP LABELS
# ==========================================================

coord_cartesian(
  
  clip = "off"
  
)


# ============================================================
# 29. DISPLAY PLOT
# ============================================================

print(p)


# ============================================================
# 30. SAVE PNG
# ============================================================

ggsave(
  
  "23_GBC_PheWAS_Volcano_without_Inflammatory_Autoimmune.png",
  
  plot = p,
  
  width = 20,
  
  height = 11,
  
  units = "in",
  
  dpi = 600,
  
  bg = "white"
  
)


# ============================================================
# 31. SAVE PDF
# ============================================================

ggsave(
  
  "23_GBC_PheWAS_Volcano_without_Inflammatory_Autoimmune.pdf",
  
  plot = p,
  
  width = 20,
  
  height = 11,
  
  units = "in",
  
  device = "pdf",
  
  bg = "white"
  
)


# ============================================================
# 32. SAVE TIFF
# ============================================================

ggsave(
  
  "23_GBC_PheWAS_Volcano_without_Inflammatory_Autoimmune.tiff",
  
  plot = p,
  
  width = 20,
  
  height = 11,
  
  units = "in",
  
  dpi = 600,
  
  compression = "lzw",
  
  bg = "white"
  
)


# ============================================================
# 33. FINAL INFORMATION
# ============================================================

cat(
  "\n\n==============================================\n"
)

cat(
  "GBC PheWAS VOLCANO PLOT COMPLETED\n"
)

cat(
  "==============================================\n"
)

cat(
  "Biological categories shown:\n"
)

cat(
  "  1. Hepatobiliary\n"
)

cat(
  "  2. Metabolic/Lipids\n"
)

cat(
  "  3. Hematological\n"
)

cat(
  "\nInflammatory/Autoimmune traits: REMOVED\n"
)

cat(
  "Other biological-category traits: REMOVED\n"
)

cat(
  "\nAssociation threshold : P < 1 × 10^-5\n"
)

cat(
  "Label threshold       : P < 1 × 10^-8\n"
)

cat(
  "Label selection       : AUTOMATIC\n"
)

cat(
  "Labels                : DISTINCT TRAITS ONLY\n"
)

cat(
  "Clinical endpoints     : K80 / J18.3\n"
)

cat(
  "Plot title             : NONE\n"
)

cat(
  "Label size             : 5.2\n"
)

cat(
  "Resolution             : 600 DPI\n"
)

cat(
  "==============================================\n\n"
)