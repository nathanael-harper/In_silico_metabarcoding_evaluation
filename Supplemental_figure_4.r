# ============================================================================
# Script: Integrated Supplemental Figure 4 Pipeline
# DESIGN: Species-level fate binning and proportional distribution plotting,
#         with gene mapping derived from the primer input file.
# ============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr required")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ----------------------------------------------------------------------------
# SECTION 1: CONFIGURATION & PATHS
# ----------------------------------------------------------------------------
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run-gr/"
EXPORT_DIR <- file.path(BASE, "results", "figures")
dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# SECTION 2: LOAD DATA
# ----------------------------------------------------------------------------
amps <- read.csv(file.path(BASE, "results", "07_amplicons.csv"), stringsAsFactors = FALSE)
prim <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE)
GR_SPP_FILE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Haplotype Analysis/2023_08_14_GR_spp.csv"
PRIMERS_TRUE_CSV <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Code/primers_input2.csv"

stopifnot(file.exists(GR_SPP_FILE), file.exists(file.path(BASE, "results", "07_amplicons.csv")))

# ----------------------------------------------------------------------------
# SECTION 3: PARSE BARCODES & MERGE PRIMER INFO
# ----------------------------------------------------------------------------
prim$assay <- sub("_[FR]$", "", prim$name)
prim$dir   <- sub("^[^_]+_", "", prim$name)
prim$fwd_len <- ifelse(prim$dir == "F", nchar(trimws(prim$seq)), 0)
prim$rev_len <- ifelse(prim$dir == "R", nchar(trimws(prim$seq)), 0)
prim_wide <- merge(prim[, c("assay", "fwd_len")], prim[, c("assay", "rev_len")], by = "assay")
amps <- merge(amps, prim_wide, by = "assay", all.x = TRUE)
amps$barcode_seq <- substr(amps$amp_seq, amps$fwd_len + 1, nchar(amps$amp_seq) - amps$rev_len)

# ----------------------------------------------------------------------------
# SECTION 4: BIN SPECIES VIA FUZZY MATCHING
# ----------------------------------------------------------------------------
gr_spp <- read.csv(GR_SPP_FILE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
gr_names <- trimws(as.character(gr_spp[["Scientific.Name"]]))

epithets <- trimws(sapply(strsplit(gr_spp$Scientific.Name, " "),
                          function(x) if(length(x) >= 2) x[2] else x[length(x)]))
epithet_to_full <- setNames(gr_spp$Scientific.Name, epithets)

matched_species <- vapply(amps$seq_id, function(id) {
  for (i in seq_along(epithets)) {
    if (length(agrep(tolower(epithets[i]), tolower(id), max.distance = 0.1, ignore.case = TRUE)) > 0) {
      return(epithet_to_full[i])
    }
  }
  NA_character_
}, character(1))

amps$species <- matched_species
cat(sprintf(">>> MATCHED: %d / %d amplicon records (%.1f%%)\n",
            sum(!is.na(matched_species)), length(matched_species),
            (sum(!is.na(matched_species)) / length(matched_species)) * 100))

# Map Taxa to standardized groups
amps$Taxa_clean <- ifelse(
  is.na(amps$Taxa) | grepl("Frogs|Toads|Anurans", amps$Taxa, ignore.case = TRUE),
  "Anurans",
  "Caudates"
)

# ----------------------------------------------------------------------------
# SECTION 5: GENE MAPPING
# ----------------------------------------------------------------------------
primers_true_df <- read.csv(PRIMERS_TRUE_CSV, stringsAsFactors = FALSE, check.names = FALSE)

assay_gene_map <- primers_true_df %>%
  mutate(Assay_Base = gsub("(_F|_R)$", "", name)) %>%
  distinct(Assay_Base, gene, .keep_all = TRUE) %>%
  # Normalize assay names to exactly match the pipeline's 07_amplicons.csv
  mutate(Assay = gsub("/", "-", Assay_Base)) %>%          # L2513/H2714 -> L2513-H2714
  mutate(Assay = sub("([A-Za-z])1$", "\\1", Assay)) %>%   # Batra1 -> Batra
  mutate(Assay = gsub("_F/extR", "-extR", Assay)) %>%     # MiAmphiL_28F/extR -> MiAmphiL_28-extR
  mutate(Assay = gsub("MiAmphiL", "MiAmphL", Assay)) %>%  # MiAmphiL -> MiAmphL
  rename(Assay = Assay)

# ----------------------------------------------------------------------------
# SECTION 6: CREATE GRID & ASSIGN GENES
# ----------------------------------------------------------------------------
assays <- unique(amps$assay)
grid <- expand.grid(assay = assays, species = gr_names, stringsAsFactors = FALSE)

# Join gene mapping
grid <- grid %>% left_join(assay_gene_map, by = c("assay" = "Assay"))

# Verify mapping success
if (any(is.na(grid$gene))) {
  stop("Gene mapping failed for: ", paste(unique(grid$assay[is.na(grid$gene)]), collapse = ", "))
}

# Rename to match downstream loop expectations
grid <- grid %>% rename(Gene = gene)

# Check for mapping errors
if (any(is.na(grid$Gene))) {
  stop("Gene mapping failed for: ", paste(unique(grid$assay[is.na(grid$Gene)]), collapse = ", "))
}

# Helper to check shared status
# Exclude NA species from comparisons to prevent logical errors
check_shared_status <- function(current_assay, current_species, current_barcode_seq) {
  other_idx <- amps$assay == current_assay & 
    !is.na(amps$species) & 
    amps$species != current_species & 
    amps$a_status == "suitable"
  if (!any(other_idx)) return(FALSE) 
  other_barcode_seqs <- amps$barcode_seq[other_idx]
  any(current_barcode_seq %in% other_barcode_seqs)
}

out_table <- data.frame(
  assay = grid$assay,
  species = grid$species,
  Gene = grid$Gene, # Gene is now baked in
  Taxa = NA_character_,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(out_table)) {
  a <- out_table$assay[i]
  s <- out_table$species[i]
  
  # Get Taxa for this species
  tax_val <- unique(amps$Taxa_clean[amps$assay == a & amps$species == s])
  out_table$Taxa[i] <- if(length(tax_val) > 0 && !is.na(tax_val[1])) tax_val[1] else "Caudates"
  
  sub_amps <- amps[amps$assay == a & amps$species == s, ]
  n <- nrow(sub_amps)
  
  if (n == 0) {
    out_table$fate[i] <- "No_coverage"
    next
  }
  
  # Priority hierarchy: Suitable > Unsuitable > Incomplete > No Coverage
  has_suitable <- any(sub_amps$a_status == "suitable")
  has_unsuitable <- any(sub_amps$a_status %in% c("unsuitable", "unmapped"))
  has_incomplete <- any(sub_amps$a_status == "incomplete_coverage")
  
  if (has_suitable) {
    suitable_barcode_seqs <- sub_amps$barcode_seq[sub_amps$a_status == "suitable"]
    if (check_shared_status(a, s, suitable_barcode_seqs)) {
      out_table$fate[i] <- "Shared_barcode"
    } else {
      out_table$fate[i] <- "Unique_barcode"
    }
  } else if (has_unsuitable) {
    out_table$fate[i] <- "No_predicted_amplification"
  } else if (has_incomplete) {
    out_table$fate[i] <- "Incomplete"
  } else {
    out_table$fate[i] <- "No_coverage" # Fallback for unmapped/no_coverage only
  }
}

# ----------------------------------------------------------------------------
# SECTION 7: TRANSFORM & EXPORT AUDIT TABLE
# ----------------------------------------------------------------------------
CATS <- c("No_coverage", "Incomplete", "No_predicted_amplification", "Shared_barcode", "Unique_barcode")
ASSAY_ORDER <- c("12Sv5","Ac12S","Am12S","AMP683","Batra","1121-1378",
                 "16S-AmTu","16sar_mod2","AM250","Amph_16S_1070-1340","BA-4445-178","L2513-H2714",
                 "MiAmphiL","MiAmphiS","MiAmphL_28-extR","Modified_16Sar","Ve16S","Vert-16S-eDNA","Ac16S")

# EXPORT: Intermediate species-level fate table for auditing
audit_file <- file.path(BASE, "results", "species_fate_audit.csv")
write.csv(out_table, audit_file, row.names = FALSE)
cat(sprintf("[AUDIT] Intermediate species-level fate table written to: %s\n", audit_file))
cat(sprintf("  Total unique species assessed: %d\n", length(unique(out_table$species))))

# Flatten grid data
df_flat <- data.frame(
  Assay = out_table$assay,
  Gene = out_table$Gene,
  Taxa = out_table$Taxa,
  Category = out_table$fate,
  stringsAsFactors = FALSE
)

# Count per Assay+Taxa+Category
df_counts <- as.data.frame(table(df_flat$Assay, df_flat$Gene, df_flat$Taxa, df_flat$Category))
names(df_counts) <- c("Assay", "Gene", "Taxa", "Category", "Count")

# Calculate Total per (Assay, Taxa, Gene) for proportion
df_totals <- aggregate(Count ~ Assay + Gene + Taxa, data = df_counts, FUN = sum)
names(df_totals)[4] <- "Total_Species"
facet_long <- merge(df_counts, df_totals, by = c("Assay", "Gene", "Taxa"))
facet_long$Prop <- facet_long$Count / facet_long$Total_Species

# Formatting
facet_long$Taxa  <- factor(facet_long$Taxa, levels = c("Anurans", "Caudates"))
facet_long$Gene  <- factor(facet_long$Gene, levels = c("mt12S", "mt16S"))
facet_long$Assay <- factor(facet_long$Assay, levels = c(ASSAY_ORDER, setdiff(unique(facet_long$Assay), ASSAY_ORDER)))
facet_long$Category <- factor(facet_long$Category, levels = CATS)

# ----------------------------------------------------------------------------
# SECTION 8: DIAGNOSTIC CHECK
# ----------------------------------------------------------------------------
cat("Gene Mapping Verification:\n")
print(table(facet_long$Assay, facet_long$Gene, useNA = "ifany"))

# ----------------------------------------------------------------------------
# SECTION 9: PLOTTING
# ----------------------------------------------------------------------------
p4 <- ggplot(facet_long[!is.na(facet_long$Prop),], aes(x = Assay, y = Prop, fill = Category)) +
  geom_col(color = "black", width = 0.7) +
  facet_grid(Taxa ~ Gene, scales = "free_x", space = "free_x",
             labeller = labeller(Taxa = c(Anurans = "Anurans", Caudates = "Caudates"))) +
  scale_fill_manual(name = "Template Fate",
                    values = c(No_coverage                = "#922b21",
                               Incomplete                 = "#ca6f1e",
                               No_predicted_amplification = "#f1c40f",
                               Shared_barcode             = "#2874a6",
                               Unique_barcode             = "#1e8449"),
                    labels = c("No Coverage",
                               "Incomplete","No Amplification", "Shared DNA Barcode", "Unique DNA Barcode")) +
  scale_y_continuous("Template Fate - Local Database", 
                     breaks = seq(0, 1, by = 0.25), 
                     limits = c(0, 1.08), 
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Candidate Metabarcoding Assay", 
       title = "Supplemental Figure 4: Species Diagnosis Distribution") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, lineheight = 0.75),
    axis.title  = element_text(size = 17, face = "bold"),
    axis.text   = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 12),
    axis.line = element_line(linewidth = 1)
  ) +
  labs(title = NULL, subtitle = NULL) +
  scale_x_discrete(labels = function(x) gsub("([_-])", "\\1\n", x)) 

# ----------------------------------------------------------------------------
# SECTION 10: OUTPUT
# ----------------------------------------------------------------------------
print(p4)
ggsave(file.path(EXPORT_DIR, "SuppFig4_species_diagnosis_barplot.png"), p4, width = 22, height = 10, dpi = 300)

cat("\n*** SUCCESS: Supplemental Figure 4 generated! ***\n")