## ============================================================================
## [5.6] SUPPLEMENTAL FIGURE 2 — Stacked Barplot of Template Fates (Non-Hardcoded)
##       RE-MADE: Follows Figure 1 logic for ordering and labeling.
## ============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# 1. CONFIGURATION & PATHS
# Using BASE from your environment/previous steps
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"
OUT  <- file.path(BASE, "results")
EXPORT_DIR <- file.path(OUT, "figures")
dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

## Visual Configuration (Preserved from legacy)
FATE_VISUALS <- list(
  colors = c(
    "No_coverage"                = "#922b21",
    "Incomplete"                 = "#ca6f1e",
    "No_predicted_amplification" = "#f1c40f",
    "Shared_barcode"             = "#2874a6",
    "Unique_barcode"             = "#1e8449"
  ),
  labels = c(
    "No_coverage"                = "No Coverage",
    "Incomplete"                 = "Incomplete Coverage",
    "No_predicted_amplification" = "No Amplification",
    "Shared_barcode"             = "Shared DNA Barcode",
    "Unique_barcode"             = "Unique DNA Barcode"
  )
)

# 2. LOAD INPUTS
primers_df <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE)
taxon_template <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE)
res08 <- read.csv(file.path(OUT, "08_barcode_resolution_summary.csv"), stringsAsFactors = FALSE)
status_df <- read.csv(file.path(OUT, "06_assay_template_status.csv"), stringsAsFactors = FALSE)

# 3. NON-HARDCODED ORDERING LOGIC (The "Figure 1 Way")

# A. Determine Assay Order: Clean primers, then sort by Gene, then Name
primers_cleaned <- primers_df %>%
  mutate(Assay_Base = gsub("(_F|_R)$", "", name)) %>%
  distinct(Assay_Base, gene, .keep_all = TRUE) %>%
  arrange(gene, Assay_Base)

ordered_assays <- primers_cleaned$Assay_Base

# B. Determine Taxa and Gene order from the template file
ordered_taxa <- sort(unique(taxon_template$Taxa))
ordered_genes <- sort(unique(taxon_template$Gene))

# 4. DATA PROCESSING

# A. Aggregate Status Counts (from 06) to get denominators
# We need: Total, NoCov, Indel, and the "Suitable" (Good) count
status_summary <- status_df %>%
  group_by(assay, Taxa, Gene) %>%
  summarise(
    A_Total      = n(),
    A_NoCov      = sum(a_status == "no_coverage"),
    A_Indel      = sum(a_status == "incomplete_coverage"),
    A_unsuitable  = sum(a_status == "unsuitable"),
    A_unmapped    = sum(a_status == "unmapped"),
    A_Good       = sum(a_status == "suitable"),
    .groups = "drop"
  )

# B. Merge with Resolution Summary (from 08)
# 08 provides: n_species_unique and n_species_in_shared
m <- merge(status_summary, res08, by = c("assay", "Taxa", "Gene"), all.x = TRUE)

# Clean up NAs from the merge
m$n_species_unique[is.na(m$n_species_unique)]       <- 0
m$n_species_in_shared[is.na(m$n_species_in_shared)] <- 0

# C. Define the 5 Fate Categories
m$No_coverage                <- m$A_NoCov
m$Incomplete                 <- m$A_Indel
m$No_predicted_amplification <- m$A_unsuitable + m$A_unmapped
m$Shared_barcode             <- m$n_species_in_shared
m$Unique_barcode             <- m$n_species_unique

# D. Calculate Proportions (for the Y-axis)
CATS <- c("No_coverage", "Incomplete", "No_predicted_amplification",
          "Shared_barcode", "Unique_barcode")

m[, CATS] <- m[, CATS] / m$A_Total

# E. Convert to Long Format for ggplot2
long <- do.call(rbind, lapply(CATS, function(f) {
  data.frame(Assay = m$assay, Taxa = m$Taxa, Gene = m$Gene,
             category = f, prop = m[[f]], stringsAsFactors = FALSE)
}))

# 5. FACTORING (Applying our non-hardcoded orders)
long$Assay    <- factor(long$Assay, levels = ordered_assays)
long$Taxa     <- factor(long$Taxa, levels = ordered_taxa)
long$Gene     <- factor(long$Gene, levels = ordered_genes)
long$category <- factor(long$category, levels = CATS)

# 6. THE PLOT
p <- ggplot(long, aes(x = Assay, y = prop, fill = category)) +
  geom_col(color = "black", linewidth = 0.2) +
  # Facet by Taxa (Rows) and Gene (Cols)
  facet_grid(Taxa ~ Gene, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = "Template Fate",
                    values = FATE_VISUALS$colors,
                    labels = FATE_VISUALS$labels) +
  scale_y_continuous("Template Fate - Global Database",
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1.05),
                     expand = expansion(mult = c(0, 0))) +
  # The Figure 1 Style Labeling (Regex gsub)
  scale_x_discrete(labels = function(x) gsub("([_-])", "\\1\n", x)) +
  labs(x = "Candidate Metabarcoding Assay", title = NULL, subtitle = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.8, size = 10),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 10),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# 7. OUTPUT
print(p)

# Save wide data for re-plotting if needed
write.csv(m[, c("assay", "Taxa", "Gene", "A_Total", CATS)],
          file.path(OUT, "14_stacked_plot_data_cleaned.csv"), row.names = FALSE)

ggsave(file.path(EXPORT_DIR, "14_stacked_barplot_cleaned.png"), p, width = 20, height = 10, dpi = 300)

cat("\n*** SUCCESS: Non-hardcoded Figure 2 generated! ***\n")
cat("Assay order derived from primers_input2.csv\n")
cat("Taxa/Gene structure derived from taxon_template_directory_2.csv\n")
cat("X-axis labels formatted via gsub (Figure 1 style)\n")