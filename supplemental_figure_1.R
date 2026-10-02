## ============================================================================
## [5.5] HYBRID FIGURE — sequence-level jitter + per-assay boxplot
## ----------------------------------------------------------------------------
## Standalone figure script (Windows R; ggplot2 + dplyr). Runs AFTER step 2
## (needs results/06_assay_template_status.csv); does NOT depend on steps 3-5.
##
## INPUTS:
##   results/06_assay_template_status.csv  (step 2) -> per-record strand scores
##   primers_input2.csv                    -> assay list + gene ordering
##   taxon_template_directory_2.csv        -> Taxa / Gene levels
##
## INCLUSION: every record with BOTH strands scored (suitable OR unsuitable).
## The SCORE_THRESHOLD gate is intentionally NOT applied to this figure.
##
## X-AXIS: bare assay names, line break after _ or -.
##
## OUTPUT: results/figures/Primer_Full_Distribution_Hybrid_Final.png
## ============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## ---- BASE — the ONE line to edit (same line as steps 0-5) ----
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

# 1. LOAD INPUTS
primers_df <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE)
taxon_template <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"),
                           stringsAsFactors = FALSE)
STATUS_IN <- file.path(BASE, "results", "06_assay_template_status.csv")
if (!file.exists(STATUS_IN))
  stop("results/06_assay_template_status.csv not found — run steps 0-2 first (step 2 on WSL).")
status_df <- read.csv(STATUS_IN, stringsAsFactors = FALSE)

# 2.. PROCESS PRIMERS (The "Input Neutral" Logic)
# A. Remove _F and _R suffixes using regex ($ matches end of string)
# B. Deduplicate so 12Sv5_F and 12Sv5_R become one single '12Sv5' entry
primers_cleaned <- primers_df %>%
  mutate(Assay_Base = gsub("(_F|_R)$", "", name)) %>%
  distinct(Assay_Base, .keep_all = TRUE) %>%
  rename(Assay = Assay_Base)

# Define the order for the X-axis (alphabetical within genes: gene first, then assay)
ordered_assay_names <- primers_cleaned %>%
  arrange(gene, Assay) %>%
  pull(Assay)

# Define the order for Taxa and Genes
groups_all <- sort(unique(taxon_template$Taxa))
ordered_genes <- sort(unique(taxon_template$Gene))

# 3. DATA PROCESSING
# ALL scored records — suitable OR unsuitable (the SCORE_THRESHOLD gate is
# deliberately not applied here). NA guards match the pipeline idiom in
# steps 4 [4.2] / 5 [5.2].
has_total <- !is.na(status_df$f_p_status) & status_df$f_p_status == "scored" &
  !is.na(status_df$r_p_status) & status_df$r_p_status == "scored" &
  !is.na(status_df$f_score + status_df$r_score)

if (sum(has_total) == 0) {
  stop("CRITICAL ERROR: No records found with 'scored' status in ",
       "results/06_assay_template_status.csv — check step 2 output.")
}

plot_data <- data.frame(
  Assay  = status_df$assay[has_total],
  Gene   = status_df$Gene[has_total],
  Group  = factor(status_df$Taxa[has_total], levels = groups_all),
  Score  = status_df$f_score[has_total] + status_df$r_score[has_total],
  seq_id = status_df$seq_id[has_total],
  stringsAsFactors = FALSE)

# Ensure Assay and Gene are factors with the intended order
plot_data$Assay <- factor(plot_data$Assay, levels = ordered_assay_names)
plot_data$Gene  <- factor(plot_data$Gene,  levels = ordered_genes)

# 4. VALIDATE (assay names in 06 must match the cleaned primer names;
# unknown assays were turned into NA factor levels above)
if (any(is.na(plot_data$Assay))) {
  cat("\n--- ERROR DIAGNOSIS ---\n")
  cat("The 'assay' names in 06 do not match the cleaned primer names.\n")
  cat("Example names in status_df$assay:", head(unique(status_df$assay), 3), "\n")
  cat("Example names in primers_input2.csv (cleaned):", head(unique(primers_cleaned$Assay), 3), "\n")
  cat("-----------------------\n")
  stop("Stop: Assay name mismatch. Check that status_df$assay carries the base name.")
}

# 5. THE PLOT
PENALTY_LIMIT <- max(150L, as.integer(ceiling(max(plot_data$Score))) + 15L)

p <- ggplot(plot_data, aes(x = Assay, y = Score)) +
  # Jitter points
  geom_jitter(aes(fill = Group),
              position = position_jitter(width = 0.3, height = 0),
              alpha = 0.6, size = 3, shape = 21, color = "black") +
  # Boxplot
  geom_boxplot(aes(group = Assay), outlier.shape = NA, alpha = 0.4, color = "black", fill = "grey90") +
  # Facet by Gene (alphabetical)
  facet_wrap(~Gene, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, PENALTY_LIMIT + 15)) +
  labs(
    title = NULL,
    subtitle = NULL,
    y = "Glotbal Database - PrimerMiner Penalty Score",
    x = "Candidate Metabarcoding Assay"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.8, size = 10),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 11),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.line = element_line(linewidth = 1)
  ) +
  # Dynamic colors for Taxa
  scale_fill_brewer(palette = "Set2", name = "Taxa") +
  # Bare assay names with a newline after _ or - for prettier X-axis labels
  scale_x_discrete(labels = function(x) gsub("([_-])", "\\1\n", x))

# 6. OUTPUT
EXPORT_DIR <- file.path(BASE, "results", "figures")
dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(EXPORT_DIR, "Primer_Full_Distribution_Hybrid_Final.png"), p,
       width = 18, height = 10, dpi = 300)

cat("\n*** SUCCESS: Hybrid plot generated with cleaned primer names! ***\n")