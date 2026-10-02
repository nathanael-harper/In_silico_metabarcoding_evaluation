# ============================================================================
# Script: Integrated Supplemental Figure 3 Pipeline (RELIABLE WORKFLOW)
# DESIGN: Maintains the original TWO-STEP approach to ensure stability.
#         Section 1 writes a CSV; Section 2 reads that CSV and plots it.
# =============================================================================

if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required")
if (!requireNamespace("tidyr", quietly = TRUE) ) stop("tidyr required")
if (!requireNamespace("ggplot2", quietly = TRUE) ) stop("ggplot2 required")
if (!requireNamespace("ggnewscale", quietly = TRUE)) stop("ggnewscale required")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggnewscale)
})

# 1. CONFIGURATION & PATHS
# ----------------------------------------------------------------------------
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run-gr"
OUT  <- file.path(BASE, "results")
EXPORT_DIR <- file.path(OUT, "figures")

AMPLICON_CSV    <- file.path(BASE, "results", "07_amplicons.csv")
STATUS_CSV      <- file.path(BASE, "results", "06_assay_template_status.csv")
GR_SPP_CSV      <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Haplotype Analysis/2023_08_14_GR_spp.csv"
PRIMERS_TRUE_CSV <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Code/primers_input2.csv"

SUMMARY_OUTPUT_CSV <- file.path(EXPORT_DIR, "SuppFig3_primer_species_summary.csv")
PLOT_OUTPUT_PNG    <- file.path(EXPORT_DIR, "SuppFig3_species_averaged_plot.png")

if (!dir.exists(EXPORT_DIR)) dir.create(EXPORT_DIR, recursive = TRUE)

# ----------------------------------------------------------------------------
# SECTION 1: DATA GENERATION & FUZZY MATCHING (UPDATED)
# ----------------------------------------------------------------------------
message(">>> STARTING STEP 1: Data Generation...")

amps_df <- read.csv(AMPLICON_CSV, stringsAsFactors = FALSE)
status_df <- read.csv(STATUS_CSV, stringsAsFactors = FALSE)
gr_spp <- read.csv(GR_SPP_CSV, stringsAsFactors = FALSE)

# Extract epithets
epithets <- trimws(sapply(strsplit(gr_spp$Scientific.Name, " "), 
                          function(x) if(length(x) >= 2) x[2] else x[length(x)]))
epithet_to_full <- setNames(gr_spp$Scientific.Name, epithets)

# Fuzzy match
matched_species <- vapply(amps_df$seq_id, function(id) {
  for (i in seq_along(epithets)) {
    if (length(agrep(tolower(epithets[i]), tolower(id), max.distance = 0.1, ignore.case = TRUE)) > 0) {
      return(epithet_to_full[i])
    }
  }
  NA_character_
}, character(1))

cat(sprintf(">>> MATCHED: %d / %d records (%.1f%%)\n", 
            sum(!is.na(matched_species)), length(matched_species),
            (sum(!is.na(matched_species)) / length(matched_species)) * 100))

parse_seq_id_details <- function(s) {
  # 1. Extract counts
  cnt_match <- regmatches(s, regexec("counts=([0-9]+)", s))[[1]]
  cnt <- if(length(cnt_match) >= 2) as.numeric(cnt_match[2]) else 1
  
  # 2. Extract species: Remove accession/range and 'counts=...'
  # Remove prefix: "NC_022696.1/1-4114 "
  s_clean <- sub("^[^/]+/[^ ]+ ", "", s)
  # Remove suffix: " counts=8"
  species <- sub("counts=[0-9]+$", "", s_clean)
  species <- trimws(species)
  
  return(data.frame(parsed_species = species, haplotype_count = cnt, stringsAsFactors = FALSE))
}


# Prepare status data
status_clean <- status_df %>% 
  select(seq_id, assay, f_score, r_score, f_p_status, r_p_status) %>% # ASSURE assay COLUMN EXISTS
  distinct(seq_id, assay, .keep_all = TRUE)

# Parse sequence info
parsed_info <- do.call(rbind, lapply(amps_df$seq_id, parse_seq_id_details))

# ✅ FIX 1: Join on BOTH seq_id AND assay to prevent broadcasting
merged <- amps_df %>%
  cbind(parsed_info) %>%
  inner_join(status_clean, by = c("seq_id", "assay")) %>%
  mutate(
    total_score = coalesce(f_score, 0) + coalesce(r_score, 0),
    weighted_val = total_score * haplotype_count,
    species = matched_species
  ) %>%
  select(-parsed_species)

# ✅ FIX 2: Correct aggregation logic
agg_data <- suppressWarnings(
  merged %>%
    group_by(assay, species, Taxa, Gene) %>%
    summarise(
      n_records     = n(),
      n_scored      = sum(f_p_status == "scored" & r_p_status == "scored", na.rm = TRUE),
      total_counts  = sum(haplotype_count, na.rm = TRUE),
      sum_weighted  = sum(ifelse(f_p_status == "scored" & r_p_status == "scored", weighted_val, 0)),
      min_v         = min(total_score[f_p_status == "scored" & r_p_status == "scored"], na.rm = TRUE), 
      max_v         = max(total_score[f_p_status == "scored" & r_p_status == "scored"], na.rm = TRUE),
      .groups       = 'drop'
    ) %>%
    mutate(
      avg_score     = ifelse(total_counts > 0, sum_weighted / total_counts, NA),
      min_score     = if_else(is.infinite(min_v) | is.na(min_v), NA_real_, min_v),
      max_score     = if_else(is.infinite(max_v) | is.na(max_v), NA_real_, max_v)
    ) %>%
    select(-min_v, -max_v)
)
# 
# # GRID EXPANSION & JOIN
# all_assays <- as.character(unique(amps_df$assay))
# all_species <- as.vector(gr_spp$Scientific.Name)
# grid <- expand.grid(assay = all_assays, species = all_species, stringsAsFactors = FALSE)
# 
# final_table <- grid %>%
#   left_join(agg_data, by = c("assay", "species")) %>%
#   mutate(
#     n_records   = replace_na(as.numeric(n_records), 0),
#     n_scored    = replace_na(as.numeric(n_scored), 0),
#     total_counts= replace_na(as.numeric(total_counts), 0),
#     min_score   = if_else(is.na(min_score), NA_real_, min_score),
#     avg_score   = if_else(is.na(avg_score), NA_real_, avg_score),
#     max_score   = if_else(is.na(max_score), NA_real_, max_score)
#   ) %>%
#   select(assay, species, Taxa, n_records, n_scored, total_counts, min_score, avg_score, max_score)
# 
# write.csv(final_table, SUMMARY_OUTPUT_CSV, row.names = FALSE)
# cat("\n*** SUCCESS: Summary table written to disk! ***\n")
# 
# 
# # ----------------------------------------------------------------------------
# # SECTION 2: PLOTTING (Step 2: Visualizing the newly created file)
# # ----------------------------------------------------------------------------
# message(">>> STARTING STEP 2: Plotting from saved CSV...")
# 
# # Re-load data fresh from disk to ensure no environment pollution
# plot_input <- read.csv(SUMMARY_OUTPUT_CSV, stringsAsFactors = FALSE)
primers_true_df <- read.csv(PRIMERS_TRUE_CSV, stringsAsFactors = FALSE, check.names = FALSE)

# Map Assay -> Gene
assay_gene_map <- primers_true_df %>%
  mutate(Assay_Base = gsub("(_F|_R)$", "", name)) %>%
  distinct(Assay_Base, gene, .keep_all = TRUE) %>%
  rename(Assay = Assay_Base)

# Prepare plot data
plot_data <- agg_data %>%
  left_join(assay_gene_map, by = c("assay" = "Assay")) %>%
  filter(n_scored > 0 & !is.na(gene))

# Calculate failure counts
failure_counts <- plot_data %>%
  group_by(assay) %>%
  summarise(fail_count = sum(avg_score > 100, na.rm = TRUE), .groups = 'drop')

plot_data <- plot_data %>%
  left_join(failure_counts, by = "assay") %>%
  mutate(fail_color = ifelse(is.na(fail_count) | fail_count == 0, "forestgreen", "firebrick"),
         fail_label = as.character(replace_na(fail_count, 0)))

# Set factor orders for clean X-axis grouping by Gene
plot_data$gene <- factor(plot_data$gene)
plot_data$assay <- factor(plot_data$assay, levels = unique(plot_data$assay[order(plot_data$gene)]))

PENALTY_LIMIT <- if (nrow(plot_data) > 0) max(150L, as.integer(ceiling(max(plot_data$avg_score, na.rm = TRUE))) + 15L) else 150L

p <- ggplot(plot_data, aes(x = assay, y = avg_score)) +
  geom_boxplot(aes(group = assay), outlier.shape = NA, alpha = 0.4, color = "black", fill="grey90") +
  geom_jitter(aes(fill = Taxa), position = position_jitter(width = 0.3, height = 0), alpha = 0.6, size = 3, shape = 21, color="black") +
  new_scale_color() +
  geom_text(aes(y = -5, label = fail_label, color = fail_color), size = 4.5, fontface = "bold") +
  scale_color_manual(values = c("forestgreen" = "forestgreen", "firebrick" = "firebrick"), guide = "none") +
  geom_hline(aes(yintercept = 100), color = "red", lty = "dashed") +
  facet_wrap(~gene, scales = "free_x", space = "free_x") + 
  scale_y_continuous(limits = c(-25, PENALTY_LIMIT + 15), expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2", name = "Taxa") +
  scale_x_discrete(labels = function(x) gsub("([_-])", "\\1\n", x)) +
  labs(y = "Local Database - PrimerMiner Penalty Score (Species-Averaged)", x = "Candidate Metabarcoding Assay") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.8, size = 10), axis.title = element_text(size = 17, face = "bold"), axis.text = element_text(size = 11), strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold", size = 12), panel.grid.major = element_blank(), panel.grid    = element_blank(), legend.position = "right", axis.line = element_line(linewidth = 1))+
  ylim(-5,100)+geom_hline(aes(yintercept =100), color = "red", lty = "dashed")

print(p)
ggsave(PLOT_OUTPUT_PNG, p, width = 18, height = 10, dpi = 300)

cat("\n============================================================================\n")
cat("PIPELINE COMPLETE.\nSummary:", SUMMARY_OUTPUT_CSV, "\nPlot:", PLOT_OUTPUT_PNG, "\n")
cat("============================================================================\n")
