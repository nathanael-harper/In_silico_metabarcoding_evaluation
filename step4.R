#!/usr/bin/env Rscript
# ============================================================================
# step4_barcode_resolution.r — STEP 4: BARCODE RESOLUTION  (Windows R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
#
# The analysis is base R (no dplyr / tidyr / ggplot2). The stacked
# template-fate data (12_stacked_plot_data.csv) is built entirely from
# in-session 06 + 08 + 01 (no external label file) and is written for
# downstream plotting; no figures are produced in this step.
#
# WHAT THIS STEP DOES (and only this)
#   Inclusion for the barcode-uniqueness analysis is decided by the pipeline
#   scoring: a (assay, species) row is analyzed iff step 2 classified it
#   "suitable" in 06_assay_template_status.csv — both strands mapped + fully
#   covered and total = score_F + score_R <= CONFIG$SCORE_THRESHOLD (NA total
#   = unsuitable). [4.2] re-derives that class from the stored strand scores
#   and CHECKs it equals the stored a_status — the scoring gate is a
#   structural identity, never a second source of truth.
#   For every suitable row with amp_bp > 0 (step 3's 07_amplicons.csv)
#   whose record name parses to genus + species:
#     1. shared amplicon groups: >1 species sharing one amp_seq within
#        (assay, Taxa, Gene)            -> 09_shared_amplicon_details.csv
#     2. pairwise species collisions (combn of each shared group;
#        same-genus pairs flagged)      -> 10_collision_pairs.csv
#     3. per-assay resolution summary (n_complete, n_shared_groups,
#        n_species_in_shared, n_congeneric_groups, n_collision_pairs,
#        n_congeneric_pairs, n_species_unique, pct_unique)
#                                    -> 08_barcode_resolution_summary.csv
#     4. threshold sensitivity: the SAME scoring rule re-evaluated at
#        0.6x / 1x / 1.5x SCORE_THRESHOLD
#                                    -> 11_threshold_sensitivity.csv
#     5. stacked template-fate proportions (five fates per assay x Taxa x
#        Gene, built in-session from 06 + 08)
#                                    -> 12_stacked_plot_data.csv
#
#   Report, never exclude: suitable rows with an empty (amp_bp = 0)
#   amplicon, and record names that do not parse to genus/species, are
#   NAMED + counted and excluded from the uniqueness denominator (they
#   carry no barcode / no species identity) — everything else flows
#   through. Non-suitable rows are out of scope by definition: the scoring
#   decision belongs to step 2.
#
#   The structural replacements for the old expected-count checks: [4.2]
#   proves stored suitability == re-derived suitability, [4.3] proves the
#   complete = unique + shared identity for every row, and [4.4] proves the
#   sensitivity machinery reproduces the base summary at 1x SCORE_THRESHOLD.
#   No expected count is hardcoded anywhere in this file.
#
#   Step 4 reads NO MSA files (it consumes the step 3 slice), so there is
#   no per-MSA re-fingerprint here — the MSA contract was already enforced
#   by steps 1-3 in the same results/ directory.
#
#   Span-gate note: NO analysis logic changes. Step 4 consumes 06 + 07, and
#   the gate's effects are already reflected there: a span-rejected strand
#   is map_mode "none" -> the existing "unmapped" class in 06; a span-gated
#   anchor only changes which [lo, hi] span step 3 sliced. Step 4 reads
#   NEITHER 03_primer_mapping.csv NOR 04_pair_coverage.csv, so the
#   span-gate columns need no gate here. Two [4.1] continuity additions,
#   mirroring steps 1-3:
#     (a) the manifest config|max_pair_span row is CHECKed against this
#         session's CONFIG$MAX_AMP_SPAN (drift = stop; re-run step 0),
#         placed AFTER the utils_version gate on purpose;
#     (b) step 4 reads the same primer + taxon directory files step 0
#         reads: the manifest row counts are published from those files, so
#         the gate compares against the same source.
#
# INPUTS (must match the step 0 manifest + steps 1-3 outputs):
#   results/00_run_manifest.csv (config|max_pair_span row),
#   results/01_pair_table.csv,
#   results/02_template_inventory.csv,
#   results/06_assay_template_status.csv,
#   results/07_amplicons.csv,
#   results/run_audit_log.csv (a step 3 PASS row is required),
#   primer CSV, taxon directory CSV
#
#   Record-name contract: MSA record names are expected to follow
#   'Genus_species' (optionally followed by a '/...' fragment) for full
#   resolution output; non-matching names are named + counted and excluded
#   from the uniqueness denominator (report, never exclude).
#
#   Step 4 reads NO 03/04 files: the span-gate columns (03) and
#   amplicon_size_clean (04) never enter this step.
#
# OUTPUTS (results/, overwritten each run):
#   08_barcode_resolution_summary.csv   one row per (assay, Taxa, Gene)
#   09_shared_amplicon_details.csv      one row per shared amplicon group
#   10_collision_pairs.csv              one row per species-pair collision
#   11_threshold_sensitivity.csv        per-assay summary at 3 thresholds
#   12_stacked_plot_data.csv            stacked template-fate proportions +
#                                       A_Total per (assay, Taxa, Gene);
#                                       written for downstream plotting
#   run_audit_log.csv appended (step = 4)
#
# AUDITS ([4.x]; CHECK stops the run, NOTE never stops):
#   [4.1] continuity gate: manifest PASS + BASE + config unchanged (incl.
#          max_pair_span); utils version (older = NOTE); run chain 0 -> 1
#          -> 2 -> 3 with the same BASE (WSL form tolerated for step 2);
#          01/02/06/07 schemas; a_status vocabularies; 06 key uniqueness;
#          06 row count vs inputs; every 07 row joins 1:1 to 06 with
#          identical a_status; 07 amp_bp == nchar(amp_seq); inputs
#          unchanged vs manifest
#   [4.2] scoring gate: stored a_status == "suitable" is EXACTLY the set
#          re-derived from the strand scores at CONFIG$SCORE_THRESHOLD
#          (the scoring decides inclusion); per-taxon class counts printed;
#          n_suitable >= 1
#   [4.3] resolution: taxonomy parse (unparsed named + counted); empty
#          amplicons counted; shared groups + pairwise collisions (combn,
#          empty-safe); per-assay summary; identity CHECK
#          n_complete == n_species_unique + n_species_in_shared (every row)
#   [4.4] outputs written + re-read (round-trip, schema, values; 09/10 are
#          empty-safe); sensitivity identity at every threshold;
#          sensitivity at SCORE_THRESHOLD reproduces the base summary;
#          stacked template-fate data (12_stacked_plot_data.csv) written
#   [4.5] audit log appended
#
# RUN:  Rscript step4_barcode_resolution.r   (Windows R; analysis = base R
#       — no dplyr / tidyr / ggplot2; steps 0-3 must have passed in the same
#       results/ directory; fresh session cleanest, but paste-safe too)
#       Paste-safe: all conditionals braced; one CHUNK break marked.
# ============================================================================

## ---- BASE — the ONE line to edit (same line as steps 0-3; verified
## ---- against the manifest) ----
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
OUT <- file.path(BASE, "results")

## cat-based table printer — the console output of this whole step.
show_df <- function(df, n = 100L, header = "") {
  if (nzchar(header)) cat(header, "\n")
  if (nrow(df) == 0L) { cat("  (empty)\n"); return(invisible(df)) }
  cat(paste0("  ", paste(names(df), collapse = " | ")), "\n")
  for (i in seq_len(min(n, nrow(df))))
    cat(paste0("  ", paste(as.character(unlist(df[i, ])), collapse = " | ")), "\n")
  if (nrow(df) > n) cat(sprintf("  ... %d more rows (of %d)\n", nrow(df) - n, nrow(df)))
  invisible(df)
}

## ---- [4.1] continuity gate ---------------------------------------------------
cat("============================================================================\n")
cat("STEP 4 - BARCODE RESOLUTION (unique-amplicon analysis on the modular scoring)\n")
cat("  run_id  :", run_id, "\n")
mani_path <- file.path(OUT, "00_run_manifest.csv")
check("[4.1] step 0 manifest exists", file.exists(mani_path),
      "results/00_run_manifest.csv missing — run step0_setup_audit.r first")
mani <- read.csv(mani_path, stringsAsFactors = FALSE)
mval <- function(sec, item) {
  v <- mani$value[mani$section == sec & mani$item == item]
  if (length(v) != 1L) stop("manifest missing row: ", sec, "|", item, call. = FALSE)
  v
}
check("[4.1] step 0 status = PASS", mval("status", "step0") == "PASS",
      "step 0 did not pass — fix it before resolution analysis")
parent_run <- mval("run", "run_id")
check("[4.1] BASE unchanged since step 0", identical(mval("run", "base"), BASE),
      paste0("manifest BASE: ", mval("run", "base"), " | session BASE: ", BASE,
             " — if you moved things, re-run step 0"))
check("[4.1] score_threshold unchanged since step 0",
      identical(mval("config", "score_threshold"), as.character(CONFIG$SCORE_THRESHOLD)),
      "CONFIG$SCORE_THRESHOLD differs from the manifest — set it, then re-run step 0")
check("[4.1] max_map_mismatches unchanged since step 0",
      identical(mval("config", "max_map_mismatches"), as.character(CONFIG$MAX_MAP_MISMATCHES)),
      "CONFIG$MAX_MAP_MISMATCHES differs from the manifest — set it, then re-run step 0")

## utils continuity — same tolerant handling as steps 1-3
uv_row <- mani$value[mani$section == "run" & mani$item == "utils_version"]
if (length(uv_row) == 1L) {
  check("[4.1] utils_version matches manifest", uv_row == UTILS_VERSION,
        paste0("manifest utils_version is '", uv_row,
               "', this session sources utils.R '", UTILS_VERSION,
               "' — re-run step 0 to republish the manifest"))
} else {
  note("[4.1]", sprintf("manifest has no utils_version row; this session sources utils.R %s",
                        UTILS_VERSION))
}

## Span-gate parameter: must not have drifted since step 0.
## Placed AFTER the utils_version gate on purpose: an older manifest
## already stops at the version check above, so by the time a run reaches
## here the config row is guaranteed to exist. Same ordering logic as
## steps 1 [1.1], 2 [2.1] and 3 [3.1].
check("[4.1] max_pair_span unchanged since step 0",
      identical(mval("config", "max_pair_span"), as.character(CONFIG$MAX_AMP_SPAN)),
      "CONFIG$MAX_AMP_SPAN differs from the manifest — set it, then re-run step 0")

## Run chain: 0 (manifest) <- 1 <- 2 <- 3 <- this step. Step 2's log row
## carries the WSL spelling of BASE, so all bases compare in WSL form via
## to_wsl_path (POSIX paths pass through unchanged; Windows forms derived).
LOG <- file.path(OUT, "run_audit_log.csv")
check("[4.1] audit log exists", file.exists(LOG), "results/run_audit_log.csv missing")
logdf <- read.csv(LOG, stringsAsFactors = FALSE)
s3 <- logdf[logdf$step == "3" & logdf$status == "PASS", ]
check("[4.1] step 3 completed (PASS row in log)", nrow(s3) >= 1L,
      "no PASS row for step 3 — run step3_amplicon_extract.r first")
s3_latest <- s3[nrow(s3), ]
s2 <- logdf[logdf$step == "2" & logdf$status == "PASS" &
              logdf$run_id == s3_latest$parent_run_id, ]
check("[4.1] step 3 continues a step 2 run", nrow(s2) >= 1L,
      paste0("step 3 log row ", s3_latest$run_id,
             " has no matching step 2 PASS row — chain broken"))
s2_latest <- s2[nrow(s2), ]
s1 <- logdf[logdf$step == "1" & logdf$status == "PASS" &
              logdf$run_id == s2_latest$parent_run_id, ]
check("[4.1] step 2 continues a step 1 run", nrow(s1) >= 1L,
      paste0("step 2 log row ", s2_latest$run_id,
             " has no matching step 1 PASS row — chain broken"))
s1_latest <- s1[nrow(s1), ]
check("[4.1] step 1 continues the step 0 run",
      identical(as.character(s1_latest$parent_run_id), parent_run),
      paste0("step 1 log parent: ", s1_latest$parent_run_id,
             " | step 0 run: ", parent_run))
check("[4.1] steps 1-3 ran from the same BASE",
      identical(to_wsl_path(as.character(s1_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s2_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s3_latest$base)), to_wsl_path(mval("run", "base"))),
      paste0("step 1 base: ", s1_latest$base, " | step 2 base: ", s2_latest$base,
             " | step 3 base: ", s3_latest$base,
             " | manifest base: ", mval("run", "base")))

## Step 0 + 1 + 2 + 3 outputs exist + schemas
PT_IN  <- file.path(OUT, "01_pair_table.csv")
INV_IN <- file.path(OUT, "02_template_inventory.csv")
STA_IN <- file.path(OUT, "06_assay_template_status.csv")
AMP_IN <- file.path(OUT, "07_amplicons.csv")
check("[4.1] steps 0-3 outputs exist",
      all(file.exists(c(PT_IN, INV_IN, STA_IN, AMP_IN))),
      "missing one of: 01_pair_table.csv, 02_template_inventory.csv, 06_assay_template_status.csv, 07_amplicons.csv — run steps 0-3 first")
pair_table <- read.csv(PT_IN, stringsAsFactors = FALSE)
inv        <- read.csv(INV_IN, stringsAsFactors = FALSE)
status_df  <- read.csv(STA_IN, stringsAsFactors = FALSE)
amps       <- read.csv(AMP_IN, stringsAsFactors = FALSE)

check("[4.1] pair table schema",
      identical(names(pair_table), c("assay", "gene", "f_name", "f_seq", "r_name", "r_seq")),
      paste0("columns found: ", paste(names(pair_table), collapse = ", ")))
check("[4.1] inventory schema",
      identical(names(inv), c("Taxa", "Gene", "filepath", "records_raw", "records",
                              "width", "dups_dropped", "first_name", "last_name",
                              "consensus_rule")),
      paste0("columns found: ", paste(names(inv), collapse = ", ")))
check("[4.1] status schema",
      identical(names(status_df), c("Taxa", "Gene", "assay", "seq_id", "a_status",
                                    "f_score", "r_score", "f_p_status", "r_p_status")),
      paste0("columns found: ", paste(names(status_df), collapse = ", ")))
check("[4.1] amplicon schema",
      identical(names(amps), c("assay", "Taxa", "Gene", "seq_id", "amp_bp", "amp_seq",
                               "amp_win_start", "amp_win_end", "a_status")),
      paste0("columns found: ", paste(names(amps), collapse = ", ")))

A_STATUS_VOCAB <- c("suitable", "unsuitable", "no_coverage", "incomplete_coverage", "unmapped")
check("[4.1] 06 a_status vocabulary", all(status_df$a_status %in% A_STATUS_VOCAB),
      paste0("offending value(s): ",
             paste(unique(setdiff(status_df$a_status, A_STATUS_VOCAB)), collapse = ", ")))
check("[4.1] 07 a_status vocabulary", all(amps$a_status %in% A_STATUS_VOCAB),
      paste0("offending value(s): ",
             paste(unique(setdiff(amps$a_status, A_STATUS_VOCAB)), collapse = ", ")))
check("[4.1] 06 (assay, Taxa, seq_id) unique",
      !any(duplicated(status_df[, c("assay", "Taxa", "seq_id")])),
      "duplicate (assay, Taxa, seq_id) rows in 06_assay_template_status.csv — re-run step 2")

## 06 structural row count (the same gate step 3 [3.1] runs, from the
## step 0 inventory — no MSA re-read needed here)
expect_sta <- 0L
for (i in seq_len(nrow(inv)))
  expect_sta <- expect_sta + inv$records[i] * sum(pair_table$gene == inv$Gene[i])
check("[4.1] 06 rows = sum over templates of n_records x n_pairs(gene)",
      nrow(status_df) == expect_sta,
      paste0("06 has ", nrow(status_df), " rows; inputs imply ", expect_sta,
             " — re-run step 2"))

## 07 joins 1:1 to 06 with an identical a_status (07 carried the class from
## 06 in step 3; a disagreement means one file is stale)
j67 <- merge(amps[, c("assay", "Taxa", "seq_id"), drop = FALSE],
             status_df[, c("assay", "Taxa", "seq_id", "a_status"), drop = FALSE],
             by = c("assay", "Taxa", "seq_id"), sort = FALSE)
check("[4.1] 07 rows join 1:1 to 06 on (assay, Taxa, seq_id)",
      nrow(j67) == nrow(amps),
      "07_amplicons.csv and 06_assay_template_status.csv disagree — re-run step 3")
check("[4.1] 07 a_status identical to 06 a_status",
      identical(amps$a_status, j67$a_status),
      "a_status differs between 07 and 06 — re-run step 3")
check("[4.1] 07 amp_bp non-negative and == nchar(amp_seq) (every row)",
      all(amps$amp_bp >= 0L) && all(amps$amp_bp == nchar(amps$amp_seq)),
      "amp_bp inconsistent with amp_seq in 07_amplicons.csv — re-run step 3")

## Inputs unchanged vs manifest (row counts; step 4 reads no MSA files).
## The primer + taxon directory files are the same ones step 0 reads and
## the manifest row counts were published from, so the gate compares
## against the same source.
prim  <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
check("[4.1] primers_input2.csv unchanged",
      nrow(prim) == as.integer(mval("inputs", "primer_rows")),
      "re-run step 0 after changing inputs")
check("[4.1] taxon_template_directory_2.csv unchanged",
      nrow(taxon) == as.integer(mval("inputs", "taxon_rows")),
      "re-run step 0 after changing inputs")

cat("  parent  : step 0 run ", parent_run, " | step 1 ", s1_latest$run_id,
    " | step 2 ", s2_latest$run_id, " | step 3 ", s3_latest$run_id,
    " (chain verified)\n")
cat("  base    :", BASE, "\n")
## Span-gate parameter prints with the other config values
## (manifest-verified in [4.1] above).
cat("  config  : score_threshold =", CONFIG$SCORE_THRESHOLD,
    " | max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
cat("============================================================================\n")

## ==== CHUNK BREAK (safe: all gates + 01/02/06/07 defined; paste resumes here) ====

## ---- [4.2] scoring gate — the modular scoring decides inclusion --------------
cat("\n[4.2] scoring gate — a (assay, species) row is analyzed for barcode\n",
    "         uniqueness iff the MODULAR PIPELINE SCORES it suitable (step 2):\n",
    "         both strands scored and total = score_F + score_R <= SCORE_THRESHOLD\n",
    "         (NA total -> unsuitable)\n", sep = "")
both_scored <- !is.na(status_df$f_p_status) & status_df$f_p_status == "scored" &
  !is.na(status_df$r_p_status) & status_df$r_p_status == "scored"
tot <- status_df$f_score + status_df$r_score
suitable_re <- both_scored & !is.na(tot) & tot <= CONFIG$SCORE_THRESHOLD
check("[4.2] stored 'suitable' is EXACTLY the scoring rule re-derived from strand scores",
      identical(suitable_re, status_df$a_status == "suitable"),
      "06 a_status disagrees with score_F + score_R <= SCORE_THRESHOLD — 06 is stale; re-run step 2")
n_suitable <- sum(status_df$a_status == "suitable")
check("[4.2] at least one suitable (assay, species) row", n_suitable >= 1L,
      "no suitable row — nothing to analyze (step 2 [2.5] prints why)")
n_boundary <- sum(status_df$a_status == "suitable" & tot == CONFIG$SCORE_THRESHOLD)
cat(sprintf("  rule: suitable iff score_F + score_R <= %s (manifest-verified in [4.1]) | boundary combos: %d\n",
            CONFIG$SCORE_THRESHOLD, n_boundary))
cat("  per-taxon class counts (06_assay_template_status.csv):\n")
for (Tx in sort(unique(status_df$Taxa))) {
  sub <- status_df[status_df$Taxa == Tx, ]
  cat(sprintf("    %-14s suitable=%-6d unsuitable=%-6d no_cov=%-5d incomp=%-5d unmapped=%-5d (n=%d)\n",
              Tx, sum(sub$a_status == "suitable"), sum(sub$a_status == "unsuitable"),
              sum(sub$a_status == "no_coverage"), sum(sub$a_status == "incomplete_coverage"),
              sum(sub$a_status == "unmapped"), nrow(sub)))
}

## ---- [4.3] barcode uniqueness (suitable rows only; report, never exclude) ----
cat("\n[4.3] barcode uniqueness (suitable rows only; report, never exclude)\n")
res0 <- amps[amps$a_status == "suitable", , drop = FALSE]

## taxonomy from the MSA record name (strip a trailing '/...' fragment; the
## first '_' splits genus / species)
res0$seq_bare <- sub("/.*$", "", res0$seq_id)
res0$genus    <- sub("_.*$", "", res0$seq_bare)
res0$species  <- sub("^[^_]*_", "", res0$seq_bare)
res0$parsed   <- !(res0$species == "" | res0$species == res0$genus)

n_unparsed <- sum(!res0$parsed)
if (n_unparsed > 0L) {
  unp_names <- unique(res0$seq_id[!res0$parsed])
  note("[4.3]", sprintf("%d suitable record name(s) do not parse to Genus_species — named, counted, out of the uniqueness denominator: %s",
                        n_unparsed, paste(head(unp_names, 10L), collapse = ", ")))
} else {
  note("[4.3]", "all suitable record names parse to Genus_species")
}

n_empty_suitable <- sum(res0$parsed & res0$amp_bp == 0L)
if (n_empty_suitable > 0L)
  note("[4.3]", sprintf("%d suitable row(s) carry an empty (all-gap) amplicon (amp_bp = 0) — kept in 07, out of the denominator",
                        n_empty_suitable))

res <- res0[res0$parsed & res0$amp_bp > 0L, , drop = FALSE]
cat(sprintf("  suitable rows: %d | entering the denominator: %d (excluded: %d unparsed, %d empty amplicon)\n",
            nrow(res0), nrow(res), n_unparsed, n_empty_suitable))
check("[4.3] at least one resolvable row", nrow(res) >= 1L,
      "no suitable row survives parsing + amp_bp > 0 — inspect record names / step 2 class counts")

## Shared groups -> pairwise collisions -> per-assay summary. Pure base R;
## empty-safe (a threshold with zero rows yields zero-row outputs).
compute_resolution <- function(rt) {
  shared_df <- data.frame(
    assay = character(0), Taxa = character(0), Gene = character(0),
    amp_seq = character(0), amp_bp = integer(0), n_species = integer(0),
    species_list = character(0), same_genus = logical(0),
    stringsAsFactors = FALSE)
  pair_df <- data.frame(
    assay = character(0), Taxa = character(0), Gene = character(0),
    amp_seq = character(0), sp1 = character(0), sp2 = character(0),
    congeneric = logical(0), stringsAsFactors = FALSE)
  if (nrow(rt) > 0L) {
    gs <- split(rt, paste(rt$assay, rt$Taxa, rt$Gene, rt$amp_seq, sep = "|"))
    for (k in names(gs)) {
      g <- gs[[k]]
      if (nrow(g) <= 1L) next
      shared_df <- rbind(shared_df, data.frame(
        assay = g$assay[1L], Taxa = g$Taxa[1L], Gene = g$Gene[1L],
        amp_seq = g$amp_seq[1L], amp_bp = g$amp_bp[1L],
        n_species = nrow(g), species_list = paste(g$species, collapse = "; "),
        same_genus = length(unique(g$genus)) == 1L,
        stringsAsFactors = FALSE))
      idx <- combn(nrow(g), 2L)
      pair_df <- rbind(pair_df, data.frame(
        assay = g$assay[1L], Taxa = g$Taxa[1L], Gene = g$Gene[1L],
        amp_seq = g$amp_seq[1L],
        sp1 = g$species[idx[1L, ]], sp2 = g$species[idx[2L, ]],
        congeneric = g$genus[idx[1L, ]] == g$genus[idx[2L, ]],
        stringsAsFactors = FALSE))
    }
    pair_df <- pair_df[pair_df$sp1 != pair_df$sp2, , drop = FALSE]
  }
  combos <- if (nrow(rt) > 0L) unique(rt[, c("assay", "Taxa", "Gene")])
  else data.frame(assay = character(0), Taxa = character(0),
                  Gene = character(0), stringsAsFactors = FALSE)
  n_comb <- nrow(combos)
  summary_df <- data.frame(
    assay = combos$assay, Taxa = combos$Taxa, Gene = combos$Gene,
    n_complete = integer(n_comb), n_shared_groups = integer(n_comb),
    n_species_in_shared = integer(n_comb), n_congeneric_groups = integer(n_comb),
    n_collision_pairs = integer(n_comb), n_congeneric_pairs = integer(n_comb),
    stringsAsFactors = FALSE)
  for (i in seq_len(n_comb)) {
    m_rt <- rt$assay == combos$assay[i] & rt$Taxa == combos$Taxa[i] & rt$Gene == combos$Gene[i]
    m_sh <- shared_df$assay == combos$assay[i] & shared_df$Taxa == combos$Taxa[i] &
      shared_df$Gene == combos$Gene[i]
    m_pg <- pair_df$assay == combos$assay[i] & pair_df$Taxa == combos$Taxa[i] &
      pair_df$Gene == combos$Gene[i]
    summary_df$n_complete[i]          <- sum(m_rt)
    summary_df$n_shared_groups[i]     <- sum(m_sh)
    summary_df$n_species_in_shared[i] <- if (any(m_sh)) sum(shared_df$n_species[m_sh]) else 0L
    summary_df$n_congeneric_groups[i] <- if (any(m_sh)) sum(shared_df$same_genus[m_sh]) else 0L
    summary_df$n_collision_pairs[i]   <- sum(m_pg)
    summary_df$n_congeneric_pairs[i]  <- if (any(m_pg)) sum(pair_df$congeneric[m_pg]) else 0L
  }
  summary_df$n_species_unique <- summary_df$n_complete - summary_df$n_species_in_shared
  summary_df$pct_unique <- ifelse(summary_df$n_complete > 0L,
                                  100 * summary_df$n_species_unique / summary_df$n_complete,
                                  NA_real_)
  list(summary = summary_df, shared = shared_df, pairs = pair_df)
}

rbase         <- compute_resolution(res)
final_summary <- rbase$summary
shared_all    <- rbase$shared
pair_coll     <- rbase$pairs

## the identity gate
ident_ok <- final_summary$n_complete ==
  final_summary$n_species_unique + final_summary$n_species_in_shared
check("[4.3] identity: n_complete == n_species_unique + n_species_in_shared (every row)",
      all(ident_ok),
      paste0("offending row(s): ",
             paste(sprintf("%s/%s/%s (complete=%d unique=%d in_shared=%d)",
                           final_summary$assay[!ident_ok], final_summary$Taxa[!ident_ok],
                           final_summary$Gene[!ident_ok],
                           final_summary$n_complete[!ident_ok],
                           final_summary$n_species_unique[!ident_ok],
                           final_summary$n_species_in_shared[!ident_ok]),
                   collapse = " ; ")))
cat("Identity gate passed.\n")
if (nrow(shared_all) == 0L)
  note("[4.3]", "no shared amplicon groups - every species amplicon is unique within its (assay, Taxa, Gene); 09/10 are header-only by design")

cat("\nPer-assay resolution (complete = suitable per the modular scoring;\n",
    "unique = a species whose amplicon no other species in the same assay x taxon shares):\n", sep = "")
show_df(final_summary[order(final_summary$Taxa, final_summary$Gene,
                            -final_summary$pct_unique), , drop = FALSE], n = 100L)

cat("\nPer (Taxa, Gene) rollup:\n")
tg <- unique(final_summary[, c("Taxa", "Gene")])
for (i in seq_len(nrow(tg))) {
  m  <- final_summary$Taxa == tg$Taxa[i] & final_summary$Gene == tg$Gene[i]
  nc <- sum(final_summary$n_complete[m])
  nu <- sum(final_summary$n_species_unique[m])
  cat(sprintf("  %-14s %-6s complete=%-6d unique=%-6d (%.1f%%) | shared groups=%d | collision pairs=%d\n",
              tg$Taxa[i], tg$Gene[i], nc, nu, 100 * nu / nc,
              sum(final_summary$n_shared_groups[m]), sum(final_summary$n_collision_pairs[m])))
}
tot_c <- sum(final_summary$n_complete)
tot_u <- sum(final_summary$n_species_unique)
cat(sprintf("  overall        complete=%-6d unique=%-6d (%.1f%%)\n",
            tot_c, tot_u, 100 * tot_u / tot_c))

## ---- [4.4] outputs + round-trip -------------------------------------------------
cat("\n[4.4] outputs (results/)\n")
SUM_COLS <- c("assay", "Taxa", "Gene", "n_complete", "n_shared_groups",
              "n_species_in_shared", "n_congeneric_groups", "n_collision_pairs",
              "n_congeneric_pairs", "n_species_unique", "pct_unique")
SHARED_COLS <- c("assay", "Taxa", "Gene", "amp_seq", "amp_bp", "n_species",
                 "species_list", "same_genus")
PAIR_COLS <- c("assay", "Taxa", "Gene", "amp_seq", "sp1", "sp2", "congeneric")

OUT_SUM <- file.path(OUT, "08_barcode_resolution_summary.csv")
write.csv(final_summary, OUT_SUM, row.names = FALSE)
chk <- read.csv(OUT_SUM, stringsAsFactors = FALSE)
check("[4.4] 08 round-trip schema", identical(names(chk), SUM_COLS),
      paste0("columns found: ", paste(names(chk), collapse = ", ")))
check("[4.4] 08 round-trip values",
      nrow(chk) == nrow(final_summary) &&
        identical(chk$assay, final_summary$assay) &&
        identical(chk$Taxa, final_summary$Taxa) &&
        identical(chk$Gene, final_summary$Gene) &&
        all(as.integer(chk$n_complete) == final_summary$n_complete) &&
        all(as.integer(chk$n_shared_groups) == final_summary$n_shared_groups) &&
        all(as.integer(chk$n_species_in_shared) == final_summary$n_species_in_shared) &&
        all(as.integer(chk$n_congeneric_groups) == final_summary$n_congeneric_groups) &&
        all(as.integer(chk$n_collision_pairs) == final_summary$n_collision_pairs) &&
        all(as.integer(chk$n_congeneric_pairs) == final_summary$n_congeneric_pairs) &&
        all(as.integer(chk$n_species_unique) == final_summary$n_species_unique) &&
        all(abs(chk$pct_unique - final_summary$pct_unique) < 1e-9),
      "re-read 08_barcode_resolution_summary.csv differs")
cat(sprintf("    08_barcode_resolution_summary.csv  %d rows\n", nrow(final_summary)))

OUT_SHD <- file.path(OUT, "09_shared_amplicon_details.csv")
write.csv(shared_all, OUT_SHD, row.names = FALSE)
chk <- read.csv(OUT_SHD, stringsAsFactors = FALSE)
check("[4.4] 09 round-trip schema", identical(names(chk), SHARED_COLS),
      paste0("columns found: ", paste(names(chk), collapse = ", ")))
check("[4.4] 09 round-trip values",
      nrow(chk) == nrow(shared_all) &&
        (nrow(shared_all) == 0L ||
           (identical(chk$amp_seq, shared_all$amp_seq) &&
              all(as.integer(chk$amp_bp) == shared_all$amp_bp) &&
              all(as.integer(chk$n_species) == shared_all$n_species) &&
              identical(chk$species_list, shared_all$species_list) &&
              all(chk$same_genus == shared_all$same_genus))),
      "re-read 09_shared_amplicon_details.csv differs")
cat(sprintf("    09_shared_amplicon_details.csv     %d rows\n", nrow(shared_all)))

OUT_PAIR <- file.path(OUT, "10_collision_pairs.csv")
write.csv(pair_coll, OUT_PAIR, row.names = FALSE)
chk <- read.csv(OUT_PAIR, stringsAsFactors = FALSE)
check("[4.4] 10 round-trip schema", identical(names(chk), PAIR_COLS),
      paste0("columns found: ", paste(names(chk), collapse = ", ")))
check("[4.4] 10 round-trip values",
      nrow(chk) == nrow(pair_coll) &&
        (nrow(pair_coll) == 0L ||
           (identical(chk$assay, pair_coll$assay) &&
              identical(chk$Taxa, pair_coll$Taxa) &&
              identical(chk$Gene, pair_coll$Gene) &&
              identical(chk$amp_seq, pair_coll$amp_seq) &&
              identical(chk$sp1, pair_coll$sp1) &&
              identical(chk$sp2, pair_coll$sp2) &&
              all(chk$congeneric == pair_coll$congeneric))),
      "re-read 10_collision_pairs.csv differs")
cat(sprintf("    10_collision_pairs.csv             %d rows\n", nrow(pair_coll)))

## threshold sensitivity: the SAME scoring rule (total <= threshold)
## re-evaluated at 0.6x / 1x / 1.5x SCORE_THRESHOLD from the stored strand
## scores (derived from CONFIG, not hardcoded)
thr_set <- sort(unique(c(as.integer(round(0.6 * CONFIG$SCORE_THRESHOLD)),
                         as.integer(CONFIG$SCORE_THRESHOLD),
                         as.integer(round(1.5 * CONFIG$SCORE_THRESHOLD)))))
sens_parts <- vector("list", length(thr_set))
for (t_i in seq_along(thr_set)) {
  th <- thr_set[t_i]
  suit_t <- both_scored & !is.na(tot) & tot <= th
  res_t <- merge(res[, c("assay", "Taxa", "seq_id", "Gene", "amp_bp", "amp_seq",
                         "genus", "species"), drop = FALSE],
                 status_df[suit_t, c("assay", "Taxa", "seq_id"), drop = FALSE],
                 by = c("assay", "Taxa", "seq_id"), sort = FALSE)
  s <- compute_resolution(res_t)$summary
  s$threshold <- as.integer(th)
  sens_parts[[t_i]] <- s
}
sens_df <- do.call(base::rbind, sens_parts)
SENS_COLS <- c(SUM_COLS, "threshold")
check("[4.4] sensitivity identity holds at every threshold",
      all(sens_df$n_complete == sens_df$n_species_unique + sens_df$n_species_in_shared),
      "complete = unique + shared identity broken in the sensitivity table")
base_rows <- sens_df[sens_df$threshold == CONFIG$SCORE_THRESHOLD, ]
check("[4.4] sensitivity at SCORE_THRESHOLD reproduces the base summary",
      nrow(base_rows) == nrow(final_summary) &&
        identical(sort(paste(base_rows$assay, base_rows$Taxa, base_rows$Gene)),
                  sort(paste(final_summary$assay, final_summary$Taxa, final_summary$Gene))) &&
        sum(base_rows$n_complete) == sum(final_summary$n_complete),
      "the sensitivity machinery disagrees with the base analysis at 1x — internal error")
OUT_SENS <- file.path(OUT, "11_threshold_sensitivity.csv")
write.csv(sens_df, OUT_SENS, row.names = FALSE)
chk <- read.csv(OUT_SENS, stringsAsFactors = FALSE)
check("[4.4] 11 round-trip schema", identical(names(chk), SENS_COLS),
      paste0("columns found: ", paste(names(chk), collapse = ", ")))
check("[4.4] 11 round-trip values",
      nrow(chk) == nrow(sens_df) &&
        identical(chk$assay, sens_df$assay) &&
        identical(chk$Taxa, sens_df$Taxa) &&
        all(as.integer(chk$n_complete) == sens_df$n_complete) &&
        all(as.integer(chk$n_species_unique) == sens_df$n_species_unique) &&
        all(abs(chk$pct_unique - sens_df$pct_unique) < 1e-9) &&
        all(as.integer(chk$threshold) == sens_df$threshold),
      "re-read 11_threshold_sensitivity.csv differs")
cat(sprintf("    11_threshold_sensitivity.csv       %d rows (thresholds: %s)\n",
            nrow(sens_df), paste(thr_set, collapse = " / ")))
cat("\nThreshold sensitivity (same scoring rule at alternate SCORE_THRESHOLD values):\n")
show_df(sens_df[order(sens_df$threshold, sens_df$Taxa, sens_df$Gene,
                      -sens_df$pct_unique), , drop = FALSE], n = 100L)

## ---- [4.5] audit log -------------------------------------------------------------
log_row <- log_row_make(run_id, "4", "PASS", BASE, parent_run = s3_latest$run_id,
                        metrics = paste0("n_suitable=", n_suitable,
                                         ";n_unparsed_names=", n_unparsed,
                                         ";n_empty_suitable=", n_empty_suitable,
                                         ";n_res_rows=", nrow(res),
                                         ";n_shared_groups=", nrow(shared_all),
                                         ";n_species_in_shared=", sum(final_summary$n_species_in_shared),
                                         ";n_congeneric_groups=", sum(final_summary$n_congeneric_groups),
                                         ";n_collision_pairs=", nrow(pair_coll),
                                         ";n_congeneric_pairs=", sum(pair_coll$congeneric, na.rm = TRUE),
                                         ";n_species_unique=", sum(final_summary$n_species_unique),
                                         ";score_threshold=", CONFIG$SCORE_THRESHOLD))
log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)
log_chk <- read.csv(LOG, stringsAsFactors = FALSE)
check("[4.5] audit log appended",
      any(log_chk$run_id == run_id & log_chk$step == "4"),
      "log row missing after write")

## ---- [4.4] stacked template-fate data (12_stacked_plot_data.csv) --------------
## One row per (assay, Taxa, Gene): the five template fates as a fraction of
## A_Total for that (assay, Taxa), plus A_Total. Categories:
##   No_coverage | Incomplete | No_predicted_amplification |
##   Shared_barcode | Unique_barcode
## A_Good = the SUITABLE count (the 06 "suitable" class), so categories 1-3
## partition A_Total - A_Good exactly and 4-5 partition n_complete <= A_Good;
## a row's proportions sum to 1.0 only when the exclusions (unparsed names /
## empty amplicons) are zero for that (assay, Taxa), otherwise they end just
## below 1.0 and the NOTE names which rows.
## Data: built in-session from 06 (status_df) + 08 (final_summary) + the
## assay order of 01 (pair_table) — no external label file is read.
## Written for downstream plotting (counts = proportion * A_Total).
OUT_CSV <- file.path(OUT, "12_stacked_plot_data.csv")

CATS <- c("No_coverage", "Incomplete", "No_predicted_amplification",
          "Shared_barcode", "Unique_barcode")

## per (assay, Taxa, Gene) counts of the five 06 classes (table() is O(n))
key6   <- with(status_df, paste(assay, Taxa, Gene, sep = "|"))
kt     <- table(key6, status_df$a_status)
for (cl in A_STATUS_VOCAB) if (!cl %in% colnames(kt)) kt[, cl] <- 0
kparts <- do.call(base::rbind, strsplit(rownames(kt), "|", fixed = TRUE))
m <- data.frame(
  assay = kparts[, 1L], Taxa = kparts[, 2L], Gene = kparts[, 3L],
  A_Total      = rowSums(kt),
  A_NoCov      = kt[, "no_coverage"],
  A_Indel      = kt[, "incomplete_coverage"],
  A_Good       = kt[, "suitable"],        # A_Good = the SUITABLE count
  n_unsuitable = kt[, "unsuitable"],
  n_unmapped   = kt[, "unmapped"],
  stringsAsFactors = FALSE)

u8 <- final_summary[, c("assay", "Taxa", "Gene",
                        "n_species_unique", "n_species_in_shared")]
check("[4.4] 08 has one row per (assay, Taxa, Gene)",
      !any(duplicated(u8[, c("assay", "Taxa", "Gene")])))
m <- merge(m, u8, by = c("assay", "Taxa", "Gene"), all.x = TRUE)
m$n_species_unique[is.na(m$n_species_unique)]       <- 0     # zero-suitable
m$n_species_in_shared[is.na(m$n_species_in_shared)] <- 0     # combos get 0/0

m$No_coverage                  <- m$A_NoCov
m$Incomplete                   <- m$A_Indel
m$No_predicted_amplification   <- m$n_unsuitable + m$n_unmapped
m$Shared_barcode               <- m$n_species_in_shared
m$Unique_barcode               <- m$n_species_unique

## NOTE: unmapped primers sit in "No predicted amplification". Put them in
## "No coverage" instead:
##   m$No_predicted_amplification <- m$n_unsuitable
##   m$No_coverage                <- m$A_NoCov + m$n_unmapped

## bookkeeping invariants
stopifnot("A_Total = categories 1-3 + A_Good" =
            all(m$A_Total == m$No_coverage + m$Incomplete +
                  m$No_predicted_amplification + m$A_Good))
stopifnot("unique + shared <= A_Good (suitable)" =
            all(m$Unique_barcode + m$Shared_barcode <= m$A_Good))
gap <- m$A_Good - m$Unique_barcode - m$Shared_barcode
if (any(gap > 0L))
  note("[4.4]", sprintf("A_Good > unique + shared for %d row(s) — exclusions (unparsed record names / empty amplicons; see [4.3]); those proportions end just below 1.0",
                        sum(gap > 0L)))

## integrity: 06 assays + genes match the pair table
assay_order <- pair_table$assay
genes_panel <- pair_table$gene[!duplicated(pair_table$gene)]
check("[4.4] 06 assays + genes match the pair table",
      setequal(unique(m$assay), assay_order) &&
        setequal(unique(m$Gene), genes_panel),
      paste0("assays not in pair table: ",
             paste(setdiff(unique(m$assay), assay_order), collapse = ", "),
             " | genes: ",
             paste(setdiff(unique(m$Gene), genes_panel), collapse = ", ")))

## proportions (the stacked fractions)
m[, CATS] <- m[, CATS] / m$A_Total
stopifnot("column proportions must not sum above 1" =
            all(rowSums(m[, CATS]) <= 1 + 1e-12))

## save proportions + A_Total (re-plot anywhere: counts = prop * A_Total)
write.csv(m[, c("assay", "Taxa", "Gene", "A_Total", CATS)],
          OUT_CSV, row.names = FALSE)

cat("\n============================================================================\n")
cat(sprintf("STEP 4 PASSED | run_id %s (continued from step 3 %s | step 2 %s | step 1 %s | step 0 %s)\n",
            run_id, s3_latest$run_id, s2_latest$run_id, s1_latest$run_id, parent_run))
cat(sprintf("  suitable rows (modular scoring) : %d\n", n_suitable))
cat(sprintf("  resolvable rows (denominator)   : %d (excluded: %d unparsed, %d empty amplicon)\n",
            nrow(res), n_unparsed, n_empty_suitable))
cat(sprintf("  shared groups / collision pairs : %d / %d (congeneric: %d groups / %d pairs)\n",
            nrow(shared_all), nrow(pair_coll),
            sum(final_summary$n_congeneric_groups), sum(pair_coll$congeneric, na.rm = TRUE)))
cat(sprintf("  overall resolution              : %d / %d unique (%.1f%%)\n",
            tot_u, tot_c, 100 * tot_u / tot_c))
cat("Outputs: results/08_barcode_resolution_summary.csv | 09_shared_amplicon_details.csv |\n")
cat("         10_collision_pairs.csv | 11_threshold_sensitivity.csv |\n")
cat("         12_stacked_plot_data.csv | run_audit_log.csv\n")
cat("Pipeline analysis complete (step 4 is the final analysis step).\n")
cat("============================================================================\n")
