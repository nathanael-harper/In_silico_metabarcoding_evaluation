#!/usr/bin/env Rscript
# ============================================================================
# step5_hybrid_labels.r — STEP 5: HYBRID-PLOT LABELS  (Windows R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
#
# Per-assay axis labels (Total / NoCov / Indel / Good per group) computed
# from the pipeline artifacts. Pure base R.
#
# WHAT THIS STEP DOES (and only this)
#   For every (assay, Taxa) combination in 06_assay_template_status.csv
#   (the step 2 per-record class table that steps 3 + 4 carry forward):
#     1. class partition:
#            A_Total = the MSA record count for (Gene, Taxa)
#                      = nrow(06[assay, Taxa]) (CHECKed against 02)
#            A_NoCov = count(a_status == "no_coverage")
#            A_Indel = count(a_status == "incomplete_coverage")
#                      (a de-gapped window of the wrong length)
#            A_Good  = count(a_status == "suitable") — the modular scoring:
#                      both strands scored and total = score_F + score_R
#                      <= SCORE_THRESHOLD (NA total = unsuitable).
#            A_Bad   = A_Total - A_Good (the exact complement: unsuitable
#                      + no_coverage + incomplete_coverage + unmapped)
#     2. the assay's Gene comes from 01_pair_table.csv (assay -> gene),
#     3. per (assay, Gene) it builds the AxisLabel string
#            "Assay\n<group> total/nocov/indel/Ggood\n<group> total/nocov/indel/Ggood"
#        (one line per group; a group line is omitted when that group is
#        absent), so the hybrid plot can join on Assay.
#
#   Input contract ("accepts the steps 3/4 inputs"): this step runs at the
#   end of the chain, AFTER steps 3 + 4, and reads their outputs:
#     - 07_amplicons.csv (step 3) and 08_barcode_resolution_summary.csv
#       (step 4) are read and cross-checked against the class table
#       (07 a_status == 06 a_status; 08 n_complete == the suitable +
#       parsable + nonempty count re-derived from 07) — a stale 07 or 08
#       stops the run;
#     - the per-record class counts themselves come from
#       06_assay_template_status.csv — the table 07 and 08 are built from
#       (07 joins a_status from it; 08's denominator is defined by it).
#
#   Report, never exclude: this step computes + writes labels only; it
#   excludes nothing.
#
#   Consensus-rule note: same posture as steps 3/4 — the per-template
#   consensus rule never enters any count here (the class partition is
#   rule-independent); the rule distribution is printed for provenance
#   only.
#
#   Span-gate note: NO label logic changes. Step 5 consumes 06 + 07 + 08,
#   and the gate's effects are already reflected there: a span-rejected
#   strand is map_mode "none" -> the existing "unmapped" class in 06 (it
#   lands in A_Bad through the 5-class partition); a span-gated anchor
#   only changes which [lo, hi] span step 3 sliced (07). Step 5 reads
#   NEITHER 03_primer_mapping.csv NOR 04_pair_coverage.csv, so the
#   span-gate columns need no gate here (same posture as step 4). Two
#   [5.1] continuity additions, mirroring steps 1-4:
#     (a) the manifest config|max_pair_span row is CHECKed against this
#         session's CONFIG$MAX_AMP_SPAN (drift = stop; re-run step 0),
#         placed AFTER the utils_version gate on purpose;
#     (b) step 5 reads the same primer + taxon directory files step 0
#         reads: the manifest row counts are published from those files,
#         so the gate compares against the same source.
#
# INPUTS (must match the step 0 manifest + steps 1-4 outputs):
#   results/00_run_manifest.csv (config|max_pair_span row),
#   results/01_pair_table.csv,
#   results/02_template_inventory.csv,
#   results/06_assay_template_status.csv,
#   results/07_amplicons.csv                        (step 3 output)
#   results/08_barcode_resolution_summary.csv       (step 4 output)
#   results/run_audit_log.csv (a step 4 PASS row is required)
#   primer CSV, taxon directory CSV
#
#   Record-name contract: same 'Genus_species' convention as step 4 — the
#   08 n_complete cross-check re-derives with step 4's exact parsing.
#
#   Step 5 reads NO 03/04 files: the span-gate columns (03) and
#   amplicon_size_clean (04) never enter this step (same as step 4).
#
# OUTPUTS (results/, overwritten each run):
#   13_assay_axis_labels.csv            one row per (assay, Taxa): assay,
#                                       Gene, Taxa, A_Total, A_NoCov,
#                                       A_Indel, A_Good, A_Bad,
#                                       n_unsuitable, n_unmapped, label_line
#   14_assay_axis_label_strings.csv     one row per (assay, Gene): Assay,
#                                       Gene, AxisLabel (join the hybrid
#                                       plot on Assay)
#   run_audit_log.csv appended (step = 5)
#
# AUDITS ([5.x]; CHECK stops the run, NOTE never stops):
#   [5.1] continuity gate: manifest PASS + BASE + config unchanged (incl.
#          max_pair_span); utils version (older = NOTE); run chain 0 -> 1
#          -> 2 -> 3 -> 4 with the same BASE (WSL form tolerated for
#          step 2); 01/02/06/07/08 schemas; a_status vocabularies; 06 key
#          uniqueness; 06 row count vs the 02 inventory; every 07 row
#          joins 1:1 to 06 with identical a_status; inputs unchanged vs
#          manifest
#   [5.2] label computation: one gene per assay (pair table); per
#          (assay, Taxa): A_Total == the inventory record count, the
#          5-class partition sums to A_Total, A_Good == the suitability
#          re-derived from the stored strand scores (the identity step 4
#          [4.2] CHECKs); 08 n_complete == the re-derivation from 07; the
#          combined AxisLabel reproduces the per-group lines exactly
#   [5.3] outputs written + re-read (round-trip, schema, values)
#   [5.4] audit log appended
#
# RUN:  Rscript step5_hybrid_labels.r   (Windows R, base R only — no
#       dplyr / tidyr; steps 0-4 must have passed in the same results/
#       directory). Only run if the hybrid figure ships.
#       Paste-safe: all conditionals braced; one CHUNK break marked.
# ============================================================================

## ---- BASE — the ONE line to edit (same line as steps 0-4; verified
## ---- against the manifest) ----
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
OUT <- file.path(BASE, "results")

## cat-based table printer (base data frames only)
show_df <- function(df, n = 100L, header = "") {
  if (nzchar(header)) cat(header, "\n")
  if (nrow(df) == 0L) { cat("  (empty)\n"); return(invisible(df)) }
  cat(paste0("  ", paste(names(df), collapse = " | ")), "\n")
  for (i in seq_len(min(n, nrow(df))))
    cat(paste0("  ", paste(as.character(unlist(df[i, ])), collapse = " | ")), "\n")
  if (nrow(df) > n) cat(sprintf("  ... %d more rows (of %d)\n", nrow(df) - n, nrow(df)))
  invisible(df)
}

## ---- [5.1] continuity gate ---------------------------------------------------
cat("============================================================================\n")
cat("STEP 5 - HYBRID-PLOT LABELS (labels from the modular pipeline artifacts)\n")
cat("  run_id  :", run_id, "\n")
mani_path <- file.path(OUT, "00_run_manifest.csv")
check("[5.1] step 0 manifest exists", file.exists(mani_path),
      "results/00_run_manifest.csv missing — run step0_setup_audit.r first")
mani <- read.csv(mani_path, stringsAsFactors = FALSE)
mval <- function(sec, item) {
  v <- mani$value[mani$section == sec & mani$item == item]
  if (length(v) != 1L) stop("manifest missing row: ", sec, "|", item, call. = FALSE)
  v
}
check("[5.1] step 0 status = PASS", mval("status", "step0") == "PASS",
      "step 0 did not pass — fix it before labels")
parent_run <- mval("run", "run_id")
check("[5.1] BASE unchanged since step 0", identical(mval("run", "base"), BASE),
      paste0("manifest BASE: ", mval("run", "base"), " | session BASE: ", BASE,
             " — if you moved things, re-run step 0"))
check("[5.1] score_threshold unchanged since step 0",
      identical(mval("config", "score_threshold"), as.character(CONFIG$SCORE_THRESHOLD)),
      "CONFIG$SCORE_THRESHOLD differs from the manifest — set it, then re-run step 0")
check("[5.1] max_map_mismatches unchanged since step 0",
      identical(mval("config", "max_map_mismatches"), as.character(CONFIG$MAX_MAP_MISMATCHES)),
      "CONFIG$MAX_MAP_MISMATCHES differs from the manifest — set it, then re-run step 0")

## utils continuity — same tolerant handling as steps 1-4
uv_row <- mani$value[mani$section == "run" & mani$item == "utils_version"]
if (length(uv_row) == 1L) {
  check("[5.1] utils_version matches manifest", uv_row == UTILS_VERSION,
        paste0("manifest utils_version is '", uv_row,
               "', this session sources utils.R '", UTILS_VERSION,
               "' — re-run step 0 to republish the manifest"))
} else {
  note("[5.1]", sprintf("manifest has no utils_version row; this session sources utils.R %s",
                        UTILS_VERSION))
}

## Span-gate parameter: must not have drifted since step 0.
## Placed AFTER the utils_version gate on purpose: an older manifest
## already stops at the version check above, so by the time a run reaches
## here the config row is guaranteed to exist. Same ordering logic as
## steps 1 [1.1], 3 [3.1] and 4 [4.1].
check("[5.1] max_pair_span unchanged since step 0",
      identical(mval("config", "max_pair_span"), as.character(CONFIG$MAX_AMP_SPAN)),
      "CONFIG$MAX_AMP_SPAN differs from the manifest — set it, then re-run step 0")

## Run chain: 0 (manifest) <- 1 <- 2 <- 3 <- 4 <- this step. Step 2's log row
## carries the WSL spelling of BASE, so all bases compare in WSL form via
## to_wsl_path (POSIX paths pass through unchanged; Windows forms derived).
LOG <- file.path(OUT, "run_audit_log.csv")
check("[5.1] audit log exists", file.exists(LOG), "results/run_audit_log.csv missing")
logdf <- read.csv(LOG, stringsAsFactors = FALSE)
s4 <- logdf[logdf$step == "4" & logdf$status == "PASS", ]
check("[5.1] step 4 completed (PASS row in log)", nrow(s4) >= 1L,
      "no PASS row for step 4 — run step4_barcode_resolution.r first")
s4_latest <- s4[nrow(s4), ]
s3 <- logdf[logdf$step == "3" & logdf$status == "PASS" &
              logdf$run_id == s4_latest$parent_run_id, ]
check("[5.1] step 4 continues a step 3 run", nrow(s3) >= 1L,
      paste0("step 4 log row ", s4_latest$run_id,
             " has no matching step 3 PASS row — chain broken"))
s3_latest <- s3[nrow(s3), ]
s2 <- logdf[logdf$step == "2" & logdf$status == "PASS" &
              logdf$run_id == s3_latest$parent_run_id, ]
check("[5.1] step 3 continues a step 2 run", nrow(s2) >= 1L,
      paste0("step 3 log row ", s3_latest$run_id,
             " has no matching step 2 PASS row — chain broken"))
s2_latest <- s2[nrow(s2), ]
s1 <- logdf[logdf$step == "1" & logdf$status == "PASS" &
              logdf$run_id == s2_latest$parent_run_id, ]
check("[5.1] step 2 continues a step 1 run", nrow(s1) >= 1L,
      paste0("step 2 log row ", s2_latest$run_id,
             " has no matching step 1 PASS row — chain broken"))
s1_latest <- s1[nrow(s1), ]
check("[5.1] step 1 continues the step 0 run",
      identical(as.character(s1_latest$parent_run_id), parent_run),
      paste0("step 1 log parent: ", s1_latest$parent_run_id,
             " | step 0 run: ", parent_run))
check("[5.1] steps 1-4 ran from the same BASE",
      identical(to_wsl_path(as.character(s1_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s2_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s3_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s4_latest$base)), to_wsl_path(mval("run", "base"))),
      paste0("step bases: ", s1_latest$base, " | ", s2_latest$base, " | ",
             s3_latest$base, " | ", s4_latest$base,
             " | manifest base: ", mval("run", "base")))

## Step 0-4 outputs exist + schemas
PT_IN  <- file.path(OUT, "01_pair_table.csv")
INV_IN <- file.path(OUT, "02_template_inventory.csv")
STA_IN <- file.path(OUT, "06_assay_template_status.csv")
AMP_IN <- file.path(OUT, "07_amplicons.csv")
SUM_IN <- file.path(OUT, "08_barcode_resolution_summary.csv")
check("[5.1] steps 0-4 outputs exist",
      all(file.exists(c(PT_IN, INV_IN, STA_IN, AMP_IN, SUM_IN))),
      "missing one of: 01_pair_table.csv, 02_template_inventory.csv, 06_assay_template_status.csv, 07_amplicons.csv, 08_barcode_resolution_summary.csv — run steps 0-4 first")
pair_table <- read.csv(PT_IN, stringsAsFactors = FALSE)
inv        <- read.csv(INV_IN, stringsAsFactors = FALSE)
status_df  <- read.csv(STA_IN, stringsAsFactors = FALSE)
amps       <- read.csv(AMP_IN, stringsAsFactors = FALSE)
res_sum    <- read.csv(SUM_IN, stringsAsFactors = FALSE)

check("[5.1] pair table schema",
      identical(names(pair_table), c("assay", "gene", "f_name", "f_seq", "r_name", "r_seq")),
      paste0("columns found: ", paste(names(pair_table), collapse = ", ")))
check("[5.1] inventory schema",
      identical(names(inv), c("Taxa", "Gene", "filepath", "records_raw", "records",
                              "width", "dups_dropped", "first_name", "last_name",
                              "consensus_rule")),
      paste0("columns found: ", paste(names(inv), collapse = ", ")))
check("[5.1] status schema",
      identical(names(status_df), c("Taxa", "Gene", "assay", "seq_id", "a_status",
                                    "f_score", "r_score", "f_p_status", "r_p_status")),
      paste0("columns found: ", paste(names(status_df), collapse = ", ")))
check("[5.1] amplicon schema",
      identical(names(amps), c("assay", "Taxa", "Gene", "seq_id", "amp_bp", "amp_seq",
                               "amp_win_start", "amp_win_end", "a_status")),
      paste0("columns found: ", paste(names(amps), collapse = ", ")))
check("[5.1] resolution summary schema",
      identical(names(res_sum), c("assay", "Taxa", "Gene", "n_complete", "n_shared_groups",
                                  "n_species_in_shared", "n_congeneric_groups",
                                  "n_collision_pairs", "n_congeneric_pairs",
                                  "n_species_unique", "pct_unique")),
      paste0("columns found: ", paste(names(res_sum), collapse = ", ")))

A_STATUS_VOCAB <- c("suitable", "unsuitable", "no_coverage", "incomplete_coverage", "unmapped")
check("[5.1] 06 a_status vocabulary", all(status_df$a_status %in% A_STATUS_VOCAB),
      paste0("offending value(s): ",
             paste(unique(setdiff(status_df$a_status, A_STATUS_VOCAB)), collapse = ", ")))
check("[5.1] 07 a_status vocabulary", all(amps$a_status %in% A_STATUS_VOCAB),
      paste0("offending value(s): ",
             paste(unique(setdiff(amps$a_status, A_STATUS_VOCAB)), collapse = ", ")))
check("[5.1] 06 (assay, Taxa, seq_id) unique",
      !any(duplicated(status_df[, c("assay", "Taxa", "seq_id")])),
      "duplicate (assay, Taxa, seq_id) rows in 06_assay_template_status.csv — re-run step 2")

## 06 structural row count (the same gate steps 3/4 run, from the 02 inventory)
expect_sta <- 0L
for (i in seq_len(nrow(inv)))
  expect_sta <- expect_sta + inv$records[i] * sum(pair_table$gene == inv$Gene[i])
check("[5.1] 06 rows = sum over templates of n_records x n_pairs(gene)",
      nrow(status_df) == expect_sta,
      paste0("06 has ", nrow(status_df), " rows; inputs imply ", expect_sta,
             " — re-run step 2"))

## 07 (step 3 output) joins 1:1 to 06 with identical a_status — a stale 07
## cannot survive here (key uniqueness in 06 was CHECKed above, so match
## finds at most one row per 07 key)
key7 <- paste(amps$assay, amps$Taxa, amps$seq_id)
key6 <- paste(status_df$assay, status_df$Taxa, status_df$seq_id)
k67  <- match(key7, key6)
check("[5.1] 07 rows join 1:1 to 06 on (assay, Taxa, seq_id) with identical a_status",
      all(!is.na(k67)) && identical(amps$a_status, status_df$a_status[k67]),
      "07_amplicons.csv and 06_assay_template_status.csv disagree — re-run step 3")

## Inputs unchanged vs manifest (row counts; step 5 reads no MSA files).
## The primer + taxon directory files are the same ones step 0 reads and
## the manifest row counts were published from, so the gate compares
## against the same source.
prim  <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
check("[5.1] primers_input2.csv unchanged",
      nrow(prim) == as.integer(mval("inputs", "primer_rows")),
      "re-run step 0 after changing inputs")
check("[5.1] taxon_template_directory_2.csv unchanged",
      nrow(taxon) == as.integer(mval("inputs", "taxon_rows")),
      "re-run step 0 after changing inputs")

cat("  parent  : step 0 run ", parent_run, " | step 1 ", s1_latest$run_id,
    " | step 2 ", s2_latest$run_id, " | step 3 ", s3_latest$run_id,
    " | step 4 ", s4_latest$run_id, " (chain verified)\n")
cat("  base    :", BASE, "\n")
## Span-gate parameter prints with the other config values
## (manifest-verified in [5.1] above).
cat("  config  : score_threshold =", CONFIG$SCORE_THRESHOLD,
    " | max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
rule_rows <- mani[mani$section == "consensus_rule", ]

if (nrow(rule_rows) > 0L)
  cat(sprintf("  consensus rules : %d gaps_can_win | %d gaps_never_win (context only — never enter the labels)\n",
              sum(rule_rows$value == "gaps_can_win"),
              sum(rule_rows$value == "gaps_never_win")))
cat("============================================================================\n")

## ==== CHUNK BREAK (safe: all gates + 01/02/06/07/08 defined; paste resumes here) ====

## ---- [5.2] label computation --------------------------------------------------
cat("\n[5.2] per (assay, Taxa) class partition (06 a_status)\n")
cat("  convention: A_Good = 'suitable' (both strands scored, total = f_score + r_score <= ",
    CONFIG$SCORE_THRESHOLD, "; NA total = unsuitable)\n",
    sep = "")
cat("  legend: per group — Total / NoCoverage / Indel(incomplete) / G(Good)\n")

## one gene per assay; 01_pair_table.csv is the Gene source
pair_u <- unique(pair_table[, c("assay", "gene"), drop = FALSE])
check("[5.2] one gene per assay (01_pair_table.csv)",
      nrow(pair_u) == length(unique(pair_u$assay)),
      paste0("assay(s) with more than one gene: ",
             paste(unique(pair_u$assay[duplicated(pair_u$assay)]), collapse = ", ")))
gene_of <- setNames(pair_u$gene, pair_u$assay)
check("[5.2] every 06 assay is in the pair table",
      all(unique(status_df$assay) %in% names(gene_of)),
      paste0("offending assay(s): ",
             paste(setdiff(unique(status_df$assay), names(gene_of)), collapse = ", ")))

## groups: every Taxa in 06 must be an inventory taxon (and vice versa)
check("[5.2] 06 Taxa set == inventory Taxa set",
      setequal(sort(unique(status_df$Taxa)), sort(unique(inv$Taxa))),
      paste0("06 Taxa: ", paste(sort(unique(status_df$Taxa)), collapse = ", "),
             " | inventory Taxa: ", paste(sort(unique(inv$Taxa)), collapse = ", ")))

combos <- unique(status_df[, c("assay", "Taxa"), drop = FALSE])
lab_parts <- vector("list", nrow(combos))
for (i in seq_len(nrow(combos))) {
  a   <- combos$assay[i]; T <- combos$Taxa[i]
  sub <- status_df[status_df$assay == a & status_df$Taxa == T, , drop = FALSE]
  n_rec <- nrow(sub)
  G <- gene_of[[a]]
  n_inv <- inv$records[inv$Taxa == T & inv$Gene == G]
  check(sprintf("[5.2] %s x %s: A_Total == inventory record count", a, T),
        length(n_inv) == 1L && n_inv[1L] == n_rec,
        paste0("06 has ", n_rec, " row(s); inventory has ",
               if (length(n_inv) == 1L) n_inv[1L] else "no (Taxa, Gene) template row"))
  n_good   <- sum(sub$a_status == "suitable")
  n_uns    <- sum(sub$a_status == "unsuitable")
  n_nocov  <- sum(sub$a_status == "no_coverage")
  n_incomp <- sum(sub$a_status == "incomplete_coverage")
  n_unmap  <- sum(sub$a_status == "unmapped")
  check(sprintf("[5.2] %s x %s: 5-class partition sums to A_Total", a, T),
        n_good + n_uns + n_nocov + n_incomp + n_unmap == n_rec,
        "class counts do not partition A_Total — a_status outside the vocabulary?")
  ## A_Good == the scoring rule re-derived from the stored strand scores
  ## (the same identity step 4 [4.2] CHECKs — re-asserted here because the
  ## labels inherit the inclusion decision)
  both_scored <- !is.na(sub$f_p_status) & sub$f_p_status == "scored" &
    !is.na(sub$r_p_status) & sub$r_p_status == "scored"
  tot <- sub$f_score + sub$r_score
  good_re <- both_scored & !is.na(tot) & tot <= CONFIG$SCORE_THRESHOLD
  check(sprintf("[5.2] %s x %s: A_Good == suitability re-derived from strand scores", a, T),
        identical(good_re, sub$a_status == "suitable"),
        "06 a_status disagrees with score_F + score_R <= SCORE_THRESHOLD — 06 is stale; re-run step 2")
  pref <- substr(T, 1L, 2L)
  lab_parts[[i]] <- data.frame(
    assay = a, Gene = G, Taxa = T,
    A_Total = n_rec, A_NoCov = n_nocov, A_Indel = n_incomp,
    A_Good = n_good, A_Bad = n_rec - n_good,
    n_unsuitable = n_uns, n_unmapped = n_unmap,
    label_line = sprintf("%s %d/%d/%d/G%d", pref, n_rec, n_nocov, n_incomp, n_good),
    stringsAsFactors = FALSE)
}
lab <- do.call(rbind, lab_parts)
lab <- lab[order(lab$assay, lab$Taxa), , drop = FALSE]

## ---- [5.2] cross-check against the step 4 artifact ----------------------------
## 08's n_complete (per assay x Taxa x Gene) = the suitable rows that parse to
## Genus_species AND have amp_bp > 0 (step 4 [4.3]). Re-derive from 07 (the
## step 3 artifact) with step 4's exact name convention; a stale 08 stops.
res07 <- amps[amps$a_status == "suitable", , drop = FALSE]
res07$seq_bare <- sub("/.*$", "", res07$seq_id)
res07$genus    <- sub("_.*$", "", res07$seq_bare)
res07$species  <- sub("^[^_]*_", "", res07$seq_bare)
res07$parsed   <- !(res07$species == "" | res07$species == res07$genus)
denom   <- res07[res07$parsed & res07$amp_bp > 0L, , drop = FALSE]
denom_u <- unique(denom[, c("assay", "Taxa", "Gene"), drop = FALSE])
re_counts <- integer(nrow(denom_u))
for (i in seq_len(nrow(denom_u)))
  re_counts[i] <- sum(denom$assay == denom_u$assay[i] &
                        denom$Taxa  == denom_u$Taxa[i] &
                        denom$Gene  == denom_u$Gene[i])
re_sum <- data.frame(assay = denom_u$assay, Taxa = denom_u$Taxa, Gene = denom_u$Gene,
                     n_complete = re_counts, stringsAsFactors = FALSE)
m8 <- merge(res_sum[, c("assay", "Taxa", "Gene", "n_complete"), drop = FALSE],
            re_sum, by = c("assay", "Taxa", "Gene"), all = TRUE, sort = FALSE,
            suffixes = c("_08", "_re"))
bad8 <- m8[is.na(m8$n_complete_08) | is.na(m8$n_complete_re) |
             m8$n_complete_08 != m8$n_complete_re, , drop = FALSE]
check("[5.2] 08 n_complete == re-derivation from 07 (suitable + parsable + amp_bp > 0)",
      nrow(bad8) == 0L,
      paste0("offender(s): ",
             paste(sprintf("%s/%s/%s (08=%s, re-derived=%s)",
                           bad8$assay, bad8$Taxa, bad8$Gene,
                           ifelse(is.na(bad8$n_complete_08), "NA", bad8$n_complete_08),
                           ifelse(is.na(bad8$n_complete_re), "NA", bad8$n_complete_re)),
                   collapse = " ; ")))

## ---- [5.2] combined AxisLabel per (assay, Gene) ------------------------------
groups_all <- sort(unique(lab$Taxa))
assays_all <- sort(unique(lab$assay))
str_parts <- vector("list", length(assays_all))
for (i in seq_along(assays_all)) {
  a   <- assays_all[i]
  sub <- lab[lab$assay == a, , drop = FALSE]
  G   <- sub$Gene[1L]
  check(sprintf("[5.2] %s: one gene across all groups", a),
        length(unique(sub$Gene)) == 1L,
        "assay spans two genes in the label table — pair table drift?")
  lines <- c(a)
  for (g in groups_all) {
    sg <- sub[sub$Taxa == g, ]
    if (nrow(sg) == 0L) next
    lines <- c(lines, sg$label_line[1L])
  }
  str_parts[[i]] <- data.frame(Assay = a, Gene = G,
                               AxisLabel = paste(lines, collapse = "\n"),
                               stringsAsFactors = FALSE)
}
lab_str <- do.call(rbind, str_parts)

show_df(lab, header = "Per (assay, group) label counts:")
disp <- lab_str
disp$AxisLabel <- gsub("\n", " | ", disp$AxisLabel, fixed = TRUE)
show_df(disp, header = "Combined AxisLabel per (assay, Gene):")

## ---- [5.3] outputs (results/) ---------------------------------------------------
cat("\n[5.3] outputs (results/)\n")
LAB_OUT <- file.path(OUT, "13_assay_axis_labels.csv")
STR_OUT <- file.path(OUT, "14_assay_axis_label_strings.csv")
write.csv(lab, LAB_OUT, row.names = FALSE)
write.csv(lab_str, STR_OUT, row.names = FALSE)
lab_chk <- read.csv(LAB_OUT, stringsAsFactors = FALSE)
str_chk <- read.csv(STR_OUT, stringsAsFactors = FALSE)
check("[5.3] label table round-trip schema",
      identical(names(lab_chk), c("assay", "Gene", "Taxa", "A_Total", "A_NoCov", "A_Indel",
                                  "A_Good", "A_Bad", "n_unsuitable", "n_unmapped", "label_line")),
      paste0("columns found: ", paste(names(lab_chk), collapse = ", ")))
check("[5.3] label table round-trip values",
      nrow(lab_chk) == nrow(lab) &&
        identical(lab_chk$assay, lab$assay) &&
        identical(lab_chk$Gene, lab$Gene) &&
        identical(lab_chk$Taxa, lab$Taxa) &&
        identical(as.integer(lab_chk$A_Total), as.integer(lab$A_Total)) &&
        identical(as.integer(lab_chk$A_NoCov), as.integer(lab$A_NoCov)) &&
        identical(as.integer(lab_chk$A_Indel), as.integer(lab$A_Indel)) &&
        identical(as.integer(lab_chk$A_Good), as.integer(lab$A_Good)) &&
        identical(as.integer(lab_chk$A_Bad), as.integer(lab$A_Bad)) &&
        identical(lab_chk$label_line, lab$label_line),
      "re-read 13_assay_axis_labels.csv differs from the in-memory table")
check("[5.3] strings round-trip schema",
      identical(names(str_chk), c("Assay", "Gene", "AxisLabel")),
      paste0("columns found: ", paste(names(str_chk), collapse = ", ")))
check("[5.3] strings round-trip values",
      nrow(str_chk) == nrow(lab_str) &&
        identical(str_chk$Assay, lab_str$Assay) &&
        identical(str_chk$Gene, lab_str$Gene) &&
        identical(str_chk$AxisLabel, lab_str$AxisLabel),
      "re-read 14_assay_axis_label_strings.csv differs from the in-memory table")
cat(sprintf("    13_assay_axis_labels.csv          %d rows (assay x group)\n", nrow(lab)))
cat(sprintf("    14_assay_axis_label_strings.csv   %d rows (assay x gene; AxisLabel for the Rmd join)\n", nrow(lab_str)))

## ---- [5.4] audit log ----------------------------------------------------------------
log_row <- log_row_make(run_id, "5", "PASS", BASE, parent_run = s4_latest$run_id,
                        metrics = paste0("n_label_rows=", nrow(lab),
                                         ";n_assays=", nrow(lab_str),
                                         ";n_groups=", length(groups_all),
                                         ";n_good=", sum(lab$A_Good),
                                         ";n_bad=", sum(lab$A_Bad),
                                         ";score_threshold=", CONFIG$SCORE_THRESHOLD))
if (file.exists(LOG))
  log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)
log_chk <- read.csv(LOG, stringsAsFactors = FALSE)
check("[5.4] audit log appended",
      any(log_chk$run_id == run_id & log_chk$step == "5"),
      "log row missing after write")

cat("\n============================================================================\n")
cat(sprintf("STEP 5 PASSED | run_id %s (continued from step 4 %s | step 3 %s | step 2 %s | step 1 %s | step 0 %s)\n",
            run_id, s4_latest$run_id, s3_latest$run_id, s2_latest$run_id, s1_latest$run_id, parent_run))
cat(sprintf("  labels    : %d (assay, group) rows over %d assay(s) x %d group(s)\n",
            nrow(lab), nrow(lab_str), length(groups_all)))
cat(sprintf("  A_Good total: %d (modular scoring: total <= %s) | A_Bad total: %d\n",
            sum(lab$A_Good), CONFIG$SCORE_THRESHOLD, sum(lab$A_Bad)))
cat("Outputs: results/13_assay_axis_labels.csv | 14_assay_axis_label_strings.csv | run_audit_log.csv\n")
cat("Join 14_assay_axis_label_strings.csv into the hybrid plot on Assay.\n")
cat("============================================================================\n")
