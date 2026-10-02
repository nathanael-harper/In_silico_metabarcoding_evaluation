#!/usr/bin/env Rscript
# ============================================================================
# step3_amplicon_extract.r — STEP 3: AMPLICON EXTRACTION  (Windows R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
#
# WHAT THIS STEP DOES (and only this)
#   For every (assay pair, taxon) where BOTH strands mapped in step 1
#   (04_pair_coverage.csv cov == "full"):
#     1. reads the two mapped primer windows from 03_primer_mapping.csv
#        (MSA gapped coordinates) and forms the amplicon span
#            lo = min(F_start, R_start) ; hi = max(F_end, R_end)
#        (inverted pairs, such as those on circular mtDNA, are handled
#        by the same min/max rule, not special-cased),
#     2. slices [lo, hi] out of EVERY MSA record (utils extract_window),
#        de-gaps (strip_gaps), and writes one row per record:
#            amp_bp  = de-gapped length (0 = record gapped across the whole
#                      span — a documented result, kept, never excluded)
#            amp_seq = the de-gapped amplicon sequence
#     3. joins each row to 06_assay_template_status.csv on (assay, Taxa,
#        seq_id) and carries a_status, so step 4 can filter to suitable
#        templates without re-deriving anything.
#
#   Report, never exclude: every full (assay, taxon) combination
#   contributes exactly n_records rows. Pairs with cov != "full" are
#   named + skipped as a documented result (step 1 [1.5] already named
#   the unmapped primers behind them). Inverted (F start > R start) and
#   overlapping (windows share columns) pairs are extracted on the union
#   span and NAMED in [3.2] — kept + noted, nothing is dropped silently.
#   Empty (all-gap) amplicons are kept with amp_bp = 0.
#
#   Consensus-rule note: step 3 never computes a consensus — it slices
#   the MSA coordinates step 1 already anchored, so the per-template rule
#   cannot affect any row here. The rule distribution is printed for
#   provenance only (same posture as step 2).
#
#   Span-gate note: NO extraction logic changes. Step 3 slices the FINAL
#   MSA coordinates step 1 anchored; the gate only changed which window
#   step 1 anchored, and its effects are already reflected in 03/04/06:
#   a span-rejected strand is map_mode "none" -> the pair is cov !=
#   "full" -> skipped + named in [3.2]; a span-gated anchor only changes
#   the [lo, hi] span that gets sliced. Two [3.1] contract checks:
#     (a) the mapping schema tolerance chain accepts the current 17-column
#         file; older 15-column and 14-column files are NOTE'd and used
#         AS-IS (step 3 reads only name / Taxa / gene / start / end; the
#         gate columns are never read); any other column set is a stop.
#     (b) the manifest config|max_pair_span row is CHECKed against this
#         session's CONFIG$MAX_AMP_SPAN (drift = stop; re-run step 0),
#         placed after the utils_version gate on purpose.
#   The span-gate distribution is printed in [3.1] for context only and
#   NEVER enters extraction (same posture as the consensus-rule line).
#
# INPUTS (must match the step 0 manifest + step 1 + step 2 outputs):
#   results/00_run_manifest.csv (config|max_pair_span row),
#   results/01_pair_table.csv,
#   results/03_primer_mapping.csv  (17 columns current; older 15- or
#       14-column files are NOTE'd and used as-is — step 3 never reads
#       the gate columns),
#   results/04_pair_coverage.csv  (amplicon_size_clean column accepted
#       by the existing subset check; not read here),
#   results/06_assay_template_status.csv,
#   results/run_audit_log.csv (a step 2 PASS row is required),
#   primer CSV, taxon directory CSV, the MSA FASTA files
#
# OUTPUTS (results/, overwritten each run):
#   07_amplicons.csv  one row per (full pair x MSA record):
#                     assay, Taxa, Gene, seq_id, amp_bp, amp_seq,
#                     amp_win_start, amp_win_end, a_status
#                     (f_score/r_score stay in 05/06, joinable on
#                     assay + Taxa + seq_id)
#   run_audit_log.csv appended (step = 3)
#
# AUDITS ([3.x]; CHECK stops the run, NOTE / printed results never stop):
#   [3.1] continuity gate: manifest PASS + BASE + config unchanged (incl.
#          max_pair_span); utils version (older = NOTE); step 1 + step 2
#          PASS rows in the log with the intact run chain and the same
#          BASE; 01/03/04/06 schemas (mapping: 17 cols current; older
#          NOTE'd — the gate columns are never read); mapping + full-pair
#          counts vs step 1 log; inputs unchanged vs manifest; per-MSA
#          re-fingerprint; 06 row count + unmapped/coverage alignment
#   [3.2] extraction loop (one printed line per pair); unmapped pairs
#          named + skipped; inverted / overlapping windows counted +
#          named; coordinates finite + in MSA range
#   [3.3] full re-derivation: every row's amp_seq + amp_bp recomputed
#          from the MSA via (seq_id, amp_win_start, amp_win_end);
#          a_status vocabulary; forward-oriented span vs 04
#          amplicon_size; per-taxon a_status counts (printed)
#   [3.4] output written + re-read (round-trip, schema, values)
#   [3.5] audit log appended
#
# RUN:  Rscript step3_amplicon_extract.r   (Windows R, base R only — no
#       external packages; steps 0-2 must have passed in the same
#       results/ directory)
#       Paste-safe: all conditionals braced; one CHUNK break marked.
# ============================================================================

## ---- BASE — the ONE line to edit (same line as steps 0 + 1; verified
## ---- against the manifest) ----
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
OUT <- file.path(BASE, "results")

## ---- [3.1] continuity gate ---------------------------------------------------
cat("============================================================================\n")
cat("STEP 3 - AMPLICON EXTRACTION (pure-R window slice of the gapped MSA)\n")
cat("  run_id  :", run_id, "\n")
mani_path <- file.path(OUT, "00_run_manifest.csv")
check("[3.1] step 0 manifest exists", file.exists(mani_path),
      "results/00_run_manifest.csv missing — run step0_setup_audit.r first")
mani <- read.csv(mani_path, stringsAsFactors = FALSE)
mval <- function(sec, item) {
  v <- mani$value[mani$section == sec & mani$item == item]
  if (length(v) != 1L) stop("manifest missing row: ", sec, "|", item, call. = FALSE)
  v
}
check("[3.1] step 0 status = PASS", mval("status", "step0") == "PASS",
      "step 0 did not pass — fix it before extraction")
parent_run <- mval("run", "run_id")
check("[3.1] BASE unchanged since step 0", identical(mval("run", "base"), BASE),
      paste0("manifest BASE: ", mval("run", "base"), " | session BASE: ", BASE,
             " — if you moved things, re-run step 0"))
check("[3.1] max_map_mismatches unchanged since step 0",
      identical(mval("config", "max_map_mismatches"), as.character(CONFIG$MAX_MAP_MISMATCHES)),
      "CONFIG$MAX_MAP_MISMATCHES differs from the manifest — set it, then re-run step 0")
check("[3.1] score_threshold unchanged since step 0",
      identical(mval("config", "score_threshold"), as.character(CONFIG$SCORE_THRESHOLD)),
      "CONFIG$SCORE_THRESHOLD differs from the manifest — set it, then re-run step 0")

## utils continuity — same tolerant handling as steps 1/2
uv_row <- mani$value[mani$section == "run" & mani$item == "utils_version"]
if (length(uv_row) == 1L) {
  check("[3.1] utils_version matches manifest", uv_row == UTILS_VERSION,
        paste0("manifest utils_version is '", uv_row,
               "', this session sources utils.R '", UTILS_VERSION,
               "' — re-run step 0 to republish the manifest"))
} else {
  note("[3.1]", sprintf("manifest has no utils_version row; this session sources utils.R %s",
                        UTILS_VERSION))
}

## Span-gate parameter: must not have drifted since step 0.
## Placed AFTER the utils_version gate on purpose: an older manifest
## already stops at the version check above, so by the time a run reaches
## here the config row is guaranteed to exist. Same ordering logic as
## steps 1 [1.1] and 2 [2.1].
check("[3.1] max_pair_span unchanged since step 0",
      identical(mval("config", "max_pair_span"), as.character(CONFIG$MAX_AMP_SPAN)),
      "CONFIG$MAX_AMP_SPAN differs from the manifest — set it, then re-run step 0")

## Run chain: step 2 (parent = the step 1 run) <- step 1 (parent = the step 0
## run). Step 2's log row carries the WSL spelling of BASE (it ran on WSL),
## so both bases are compared in WSL form via to_wsl_path (POSIX paths pass
## through unchanged; Windows forms are derived).
LOG <- file.path(OUT, "run_audit_log.csv")
check("[3.1] audit log exists", file.exists(LOG), "results/run_audit_log.csv missing")
logdf <- read.csv(LOG, stringsAsFactors = FALSE)
s2 <- logdf[logdf$step == "2" & logdf$status == "PASS", ]
check("[3.1] step 2 completed (PASS row in log)", nrow(s2) >= 1L,
      "no PASS row for step 2 — run step2_primerminer.r (WSL) first")
s2_latest <- s2[nrow(s2), ]
s1 <- logdf[logdf$step == "1" & logdf$status == "PASS" &
              logdf$run_id == s2_latest$parent_run_id, ]
check("[3.1] step 2 continues a step 1 run", nrow(s1) >= 1L,
      paste0("step 2 log row ", s2_latest$run_id,
             " has no matching step 1 PASS row — chain broken"))
s1_latest <- s1[nrow(s1), ]
check("[3.1] step 1 continues the step 0 run",
      identical(as.character(s1_latest$parent_run_id), parent_run),
      paste0("step 1 log parent: ", s1_latest$parent_run_id,
             " | step 0 run: ", parent_run))
check("[3.1] steps 1 + 2 ran from the same BASE",
      identical(to_wsl_path(as.character(s1_latest$base)), to_wsl_path(mval("run", "base"))) &&
        identical(to_wsl_path(as.character(s2_latest$base)), to_wsl_path(mval("run", "base"))),
      paste0("step 1 base: ", s1_latest$base, " | step 2 base: ", s2_latest$base,
             " | manifest base: ", mval("run", "base")))

## Step 1 + step 2 outputs exist + schemas
PT_IN  <- file.path(OUT, "01_pair_table.csv")
MAP_IN <- file.path(OUT, "03_primer_mapping.csv")
COV_IN <- file.path(OUT, "04_pair_coverage.csv")
STA_IN <- file.path(OUT, "06_assay_template_status.csv")
check("[3.1] step 1 + step 2 outputs exist",
      all(file.exists(c(PT_IN, MAP_IN, COV_IN, STA_IN))),
      "missing one of: 01_pair_table.csv, 03_primer_mapping.csv, 04_pair_coverage.csv, 06_assay_template_status.csv")
pair_table <- read.csv(PT_IN, stringsAsFactors = FALSE)
map_df     <- read.csv(MAP_IN, stringsAsFactors = FALSE)
cov_df     <- read.csv(COV_IN, stringsAsFactors = FALSE)
status_df  <- read.csv(STA_IN, stringsAsFactors = FALSE)

## Mapping schema: the current 17 columns = the 15 base columns +
## span_gate + span_gate_mm. Step 3 reads only name / Taxa / gene /
## start / end from the mapping file (it never reads consensus_rule and
## never reads the gate columns), so older 15-column and 14-column files
## are NOTE'd and used AS-IS — no backfill needed; any OTHER column set
## is a stop.
MAP_COLS_BASE <- c("name", "seq", "gene", "Taxa", "start", "end",
                  "target_with_gaps", "n_mismatches", "min_mismatches",
                  "map_orientation", "other_orientation_mm", "map_mode",
                  "win_min_occ", "win_n_low_occ", "consensus_rule")
MAP_COLS <- c(MAP_COLS_BASE, "span_gate", "span_gate_mm")
map_current    <- identical(names(map_df), MAP_COLS)
map_prev    <- identical(names(map_df), MAP_COLS_BASE)
map_legacy <- identical(names(map_df), MAP_COLS_BASE[-length(MAP_COLS_BASE)])
check("[3.1] mapping schema (17 cols current; 15 or 14 col files used as-is)",
      map_current || map_prev || map_legacy,
      paste0("columns found: ", paste(names(map_df), collapse = ", ")))
if (!map_current) {
  note("[3.1]", sprintf("03_primer_mapping.csv has no span-gate columns (%d cols) — used as-is (step 3 reads only name / Taxa / gene / start / end; the gate columns never enter extraction); re-run steps 0 + 1 to publish the span-gate columns",
                        ncol(map_df)))
}
check("[3.1] coverage schema",
      all(c("assay", "gene", "Taxa", "cov", "amplicon_size", "orientation_ok") %in% names(cov_df)),
      paste0("columns found: ", paste(names(cov_df), collapse = ", ")))
check("[3.1] status schema",
      identical(names(status_df), c("Taxa", "Gene", "assay", "seq_id", "a_status",
                                    "f_score", "r_score", "f_p_status", "r_p_status")),
      paste0("columns found: ", paste(names(status_df), collapse = ", ")))

## Count cross-checks against the step 1 log (the same gate step 2 [2.1] runs)
mm_get <- function(row, key) {
  parts <- strsplit(as.character(row$metrics), ";", fixed = TRUE)[[1L]]
  hit <- parts[startsWith(parts, paste0(key, "="))]
  if (length(hit) != 1L) stop("step 1 log metrics missing: ", key, call. = FALSE)
  as.integer(sub("^.*=", "", hit))
}
check("[3.1] mapping row count matches step 1 log",
      nrow(map_df) == mm_get(s1_latest, "n_mapping_rows"),
      paste0("03 has ", nrow(map_df), " rows; step 1 log has ",
             mm_get(s1_latest, "n_mapping_rows"), " — re-run step 1"))
check("[3.1] full pair count matches step 1 log",
      sum(cov_df$cov == "full") == mm_get(s1_latest, "n_full_pairs"),
      paste0("04 has ", sum(cov_df$cov == "full"), " full combos; step 1 log has ",
             mm_get(s1_latest, "n_full_pairs"), " — re-run step 1"))
check("[3.1] full pair count >= 1", sum(cov_df$cov == "full") >= 1L,
      "no full (assay, taxon) combination — step 1 [1.6] should have caught this")

## Inputs unchanged vs manifest (same gate as steps 1/2)
prim <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon[] <- lapply(taxon, trimws)
taxon$full_path <- file.path(BASE, taxon$filepath)
check("[3.1] primers_input.csv unchanged",
      nrow(prim) == as.integer(mval("inputs", "primer_rows")),
      "re-run step 0 after changing inputs")
check("[3.1] taxon_template_directory.csv unchanged",
      nrow(taxon) == as.integer(mval("inputs", "taxon_rows")),
      "re-run step 0 after changing inputs")

## Per-MSA reload + fingerprint vs manifest (files must be unchanged since
## step 2 sliced them)
TEMPLATES <- vector("list", nrow(taxon))
fp_get <- function(v, tag) {
  parts <- strsplit(v, ";", fixed = TRUE)[[1L]]
  pre <- paste0(tag, "=")
  hit <- parts[startsWith(parts, pre)]
  if (length(hit) != 1L)
    stop("manifest fingerprint: expected exactly one '", pre, "', found ", length(hit), call. = FALSE)
  substring(hit, nchar(pre) + 1L)
}
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  key <- paste0(T, "|", G)
  f <- read_fasta(taxon$full_path[i])
  f2 <- fasta_dedup(f)
  fp <- mval("templates", key)
  last <- f2$seqs[f2$records]
  diffs <- character(0)
  if (f$ragged) diffs <- c(diffs, "file is now ragged (uniform width required)")
  if (f2$records != as.integer(fp_get(fp, "records")))
    diffs <- c(diffs, sprintf("records %d vs manifest %s", f2$records, fp_get(fp, "records")))
  if (f2$width[1L] != as.integer(fp_get(fp, "width")))
    diffs <- c(diffs, sprintf("width %d vs manifest %s", f2$width[1L], fp_get(fp, "width")))
  if (f2$names[1L] != fp_get(fp, "first_name"))
    diffs <- c(diffs, sprintf("first name %s vs manifest %s", f2$names[1L], fp_get(fp, "first_name")))
  tail_now <- substr(last, pmax(1L, nchar(last) - 39L), nchar(last))
  if (tail_now != fp_get(fp, "last_tail40"))
    diffs <- c(diffs, "last record tail-40 differs from manifest")
  check(sprintf("[3.1] %s %s unchanged", T, G), length(diffs) == 0L, paste(diffs, collapse = " ; "))
  TEMPLATES[[i]] <- f2
}

## 06 structural gates: row count per (gene) + unmapped rows align with
## coverage (a pair is cov "full" iff step 2 wrote no "unmapped" rows for it)
expect_sta <- 0L
for (i in seq_len(nrow(taxon)))
  expect_sta <- expect_sta + TEMPLATES[[i]]$records *
    sum(pair_table$gene == taxon$Gene[i])
check("[3.1] 06 rows = sum over templates of n_records x n_pairs(gene)",
      nrow(status_df) == expect_sta,
      paste0("06 has ", nrow(status_df), " rows; inputs imply ", expect_sta,
             " — re-run step 2"))
A_STATUS_VOCAB <- c("suitable", "unsuitable", "no_coverage", "incomplete_coverage", "unmapped")
check("[3.1] 06 a_status vocabulary",
      all(status_df$a_status %in% A_STATUS_VOCAB),
      paste0("offending value(s): ",
             paste(unique(setdiff(status_df$a_status, A_STATUS_VOCAB)), collapse = ", ")))
bad_unm <- character(0)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  nT <- TEMPLATES[[i]]$records
  pg <- pair_table[pair_table$gene == G, ]
  for (k in seq_len(nrow(pg))) {
    p    <- pg$assay[k] 
    covk <- cov_df$cov[cov_df$assay == p & cov_df$Taxa == T]
    if (length(covk) != 1L) {
      bad_unm <- c(bad_unm, sprintf("%s/%s (%d coverage rows)", T, p, length(covk)))
    } else {
      stsub <- status_df[status_df$assay == p & status_df$Taxa == T, ]
      u <- sum(stsub$a_status == "unmapped")
      if (covk == "full" && u != 0L)
        bad_unm <- c(bad_unm, sprintf("%s/%s (cov full but %d unmapped rows)", T, p, u))
      if (covk != "full" && u != nT)
        bad_unm <- c(bad_unm, sprintf("%s/%s (cov %s but only %d of %d rows unmapped)",
                                      T, p, covk, u, nT))
    }
  }
}
check("[3.1] 06 unmapped rows align with 04 coverage", length(bad_unm) == 0L,
      paste(bad_unm, collapse = " ; "))

cat("  parent  : step 0 run ", parent_run, " | step 1 ", s1_latest$run_id,
    " | step 2 ", s2_latest$run_id, " (chain verified)\n")
cat("  base    :", BASE, "\n")
## Span-gate parameter prints with the other config values
## (manifest-verified in [3.1] above).
cat("  config  : max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | score_threshold =", CONFIG$SCORE_THRESHOLD,
    " | max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
rule_rows <- mani[mani$section == "consensus_rule", ]
if (nrow(rule_rows) > 0L)
  cat(sprintf("  consensus rules : %d gaps_can_win | %d gaps_never_win (context only — never enter extraction)\n",
              sum(rule_rows$value == "gaps_can_win"),
              sum(rule_rows$value == "gaps_never_win")))
## Span-gate distribution, context only (current 17-column file only;
## same posture as the consensus-rule line — the gate never enters
## extraction).
if (map_current) {
  cat(sprintf("  span gate     : MAX_AMP_SPAN = %d bp | gated=%d | no_op=%d | rejected=%d | ungated=%d (context only — the gate never enters extraction)\n",
              CONFIG$MAX_AMP_SPAN,
              sum(map_df$span_gate == "gated"),
              sum(map_df$span_gate == "no_op"),
              sum(map_df$span_gate == "rejected"),
              sum(map_df$span_gate == "ungated")))
}
cat("============================================================================\n")

## ==== CHUNK BREAK (safe: all gates + TEMPLATES defined; paste resumes here) ====

## ---- [3.2] extraction loop ------------------------------------------------------
cat("\n[3.2] extraction (amplicon span = union of the two mapped primer windows;\n",
    "         every MSA record yields one row; de-gapped; a_status carried from step 2;\n",
    "         inverted / overlapping windows kept + named — nothing dropped)\n", sep = "")

amp_all <- data.frame()     # includes helper columns (tpl_idx, seq_pos) for [3.3]
n_extracted_pairs <- 0L
n_skipped_pairs   <- 0L
n_inverted        <- 0L
n_overlap         <- 0L
n_empty           <- 0L
skipped_names     <- character(0)

for (i in seq_len(nrow(taxon))) {
  T  <- taxon$Taxa[i]; G <- taxon$Gene[i]
  f2 <- TEMPLATES[[i]]
  nT <- f2$records
  wM <- f2$width[1L]
  pg <- pair_table[pair_table$gene == G, ]
  for (k in seq_len(nrow(pg))) {
    p    <- pg$assay[k]
    frow <- map_df[map_df$name == paste0(p, "_F") & map_df$Taxa == T, ]
    rrow <- map_df[map_df$name == paste0(p, "_R") & map_df$Taxa == T, ]
    covk <- cov_df$cov[cov_df$assay == p & cov_df$Taxa == T & cov_df$gene == G] 
    check(sprintf("[3.2] %s %s %s has one coverage row", T, G, p), length(covk) == 1L,
          paste0("expected 1 coverage row, found ", length(covk), " — re-run step 1"))
    if (covk != "full") {
      n_skipped_pairs <- n_skipped_pairs + 1L
      skipped_names <- c(skipped_names, sprintf("%s/%s (cov=%s)", T, p, covk))
      next
    }
    check(sprintf("[3.2] %s %s %s both strands have one mapping row", T, G, p),
          nrow(frow) == 1L && nrow(rrow) == 1L,
          "a full pair without exactly one mapping row per strand — re-run step 1")
    fs <- frow$start[1L]; fe <- frow$end[1L]
    rs <- rrow$start[1L]; re <- rrow$end[1L]
    check(sprintf("[3.2] %s %s %s coordinates finite + in MSA range", T, G, p),
          is.finite(fs) && is.finite(fe) && is.finite(rs) && is.finite(re) &&
            fs >= 1L && fe >= fs && rs >= 1L && re >= rs && fe <= wM && re <= wM,
          "mapping coordinates invalid for this MSA — re-run step 1")

    inverted <- fs > rs
    overlap  <- pmin(fe, re) >= pmax(fs, rs)
    if (inverted) n_inverted <- n_inverted + 1L
    if (overlap)  n_overlap  <- n_overlap + 1L
    lo <- pmin(fs, rs)
    hi <- pmax(fe, re)

    ## per-record a_status from step 2 (join on seq_id, MSA record order)
    stsub <- status_df[status_df$assay == p & status_df$Taxa == T, ]
    check(sprintf("[3.2] %s %s %s status rows = n_records in MSA order", T, G, p),
          nrow(stsub) == nT && identical(stsub$seq_id, f2$names),
          "06 rows out of step for this pair — re-run step 2")
    st <- setNames(stsub$a_status, stsub$seq_id)[f2$names]

    wins  <- extract_window(f2$seqs, lo, hi)
    ampsq <- strip_gaps(wins)
    bp    <- as.integer(nchar(ampsq))
    n_empty <- n_empty + sum(bp == 0L)

    chunk <- data.frame(
      assay = p, Taxa = T, Gene = G, seq_id = f2$names,
      amp_bp = bp, amp_seq = ampsq,
      amp_win_start = as.integer(lo), amp_win_end = as.integer(hi),
      a_status = st,
      tpl_idx = i, seq_pos = seq_len(nT),
      stringsAsFactors = FALSE)
    amp_all <- rbind(amp_all, chunk)
    n_extracted_pairs <- n_extracted_pairs + 1L

    flag <- ""
    if (inverted) flag <- paste0(flag, " [inverted F>R]")
    if (overlap)  flag <- paste0(flag, " [windows overlap]")
    med <- if (sum(bp > 0L) > 0L) sprintf("%.0f", median(bp[bp > 0L])) else "NA"
    cat(sprintf("  %-14s %-4s %-22s window %d-%d (span %d) | %d empty, %d with amplicon | med %s bp (max %d)%s\n",
                T, G, p, lo, hi, hi - lo + 1L, sum(bp == 0L), nT - sum(bp == 0L), med, max(bp), flag))
  }
}
check("[3.2] at least one full pair extracted", n_extracted_pairs >= 1L,
      "no full (assay, taxon) combination extracted — step 1 [1.6] should have caught this")
if (n_skipped_pairs > 0L)
  note("[3.2]", sprintf("%d pair(s) skipped (cov != full; named): %s",
                        n_skipped_pairs, paste(skipped_names, collapse = "; ")))
if (n_inverted > 0L)
  note("[3.2]", sprintf("%d pair(s) inverted (F start > R start) — union span used (step 1 [1.6] flagged them)", n_inverted))
if (n_overlap > 0L)
  note("[3.2]", sprintf("%d pair(s) with overlapping primer windows — KEPT on the union span + named", n_overlap))

## ---- [3.3] row invariants (every row re-derives from MSA + own coordinates) ----
cat("\n[3.3] row invariants (per-row re-derivation from the MSA; a_status vocabulary; 04 cross-check)\n")
PUB_COLS <- c("assay", "Taxa", "Gene", "seq_id", "amp_bp", "amp_seq",
              "amp_win_start", "amp_win_end", "a_status")
amps_pub <- amp_all[, PUB_COLS, drop = FALSE]
bad <- character(0)
for (i in seq_len(nrow(amp_all))) {
  r  <- amp_all[i, ]
  f2 <- TEMPLATES[[r$tpl_idx]]
  reseq <- strip_gaps(extract_window(f2$seqs[r$seq_pos], r$amp_win_start, r$amp_win_end))
  if (!identical(reseq, r$amp_seq) || nchar(reseq) != r$amp_bp)
    bad <- c(bad, sprintf("%s/%s/%s: re-derived amplicon differs from the stored one",
                          r$Taxa, r$assay, r$seq_id))
  if (!r$a_status %in% A_STATUS_VOCAB)
    bad <- c(bad, sprintf("%s/%s/%s: a_status '%s' not in vocabulary",
                          r$Taxa, r$assay, r$seq_id, r$a_status))
}
check("[3.3] all rows re-derive (amp_seq + amp_bp + a_status)", length(bad) == 0L,
      paste(head(bad, 10L), collapse = " ; "))
note("[3.3]", sprintf("%d of %d rows carry an empty (all-gap) amplicon (amp_bp = 0) — kept, never excluded",
                      n_empty, nrow(amps_pub)))

## Forward-oriented pairs: the extracted union span must equal step 1's
## amplicon_size (04). Inverted pairs carry the F-left-of-R value from step
## 1 [1.6], which is NOT the union span — NOTE'd, not CHECKed.
win_uniq <- amps_pub[!duplicated(amps_pub[, c("assay", "Taxa"), drop = FALSE]),
                     c("assay", "Taxa", "Gene", "amp_win_start", "amp_win_end"), drop = FALSE]
m <- merge(win_uniq,
           cov_df[, c("assay", "Taxa", "gene", "cov", "amplicon_size", "orientation_ok")], # <--- Change pair to assay
           by.x = c("assay", "Taxa", "Gene"), 
           by.y = c("assay", "Taxa", "gene"))                                              # <--- Change pair to assay
check("[3.3] coverage rows align 1:1 with extracted pairs, all full",
      nrow(m) == nrow(win_uniq) && all(m$cov == "full"),
      "extracted pairs and 04_pair_coverage.csv disagree — re-run step 1")
spanm    <- m$amp_win_end - m$amp_win_start + 1L
fwd      <- m$orientation_ok == TRUE
fwdbad   <- m[fwd & spanm != m$amplicon_size, ]
check("[3.3] forward pairs: extracted span == 04 amplicon_size",
      nrow(fwdbad) == 0L,
      paste0("offender(s): ",
             paste(sprintf("%s/%s (span %d vs %d)", fwdbad$assay, fwdbad$Taxa,
                           spanm[fwd & spanm != m$amplicon_size], fwdbad$amplicon_size),
                   collapse = " ; ")))
inv <- m$orientation_ok == FALSE
if (any(inv))
  note("[3.3]", sprintf("%d inverted pair(s): 04 amplicon_size is computed F-left-of-R (step 1 [1.6]) and does NOT equal the extracted union span: %s",
                        sum(inv),
                        paste(sprintf("%s/%s", m$assay[inv], m$Taxa[inv]), collapse = ", ")))

cat("  per-taxon a_status counts (rows = record x pair):\n")
for (Tx in sort(unique(amps_pub$Taxa))) {
  sub <- amps_pub[amps_pub$Taxa == Tx, ]
  cat(sprintf("    %-14s rows=%-7d suitable=%-7d unsuitable=%-7d no_cov=%-5d incomp=%-5d unmapped=%-5d | suitable amplicon bp total = %s\n",
              Tx, nrow(sub), sum(sub$a_status == "suitable"), sum(sub$a_status == "unsuitable"),
              sum(sub$a_status == "no_coverage"), sum(sub$a_status == "incomplete_coverage"),
              sum(sub$a_status == "unmapped"),
              format(sum(sub$amp_bp[sub$a_status == "suitable"]), big.mark = ",")))
}

## ---- [3.4] output ------------------------------------------------------------------
cat("\n[3.4] output (results/)\n")
AMP_OUT <- file.path(OUT, "07_amplicons.csv")
write.csv(amps_pub, AMP_OUT, row.names = FALSE)
amp_chk <- read.csv(AMP_OUT, stringsAsFactors = FALSE)
check("[3.4] amplicons round-trip schema", identical(names(amp_chk), PUB_COLS),
      paste0("columns found: ", paste(names(amp_chk), collapse = ", ")))
check("[3.4] amplicons round-trip values",
      nrow(amp_chk) == nrow(amps_pub) &&
        identical(amp_chk$assay, amps_pub$assay) &&
        identical(amp_chk$Taxa, amps_pub$Taxa) &&
        identical(amp_chk$Gene, amps_pub$Gene) &&
        identical(amp_chk$seq_id, amps_pub$seq_id) &&
        identical(amp_chk$amp_bp, amps_pub$amp_bp) &&
        identical(amp_chk$amp_seq, amps_pub$amp_seq) &&
        identical(amp_chk$amp_win_start, amps_pub$amp_win_start) &&
        identical(amp_chk$amp_win_end, amps_pub$amp_win_end) &&
        identical(amp_chk$a_status, amps_pub$a_status),
      "re-read amplicon table differs from the in-memory table")
cat(sprintf("    07_amplicons.csv  %d rows (empty amp_bp=0: %d)\n", nrow(amps_pub), n_empty))

## ---- [3.5] audit log -------------------------------------------------------------------
log_row <- log_row_make(run_id, "3", "PASS", BASE, parent_run = s2_latest$run_id,
                        metrics = paste0("n_amp_rows=", nrow(amps_pub),
                                         ";n_pairs_extracted=", n_extracted_pairs,
                                         ";n_pairs_skipped=", n_skipped_pairs,
                                         ";n_inverted_pairs=", n_inverted,
                                         ";n_overlap_pairs=", n_overlap,
                                         ";n_empty_amplicons=", n_empty,
                                         ";n_suitable_rows=", sum(amps_pub$a_status == "suitable")))
log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)
log_chk <- read.csv(LOG, stringsAsFactors = FALSE)
check("[3.5] audit log appended",
      any(log_chk$run_id == run_id & log_chk$step == "3"),
      "log row missing after write")

cat("\n============================================================================\n")
cat(sprintf("STEP 3 PASSED | run_id %s (continued from step 2 %s | step 1 %s | step 0 %s)\n",
            run_id, s2_latest$run_id, s1_latest$run_id, parent_run))
cat(sprintf("  extracted pairs : %d full (assay, taxon) combinations (%d skipped, cov != full)\n",
            n_extracted_pairs, n_skipped_pairs))
cat(sprintf("  inverted / overlap windows: %d / %d (kept + named)\n", n_inverted, n_overlap))
cat(sprintf("  amplicon rows   : %d (empty amp_bp=0: %d) | suitable rows: %d\n",
            nrow(amps_pub), n_empty, sum(amps_pub$a_status == "suitable")))
cat("Outputs: results/07_amplicons.csv | run_audit_log.csv\n")
cat("Next: step 4 (Windows R)\n")
cat("============================================================================\n")
