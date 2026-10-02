#!/usr/bin/env Rscript
# ============================================================================
# step2_primerminer.r — STEP 2: PRIMERMINER SCORING  (Ubuntu/WSL R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
#
# WHAT THIS STEP DOES (and only this)
#   For every complete assay pair x taxon (step 1's pair table + mapping):
#     1. for each strand step 1 mapped, slices the gapped MSA window
#        [start, end] out of every template record and de-gaps it,
#     2. classifies each primer x template record (05_primer_scores.csv):
#            de-gapped length 0            -> no_coverage
#            de-gapped length != primer    -> incomplete_coverage
#            else                          -> scored
#     3. scores the scored records: one temporary FASTA per (Taxa, Gene,
#        primer), and evaluate_primer called with the STRAND-CORRECT
#        orientation: the primer is ALWAYS passed as written; forward = T
#        for F primers, forward = F for R primers. For R primers
#        PrimerMiner reverse-complements the template window internally,
#        so the position penalty table (high at the 3' end, low at the
#        5' end) is applied identically to BOTH strands — a 3' mismatch
#        costs the same for F and R. A second diagnostic call —
#        RC(primer) with forward flipped — is recorded as score_rc: it
#        carries the same mismatch set but with the 3'/5' position
#        weights swapped, so a large divergence from score_written flags
#        an orientation anomaly; score_rc is NEVER used for
#        classification.
#        PrimerMiner is called per mapped primer x taxon (two calls:
#        strand-correct + diagnostic), NOT per template.
#     4. joins F and R per template (06_assay_template_status.csv) with the
#        five-class precedence:
#            either strand unmapped (step 1 NA) -> unmapped
#            any primer no_coverage             -> no_coverage
#            any primer incomplete_coverage     -> incomplete_coverage
#            both scored:
#                total = score_F + score_R
#                total <= SCORE_THRESHOLD      -> suitable
#                total >  SCORE_THRESHOLD      -> unsuitable
#                (NA score -> unsuitable; counted in [2.6])
#
#   FAULT-TOLERANCE CONTRACT: the [2.3] scoring loop NEVER halts on a
#   per-primer problem:
#     - an UNMAPPED strand (map_mode = "none", incl. span-rejected) is
#       skipped with `next` — the loop continues;
#     - a PM CALL FAILURE (strand-correct OR diagnostic) is SOFT: the full
#       error (condition class, call, message, + traceback when rlang is
#       available) is emitted as a NOTE, every template of that primer
#       gets NA scores (-> unsuitable, named in [2.6]), and the loop
#       continues to the next primer;
#     - only TOTAL failure of the strand-correct calls (every call failed)
#       is a stop: the [2.2] canary already proved PM runs on synthetic
#       data, so total failure is an environment problem, not data.
#   The [2.2] canary remains a HARD environment gate (fatal = TRUE): if
#   evaluate_primer cannot run on synthetic data in either orientation,
#   nothing in this step can be trusted and the run stops before touching
#   any project file.
#
#   Consensus-rule contract: step 1 anchors windows under a PER-TEMPLATE
#   consensus engine (gaps_never_win / gaps_can_win, selected from the
#   manifest consensus_rule section step 0 publishes) and writes the rule
#   into every row of 03_primer_mapping.csv plus rule counters in its log.
#   Step 2 does NOT compute any consensus — it slices MSA coordinates step
#   1 already anchored — so the rule NEVER enters any classification below
#   (same posture as the occupancy columns). But the rule is part of the
#   continuity contract, so [2.1] verifies it end-to-end:
#     - CONFIG consensus default + overrides in the rule vocabulary;
#     - manifest consensus_rule section present -> its config rows must
#       equal this session's CONFIG (drift = stop, re-run step 0); section
#       absent (older manifest) -> NOTE + CONFIG fallback (tolerant);
#     - 03_primer_mapping.csv: current 17-column schema; older 15-column
#       or 14-column files = NOTE + backfill (tolerant); every row's
#       consensus_rule in the vocabulary and, when the manifest section is
#       present, equal to the manifest row for that template;
#     - rule-aware win_min_occ invariant (mirrors step 1 [1.7]): under
#       gaps_never_win a mapped window cannot cross a zero-coverage column
#       -> win_min_occ >= 1 is CHECKed; under gaps_can_win holes are
#       legitimate -> win_min_occ = 0 is allowed, NOTE-level only;
#     - step 1 log rule counters (n_templates_gaps_can_win,
#       n_win_min_occ_zero_rows) cross-checked against the mapping file.
#   Under gaps_can_win a mapped window span may include consensus-hole
#   columns: template records gapped at a hole column de-gap to a shorter
#   window and are classified incomplete_coverage by the existing [2.3]
#   partition — no new classification class is needed, and nothing is
#   excluded on rule grounds.
#
#   Span-gate contract: step 1 also applies the PAIR-AWARE SPAN GATE
#   (CONFIG$MAX_AMP_SPAN) and writes span_gate / span_gate_mm into every
#   row of 03_primer_mapping.csv plus span-gate counters in its log. Step
#   2 does NOT apply the gate — it slices MSA coordinates step 1 already
#   anchored — so the gate NEVER enters any classification or score below
#   (same posture as the consensus rule and the occupancy columns): a
#   span-rejected primer is an unmapped strand (the existing unmapped
#   class in [2.4]) and a span-gated anchor only changes which window gets
#   sliced. But the gate IS part of the continuity contract, so [2.1]
#   verifies it end-to-end:
#     - manifest config|max_pair_span == this session's CONFIG$MAX_AMP_SPAN
#       (drift = stop, re-run step 0; placed AFTER the utils_version gate,
#       which stops on older manifests, so the row is guaranteed to exist
#       here);
#     - 03_primer_mapping.csv schema chain: 17 cols (current) = clean;
#       15 cols = NOTE + backfill span_gate = "ungated", span_gate_mm = NA;
#       14 cols = NOTE + backfill consensus_rule + the two span-gate
#       columns; anything else = stop;
#     - span-gate structural contract on the mapping file: span_gate in
#       the vocabulary {ungated, no_op, gated, rejected}; span_gate_mm NA
#       iff span_gate in {ungated, no_op}; rejected implies map_mode none;
#     - step 1 log span-gate counters (n_span_gated, n_span_noop,
#       n_span_rejected) cross-checked against the mapping file (older log
#       = NOTE, tolerant).
#
#   Nothing is hardcoded: every acceptance check is a structural identity
#   (partition sums, row counts, round-trips), so any primer set, any
#   genes, any MSA dimensions run through.
#
# WHY WSL
#   PrimerMiner's dependency tree installs reliably only on Ubuntu/WSL R.
#   BASE below is the WSL spelling of the SAME folder step 0 wrote the
#   manifest from. [2.1] verifies BASE against the manifest by deriving the
#   WSL form at runtime (utils to_wsl_path) — strict gate, no second
#   hardcoded path, and a moved folder still stops the run.
#
# INPUTS (all must match the step 0 manifest + step 1 outputs):
#   results/00_run_manifest.csv (consensus_rule section + config|max_pair_span),
#   results/01_pair_table.csv,
#   results/03_primer_mapping.csv  (17 columns; older 15- or 14-column
#       files tolerated with NOTE + backfill),
#   results/04_pair_coverage.csv  (amplicon_size_clean column accepted
#       automatically; the check below is a subset check),
#   primer CSV, taxon directory CSV, the MSA FASTA files
#
# OUTPUTS (results/, overwritten each run):
#   05_primer_scores.csv            one row per (mapped primer x template):
#                                   p_status, degapped_len, score_written
#                                   (THE strand-correct score: F ->
#                                   forward=T, R -> forward=F; primer as
#                                   written in both cases), score_rc
#                                   (diagnostic alt call: RC primer,
#                                   forward flipped — same mismatches,
#                                   3'/5' weights swapped; never
#                                   classified), score = score_written
#   06_assay_template_status.csv    one row per (assay x template):
#                                   a_status + f_score/r_score + per-strand
#                                   p_status
#   run_audit_log.csv               appended (step = 2)
#
# AUDITS ([2.x]; CHECK stops the run, NOTE never stops):
#   [2.1] continuity gate: manifest PASS, BASE == WSL form of manifest base,
#          config unchanged (incl. max_pair_span), consensus-rule contract
#          (vocabulary, manifest config rows == CONFIG, per-template mapping
#          rule == manifest row; older = NOTE + fallback), utils version,
#          mapping schema chain (17 cols clean; older = NOTE + backfill) +
#          span-gate structural contract (vocabulary, NA pattern, rejected
#          implies unmapped), step 1 outputs + log consistent (incl. rule
#          and span-gate counters), inputs unchanged, per-MSA re-fingerprint,
#          rule-aware occupancy contract on the mapping file (NA iff
#          unmapped; non-negative; gaps_never_win win_min_occ >= 1;
#          gaps_can_win win_min_occ = 0 allowed + named; low-occupancy
#          windows named — never excluded)
#   [2.2] PrimerMiner installed; schema canary on synthetic data covering
#          BOTH orientations (forward=T and forward=F; FATAL — environment
#          gate); every real PM output has Template + sum, numeric, exact
#          input template set
#   [2.3] per (Taxa, Gene, primer): scored + no_coverage + incomplete
#          = n_templates; unmapped strands skipped (loop continues); PM
#          call failures SOFT (NOTE + NA + loop continues); total
#          strand-correct failure = stop (environment problem)
#   [2.4] per (Taxa, Gene): assay-status rows = n_templates x n_pairs
#   [2.5] per taxon x assay class counts (printed)
#   [2.6] NA-score count + no scored-NA row is suitable + PM call-failure
#          counts
#   [2.7] threshold + adjacent provenance lines, strand-correct scoring
#          rule, PM call-failure counts, MAX_AMP_SPAN, consensus rule
#          distribution, span-gate distribution, boundary-combo count
#   [2.8] outputs written + re-read (round-trip, schema, values; score ==
#          score_written identity)
#   [2.9] audit log appended
#
# RUN:  Rscript step2_primerminer.r   (Ubuntu/WSL R, FRESH session — do not
#       load an old workspace; PrimerMiner installed ONCE via
#       BiocManager::install("PrimerMiner"); step 1 must have passed in the
#       same results/ directory)
#       Paste-safe: all conditionals braced; one CHUNK break marked.
# ============================================================================

## ---- BASE — the ONE line to edit (WSL spelling of the folder; verified
## ---- against the manifest in [2.1] via to_wsl_path) ----
BASE <- "/mnt/c/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"



source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
OUT <- file.path(BASE, "results")

## ---- [2.1] continuity gate: manifest -----------------------------------------
cat("============================================================================\n")
cat("STEP 2 - PRIMERMINER SCORING (Ubuntu/WSL R; consensus-rule + span-gate continuity checks; strand-correct scoring; fault-tolerant loop)\n")
cat("  run_id  :", run_id, "\n")
mani_path <- file.path(OUT, "00_run_manifest.csv")
check("[2.1] step 0 manifest exists", file.exists(mani_path),
      "results/00_run_manifest.csv missing — run step0_setup_audit.r first")
mani <- read.csv(mani_path, stringsAsFactors = FALSE)
mval <- function(sec, item) {
  v <- mani$value[mani$section == sec & mani$item == item]
  if (length(v) != 1L) stop("manifest missing row: ", sec, "|", item, call. = FALSE)
  v
}
check("[2.1] step 0 status = PASS", mval("status", "step0") == "PASS",
      "step 0 did not pass — fix it before scoring")
parent_run <- mval("run", "run_id")
check("[2.1] BASE matches manifest (WSL form)",
      identical(BASE, to_wsl_path(mval("run", "base"))),
      paste0("manifest base: ", mval("run", "base"),
             " (WSL: ", to_wsl_path(mval("run", "base")),
             ") | session BASE: ", BASE, " — re-run step 0 if you moved things"))
check("[2.1] score_threshold unchanged since step 0",
      identical(mval("config", "score_threshold"), as.character(CONFIG$SCORE_THRESHOLD)),
      "CONFIG$SCORE_THRESHOLD differs from the manifest — set it, then re-run step 0")
check("[2.1] max_map_mismatches unchanged since step 0",
      identical(mval("config", "max_map_mismatches"), as.character(CONFIG$MAX_MAP_MISMATCHES)),
      "CONFIG$MAX_MAP_MISMATCHES differs from the manifest — set it, then re-run step 0")

## utils continuity: the utils_version published by step 0 must equal the
## utils.R sourced in THIS session. Older manifests lack the row -> NOTE
## instead of FAIL (tolerant by design).
if ("utils_version" %in% mani$item[mani$section == "run"]) {
  check("[2.1] utils version matches manifest",
        identical(mval("run", "utils_version"), UTILS_VERSION),
        paste0("manifest utils_version: ", mval("run", "utils_version"),
               " | session utils: ", UTILS_VERSION,
               " — utils.R changed after step 0; re-run step 0"))
} else {
  note("[2.1]", "manifest has no utils_version row — re-run step 0 to enable this gate")
}

## Span-gate parameter: the manifest value must equal this session's
## CONFIG. Placed AFTER the utils_version gate on purpose: an older
## manifest stops at the version check above, and a current manifest
## always carries config|max_pair_span (step 0 publishes it), so mval is
## guaranteed to find the row here.
check("[2.1] max_pair_span unchanged since step 0",
      identical(mval("config", "max_pair_span"), as.character(CONFIG$MAX_AMP_SPAN)),
      "CONFIG$MAX_AMP_SPAN differs from the manifest — set it, then re-run step 0")

## ---- [2.1] consensus-rule contract ------------------------------------------------
## Step 2 never computes a consensus (it slices MSA coordinates step 1
## anchored), but the rule IS part of the continuity contract: if the
## manifest rules drift from this session's CONFIG, the mapping file's
## coordinates + occupancy columns were derived under a different engine
## than CONFIG now implies -> stop (re-run step 0). Older manifest (no
## section) = NOTE + CONFIG fallback (tolerant, same posture as the
## utils_version handling above).
CFG_DEFAULT <- CONFIG$CONSENSUS_RULE
CFG_OVERR   <- CONFIG$CONSENSUS_RULE_OVERRIDES
check("[2.1] CONFIG consensus default in vocabulary",
      CFG_DEFAULT %in% CONSENSUS_RULES_VALID,
      paste0("CONFIG$CONSENSUS_RULE '", CFG_DEFAULT, "' is not a valid rule (valid: ",
             paste(CONSENSUS_RULES_VALID, collapse = "/"), ")"))
check("[2.1] CONFIG override values in vocabulary",
      is.list(CFG_OVERR) &&
        all(vapply(CFG_OVERR, function(v) v %in% CONSENSUS_RULES_VALID, logical(1))),
      "CONFIG$CONSENSUS_RULE_OVERRIDES must be a list whose values are valid rules")

rule_rows            <- mani[mani$section == "consensus_rule", ]
rule_section_present <- nrow(rule_rows) > 0L
mani_rule_default    <- mani$value[mani$section == "config" & mani$item == "consensus_rule"]
mani_rule_ov         <- mani$value[mani$section == "config" & mani$item == "consensus_rule_overrides"]

parse_packed_overrides <- function(s) {
  s <- trimws(s)
  if (s == "" || s == "(none)") {
    return(list())
  }
  parts <- strsplit(s, ";", fixed = TRUE)[[1L]]
  out <- list()
  for (p in parts) {
    kv <- strsplit(trimws(p), "=", fixed = TRUE)[[1L]]
    if (length(kv) != 2L) {
      stop("manifest consensus_rule_overrides pack malformed: '", p, "'", call. = FALSE)
    }
    out[[trimws(kv[1L])]] <- trimws(kv[2L])
  }
  out
}

if (rule_section_present) {
  check("[2.1] manifest consensus_rule default == CONFIG",
        length(mani_rule_default) == 1L && mani_rule_default[1L] == CFG_DEFAULT,
        paste0("manifest config|consensus_rule = '",
               if (length(mani_rule_default) == 1L) mani_rule_default[1L] else "(missing)",
               "' vs CONFIG '", CFG_DEFAULT, "' — re-run step 0"))
  mani_ov <- parse_packed_overrides(
    if (length(mani_rule_ov) == 1L) mani_rule_ov[1L] else "(none)")
  check("[2.1] manifest overrides == CONFIG (key = value, set-equal)",
        setequal(names(mani_ov), names(CFG_OVERR)) &&
          all(vapply(names(mani_ov),
                     function(k) isTRUE(mani_ov[[k]] == CFG_OVERR[[k]]),
                     logical(1))),
        paste0("manifest overrides: '",
               if (length(mani_rule_ov) == 1L) mani_rule_ov[1L] else "(missing)",
               "' vs CONFIG: '",
               if (length(CFG_OVERR) == 0L) "(none)"
               else paste(sprintf("%s=%s", names(CFG_OVERR), unlist(CFG_OVERR)), collapse = "; "),
               "' — re-run step 0"))
} else {
  note("[2.1]", "manifest has no consensus_rule section — falling back to CONFIG-derived rules (consensus_rule_for) for every template")
}

cat("  parent  : step 0 run ", parent_run, " (manifest)\n")
cat("  base    :", BASE, "\n")
cat("  config  : score_threshold =", CONFIG$SCORE_THRESHOLD,
    " | max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
cat("           pm_adjacent =", CONFIG$PRIMERMINER_ADJACENT, "\n")
cat("  scoring : strand-correct (F -> forward=T, R -> forward=F; primer as written)\n")
cat("  loop    : fault-tolerant (unmapped -> skip; PM call failure -> NOTE + NA + continue)\n")
cat("  consensus_rule : default =", CFG_DEFAULT,
    "| overrides =",
    if (length(CFG_OVERR) == 0L) "(none)"
    else paste(sprintf("%s=%s", names(CFG_OVERR), unlist(CFG_OVERR)), collapse = "; "), "\n")
cat("  manifest   : consensus_rule section",
    if (rule_section_present) "present" else "ABSENT (older step 0 manifest; CONFIG fallback)", "\n")
cat("============================================================================\n")

## ---- [2.1] step 1 outputs + audit log consistent ---------------------------------
PT_IN  <- file.path(OUT, "01_pair_table.csv")
MAP_IN <- file.path(OUT, "03_primer_mapping.csv")
COV_IN <- file.path(OUT, "04_pair_coverage.csv")
check("[2.1] step 0 + step 1 outputs exist",
      file.exists(PT_IN) && file.exists(MAP_IN) && file.exists(COV_IN),
      "missing one of: 01_pair_table.csv, 03_primer_mapping.csv, 04_pair_coverage.csv — run steps 0 and 1 first")
pair_table <- read.csv(PT_IN, stringsAsFactors = FALSE)
map_df     <- read.csv(MAP_IN, stringsAsFactors = FALSE)
cov_df     <- read.csv(COV_IN, stringsAsFactors = FALSE)
n_cols_read <- ncol(map_df)

## Mapping schema: the current 17 columns = the 15 base columns +
## span_gate + span_gate_mm (16th/17th). Tolerance chain, same posture as
## the consensus_rule handling: a 15-column file = NOTE + backfill
## span_gate = "ungated" / span_gate_mm = NA; a 14-column file = NOTE +
## backfill consensus_rule (CONFIG-derived) + the two span-gate columns.
## Any OTHER column set is a stop. Backfilled columns are appended in
## schema order, so downstream code (incl. the checks below) sees the full
## 17-column frame either way.
MAP_COLS <- c("name", "seq", "gene", "Taxa", "start", "end",
              "target_with_gaps", "n_mismatches", "min_mismatches",
              "map_orientation", "other_orientation_mm", "map_mode",
              "win_min_occ", "win_n_low_occ", "consensus_rule",
              "span_gate", "span_gate_mm")
map_current <- identical(names(map_df), MAP_COLS)
map_prev    <- identical(names(map_df), MAP_COLS[-c(16L, 17L)])
map_legacy  <- identical(names(map_df), MAP_COLS[-c(15L, 16L, 17L)])
check("[2.1] mapping schema (17 cols current; 15 or 14 col files backfilled)",
      map_current || map_prev || map_legacy,
      paste0("columns found: ", paste(names(map_df), collapse = ", ")))
if (!map_current) {
  if (map_legacy) {
    note("[2.1]", "03_primer_mapping.csv has no consensus_rule column — falling back to CONFIG-derived rules (consensus_rule_for) per row; re-run steps 0 + 1 to publish the rule contract")
    map_df$consensus_rule <- vapply(seq_len(nrow(map_df)), function(k)
      consensus_rule_for(paste0(map_df$Taxa[k], "|", map_df$gene[k]), CFG_OVERR, CFG_DEFAULT),
      character(1))
  }
  note("[2.1]", sprintf("03_primer_mapping.csv has no span-gate columns (%d cols) — backfilling span_gate = 'ungated', span_gate_mm = NA; re-run steps 0 + 1 to publish the span-gate columns",
                        n_cols_read))
  map_df$span_gate    <- "ungated"
  map_df$span_gate_mm <- NA_integer_
}
check("[2.1] mapping consensus_rule values in vocabulary",
      all(map_df$consensus_rule %in% CONSENSUS_RULES_VALID),
      paste0("offending value(s): ",
             paste(unique(setdiff(map_df$consensus_rule, CONSENSUS_RULES_VALID)), collapse = ", ")))
if (rule_section_present) {
  mani_rules <- setNames(rule_rows$value, rule_rows$item)
  map_keys   <- paste0(map_df$Taxa, "|", map_df$gene)
  ok_rule    <- vapply(seq_len(nrow(map_df)), function(k) {
    r <- mani_rules[[map_keys[k]]]
    isTRUE(!is.null(r) && r == map_df$consensus_rule[k])
  }, logical(1))
  bad_keys   <- unique(map_keys[!ok_rule])
  check("[2.1] every mapping row carries the manifest rule for its template",
        all(ok_rule),
        paste0("template(s) where the mapping row rule differs from the manifest: ",
               paste(head(bad_keys, 10L), collapse = ", "),
               " — re-run step 1 (or step 0 if the manifest itself drifted)"))
}

## Span-gate structural contract (same posture as the consensus_rule
## vocabulary check — pure structural identities that can only fire on a
## stale or hand-edited 03, never on a clean current run):
##   - span_gate in the four-value vocabulary;
##   - span_gate_mm NA iff span_gate in {ungated, no_op} (gated/rejected
##     rows carry the displaced pre-gate anchor's mm — reported, never
##     blank);
##   - a span-rejected row is an unmapped row (map_mode = "none"): the
##     gate's rejection IS the existing unmapped class, so [2.4] counts it
##     with the unmapped precedence branch — no new class.
SPAN_GATE_VALID <- c("ungated", "no_op", "gated", "rejected")
check("[2.1] mapping span_gate values in vocabulary",
      all(map_df$span_gate %in% SPAN_GATE_VALID),
      paste0("offending value(s): ",
             paste(unique(setdiff(map_df$span_gate, SPAN_GATE_VALID)), collapse = ", ")))
check("[2.1] span_gate_mm NA pattern (NA iff ungated/no_op)",
      identical(is.na(map_df$span_gate_mm),
                map_df$span_gate %in% c("ungated", "no_op")),
      "span_gate_mm present on an ungated/no_op row or missing on a gated/rejected one")
check("[2.1] span-rejected rows are unmapped (map_mode = none)",
      all(map_df$map_mode[map_df$span_gate == "rejected"] == "none"),
      "a span-rejected row carries map_mode != 'none' — re-run step 1")

## Rule-aware occupancy contract (report, never exclude).
## CHECK the structure; NOTE the values.
##   - NA iff unmapped (every row, both rules).
##   - non-negative (every row, both rules).
##   - rule-aware win_min_occ (mirrors step 1 [1.7]; this is the step 2
##     side of the same contract):
##       gaps_never_win: a consensus '-' column exists only at occ = 0
##       (zero coverage), so a mapped window span carrying one is a
##       structural anomaly -> win_min_occ >= 1 is CHECKed.
##       gaps_can_win: gap-majority columns are legitimate holes; a mapped
##       span may cross one -> win_min_occ = 0 is ALLOWED, NOTE-level only.
##   (fallback rows carry their CONFIG-derived rule here too, so the same
##   partition applies.)
mmap <- map_df$map_mode != "none"
check("[2.1] occupancy columns: NA iff unmapped",
      identical(is.na(map_df$win_min_occ), !mmap) &&
        identical(is.na(map_df$win_n_low_occ), !mmap),
      "win_min_occ / win_n_low_occ missing on a mapped row or present on an unmapped one")
check("[2.1] occupancy values non-negative",
      all(map_df$win_min_occ[mmap] >= 0) && all(map_df$win_n_low_occ[mmap] >= 0),
      "negative occupancy count in 03_primer_mapping.csv")
gnw <- mmap & map_df$consensus_rule == "gaps_never_win"
bad_gnw <- map_df[gnw & map_df$win_min_occ < 1L, ]
check("[2.1] gaps_never_win mapped rows: win_min_occ >= 1 (no zero-coverage column in span)",
      nrow(bad_gnw) == 0L,
      paste0("offender(s): ",
             paste(sprintf("%s/%s (win_min_occ=%d)", bad_gnw$Taxa, bad_gnw$name, bad_gnw$win_min_occ),
                   collapse = ", "),
             " — step 1 [1.7] gates the same invariant; re-run step 1"))
z0 <- map_df[mmap & map_df$win_min_occ == 0L, ]
if (nrow(z0) > 0L)
  note("[2.1]", sprintf("%d mapped window(s) with win_min_occ = 0: %d under gaps_can_win (legitimate holes; never excluded) | %d under gaps_never_win (CHECKed above): %s",
                        nrow(z0), sum(z0$consensus_rule == "gaps_can_win"),
                        sum(z0$consensus_rule == "gaps_never_win"),
                        paste(sprintf("%s/%s [%s]", z0$Taxa, z0$name, z0$consensus_rule), collapse = ", ")))
lowwin <- map_df[mmap & map_df$win_n_low_occ > 0, ]
note("[2.1]", sprintf("%d of %d mapped window(s) enter step 2 with >=1 column below MIN_OCC=%d (report only; never excluded): %s",
                      nrow(lowwin), sum(mmap), CONFIG$MIN_OCC,
                      if (nrow(lowwin) > 0L) paste(sprintf("%s/%s", lowwin$Taxa, lowwin$name), collapse = ", ")
                      else "none"))

check("[2.1] coverage schema",
      all(c("assay", "gene", "Taxa", "cov") %in% names(cov_df)),
      paste0("columns found: ", paste(names(cov_df), collapse = ", ")))
LOG <- file.path(OUT, "run_audit_log.csv")
check("[2.1] audit log exists", file.exists(LOG), "results/run_audit_log.csv missing")
logdf <- read.csv(LOG, stringsAsFactors = FALSE)
s1 <- logdf[logdf$step == "1" & logdf$status == "PASS" &
              logdf$parent_run_id == parent_run, ]
check("[2.1] step 1 completed (PASS row in log)", nrow(s1) >= 1L,
      paste0("no PASS row for step 1 attached to step 0 run ", parent_run,
             " — run step1_primer_mapping.r first"))
s1_latest <- s1[nrow(s1), ]                       # log appends newest last
check("[2.1] step 1 ran from the same BASE as step 0",
      identical(as.character(s1_latest$base), mval("run", "base")),
      paste0("step 1 log base: ", s1_latest$base, " | manifest base: ", mval("run", "base")))
mm_get <- function(row, key) {
  parts <- strsplit(as.character(row$metrics), ";", fixed = TRUE)[[1L]]
  hit <- parts[startsWith(parts, paste0(key, "="))]
  if (length(hit) != 1L) stop("step 1 log metrics missing: ", key, call. = FALSE)
  as.integer(sub("^.*=", "", hit))
}
check("[2.1] mapping row count matches step 1 log",
      nrow(map_df) == mm_get(s1_latest, "n_mapping_rows"),
      paste0("03_primer_mapping.csv has ", nrow(map_df), " rows; step 1 log has ",
             mm_get(s1_latest, "n_mapping_rows"), " — outputs out of sync, re-run step 1"))
check("[2.1] full (assay, taxon) count matches step 1 log",
      sum(cov_df$cov == "full") == mm_get(s1_latest, "n_full_pairs"),
      paste0("04_pair_coverage.csv has ", sum(cov_df$cov == "full"), " full combos; step 1 log has ",
             mm_get(s1_latest, "n_full_pairs"), " — outputs out of sync, re-run step 1"))

## Rule counters in the step 1 log (present in current logs): cross-check
## them against the mapping file so a stale 03 cannot survive a step 1
## re-run. Absent (older log) = NOTE (tolerant).
tpl_keys  <- paste0(map_df$Taxa, "|", map_df$gene)
uq_tpl    <- !duplicated(tpl_keys)
n_tpl_gcw <- sum(map_df$consensus_rule[uq_tpl] == "gaps_can_win")
n_tpl_gnw <- sum(map_df$consensus_rule[uq_tpl] == "gaps_never_win")
n_win0    <- sum(mmap & map_df$win_min_occ == 0L)
m1s <- as.character(s1_latest$metrics)
if (grepl("n_templates_gaps_can_win=", m1s, fixed = TRUE)) {
  check("[2.1] gaps_can_win template count matches step 1 log",
        n_tpl_gcw == mm_get(s1_latest, "n_templates_gaps_can_win"),
        paste0("mapping implies ", n_tpl_gcw, " gaps_can_win template(s); step 1 log has ",
               mm_get(s1_latest, "n_templates_gaps_can_win"), " — outputs out of sync, re-run step 1"))
  check("[2.1] win_min_occ=0 row count matches step 1 log",
        n_win0 == mm_get(s1_latest, "n_win_min_occ_zero_rows"),
        paste0("mapping has ", n_win0, " win_min_occ=0 row(s); step 1 log has ",
               mm_get(s1_latest, "n_win_min_occ_zero_rows"), " — outputs out of sync, re-run step 1"))
} else {
  note("[2.1]", "step 1 log has no consensus-rule counters — skipping rule sync checks (tolerant)")
}

## Span-gate counters in the step 1 log (present in current logs):
## cross-check them against the mapping file so a stale 03 cannot survive a
## step 1 re-run (a backfilled 03 would recount 0/0/0 and mismatch a
## real current log). Absent (older log) = NOTE (tolerant, same posture
## as the rule counters).
n_span_gated    <- sum(map_df$span_gate == "gated")
n_span_noop     <- sum(map_df$span_gate == "no_op")
n_span_rejected <- sum(map_df$span_gate == "rejected")
if (grepl("n_span_gated=", m1s, fixed = TRUE)) {
  check("[2.1] span-gate gated count matches step 1 log",
        n_span_gated == mm_get(s1_latest, "n_span_gated"),
        paste0("mapping implies ", n_span_gated, " gated row(s); step 1 log has ",
               mm_get(s1_latest, "n_span_gated"), " — outputs out of sync, re-run step 1"))
  check("[2.1] span-gate no_op count matches step 1 log",
        n_span_noop == mm_get(s1_latest, "n_span_noop"),
        paste0("mapping implies ", n_span_noop, " no_op row(s); step 1 log has ",
               mm_get(s1_latest, "n_span_noop"), " — outputs out of sync, re-run step 1"))
  check("[2.1] span-gate rejected count matches step 1 log",
        n_span_rejected == mm_get(s1_latest, "n_span_rejected"),
        paste0("mapping implies ", n_span_rejected, " rejected row(s); step 1 log has ",
               mm_get(s1_latest, "n_span_rejected"), " — outputs out of sync, re-run step 1"))
} else {
  note("[2.1]", "step 1 log has no span-gate counters — skipping span sync checks (tolerant)")
}

## ---- [2.1] inputs unchanged + per-MSA re-fingerprint -------------------------------
prim <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
prim[] <- lapply(prim, trimws)
taxon <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon[] <- lapply(taxon, trimws)
taxon$full_path <- file.path(BASE, taxon$filepath)
check("[2.1] primers_input2.csv unchanged",
      nrow(prim) == as.integer(mval("inputs", "primer_rows")),
      paste0("manifest has ", mval("inputs", "primer_rows"), " rows, file has ", nrow(prim),
             " — re-run step 0 after changing inputs"))
check("[2.1] taxon_template_directory_2.csv unchanged",
      nrow(taxon) == as.integer(mval("inputs", "taxon_rows")),
      "directory row count differs from manifest — re-run step 0")
check("[2.1] pair table consistent",
      nrow(pair_table) == as.integer(mval("inputs", "pairs_complete")),
      "results/01_pair_table.csv row count differs from manifest — re-run step 0")
expect_map_rows <- sum(vapply(seq_len(nrow(taxon)),
                              function(i) sum(prim$gene == taxon$Gene[i]), integer(1)))
check("[2.1] mapping covers every directory row",
      nrow(map_df) == expect_map_rows,
      paste0("03 has ", nrow(map_df), " rows; inputs imply ", expect_map_rows, " — re-run step 1"))

cat("\n[2.1] MSA reload + fingerprint vs manifest (files must be unchanged since step 1)\n")
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
  check(sprintf("[2.1] %s %s unchanged", T, G), length(diffs) == 0L, paste(diffs, collapse = " ; "))
  TEMPLATES[[i]] <- f2
}

## ==== CHUNK BREAK (safe: inputs, TEMPLATES, all gate data defined; paste resumes here) ====

## ---- [2.2] PrimerMiner environment + schema canary ---------------------------------
cat("\n[2.2] PrimerMiner environment + schema canary (synthetic data only)\n")
check("[2.2] PrimerMiner installed",
      requireNamespace("PrimerMiner", quietly = TRUE),
      paste0("install ONCE on this WSL R: ",
             "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); ",
             "BiocManager::install('PrimerMiner')"))
suppressMessages(library(PrimerMiner))

## pm_evaluate: call evaluate_primer on a temporary FASTA of de-gapped
## windows and return a named numeric vector (template name -> sum).
## forward: the ORIENTATION of the call. forward = T scores the primer
## (as written) 5'->3' against the window as sliced; forward = F makes
## PrimerMiner reverse-complement the window internally (use for R
## primers: the primer's 3' end then lands on the high-penalty position
## table row, exactly as for F primers).
## fatal: TRUE (default) = a PM error stops the run via CHECK (the [2.2]
## canary — the environment gate: if PM cannot run on synthetic data in
## either orientation, nothing in this step can be trusted). FALSE
## (fault-tolerance contract) = a PM error is SOFT: the full error
## (class, call, message, + traceback when rlang is available) is emitted
## as a NOTE, NA is returned for every template (-> unsuitable, named in
## [2.6]), and the [2.3] scoring loop KEEPS RUNNING.
## Guards (fatal mode): uniform primer-length input, output CSV written,
## Template + sum columns, numeric sum, and the EXACT input template set
## returned (order-independent). If your PrimerMiner version names the
## columns differently, the schema CHECK names them — edit the two
## literals below.
pm_evaluate <- function(primer_seq, windows, tmpdir, tag, forward = TRUE, fatal = TRUE) {
  Lp <- nchar(primer_seq)
  check(sprintf("%s PM input uniform primer-length", tag),
        length(windows) >= 1L && all(nchar(windows) == Lp),
        "scored windows are not all exactly primer-length — coverage classification drifted")
  nW <- length(windows)
  ## PrimerMiner's output 'Template' column does NOT reliably echo the full
  ## FASTA header: it truncates to the first whitespace-delimited token
  ## ('NC_045917.1/1-18800 SomeGenus sp. ...' -> 'NC_045917.1/1-18800'), so a
  ## full-header comparison fails even when PM returned every template.
  ## The temp FASTA therefore uses synthetic space-free IDs (pm_0001, ...);
  ## sums are mapped back to the REAL names by ID match, no order assumed.
  ids <- paste0("pm_", sprintf("%04d", seq_len(nW)))
  ln <- character(2L * nW)
  for (j in seq_len(nW)) {
    ln[2L * j - 1L] <- paste0(">", ids[j])
    ln[2L * j]      <- windows[j]
  }
  fa <- tempfile(fileext = ".fasta", tmpdir = tmpdir)
  writeLines(ln, fa)
  cs <- tempfile(fileext = ".csv", tmpdir = tmpdir)
  pm_err <- ""
  ok <- tryCatch({
    evaluate_primer(alignment_imp = fa, primer_sequ = primer_seq,
                    start = 1L, stop = Lp, forward = forward,
                    adjacent = CONFIG$PRIMERMINER_ADJACENT, save = cs)
    TRUE
  }, error = function(e) {
    ## rich capture: conditionMessage(e) alone can be EMPTY (as seen with
    ## the forward=F path on real data) — record class, call, message and
    ## a traceback so the NOTE is actionable. `<<-` assigns into pm_evaluate's
    ## environment (the handler's plain <- would create a handler-local).
    pm_err <<- paste0(
      "[", paste(class(e), collapse = "/"), "] ",
      if (is.null(conditionCall(e))) ""
      else paste0("call: ", paste(deparse(conditionCall(e), width.cutoff = 200L), collapse = " "), " | "),
      "msg: ", conditionMessage(e),
      if (nzchar(conditionMessage(e))) "" else " (EMPTY MESSAGE — see traceback)",
      tryCatch({
        if (requireNamespace("rlang", quietly = TRUE)) {
          tb <- rlang::trace_back(max_bullets = 10L)
          paste0("\n  traceback:\n  ", paste(format(tb), collapse = "\n  "))
        } else ""
      }, error = function(x) ""))
    FALSE
  })
  if (!isTRUE(ok)) {
    if (fatal) {
      check(sprintf("%s evaluate_primer ran", tag), FALSE,
            paste0("evaluate_primer errored (environment gate — fix before resuming):\n", pm_err))
    } else {
      note(sprintf("%s PM call failed (SOFT — loop continues)", tag),
           paste0("evaluate_primer errored -> ALL templates of this primer get NA scores (unsuitable, named in [2.6]):\n", pm_err))
      out <- setNames(rep(NA_real_, nW), names(windows))
      attr(out, "pm_call_error") <- TRUE
      return(out)
    }
  }
  check(sprintf("%s PM wrote output", tag),
        file.exists(cs) && file.info(cs)$size > 0,
        "evaluate_primer produced no non-empty output CSV at the save path")
  out <- read.csv(cs, stringsAsFactors = FALSE, check.names = FALSE)
  check(sprintf("%s PM output schema (Template + sum)", tag),
        "Template" %in% names(out) && "sum" %in% names(out),
        paste0("columns found: ", paste(names(out), collapse = ", ")))
  sum_ok <- is.numeric(out$sum) || (is.logical(out$sum) && all(is.na(out$sum)))
  check(sprintf("%s PM sum is numeric or all-NA", tag), sum_ok,
        paste0("sum column type: ", class(out$sum)[1L]))
  got   <- as.character(out$Template)
  miss  <- setdiff(ids, got)
  extra <- setdiff(got, ids)
  check(sprintf("%s PM returned exactly the input template set", tag),
        nrow(out) == nW && length(miss) == 0L && length(extra) == 0L,
        sprintf("expected %d | got %d | missing: %s | unexpected: %s",
                nW, nrow(out),
                if (length(miss)) paste(head(miss, 5L), collapse = ", ") else "none",
                if (length(extra)) paste(head(extra, 5L), collapse = ", ") else "none"))
  ## map synthetic IDs back to input positions (PM row order not assumed)
  pos <- match(ids, got)
  if (any(is.na(pos)))
    stop("pm_evaluate: internal ID alignment failed", call. = FALSE)
  sums <- as.numeric(out$sum)[pos]
  n_unsc <- sum(is.na(sums))
  if (n_unsc > 0L)
    note(sprintf("%s PM unscorable windows", tag),
         sprintf("%d of %d scored template(s) returned an NA penalty (window content outside the PM penalty table) -> unsuitable, named in [2.6]",
                 n_unsc, nW))
  setNames(sums, names(windows))
}

## canary: run the real call shape once per ORIENTATION on synthetic data
## before touching the project files — verifies install, API shape (incl.
## the forward argument), output schema, template echo. FATAL by design:
## this is the environment gate (a PM that cannot run here at all cannot
## be run on the project data).
canary_dir <- tempfile(pattern = "pm_canary_"); dir.create(canary_dir)
canary_primer <- "ACGTACGTACGTACGTACGT"
canary_wins <- c(ct1 = "ACGTACGTACGTACGTACGT",
                 ct2 = "ACGTACGTACGTACGTACGT", 
                 ct3 = "TTTTTTTTTTTTTTTTTTTT",
                 "NC_000001.1/1-20 Test species, complete genome" = "TTTTTTTTTTTTTTTTTTTA")
can_out <- pm_evaluate(canary_primer, canary_wins, canary_dir,
                       "[2.2] canary fwd", forward = TRUE)
check("[2.2] canary (forward=T) keyed by full header (space/comma safe)",
      identical(names(can_out), names(canary_wins)),
      paste0("returned names: ", paste(names(can_out), collapse = " | ")))
can_out_rc <- pm_evaluate(reverse_complement(canary_primer), canary_wins,
                          canary_dir, "[2.2] canary rev", forward = FALSE)
check("[2.2] canary (forward=F) keyed by full header (space/comma safe)",
      identical(names(can_out_rc), names(canary_wins)),
      paste0("returned names: ", paste(names(can_out_rc), collapse = " | ")))

unlink(canary_dir, recursive = TRUE)
cat("    canary passed: evaluate_primer output has Template + sum (numeric), full template set, both orientations\n")

## ---- [2.3] scoring loop ---------------------------------------------------------------
cat("\n[2.3] scoring (one PM call-pair per mapped primer x taxon, not per template;\n",
    "         strand-correct: F -> forward=T, R -> forward=F, primer ALWAYS as written\n",
    "         (PM rev-comps the window internally for R -> 3' end carries the high\n",
    "         penalty on BOTH strands); fault-tolerant: unmapped strands are\n",
    "         skipped (loop continues) and PM call failures are SOFT (NOTE + NA\n",
    "         scores -> unsuitable, loop continues);\n",
    "         rule-aware: gaps_can_win windows may include hole columns — a record\n",
    "         gapped at a hole de-gaps shorter and is classified incomplete_coverage)\n", sep = "")
score_rows <- vector("list", 0)
## Fault-tolerance counters (reported in [2.6]/[2.7]/[2.9] + final summary)
n_pm_calls_primary <- 0L
n_pm_fail          <- 0L   # strand-correct call failures (soft)
n_pm_rc_fail       <- 0L   # diagnostic rc call failures (soft)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  f2 <- TEMPLATES[[i]]
  nT <- f2$records
  subprim <- prim[prim$gene == G, ]
  for (j in seq_len(nrow(subprim))) {
    nm <- subprim$name[j]
    ps <- subprim$seq[j]
    Lp <- nchar(ps)
    rmap <- map_df[map_df$name == nm & map_df$Taxa == T & map_df$gene == G, ]
    check(sprintf("[2.3] %s %s %s has one mapping row", T, G, nm),
          nrow(rmap) == 1L,
          sprintf("expected 1 mapping row, found %d — re-run step 1", nrow(rmap)))
    if (rmap$map_mode[1L] == "none") next        # unmapped strand: no primer-level rows;
    # the loop CONTINUES (never a halt).
    # (span-rejected rows land here too —
    # the existing unmapped class)
    start <- rmap$start[1L]; end <- rmap$end[1L]
    check(sprintf("[2.3] %s %s %s window in MSA range", T, G, nm),
          is.finite(start) && is.finite(end) && end >= start && end <= f2$width[1L],
          "mapping coordinates invalid for this MSA — re-run step 1")
    
    ## Strand-correct orientation: the strand comes from the name suffix;
    ## it selects the forward flag. The primer is ALWAYS passed as written
    ## — forward=F makes PrimerMiner reverse-complement the template
    ## window internally, so the position penalty table (high at the 3'
    ## end) applies identically to F and R.
    suffix <- substr(nm, nchar(nm), nchar(nm))
    if (!suffix %in% c("F", "R"))
      stop(sprintf("%s: primer name lacks an _F/_R suffix (name='%s') — cannot determine strand for strand-correct scoring", T, nm),
           call. = FALSE)
    fwd <- identical(suffix, "F")
    ## report-only cross-check against step 1's empirical orientation:
    ## if it disagrees with the name suffix, the primer may actually bind
    ## the opposite way around — flagged, never switched (report only).
    mo <- rmap$map_orientation[1L]
    mo_fwd <- switch(mo, forward = TRUE, reverse = FALSE, NA)
    if (!is.na(mo_fwd) && !identical(mo_fwd, fwd))
      note(sprintf("[2.3] %s %s %s orientation", T, G, nm),
           sprintf("name suffix implies forward=%s but step 1 map_orientation='%s' — scoring follows the name suffix (report only; inspect step 1 if unexpected)",
                   fwd, mo))
    
    ## Individual support audit — templates that actually 0-mm-match the
    ## primer in the mapped (de-gapped) window. A mapping anchored on a
    ## degenerate consensus block (N-runs) scores 0 here = VACUOUS.
    pr_ch <- strsplit(toupper(ps), "", fixed = TRUE)[[1L]]
    win_d <- strip_gaps(extract_window(f2$seqs, start, end))
    supp_n <- sum(vapply(win_d, function(x) {
      if (nchar(x) != Lp) return(FALSE)
      wch <- strsplit(toupper(x), "", fixed = TRUE)[[1L]]
      all(vapply(seq_len(Lp), function(jj)
        length(intersect(base_set(wch[jj]), base_set(pr_ch[jj]))) > 0L,
        logical(1)))
    }, logical(1)))
    note(sprintf("[2.3] %s %s %s individual support", T, G, nm),
         sprintf("%d of %d templates 0-mm-match the primer in the mapped window%s",
                 supp_n, nT,
                 if (supp_n == 0L) "  -- VACUOUS (consensus anchor only; no individual binds)" else ""))
    
    
    dw  <- strip_gaps(extract_window(f2$seqs, start, end))
    len <- nchar(dw)
    pstat <- ifelse(len == 0L, "no_coverage",
                    ifelse(len != Lp, "incomplete_coverage", "scored"))
    iscored <- pstat == "scored"
    sw   <- rep(NA_real_, nT)
    sw_w <- sw; sw_r <- sw
    if (any(iscored)) {
      sname <- f2$names[iscored]
      swin  <- dw[iscored]; names(swin) <- sname
      tmpdir <- tempfile(pattern = sprintf("pm_%s_%s_%s_", G, nm, T))
      dir.create(tmpdir, recursive = TRUE)
      ## THE strand-correct call — this is the score. SOFT on failure:
      ## NOTE with full error text + traceback, NA scores for all
      ## templates (-> unsuitable), loop KEEPS RUNNING.
      n_pm_calls_primary <- n_pm_calls_primary + 1L
      sumw <- pm_evaluate(ps, swin, tmpdir,
                          sprintf("[2.2] %s %s %s written (forward=%s)", T, G, nm, fwd),
                          forward = fwd, fatal = FALSE)
      n_pm_fail <- n_pm_fail + as.integer(isTRUE(attr(sumw, "pm_call_error")))
      ## diagnostic alt call: RC(primer), forward flipped. Carries the same
      ## mismatch set but with the 3'/5' position weights SWAPPED — a large
      ## divergence from sumw flags an orientation anomaly. NEVER classified.
      ## Also SOFT: a diagnostic failure can never halt the run.
      sumr <- pm_evaluate(reverse_complement(ps), swin, tmpdir,
                          sprintf("[2.2] %s %s %s rc (forward=%s)", T, G, nm, !fwd),
                          forward = !fwd, fatal = FALSE)
      n_pm_rc_fail <- n_pm_rc_fail + as.integer(isTRUE(attr(sumr, "pm_call_error")))
      unlink(tmpdir, recursive = TRUE)
      sw_w[iscored] <- sumw[sname]
      sw_r[iscored] <- sumr[sname]
      sw[iscored]   <- sw_w[iscored]   # score = the strand-correct call
      # (the min-of-both-calls rule is not used: it flipped 3'/5' weights
      # for R and could mask true scores with orientation noise)
    }
    score_rows[[length(score_rows) + 1L]] <- data.frame(
      Taxa = T, Gene = G, assay = sub("_(F|R)$", "", nm), name = nm,
      strand = suffix, seq = ps, seq_id = f2$names,
      primer_len = Lp, degapped_len = len, p_status = pstat,
      score_written = sw_w, score_rc = sw_r, score = sw,
      stringsAsFactors = FALSE)
    check(sprintf("[2.3] %s %s %s partition = n_templates", T, G, nm),
          sum(pstat == "scored") + sum(pstat == "no_coverage") +
            sum(pstat == "incomplete_coverage") == nT,
          "primer x template partition does not cover every record")
    ## Occupancy suffix (context only; never excludes): prints only when
    ## the mapped window carries >=1 column below CONFIG$MIN_OCC.
    occ_str <- ""
    if (!is.na(rmap$win_n_low_occ[1L]) && rmap$win_n_low_occ[1L] > 0L)
      occ_str <- sprintf("  [occ min=%d, %d col(s) < MIN_OCC]",
                         rmap$win_min_occ[1L], rmap$win_n_low_occ[1L])
    cat(sprintf("  %-14s %-4s %-22s [%-14s | fwd=%s] n=%-5d scored=%-5d no_cov=%-3d incomp=%-3d mean_pm=%s%s\n",
                T, G, nm, rmap$consensus_rule[1L], fwd, nT, sum(iscored),
                sum(pstat == "no_coverage"), sum(pstat == "incomplete_coverage"),
                if (any(iscored & !is.na(sw))) sprintf("%.2f", mean(sw[iscored & !is.na(sw)]))
                else "NA", occ_str))
  }
}
if (!length(score_rows))
  stop("no mapped primer rows to score — step 1 [1.6] should have caught this", call. = FALSE)
score_df <- do.call(rbind, score_rows)

## Fault-tolerance contract: PARTIAL failure = continue (NA scores),
## TOTAL failure of the strand-correct calls = environment problem -> stop.
## The [2.2] canary already proved PM runs on synthetic data in both
## orientations, so total failure means something is wrong at scale —
## inspect the SOFT NOTEs above (full error text + traceback).
if (n_pm_calls_primary > 0L && n_pm_fail >= n_pm_calls_primary)
  stop("every strand-correct PM call failed — environment/install problem, not data; see the SOFT NOTEs above (error text + traceback) and fix PrimerMiner before re-running",
       call. = FALSE)

## ---- [2.4] assay x template status (five-class precedence) ----------------------------
cat("\n[2.4] assay x template status\n")
status_rows <- vector("list", 0)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  f2 <- TEMPLATES[[i]]
  nT <- f2$records
  pg <- pair_table[pair_table$gene == G, ]
  for (k in seq_len(nrow(pg))) {
    p <- pg$assay[k]
    fname <- paste0(p, "_F"); rname <- paste0(p, "_R")
    fsub <- score_df[score_df$name == fname & score_df$Taxa == T, ]
    rsub <- score_df[score_df$name == rname & score_df$Taxa == T, ]
    f_mapped <- nrow(fsub) > 0L
    r_mapped <- nrow(rsub) > 0L
    if (f_mapped)
      check(sprintf("[2.4] %s %s %s_F record order", T, G, p),
            identical(fsub$seq_id, f2$names), "primer score rows out of MSA record order")
    if (r_mapped)
      check(sprintf("[2.4] %s %s %s_R record order", T, G, p),
            identical(rsub$seq_id, f2$names), "primer score rows out of MSA record order")
    fscore <- rep(NA_real_, nT); rscore <- rep(NA_real_, nT)
    fstat  <- rep(NA_character_, nT); rstat <- rep(NA_character_, nT)
    if (f_mapped) { fscore[] <- fsub$score; fstat[] <- fsub$p_status }
    if (r_mapped) { rscore[] <- rsub$score; rstat[] <- rsub$p_status }
    a_status <- character(nT)
    if (!f_mapped || !r_mapped) {
      a_status[] <- "unmapped"
    } else {
      a_status[] <- "unsuitable"
      for (j in seq_len(nT)) {
        if (fstat[j] == "no_coverage" || rstat[j] == "no_coverage")
          a_status[j] <- "no_coverage"
        else if (fstat[j] == "incomplete_coverage" || rstat[j] == "incomplete_coverage")
          a_status[j] <- "incomplete_coverage"
        else {
          tot <- fscore[j] + rscore[j]
          if (!is.na(tot) && tot <= CONFIG$SCORE_THRESHOLD) a_status[j] <- "suitable"
          # else stays unsuitable (NA total -> unsuitable)
        }
      }
    }
    status_rows[[length(status_rows) + 1L]] <- data.frame(
      Taxa = T, Gene = G, assay = p, seq_id = f2$names, a_status = a_status,
      f_score = fscore, r_score = rscore, f_p_status = fstat, r_p_status = rstat,
      stringsAsFactors = FALSE)
  }
  n_block <- sum(vapply(status_rows,
                        function(d) sum(d$Taxa == T & d$Gene == G), integer(1)))
  check(sprintf("[2.4] %s %s assay rows = n_templates x n_pairs", T, G),
        n_block == nT * nrow(pg),
        paste0("rows ", n_block, " vs expected ", nT * nrow(pg)))
}
status_df <- do.call(rbind, status_rows)

## ---- [2.5] per taxon x assay class counts ----------------------------------------
cat("\n[2.5] per taxon x assay class counts\n")
for (Tx in unique(taxon$Taxa)) {
  sub <- status_df[status_df$Taxa == Tx, ]
  cat(sprintf("  %s:\n", Tx))
  for (u in sort(unique(sub$Gene))) {
    gu <- sub[sub$Gene == u, ]
    for (p in unique(gu$assay)) {
      a <- gu$a_status[gu$assay == p]
      cat(sprintf("    %-4s %-18s suitable=%-4d unsuitable=%-4d no_cov=%-4d incomp=%-4d unmapped=%-4d (n=%d)\n",
                  u, p, sum(a == "suitable"), sum(a == "unsuitable"),
                  sum(a == "no_coverage"), sum(a == "incomplete_coverage"),
                  sum(a == "unmapped"), length(a)))
    }
  }
}

## ---- [2.6] NA scores ---------------------------------------------------------------
cat("\n[2.6] NA scores (NA = unsuitable, counted + named)\n")
na_rows <- score_df[score_df$p_status == "scored" & is.na(score_df$score), ]
n_na <- nrow(na_rows)
if (n_na > 0L)
  note("[2.6]", sprintf("%d NA primer score(s), treated as unsuitable: %s",
                        n_na,
                        paste(sprintf("%s/%s/%s", na_rows$Taxa, na_rows$name, na_rows$seq_id),
                              collapse = ", "))) else
                                note("[2.6]", "0 NA scores")

## PM call failures (soft) — separate from legitimate PM NA penalties.
if (n_pm_fail > 0L)
  note("[2.6]", sprintf("%d strand-correct PM call failure(s) (SOFT; affected templates scored NA -> unsuitable): see the SOFT NOTEs above for the full error text + traceback",
                        n_pm_fail))
if (n_pm_rc_fail > 0L)
  note("[2.6]", sprintf("%d diagnostic rc PM call failure(s) (SOFT; score_rc NA only — never classified)", n_pm_rc_fail))

## Unmapped strands carry NA_character_ p_status; guard with !is.na() first
## or NA propagates into `scd` and any() -> NA, which check() (isTRUE)
## reads as FAIL. Unmapped strands are outside the NA-score rule by
## definition (the rule applies to a strand that is present).
f_na <- !is.na(status_df$f_p_status) & status_df$f_p_status == "scored" & is.na(status_df$f_score)
r_na <- !is.na(status_df$r_p_status) & status_df$r_p_status == "scored" & is.na(status_df$r_score)
scd  <- f_na | r_na


check("[2.6] no scored-NA row is suitable", !any(status_df$a_status[scd] == "suitable"),
      "an NA PM score produced a suitable classification")

## ---- [2.7] provenance --------------------------------------------------------------------
cat("\n[2.7] provenance\n")
suit <- status_df[status_df$a_status == "suitable", ]
n_boundary <- sum(suit$f_score + suit$r_score == CONFIG$SCORE_THRESHOLD)
cat(sprintf("  SCORE_THRESHOLD      = %s  (CONFIG; manifest-verified in [2.1])\n",
            CONFIG$SCORE_THRESHOLD))
cat(sprintf("  PRIMERMINER_ADJACENT = %s  (CONFIG)\n",
            CONFIG$PRIMERMINER_ADJACENT))
cat(sprintf("  MAX_AMP_SPAN         = %s  (CONFIG; manifest-verified in [2.1])\n",
            CONFIG$MAX_AMP_SPAN))
cat("  scoring rule         : F -> forward=T | R -> forward=F (primer as written; PM rev-comps the window internally for R)\n")
cat("                         -> 3' end carries the high position penalty on BOTH strands\n")
cat("                         score_rc = diagnostic alt call (RC primer, forward flipped; same mismatches, 3'/5' weights swapped) — reported, NEVER classified\n")
cat(sprintf("  fault tolerance      : %d of %d strand-correct PM calls failed (SOFT; NA -> unsuitable) | %d diagnostic rc failures (loop continued)\n",
            n_pm_fail, n_pm_calls_primary, n_pm_rc_fail))
cat(sprintf("  rule: total = score_F + score_R; suitable iff total <= %s\n",
            CONFIG$SCORE_THRESHOLD))
cat(sprintf("  consensus rules (mapping): %d template(s) gaps_can_win | %d gaps_never_win (never enter classification)\n",
            n_tpl_gcw, n_tpl_gnw))
cat(sprintf("  span gate (mapping): %d gated | %d no_op | %d rejected (never enter classification)\n",
            n_span_gated, n_span_noop, n_span_rejected))
cat(sprintf("  win_min_occ = 0 mapped rows: %d (gaps_can_win holes; rule-aware invariant CHECKed in [2.1])\n",
            n_win0))
cat(sprintf("  suitable combos exactly on the boundary (total == %s): %d\n",
            CONFIG$SCORE_THRESHOLD, n_boundary))

## ---- [2.8] outputs --------------------------------------------------------------------------
cat("\n[2.8] outputs (results/)\n")
SCORE_OUT  <- file.path(OUT, "05_primer_scores.csv")
STATUS_OUT <- file.path(OUT, "06_assay_template_status.csv")
write.csv(score_df, SCORE_OUT, row.names = FALSE)
write.csv(status_df, STATUS_OUT, row.names = FALSE)
score_chk  <- read.csv(SCORE_OUT, stringsAsFactors = FALSE)
status_chk <- read.csv(STATUS_OUT, stringsAsFactors = FALSE)
check("[2.8] primer scores round-trip schema",
      identical(names(score_chk), c("Taxa", "Gene", "assay", "name", "strand", "seq",
                                    "seq_id", "primer_len", "degapped_len", "p_status",
                                    "score_written", "score_rc", "score")),
      paste0("columns found: ", paste(names(score_chk), collapse = ", ")))
check("[2.8] primer scores round-trip values",
      nrow(score_chk) == nrow(score_df) &&
        identical(score_chk$seq_id, score_df$seq_id) &&
        identical(score_chk$p_status, score_df$p_status) &&
        identical(is.na(score_chk$score), is.na(score_df$score)) &&
        identical(as.numeric(score_chk$degapped_len), as.numeric(score_df$degapped_len)) &&
        all(abs(score_chk$score[!is.na(score_chk$score)] -
                  score_df$score[!is.na(score_df$score)]) < 1e-9),
      "re-read primer scores differ from the in-memory table")
## Invariant: score IS the strand-correct written score (score_written);
## any other combination rule would be caught here.
idx_w <- !(is.na(score_chk$score) | is.na(score_chk$score_written))
check("[2.8] score equals strand-correct score_written",
      all(abs(score_chk$score[idx_w] - score_chk$score_written[idx_w]) < 1e-9),
      "score differs from score_written — the retired min(written, RC) rule may have re-entered")
check("[2.8] status round-trip schema",
      identical(names(status_chk), c("Taxa", "Gene", "assay", "seq_id", "a_status",
                                     "f_score", "r_score", "f_p_status", "r_p_status")),
      paste0("columns found: ", paste(names(status_chk), collapse = ", ")))
check("[2.8] status round-trip values",
      nrow(status_chk) == nrow(status_df) &&
        identical(status_chk$a_status, status_df$a_status) &&
        identical(is.na(status_chk$f_score), is.na(status_df$f_score)),
      "re-read assay status differs")
cat(sprintf("    05_primer_scores.csv          %d rows (scored=%d, no_coverage=%d, incomplete=%d)\n",
            nrow(score_df), sum(score_df$p_status == "scored"),
            sum(score_df$p_status == "no_coverage"), sum(score_df$p_status == "incomplete_coverage")))
cat(sprintf("    06_assay_template_status.csv  %d rows (suitable=%d, unsuitable=%d, no_cov=%d, incomp=%d, unmapped=%d)\n",
            nrow(status_df), sum(status_df$a_status == "suitable"),
            sum(status_df$a_status == "unsuitable"), sum(status_df$a_status == "no_coverage"),
            sum(status_df$a_status == "incomplete_coverage"), sum(status_df$a_status == "unmapped")))

## ---- [2.9] audit log ----------------------------------------------------------------------------
log_row <- log_row_make(run_id, "2", "PASS", BASE, parent_run = s1_latest$run_id,
                        metrics = paste0("n_score_rows=", nrow(score_df),
                                         ";n_status_rows=", nrow(status_df),
                                         ";n_suitable=", sum(status_df$a_status == "suitable"),
                                         ";n_unsuitable=", sum(status_df$a_status == "unsuitable"),
                                         ";n_unmapped=", sum(status_df$a_status == "unmapped"),
                                         ";n_na_scores=", n_na,
                                         ";score_threshold=", CONFIG$SCORE_THRESHOLD,
                                         ";pm_adjacent=", CONFIG$PRIMERMINER_ADJACENT,
                                         ";max_pair_span=", CONFIG$MAX_AMP_SPAN,
                                         ";scoring_rule=strand_correct",
                                         ";n_pm_calls=", n_pm_calls_primary,
                                         ";n_pm_call_failures=", n_pm_fail,
                                         ";n_pm_rc_failures=", n_pm_rc_fail,
                                         ";n_templates_gaps_can_win=", n_tpl_gcw,
                                         ";n_win_min_occ_zero_rows=", n_win0))
log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)
log_chk <- read.csv(LOG, stringsAsFactors = FALSE)
check("[2.9] audit log appended",
      any(log_chk$run_id == run_id & log_chk$step == "2"),
      "log row missing after write")

cat("\n============================================================================\n")
cat(sprintf("STEP 2 PASSED | run_id %s (continued from step 1 %s | step 0 %s)\n",
            run_id, s1_latest$run_id, parent_run))
cat(sprintf("  scoring     : strand-correct (F -> forward=T, R -> forward=F; primer as written) — 3' mismatches penalised identically on both strands\n"))
cat(sprintf("  loop        : fault-tolerant — PM call failures: %d strand-correct | %d diagnostic (SOFT; NA -> unsuitable; loop continued)\n",
            n_pm_fail, n_pm_rc_fail))
cat(sprintf("  rules     : %d template(s) gaps_can_win | %d gaps_never_win (manifest/CONFIG; never enter classification)\n",
            n_tpl_gcw, n_tpl_gnw))
cat(sprintf("  span gate : MAX_AMP_SPAN = %d | gated=%d | no_op=%d | rejected=%d (never enter classification)\n",
            CONFIG$MAX_AMP_SPAN, n_span_gated, n_span_noop, n_span_rejected))
cat(sprintf("  primer-level rows : %d (scored=%d, no_coverage=%d, incomplete=%d, NA scores=%d)\n",
            nrow(score_df), sum(score_df$p_status == "scored"),
            sum(score_df$p_status == "no_coverage"), sum(score_df$p_status == "incomplete_coverage"),
            n_na))
cat(sprintf("  assay-level rows  : %d (suitable=%d, unsuitable=%d, no_cov=%d, incomp=%d, unmapped=%d)\n",
            nrow(status_df), sum(status_df$a_status == "suitable"),
            sum(status_df$a_status == "unsuitable"), sum(status_df$a_status == "no_coverage"),
            sum(status_df$a_status == "incomplete_coverage"), sum(status_df$a_status == "unmapped")))
cat("Outputs: results/05_primer_scores.csv | 06_assay_template_status.csv | run_audit_log.csv\n")
cat("Next: step3_amplicon_extract.r (Windows R — fresh session fine)\n")
cat("============================================================================\n")
