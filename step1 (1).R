#!/usr/bin/env Rscript
# ============================================================================
# step1_primer_mapping.r — STEP 1: PRIMER MAPPING  (Windows R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
#
# WHAT THIS STEP DOES (and only this)
#   For every (Taxa, Gene) MSA from the step 0 directory:
#     1. computes the population-level consensus with the consensus engine
#        selected PER TEMPLATE (from the manifest published by step 0):
#          gaps_never_win (base_consensus): the most frequent NON-GAP symbol
#          per column wins; gaps never win; a tie among non-gap symbols ->
#          IUPAC code of their union; zero-coverage column -> '-';
#          gaps_can_win (gap_inclusive_consensus): the most frequent symbol
#          overall (bases + IUPAC + gap) wins; a gap-majority or gap-tied
#          column -> '-' (a hole); a base-only tie -> IUPAC union;
#          zero-coverage column -> '-'.
#          Per-column occupancy (base count) is retained for window reporting.
#     2. matches each primer of that gene against the gap-free consensus
#        with set-based IUPAC overlap (a position mismatches ONLY when the
#        two base sets are disjoint — degenerate bases never mismatch when
#        they overlap), tolerance = CONFIG$MAX_MAP_MISMATCHES,
#     3. maps the chosen window back to GAPPED MSA column coordinates
#        (start/end are MSA positions, so step 2 slices [start, end]
#        directly out of the alignment),
#     4. reports every unmapped primer by name with its min_mismatches
#        (how close it came) — a documented result, never an error.
#
#   Occupancy reporting: every mapped row carries win_min_occ (minimum
#   per-column base count across the gapped window span [start, end]),
#   win_n_low_occ (window-span columns with occ < CONFIG$MIN_OCC), and the
#   consensus_rule that produced the row. NOTE-level only — it never excludes
#   a column, window, or mapping result. Under gaps_never_win a consensus '-'
#   column exists only at occ = 0 (zero coverage), so a mapped window with
#   win_min_occ < 1 is a stop-worthy structural anomaly (checked in [1.7]);
#   under gaps_can_win holes are legitimate, so win_min_occ = 0 is allowed
#   and NOTE-reported.
#
#   Nothing in this file refers to any previous mapping artifact, and no
#   expected count is hardcoded — the audits are structural invariants, so
#   any primer set, any genes, any MSA dimensions run through.
#
#   Pair-aware span gate: after the per-primer matching (pass 1), a second
#   pass gates each complete pair: when one strand's anchor is exact
#   (0 mm = high-confidence), the partner primer's candidate windows are
#   restricted to windows that FACE that anchor, are non-overlapping with it,
#   and sit at a de-gapped (clean) consensus span of at most
#   CONFIG$MAX_AMP_SPAN bp. Survivors are re-ranked under the existing
#   tie-break (mm -> concrete support -> leftmost). An empty kept set is a
#   documented span-rejection (an unmapped row; the distant window is never
#   scored). Both-exact pairs out of span are gated one-shot against each
#   other's exact windows. The gate is orthogonal to the consensus rule (all
#   gate math on clean coordinates) and never loosens MAX_MAP_MISMATCHES.
#   Every row carries span_gate (ungated/no_op/gated/rejected) + span_gate_mm
#   (mm of the displaced pre-gate anchor; NA unless gated/rejected).
#   04_pair_coverage.csv includes amplicon_size_clean. Steps 3-5 are
#   unaffected by the gate.
#
# INPUTS (all must match the step 0 manifest — a file change stops the run):
#   results/00_run_manifest.csv (must contain the consensus_rule section and
#   config|max_pair_span; an older manifest without these sections triggers a
#   NOTE + CONFIG fallback or stops at the version gate),
#   results/01_pair_table.csv,
#   primer CSV, taxon directory CSV, the MSA FASTA files
#
# OUTPUTS (results/, overwritten each run):
#   03_primer_mapping.csv   one row per (primer, Taxa, Gene): start, end
#                           (NA = unmapped), target_with_gaps, n_mismatches,
#                           min_mismatches, map_orientation,
#                           other_orientation_mm, map_mode, win_min_occ,
#                           win_n_low_occ, consensus_rule, span_gate,
#                           span_gate_mm
#   04_pair_coverage.csv    one row per (assay, Taxa): F/R mapped, cov,
#                           amplicon_size, orientation_ok,
#                           amplicon_size_clean
#   run_audit_log.csv       appended (step = 1)
#
# AUDITS ([1.x]; CHECK stops the run, NOTE / printed results never stop):
#   [1.1] continuity gate: manifest exists + PASS, BASE + config unchanged,
#          utils_version vs this session; max_pair_span unchanged; consensus-
#          rule contract (CONFIG default + overrides in vocabulary; manifest
#          config rows must equal CONFIG if the section is present)
#   [1.2] input files unchanged vs manifest (row counts)
#   [1.3] per MSA: rule selection (manifest row unique, in vocabulary, ==
#          CONFIG-derived) + reload + fingerprint vs manifest + ragged guard
#          + consensus via the selected rule
#   [1.4] mapping loop, TWO PASS: pass 1 = per-primer matching retaining
#          candidate frames; pass 2 = span gate per complete pair (exact
#          anchor gates the partner; no-op / re-rank / documented rejection).
#          One printed line per primer with the FINAL anchor. Occupancy and
#          hole counts NOTE-reported, never excluded.
#   [1.5] unmapped primers named per taxon x gene (min_mm reported)
#   [1.6] pair coverage per taxon; CHECK >= 1 full combo; amplicon_size_clean
#          column; span-gate backstop CHECK + both-tolerant long-span NOTE
#   [1.7] mapped-row invariants: every row re-derives (span, stripped length,
#          round-trip mm) + rule-aware win_min_occ + span-gate re-derivations
#   [1.8] outputs written + re-read (round-trip; 17-column mapping schema +
#          9-column coverage schema)
#   [1.9] audit log appended
#
# RUN:  Rscript step1_primer_mapping.r   (Windows R, fresh or same session;
#       step 0 must have passed in the same results/ directory)
#       Paste-safe: all conditionals braced; one CHUNK break marked.
# ============================================================================

## ---- BASE — the ONE line to edit (same line as step 0; verified vs manifest)
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
OUT <- file.path(BASE, "results")

## ---- [1.1] continuity gate: step 0 manifest ---------------------------------
cat("============================================================================\n")
cat("STEP 1 - PRIMER MAPPING (set-based IUPAC overlap; consensus rule per template; pair-aware span gate)\n")
cat("  run_id  :", run_id, "\n")
mani_path <- file.path(OUT, "00_run_manifest.csv")
check("[1.1] step 0 manifest exists", file.exists(mani_path),
      "results/00_run_manifest.csv missing — run step0_setup_audit.r first")
mani <- read.csv(mani_path, stringsAsFactors = FALSE)
mval <- function(sec, item) {
  v <- mani$value[mani$section == sec & mani$item == item]
  if (length(v) != 1L) stop("manifest missing row: ", sec, "|", item, call. = FALSE)
  v
}
check("[1.1] step 0 status = PASS", mval("status", "step0") == "PASS",
      "step 0 did not pass — fix it before mapping")
parent_run <- mval("run", "run_id")
check("[1.1] BASE unchanged since step 0", identical(mval("run", "base"), BASE),
      paste0("manifest BASE: ", mval("run", "base"), " | session BASE: ", BASE,
             " — if you moved things, re-run step 0"))
check("[1.1] max_map_mismatches unchanged since step 0",
      identical(mval("config", "max_map_mismatches"), as.character(CONFIG$MAX_MAP_MISMATCHES)),
      "CONFIG$MAX_MAP_MISMATCHES differs from the manifest — set it, then re-run step 0")
check("[1.1] score_threshold unchanged since step 0",
      identical(mval("config", "score_threshold"), as.character(CONFIG$SCORE_THRESHOLD)),
      "CONFIG$SCORE_THRESHOLD differs from the manifest — set it, then re-run step 0")

## utils continuity: the utils_version published by step 0 must equal the
## utils.R sourced in THIS session. A mismatch means the core behavior
## changed after the manifest was written; stop and re-run step 0.
## Older manifests without the row -> NOTE instead of FAIL (tolerant).
uv_row <- mani$value[mani$section == "run" & mani$item == "utils_version"]
if (length(uv_row) == 1L) {
  check("[1.1] utils_version matches manifest", uv_row == UTILS_VERSION,
        paste0("manifest utils_version is '", uv_row,
               "', this session sources utils.R '", UTILS_VERSION,
               "' — re-run step 0 to republish the manifest"))
} else {
  note("[1.1]", sprintf("manifest has no utils_version row; this session sources utils.R %s",
                        UTILS_VERSION))
}

## ---- [1.1] span-gate parameter ----------------------------------------------
## The utils_version gate above runs first and stops on a manifest written
## before the span-gate feature existed, so mval will always find the row here.
check("[1.1] max_pair_span unchanged since step 0",
      identical(mval("config", "max_pair_span"), as.character(CONFIG$MAX_AMP_SPAN)),
      "CONFIG$MAX_AMP_SPAN differs from the manifest — set it, then re-run step 0")

## ---- [1.1] consensus-rule contract -------------------------------------------
## Per-template rule source: the consensus_rule manifest section published by
## step 0:
##   config       | consensus_rule            -> global default rule
##   config       | consensus_rule_overrides  -> packed "Taxa|Gene=rule;..."
##   consensus_rule| "Taxa|Gene"              -> the rule for that template
## Section absent (older manifest) -> NOTE + the CONFIG-derived rule per
## template (consensus_rule_for) — tolerant. Section present -> the manifest
## config rows must equal this session's CONFIG (drift = stop, re-run step 0);
## the per-template rows are gated in [1.3] where the templates are enumerated.
CFG_DEFAULT <- CONFIG$CONSENSUS_RULE
CFG_OVERR   <- CONFIG$CONSENSUS_RULE_OVERRIDES
check("[1.1] CONFIG consensus default in vocabulary",
      CFG_DEFAULT %in% CONSENSUS_RULES_VALID,
      paste0("CONFIG$CONSENSUS_RULE '", CFG_DEFAULT, "' is not a valid rule (valid: ",
             paste(CONSENSUS_RULES_VALID, collapse = "/"), ")"))
check("[1.1] CONFIG override values in vocabulary",
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
  check("[1.1] manifest consensus_rule default == CONFIG",
        length(mani_rule_default) == 1L && mani_rule_default[1L] == CFG_DEFAULT,
        paste0("manifest config|consensus_rule = '",
               if (length(mani_rule_default) == 1L) mani_rule_default[1L] else "(missing)",
               "' vs CONFIG '", CFG_DEFAULT, "' — re-run step 0"))
  mani_ov <- parse_packed_overrides(
    if (length(mani_rule_ov) == 1L) mani_rule_ov[1L] else "(none)")
  check("[1.1] manifest overrides == CONFIG (key = value, set-equal)",
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
  note("[1.1]", "manifest has no consensus_rule section — falling back to CONFIG-derived rules (consensus_rule_for) for every template")
}

cat("  parent  : step 0 run ", parent_run, " (manifest)\n")
cat("  base    :", BASE, "\n")
cat("  config  : max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | score_threshold =", CONFIG$SCORE_THRESHOLD, "\n")
cat("  span gate  : max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
cat("  consensus_rule : default =", CFG_DEFAULT,
    "| overrides =",
    if (length(CFG_OVERR) == 0L) "(none)"
    else paste(sprintf("%s=%s", names(CFG_OVERR), unlist(CFG_OVERR)), collapse = "; "), "\n")
cat("  manifest   : consensus_rule section",
    if (rule_section_present) "present" else "ABSENT (older step 0 manifest; CONFIG fallback)", "\n")
cat("============================================================================\n")

## ---- [1.2] input files unchanged vs manifest ---------------------------------
cat("\n[1.2] input files vs manifest\n")
prim <- read.csv(file.path(BASE, "primers_input2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
prim[] <- lapply(prim, trimws)
taxon <- read.csv(file.path(BASE, "taxon_template_directory_2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
taxon[] <- lapply(taxon, trimws)
taxon$full_path <- file.path(BASE, taxon$filepath)
pair_table <- read.csv(file.path(OUT, "01_pair_table.csv"), stringsAsFactors = FALSE)
check("[1.2] primers_input2.csv unchanged",
      nrow(prim) == as.integer(mval("inputs", "primer_rows")),
      paste0("manifest has ", mval("inputs", "primer_rows"), " rows, file has ", nrow(prim),
             " — re-run step 0 after changing inputs"))
check("[1.2] taxon_template_directory_2.csv unchanged",
      nrow(taxon) == as.integer(mval("inputs", "taxon_rows")),
      "directory row count differs from manifest — re-run step 0")
check("[1.2] pair table consistent",
      nrow(pair_table) == as.integer(mval("inputs", "pairs_complete")),
      "results/01_pair_table.csv row count differs from manifest — re-run step 0")
prim$strand <- substr(prim$name, nchar(prim$name), nchar(prim$name))
genes <- sort(unique(prim$gene))
IUPAC_CODES <- setdiff(names(BASE_SET), c("A", "C", "G", "T"))

## ---- [1.3] per-MSA reload, fingerprint check, rule selection, consensus ------
cat("\n[1.3] MSA reload + fingerprint vs manifest + consensus (rule per template)\n")
CONS <- list(); CONS_CHARS <- list(); NG <- list(); CLEAN <- list(); OCC <- list(); RULE <- list()
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
  key_t <- paste0(T, "_", G)

  ## ---- rule selection: manifest row vs CONFIG-derived -----------------------
  if (rule_section_present) {
    rr <- rule_rows$value[rule_rows$item == key]
    check(sprintf("[1.3] %s has exactly one manifest rule row", key), length(rr) == 1L,
          sprintf("expected exactly one consensus_rule row for '%s' in the manifest, found %d — re-run step 0", key, length(rr)))
    check(sprintf("[1.3] %s manifest rule in vocabulary", key), rr[1L] %in% CONSENSUS_RULES_VALID,
          paste0("manifest rule '", rr[1L], "' for '", key, "' is not a valid rule — re-run step 0"))
    rule_expected <- tryCatch(consensus_rule_for(key, CFG_OVERR, CFG_DEFAULT),
                              error = function(e) stop("consensus_rule_for failed for '", key, "': ",
                                                       conditionMessage(e), call. = FALSE))
    check(sprintf("[1.3] %s manifest rule == CONFIG-derived rule", key), rr[1L] == rule_expected,
          paste0("manifest says '", rr[1L], "', CONFIG derives '", rule_expected,
                 "' for '", key, "' — drift; re-run step 0"))
    rule <- rr[1L]
  } else {
    rule <- tryCatch(consensus_rule_for(key, CFG_OVERR, CFG_DEFAULT),
                     error = function(e) stop("consensus_rule_for failed for '", key, "': ",
                                              conditionMessage(e), call. = FALSE))
  }

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
  check(sprintf("[1.3] %s %s unchanged", T, G), length(diffs) == 0L,
        paste(diffs, collapse = " ; "))
  bc    <- consensus(f2$seqs, rule)   # dispatcher: rule selects the engine
  cons  <- bc$cons
  chrs  <- strsplit(cons, "", fixed = TRUE)[[1L]]
  ng    <- gap_positions(cons)
  CONS[[key_t]]       <- cons
  OCC[[key_t]]        <- bc$occ      # per-column base count (occupancy)
  CONS_CHARS[[key_t]] <- chrs
  NG[[key_t]]         <- ng
  RULE[[key_t]]       <- rule
  CLEAN[[key_t]]      <- paste(chrs[ng], collapse = "")
  cat(sprintf("    %-14s %-4s records=%-5d width=%-5d rule=%-14s | consensus gaps=%-5d iupac_cols=%-5d clean=%-5d\n",
              T, G, f2$records, f2$width[1L], rule,
              length(chrs) - length(ng), sum(chrs %in% IUPAC_CODES), length(ng)))
}

## ---- [1.4] span-gate machinery -------------------------------------------------
## gate_strand: gate ONE strand's candidate frames against the partner's exact
## (0 mm) anchor window [a_start, a_end] (clean coords). `rec` is the primer's
## pass-1 record and is MAPPED (the gate never touches a tolerance rejection);
## side = "F" (rec is the F strand; the anchor is the partner R window) or
## "R" (rec is the R strand; the anchor is the partner F window). Returns the
## record with span_gate set (+ span_gate_mm, and a rebuilt row when the
## anchor moved or the primer was rejected); prints the canonical [1.4] NOTE
## on gated/rejected events.
##
## No-op fast path: if the pre-gate winner already faces the anchor
## non-overlapping within the span bound (span_keep TRUE), the in-span kept
## set contains the global winner of the TOTAL tie-break order (mm ->
## concrete support -> leftmost), so the re-rank is provably identical and
## the pre-gate row stands unchanged with span_gate = "no_op".
##
## Gated: the gated winner differs -> the row is rebuilt from the gated
## winner's clean start (st <- ng[cs], en <- ng[cs + Lp - 1], re-sliced
## target, re-derived occupancy flags); span_gate_mm records the mm of the
## displaced pre-gate anchor.
##
## Rejected: both orientation frames empty after gating -> the standard
## unmapped row (NA coords, min_mismatches = the global minimum,
## map_mode = "none"); the distant window is never scored.
gate_strand <- function(rec, side, a_start, a_end, p, partner_name, scope) {
  T <- scope$T; G <- scope$G
  ng <- scope$ng; clean <- scope$clean; chrs <- scope$chrs
  occ <- scope$occ; rule <- scope$rule
  Lp <- rec$Lp
  ind <- strrep(" ", 18L)   # NOTE continuation indent: aligns under note()'s tag prefix

  ## no-op fast path (provably identical result; see the block above)
  if (!is.na(rec$clean_s) &&
      isTRUE(span_keep(rec$clean_s, Lp, a_start, a_end, side, CONFIG$MAX_AMP_SPAN))) {
    rec$sg <- "no_op"
    return(rec)
  }
  

  ## gate each orientation frame (F: one frame; R: as-written + RC)
  ga <- gated_best_match(rec$cand_as, rec$seq, clean, a_start, a_end, side, CONFIG$MAX_AMP_SPAN)
  gr <- if (side == "R")
    gated_best_match(rec$cand_rc, rec$ps_rc, clean, a_start, a_end, side, CONFIG$MAX_AMP_SPAN)
  else NULL

  ## orientation choice on the GATED per-orientation winners — the existing
  ## rule: fewer mm; exact tie -> prefer RC; both empty -> unmapped
  if (side == "F") {
    orient <- "as_written"; winner <- ga$best; other <- NA_integer_
  } else {
    as_score <- if (nrow(ga$best) == 1L) ga$best$n_mismatches[1L] else Inf
    rc_score <- if (nrow(gr$best) == 1L) gr$best$n_mismatches[1L] else Inf
    if (rc_score < as_score) {
      use_rc <- TRUE
    } else if (as_score < rc_score) {
      use_rc <- FALSE
    } else if (is.finite(rc_score)) {
      use_rc <- TRUE          # tie on a real match: prefer RC
    } else {
      use_rc <- FALSE         # both empty (rejection)
    }
    orient <- if (use_rc) "rc" else "as_written"
    winner <- if (use_rc) gr$best else ga$best
    other  <- if (use_rc)
      (if (nrow(ga$best) == 1L) ga$best$n_mismatches[1L] else NA_integer_)
      else (if (nrow(gr$best) == 1L) gr$best$n_mismatches[1L] else NA_integer_)
  }

  if (nrow(winner) == 0L) {
    ## span rejection: standard unmapped row; the distant window is never scored
    cand <- if (side == "R") c(rec$mm_as, rec$mm_rc) else c(rec$mm_as)
    mm_report <- if (all(is.na(cand))) NA_integer_ else min(cand, na.rm = TRUE)
    old_span <- if (side == "F") a_end - rec$clean_s + 1L
                else (rec$clean_s + Lp - 1L) - a_start + 1L
    note("[1.4]",
      paste0(sprintf("%s (%s x %s): %s span-REJECTED vs exact partner\n",
                     p, T, G, rec$name),
             ind, sprintf("%s (0 mm @ clean %d-%d): best passing window %d mm @ clean %d\n",
                          partner_name, a_start, a_end, rec$mm_u, rec$clean_s),
             ind, sprintf("(facing span %d bp) has no in-span alternative (MAX_AMP_SPAN=%d)",
                          old_span, CONFIG$MAX_AMP_SPAN)))
    rec$sg   <- "rejected"
    rec$sgmm <- rec$mm_u
    rec$row <- data.frame(
      name = rec$name, seq = rec$seq, gene = G, Taxa = T,
      start = NA_integer_, end = NA_integer_, target_with_gaps = NA_character_,
      n_mismatches = NA_integer_, min_mismatches = mm_report,
      map_orientation = "as_written", other_orientation_mm = NA_integer_,
      map_mode = "none", win_min_occ = NA_integer_, win_n_low_occ = NA_integer_,
      consensus_rule = rule, span_gate = "rejected", span_gate_mm = rec$mm_u,
      stringsAsFactors = FALSE)
    return(rec)
  }

  ## gated: rebuild the row from the gated winner's clean start
  gcs <- winner$start[1L]; gmm <- winner$n_mismatches[1L]
  st  <- ng[gcs]; en <- ng[gcs + Lp - 1L]
  tgt <- paste(chrs[st:en], collapse = "")
  wocc <- occ[st:en]
  old_span <- if (side == "F") a_end - rec$clean_s + 1L
              else (rec$clean_s + Lp - 1L) - a_start + 1L
  new_span <- if (side == "F") a_end - gcs + 1L
              else (gcs + Lp - 1L) - a_start + 1L
  span_txt <- sprintf("(facing span %d bp > MAX_AMP_SPAN=%d)", old_span, CONFIG$MAX_AMP_SPAN)
  if (old_span <= CONFIG$MAX_AMP_SPAN)
    span_txt <- sprintf("(facing span %d bp; window overlaps/anti-faces the exact anchor; MAX_AMP_SPAN=%d)",
                        old_span, CONFIG$MAX_AMP_SPAN)
  note("[1.4]",
    paste0(sprintf("%s (%s x %s): %s span-GATED vs exact partner\n",
                   p, T, G, rec$name),
           ind, sprintf("%s (0 mm @ clean %d-%d): pre-gate anchor %d mm @ clean %d\n",
                        partner_name, a_start, a_end, rec$mm_u, rec$clean_s),
           ind, paste0(span_txt, " displaced -> ", gmm, " mm @ clean ", gcs, "\n"),
           ind, sprintf("(facing span %d bp)", new_span)))
  rec$sg   <- "gated"
  rec$sgmm <- rec$mm_u
  rec$disp_clean_s <- rec$clean_s
  rec$row <- data.frame(
    name = rec$name, seq = rec$seq, gene = G, Taxa = T,
    start = st, end = en, target_with_gaps = tgt,
    n_mismatches = gmm, min_mismatches = gmm,
    map_orientation = orient, other_orientation_mm = other,
    map_mode = if (gmm == 0L) "exact" else "tolerant",
    win_min_occ = as.integer(min(wocc)), win_n_low_occ = as.integer(sum(wocc < CONFIG$MIN_OCC)),
    consensus_rule = rule, span_gate = "gated", span_gate_mm = rec$mm_u,
    stringsAsFactors = FALSE)
  rec
}

## ---- [1.4] mapping loop (two-pass) ---------------------------------------------
cat("\n[1.4] mapping (two-pass: per-primer matching + pair-aware span gate)\n")
cat("      pass 1: per-primer matching (tolerance = ", CONFIG$MAX_MAP_MISMATCHES,
    "; set-based IUPAC overlap; R primers tested as-written AND reverse-complemented;\n",
    "         best of the two wins — fewer mm, ties prefer RC orientation;\n",
    "         consensus rule per template from the manifest; win_min_occ /\n",
    "         win_n_low_occ vs MIN_OCC = ", CONFIG$MIN_OCC,
    " reported per mapped window, never enforced);\n",
    "      pass 2: per complete pair — when one strand is exact (0 mm), the partner's\n",
    "         candidate windows must FACE that anchor, not overlap it, and sit within\n",
    "         MAX_AMP_SPAN = ", CONFIG$MAX_AMP_SPAN,
    " bp clean span; survivors re-ranked under the existing tie-break;\n",
    "         empty kept set = documented span-rejection (an unmapped row)\n", sep = "")

## ONE global accumulator, appended exactly once per (taxon, primer) in the
## loop below — the rbind after the loop can never see a partial/empty list
## (the [1.4] pre-rbind CHECK below names it if that ever happens).
final_rows <- vector("list", 0)

n_expect <- 0L
for (i in seq_len(nrow(taxon))) n_expect <- n_expect + sum(prim$gene == taxon$Gene[i])

for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  key_t <- paste0(T, "_", G)
  subprim <- prim[prim$gene == G, ]
  check(sprintf("[1.4] %s %s has primers", T, G), nrow(subprim) >= 1L,
        "directory row with a gene that has no primers — step 0 [0.8] should have caught this")
  ng <- NG[[key_t]]; clean <- CLEAN[[key_t]]; chrs <- CONS_CHARS[[key_t]]
  cat(sprintf("  %s x %s (clean width %d, rule=%s):\n", T, G, nchar(clean), RULE[[key_t]]))

  ## ===== pass 1: per-primer matching, candidate frames retained =============
  n_p <- nrow(subprim)
  recs <- vector("list", n_p)
  for (j in seq_len(n_p)) {
    ps     <- subprim$seq[j]
    Lp     <- nchar(ps)
    strand <- subprim$strand[j]

    ## always map as-written
    bm_as  <- best_match(ps, clean, CONFIG$MAX_MAP_MISMATCHES)
    cand_as <- primer_match(ps, clean, CONFIG$MAX_MAP_MISMATCHES)   # retained frame for pass 2
    mm_as  <- primer_min_mismatches(ps, clean)

    ps_rc  <- NA_character_; cand_rc <- NULL; mm_rc <- NA_integer_

    ## R primers: also map the reverse complement
    if (strand == "R") {
      ps_rc  <- reverse_complement(ps)
      bm_rc  <- best_match(ps_rc, clean, CONFIG$MAX_MAP_MISMATCHES)
      cand_rc <- primer_match(ps_rc, clean, CONFIG$MAX_MAP_MISMATCHES)  # retained frame for pass 2
      mm_rc  <- primer_min_mismatches(ps_rc, clean)

      ## score: mismatch count of the chosen window (Inf = unmapped)
      as_score <- if (nrow(bm_as) == 1L) bm_as$n_mismatches[1L] else Inf
      rc_score <- if (nrow(bm_rc) == 1L) bm_rc$n_mismatches[1L] else Inf

      ## pick the winner: fewer mm wins; exact tie -> prefer RC (expected
      ## orientation for a reverse primer); both unmapped -> as_written
      if (rc_score < as_score) {
        use_rc <- TRUE
      } else if (as_score < rc_score) {
        use_rc <- FALSE
      } else if (is.finite(rc_score)) {
        use_rc <- TRUE          # tie on a real match: prefer RC
      } else {
        use_rc <- FALSE         # both Inf (unmapped)
      }

      if (use_rc) {
        bm   <- bm_rc
        mm_w <- if (nrow(bm_rc) == 1L) bm_rc$n_mismatches[1L] else NA_integer_
        mm_other <- if (nrow(bm_as) == 1L) bm_as$n_mismatches[1L] else NA_integer_
        orient <- "rc"
      } else {
        bm   <- bm_as
        mm_w <- if (nrow(bm_as) == 1L) bm_as$n_mismatches[1L] else NA_integer_
        mm_other <- if (nrow(bm_rc) == 1L) bm_rc$n_mismatches[1L] else NA_integer_
        orient <- "as_written"
      }
    } else {
      bm       <- bm_as
      mm_w     <- if (nrow(bm_as) == 1L) bm_as$n_mismatches[1L] else NA_integer_
      mm_other <- NA_integer_
      orient   <- "as_written"
    }

    ## pre-gate winner coordinates (clean space; NA if unmapped)
    mm_u    <- if (nrow(bm) == 1L) bm$n_mismatches[1L] else NA_integer_
    clean_s <- if (nrow(bm) == 1L) bm$start[1L] else NA_integer_

    ## build the row — 17 columns including span_gate / span_gate_mm
    if (nrow(bm) == 1L) {
      st  <- ng[bm$start[1L]]
      en  <- ng[bm$start[1L] + Lp - 1L]
      tgt <- paste(chrs[st:en], collapse = "")
      mode <- if (mm_w == 0L) "exact" else "tolerant"
      ## occupancy flags (NOTE-level only; never excludes):
      ## per-column base count across the gapped window span [st, en]
      wocc          <- OCC[[key_t]][st:en]
      win_min_occ   <- as.integer(min(wocc))
      win_n_low_occ <- as.integer(sum(wocc < CONFIG$MIN_OCC))
      row <- data.frame(
        name = subprim$name[j], seq = ps, gene = G, Taxa = T,
        start = st, end = en, target_with_gaps = tgt,
        n_mismatches = mm_w, min_mismatches = mm_w,
        map_orientation = orient, other_orientation_mm = mm_other,
        map_mode = mode,
        win_min_occ = win_min_occ, win_n_low_occ = win_n_low_occ,
        consensus_rule = RULE[[key_t]],
        span_gate = "ungated", span_gate_mm = NA_integer_,
        stringsAsFactors = FALSE)
    } else {
      ## primer failed to map in the orientation(s) evaluated
      cand <- if (strand == "R") c(mm_as, mm_rc) else c(mm_as)
      mm_report <- if (all(is.na(cand))) NA_integer_ else min(cand, na.rm = TRUE)
      row <- data.frame(
        name = subprim$name[j], seq = ps, gene = G, Taxa = T,
        start = NA_integer_, end = NA_integer_,
        target_with_gaps = NA_character_,
        n_mismatches = NA_integer_, min_mismatches = mm_report,
        map_orientation = orient, other_orientation_mm = NA_integer_,
        map_mode = "none",
        win_min_occ = NA_integer_, win_n_low_occ = NA_integer_,
        consensus_rule = RULE[[key_t]],
        span_gate = "ungated", span_gate_mm = NA_integer_,
        stringsAsFactors = FALSE)
    }

    recs[[j]] <- list(
      name = subprim$name[j], seq = ps, Lp = Lp, strand = strand,
      cand_as = cand_as, cand_rc = cand_rc, ps_rc = ps_rc,
      mm_as = mm_as, mm_rc = mm_rc,
      mm_u = mm_u, clean_s = clean_s,
      sg = "ungated", sgmm = NA_integer_,
      row = row)
  }

  ## ===== pass 2: span gate per complete pair — NOTEs print here ==============
  scope <- list(ng = ng, clean = clean, chrs = chrs, occ = OCC[[key_t]],
                rule = RULE[[key_t]], T = T, G = G)
  for (k in seq_len(nrow(pair_table))) {
    if (pair_table$gene[k] != G) next
    p  <- pair_table$assay[k]
    jf <- which(subprim$name == paste0(p, "_F"))
    jr <- which(subprim$name == paste0(p, "_R"))
    if (length(jf) != 1L || length(jr) != 1L) next   # defensive: step 0 guarantees completeness
    frec <- recs[[jf]]; rrec <- recs[[jr]]
    f_mapped <- !is.na(frec$mm_u); r_mapped <- !is.na(rrec$mm_u)
    f_exact  <- f_mapped && frec$mm_u == 0L
    r_exact  <- r_mapped && rrec$mm_u == 0L

    if (f_exact && r_exact) {
      ## both exact: facing, non-overlapping and within the bound -> no_op for
      ## both; else NOTE + one-shot mutual gate against each other's exact
      ## windows (captured BEFORE either row is replaced)
      fs <- frec$clean_s; fe <- fs + frec$Lp - 1L
      rs <- rrec$clean_s; re <- rs + rrec$Lp - 1L
      sp <- pair_span(fs, fe, rs, re)
      if (fe < rs && sp <= CONFIG$MAX_AMP_SPAN) {
        frec$sg <- "no_op"; rrec$sg <- "no_op"
      } else {
        note("[1.4]",
          sprintf("%s (%s x %s): both strands exact — F 0 mm @ clean %d-%d, R 0 mm @ clean %d-%d, span %d bp (MAX_AMP_SPAN=%d) — gating both strands one-shot",
                  p, T, G, fs, fe, rs, re, sp, CONFIG$MAX_AMP_SPAN))
        frec <- gate_strand(frec, "F", rs, re, p, rrec$name, scope)
        rrec <- gate_strand(rrec, "R", fs, fe, p, frec$name, scope)
      }
    } else if (f_exact) {
      ## F exact: gate R against F's exact window (orientation-independent).
      ## R unmapped for tolerance has an empty candidate set: nothing to
      ## gate — it keeps its pre-gate unmapped row with span_gate = "ungated".
      if (r_mapped)
        rrec <- gate_strand(rrec, "R", frec$clean_s, frec$clean_s + frec$Lp - 1L,
                            p, frec$name, scope)
    } else if (r_exact) {
      ## R exact: gate F against R's exact window (same tolerance rule)
      if (f_mapped)
        frec <- gate_strand(frec, "F", rrec$clean_s, rrec$clean_s + rrec$Lp - 1L,
                            p, rrec$name, scope)
    }
    ## neither exact: no gate; both stay "ungated"
    recs[[jf]] <- frec; recs[[jr]] <- rrec
  }

  ## ===== per-primer print lines (FINAL anchor state after the span gate) =====
  for (j in seq_len(n_p)) {
    r   <- recs[[j]]
    row <- r$row
    if (is.na(row$n_mismatches)) {
      rej_str <- if (r$sg == "rejected") "  [span-rejected: no in-span alternative]" else ""
      if (r$strand == "R")
        cat(sprintf("      %-22s UNMAPPED  min_mm(as=%s, rc=%s)%s\n",
                    row$name,
                    if (is.na(r$mm_as)) "NA" else r$mm_as,
                    if (is.na(r$mm_rc)) "NA" else r$mm_rc, rej_str))
      else
        cat(sprintf("      %-22s UNMAPPED  min_mm(as=%s)%s\n",
                    row$name,
                    if (is.na(r$mm_as)) "NA" else r$mm_as, rej_str))
    } else {
      st <- row$start; en <- row$end
      other_str <- ""
      if (r$strand == "R" && !is.na(row$other_orientation_mm))
        other_str <- sprintf("  [as-written mm=%d]", row$other_orientation_mm)
      occ_str <- ""
      if (row$win_n_low_occ > 0L)
        occ_str <- sprintf("  [occ min=%d, %d col(s) < MIN_OCC=%d]",
                           row$win_min_occ, row$win_n_low_occ, CONFIG$MIN_OCC)
      hole_str <- ""
      if (RULE[[key_t]] == "gaps_can_win") {
        n_holes <- sum(chrs[st:en] == "-")
        if (n_holes > 0L)
          hole_str <- sprintf("  [%d hole col(s) in window]", n_holes)
      }
      gate_str <- ""
      if (r$sg == "gated")
        gate_str <- sprintf("  [span-gated: displaced %d mm @ clean %d]",
                            r$sgmm, r$disp_clean_s)
      cat(sprintf("      %-22s msa %d-%d (span %d)  %d mm  [%s]%s%s%s%s\n",
                  row$name, st, en, en - st + 1L,
                  row$n_mismatches, row$map_orientation, other_str, occ_str, hole_str, gate_str))
    }
    final_rows[[length(final_rows) + 1L]] <- r
  }
}

check("[1.4] all primer rows collected", length(final_rows) == n_expect,
      paste0("collected ", length(final_rows), " of ", n_expect,
             " primer rows — a row was not appended in the [1.4] loop; stopping before rbind"))
map_df <- do.call(rbind, lapply(final_rows, function(p) p$row))
check("[1.4] mapping table complete", nrow(map_df) == n_expect,
      paste0("rows ", nrow(map_df), " vs expected ", n_expect))

## span-gate counters (log metrics in [1.9] + closing banner)
n_span_gated    <- sum(map_df$span_gate == "gated")
n_span_noop     <- sum(map_df$span_gate == "no_op")
n_span_rejected <- sum(map_df$span_gate == "rejected")

## ---- [1.5] unmapped (documented results — excluded for that taxon) -------------
cat("\n[1.5] unmapped primers (named; excluded downstream for that taxon)\n")
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  u <- map_df[map_df$Taxa == T & map_df$gene == G & map_df$map_mode == "none", ]
  if (nrow(u) == 0L) {
    note("[1.5]", sprintf("%s %s: all %d primers mapped", T, G, sum(prim$gene == G)))
  } else {
    labs <- sprintf("%s (min_mm=%s%s)", u$name,
                    ifelse(is.na(u$min_mismatches), "NA", u$min_mismatches),
                    ifelse(u$span_gate == "rejected", ", span-rejected", ""))
    note("[1.5]", sprintf("%s %s: %d unmapped: %s", T, G, nrow(u),
                          paste(labs, collapse = ", ")))
  }
}

## ---- [1.6] pair coverage per taxon + amplicon sizes (+ clean spans) -------------
cat("\n[1.6] pair coverage per taxon (both strands mapped -> full; amplicon size = R_end - F_start + 1; + amplicon_size_clean + span-gate audits)\n")

## MSA [st, en] -> clean start/end (1-based) — the gap-free projection of the
## amplicon span (CONS_CHARS is in scope from [1.3]).
clean_proj <- function(key_t, st, en) {
  cc <- CONS_CHARS[[key_t]]
  c(sum(cc[seq_len(st - 1L)] != "-") + 1L, sum(cc[seq_len(en)] != "-"))
}

cov_rows <- vector("list", 0)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  key_t <- paste0(T, "_", G)
  mapped_names <- map_df$name[map_df$Taxa == T & !is.na(map_df$start)]
  for (k in seq_len(nrow(pair_table))) {
    if (pair_table$gene[k] != G) next
    p  <- pair_table$assay[k]
    fm <- paste0(p, "_F") %in% mapped_names
    rm_ <- paste0(p, "_R") %in% mapped_names
    cov <- if (fm && rm_) "full" else if (fm || rm_) "one_strand" else "none"

    ## amplicon size (MSA coordinates, gap-inclusive span) + clean span
    amp_size  <- NA_integer_
    amp_clean <- NA_integer_
    orientation_ok <- NA
    if (cov == "full") {
      frow <- map_df[map_df$Taxa == T & map_df$name == paste0(p, "_F"), ]
      rrow <- map_df[map_df$Taxa == T & map_df$name == paste0(p, "_R"), ]
      f_start <- frow$start[1L]; r_end <- rrow$end[1L]
      r_start <- rrow$start[1L]; f_end <- frow$end[1L]
      amp_size <- r_end - f_start + 1L
      ## flag if primers are in the wrong order (F should be left of R on fwd strand)
      orientation_ok <- (f_start <= r_start)
      if (!orientation_ok)
        note("[1.6]", sprintf("%s x %s: F start (%d) > R start (%d) — primers may be mis-mapped or non-nested",
                              T, p, f_start, r_start))
      ## de-gapped (clean) amplicon length; negative if anti-facing
      cp <- clean_proj(key_t, f_start, r_end)
      amp_clean <- cp[2L] - cp[1L] + 1L
    }

    cov_rows[[length(cov_rows) + 1L]] <- data.frame(
      assay = p, gene = pair_table$gene[k], Taxa = T,
      F_mapped = fm, R_mapped = rm_, cov = cov,
      amplicon_size = amp_size,
      orientation_ok = orientation_ok,
      amplicon_size_clean = amp_clean,
      stringsAsFactors = FALSE)
  }
}
pair_cov <- do.call(rbind, cov_rows)

for (i in seq_len(nrow(taxon))) {
  T   <- taxon$Taxa[i]
  sub <- pair_cov[pair_cov$Taxa == T, ]
  cat(sprintf("  %-14s full=%-3d one_strand=%-3d none=%-3d\n", T,
              sum(sub$cov == "full"), sum(sub$cov == "one_strand"), sum(sub$cov == "none")))
  ones  <- sub$assay[sub$cov == "one_strand"]
  nones <- sub$assay[sub$cov == "none"]
  if (length(ones) > 0L)  cat("      one_strand:", paste(ones, collapse = ", "), "\n")
  if (length(nones) > 0L) cat("      none      :", paste(nones, collapse = ", "), "\n")
  ## print amplicon sizes for full pairs (gap-inclusive + clean)
  fulls <- sub[sub$cov == "full", ]
  if (nrow(fulls) > 0L) {
    cat("      amplicons (bp):\n")
    for (k in seq_len(nrow(fulls)))
      cat(sprintf("        %-22s %6d bp (clean %d)%s\n",
                  fulls$assay[k], fulls$amplicon_size[k], fulls$amplicon_size_clean[k],
                  if (isTRUE(fulls$orientation_ok[k]) == FALSE) "  ** F>R — CHECK ** " else ""))
  }
}

## Span-gate backstop audits:
##   CHECK — every FULL pair with at least one non-ungated strand must face
##   within the clean span bound (the invariant the gate exists to enforce;
##   a rejected strand makes the pair one_strand, so it is not checked —
##   this CHECK is a backstop that should never fire).
##   NOTE  — full pair, both strands ungated (both tolerant; no exact anchor
##   to gate on), clean span over the bound: documented known limitation;
##   the pair stays full and is named for inspection.
span_viol <- character(0)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  sub <- pair_cov[pair_cov$Taxa == T & pair_cov$gene == G & pair_cov$cov == "full", ]
  for (k in seq_len(nrow(sub))) {
    p  <- sub$assay[k]
    fsg <- map_df$span_gate[map_df$Taxa == T & map_df$name == paste0(p, "_F")][1L]
    rsg <- map_df$span_gate[map_df$Taxa == T & map_df$name == paste0(p, "_R")][1L]
    ac  <- sub$amplicon_size_clean[k]
    if (fsg != "ungated" || rsg != "ungated") {
      if (!(ac >= 1L && ac <= CONFIG$MAX_AMP_SPAN))
        span_viol <- c(span_viol,
          sprintf("%s (%s x %s): clean span %d bp outside [1, %d] (span_gate F=%s, R=%s)",
                  p, T, G, ac, CONFIG$MAX_AMP_SPAN, fsg, rsg))
    } else if (ac > CONFIG$MAX_AMP_SPAN) {
      note("[1.6]",
        sprintf("%s (%s): clean span %d bp > MAX_AMP_SPAN=%d, no exact anchor to gate on (both strands tolerant) — inspect",
                p, T, ac, CONFIG$MAX_AMP_SPAN))
    }
  }
}
check("[1.6] every full pair with a gated strand faces within the clean span bound",
      length(span_viol) == 0L,
      paste0(paste(span_viol, collapse = " ; "),
             " — a paralogous exact hit must never produce a multi-kbp amplicon"))

n_full <- sum(pair_cov$cov == "full")
check("[1.6] at least one full (assay, taxon) combination", n_full >= 1L,
      "no assay has both strands mapped on any taxon — step 2 would have nothing to score")


## ==== CHUNK BREAK (safe: everything above is defined; paste resumes here) ====

## ---- [1.7] mapped-row invariants (every row re-derives from its own values) ----
cat("\n[1.7] mapped-row invariants (per-row re-derivation; RC-aware; rule-aware win_min_occ) + span-gate re-derivations\n")
m <- map_df[!is.na(map_df$start), ]
bad <- character(0)
n_occ_zero <- 0L; n_occ_low <- 0L
for (i in seq_len(nrow(m))) {
  r   <- m[i, ]
  tag <- sprintf("%s/%s", r$Taxa, r$name)
  if (r$end < r$start || r$end - r$start + 1L < nchar(r$seq))
    bad <- c(bad, paste0(tag, ": gapped span ", r$end - r$start + 1L,
                         " < primer length ", nchar(r$seq)))
  sw <- strip_gaps(r$target_with_gaps)
  if (nchar(sw) != nchar(r$seq))
    bad <- c(bad, paste0(tag, ": stripped window length ", nchar(sw),
                         " != primer length ", nchar(r$seq)))
  ## round-trip mismatch: use the orientation that actually won
  check_seq <- if (r$map_orientation == "rc") reverse_complement(r$seq) else r$seq
  w <- window_mismatch_counts(check_seq, sw)
  if (nrow(w) != 1L || w$n_mismatches[1L] != r$n_mismatches)
    bad <- c(bad, paste0(tag, ": round-trip mm (", r$map_orientation, ") ",
                         if (nrow(w) == 1L) w$n_mismatches[1L] else "ERR",
                         " != recorded ", r$n_mismatches))
  ## rule-aware win_min_occ invariant:
  ##   gaps_never_win: a consensus '-' column exists only at occ = 0 (zero
  ##   coverage), so a mapped window span carrying one is a stop-worthy
  ##   structural anomaly.
  ##   gaps_can_win: a consensus '-' column is a legitimate hole (occ may be
  ##   0); 0 is allowed — NOTE-level only.
  if (r$consensus_rule == "gaps_never_win") {
    if (r$win_min_occ < 1L)
      bad <- c(bad, paste0(tag, ": win_min_occ=", r$win_min_occ,
                           " < 1 under gaps_never_win — window span crossed a zero-coverage column (inspect the MSA)"))
  } else {
    if (r$win_min_occ == 0L) n_occ_zero <- n_occ_zero + 1L
    if (r$win_min_occ < CONFIG$MIN_OCC) n_occ_low <- n_occ_low + 1L
  }
}
check("[1.7] all mapped rows re-derive", length(bad) == 0L, paste(bad, collapse = " ; "))
if (n_occ_zero > 0L || n_occ_low > 0L)
  note("[1.7]", sprintf("gaps_can_win rows: %d with win_min_occ = 0 (span crosses an all-gap hole), %d with win_min_occ < MIN_OCC=%d (includes the zero rows) — NOTE only",
                        n_occ_zero, n_occ_low, CONFIG$MIN_OCC))
cat(sprintf("    %d mapped rows checked: span >= length | stripped length == primer length | round-trip mm (RC-aware) | win_min_occ rule-aware\n",
            nrow(m)))

## Span-gate re-derivations: the span-gate columns re-derive from the table
## itself (no new columns needed).
check("[1.7] span_gate vocabulary {ungated, no_op, gated, rejected}",
      all(map_df$span_gate %in% c("ungated", "no_op", "gated", "rejected")),
      paste0("unexpected value(s): ",
             paste(unique(map_df$span_gate[!map_df$span_gate %in% c("ungated", "no_op", "gated", "rejected")]),
                   collapse = ", ")))
check("[1.7] span_gate_mm is NA iff span_gate in {ungated, no_op}",
      identical(is.na(map_df$span_gate_mm), map_df$span_gate %in% c("ungated", "no_op")),
      "span_gate_mm NA-pattern disagrees with span_gate")
rej_rows <- map_df$span_gate == "rejected"
check("[1.7] every rejected row is unmapped (map_mode = none)",
      all(map_df$map_mode[rej_rows] == "none"),
      "a span-rejected row is not an unmapped row")
## recompute amplicon_size_clean for every full pair from the rows' MSA starts
clean_reviol <- character(0)
for (i in seq_len(nrow(taxon))) {
  T <- taxon$Taxa[i]; G <- taxon$Gene[i]
  key_t <- paste0(T, "_", G)
  sub <- pair_cov[pair_cov$Taxa == T & pair_cov$gene == G & pair_cov$cov == "full", ]
  for (k in seq_len(nrow(sub))) {
    p <- sub$assay[k]
    frow <- map_df[map_df$Taxa == T & map_df$name == paste0(p, "_F"), ]
    rrow <- map_df[map_df$Taxa == T & map_df$name == paste0(p, "_R"), ]
    cp <- clean_proj(key_t, frow$start[1L], rrow$end[1L])
    ac2 <- cp[2L] - cp[1L] + 1L
    if (is.na(sub$amplicon_size_clean[k]) || ac2 != as.numeric(sub$amplicon_size_clean[k]))
      clean_reviol <- c(clean_reviol,
        sprintf("%s (%s x %s): clean span recomputes to %d, table has %s",
                p, T, G, ac2,
                if (is.na(sub$amplicon_size_clean[k])) "NA" else sub$amplicon_size_clean[k]))
    fsg <- frow$span_gate[1L]; rsg <- rrow$span_gate[1L]
    if ((fsg != "ungated" || rsg != "ungated") && !(ac2 >= 1L && ac2 <= CONFIG$MAX_AMP_SPAN))
      clean_reviol <- c(clean_reviol,
        sprintf("%s (%s x %s): recomputed clean span %d bp outside [1, %d]",
                p, T, G, ac2, CONFIG$MAX_AMP_SPAN))
  }
}
check("[1.7] amplicon_size_clean re-derives from the rows' MSA starts (bound re-verified)",
      length(clean_reviol) == 0L, paste(clean_reviol, collapse = " ; "))

## ---- [1.8] outputs ----------------------------------------------------------------
cat("\n[1.8] outputs (results/)\n")
MAP_OUT <- file.path(OUT, "03_primer_mapping.csv")
COV_OUT <- file.path(OUT, "04_pair_coverage.csv")
write.csv(map_df, MAP_OUT, row.names = FALSE)
write.csv(pair_cov, COV_OUT, row.names = FALSE)
map_chk <- read.csv(MAP_OUT, stringsAsFactors = FALSE)
cov_chk <- read.csv(COV_OUT, stringsAsFactors = FALSE)

## 17-column mapping schema and 9-column coverage schema.
MAP_COLS <- c("name", "seq", "gene", "Taxa", "start", "end",
              "target_with_gaps", "n_mismatches", "min_mismatches",
              "map_orientation", "other_orientation_mm", "map_mode",
              "win_min_occ", "win_n_low_occ", "consensus_rule",
              "span_gate", "span_gate_mm")
check("[1.8] mapping round-trip schema (17 columns)",
      identical(names(map_chk), MAP_COLS),
      paste0("columns found: ", paste(names(map_chk), collapse = ", ")))
check("[1.8] mapping round-trip values",
      nrow(map_chk) == nrow(map_df) &&
        identical(map_chk$name, map_df$name) &&
        identical(map_chk$map_mode, map_df$map_mode) &&
        identical(map_chk$consensus_rule, map_df$consensus_rule) &&
        identical(is.na(map_chk$start), is.na(map_df$start)) &&
        identical(as.character(map_chk$target_with_gaps), as.character(map_df$target_with_gaps)) &&
        identical(map_chk$n_mismatches, map_df$n_mismatches) &&
        identical(map_chk$min_mismatches, map_df$min_mismatches) &&
        all(map_chk$start[!is.na(map_chk$start)] == map_df$start[!is.na(map_df$start)]) &&
        all(map_chk$end[!is.na(map_chk$end)]   == map_df$end[!is.na(map_df$end)]) &&
        all(map_chk$win_min_occ[!is.na(map_chk$win_min_occ)] ==
            map_df$win_min_occ[!is.na(map_df$win_min_occ)]) &&
        all(map_chk$win_n_low_occ[!is.na(map_chk$win_n_low_occ)] ==
            map_df$win_n_low_occ[!is.na(map_df$win_n_low_occ)]) &&
        identical(map_chk$span_gate, map_df$span_gate) &&
        identical(is.na(map_chk$span_gate_mm), is.na(map_df$span_gate_mm)) &&
        (all(is.na(map_chk$span_gate_mm)) ||
          all(map_chk$span_gate_mm[!is.na(map_chk$span_gate_mm)] ==
              map_df$span_gate_mm[!is.na(map_df$span_gate_mm)])),
      "re-read mapping differs from the in-memory table")
COV_COLS <- c("assay", "gene", "Taxa", "F_mapped", "R_mapped", "cov",
              "amplicon_size", "orientation_ok", "amplicon_size_clean")
check("[1.8] coverage round-trip schema (9 columns)",
      identical(names(cov_chk), COV_COLS),
      paste0("columns found: ", paste(names(cov_chk), collapse = ", ")))
check("[1.8] coverage round-trip values",
      nrow(cov_chk) == nrow(pair_cov) &&
        identical(cov_chk$cov, pair_cov$cov) &&
        identical(is.na(cov_chk$amplicon_size_clean), is.na(pair_cov$amplicon_size_clean)) &&
        (all(is.na(cov_chk$amplicon_size_clean)) ||
          all(cov_chk$amplicon_size_clean[!is.na(cov_chk$amplicon_size_clean)] ==
              pair_cov$amplicon_size_clean[!is.na(pair_cov$amplicon_size_clean)])),
      "re-read pair coverage differs from the in-memory table")
cat(sprintf("    03_primer_mapping.csv  %d rows (%d mapped / %d unmapped) + span_gate / span_gate_mm\n",
            nrow(map_df), nrow(m), nrow(map_df) - nrow(m)))
cat(sprintf("    04_pair_coverage.csv   %d rows (%d full) + amplicon_size_clean\n", nrow(pair_cov), n_full))

## ---- [1.9] audit log -----------------------------------------------------------------
LOG <- file.path(OUT, "run_audit_log.csv")
n_gcw_templates <- sum(vapply(RULE, function(r) r == "gaps_can_win", logical(1)))
log_row <- log_row_make(run_id, "1", "PASS", BASE, parent_run = parent_run,
                        metrics = paste0("n_mapping_rows=", nrow(map_df),
                                         ";n_mapped=", nrow(m),
                                         ";n_unmapped=", nrow(map_df) - nrow(m),
                                         ";n_full_pairs=", n_full,
                                         ";max_map_mismatches=", CONFIG$MAX_MAP_MISMATCHES,
                                         ";n_templates_gaps_can_win=", n_gcw_templates,
                                         ";n_win_min_occ_zero_rows=", n_occ_zero,
                                         ";n_span_gated=", n_span_gated,
                                         ";n_span_noop=", n_span_noop,
                                         ";n_span_rejected=", n_span_rejected,
                                         ";max_pair_span=", CONFIG$MAX_AMP_SPAN))
if (file.exists(LOG))
  log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)
log_chk <- read.csv(LOG, stringsAsFactors = FALSE)
check("[1.9] audit log appended",
      any(log_chk$run_id == run_id & log_chk$step == "1"),
      "log row missing after write")

cat("\n============================================================================\n")
cat(sprintf("STEP 1 PASSED | run_id %s (continued from step 0 %s)\n", run_id, parent_run))
cat(sprintf("  rules   : %d template(s) gaps_can_win | %d gaps_never_win (manifest/CONFIG)\n",
            sum(vapply(RULE, function(r) r == "gaps_can_win", logical(1))),  # ✅ Correct
            sum(vapply(RULE, function(r) r == "gaps_never_win", logical(1)))))
cat(sprintf("  mapped: %d of %d primer x template rows (exact=%d, tolerant=%d)\n",
            nrow(m), nrow(map_df),
            sum(map_df$map_mode == "exact"), sum(map_df$map_mode == "tolerant")))
cat(sprintf("  unmapped: %d (named in [1.5]; carried as NA rows in 03_primer_mapping.csv)\n",
            nrow(map_df) - nrow(m)))
cat(sprintf("  span gate : MAX_AMP_SPAN = %d | gated=%d | noop=%d | rejected=%d\n",
            CONFIG$MAX_AMP_SPAN, n_span_gated, n_span_noop, n_span_rejected))
cat(sprintf("  full (assay, taxon) combinations: %d -> these go to step 2\n", n_full))
cat("Outputs: results/03_primer_mapping.csv | 04_pair_coverage.csv | run_audit_log.csv\n")
cat("Next: step2_primerminer.r (Ubuntu/WSL R — PrimerMiner dependency)\n")
cat("============================================================================\n")