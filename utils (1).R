# ============================================================================
# utils.R — shared utilities for a primer mapping / DNA metabarcoding pipeline
# ----------------------------------------------------------------------------
# Version   : 2.0
#             v2.0: Generalized for any taxa and any primer set.
#               - Two consensus rules (gaps_never_win, gaps_can_win) selected
#                 per template via CONFIG.
#               - Pair-aware span gate limits amplicon size during mapping.
#               - Pure base R; no external packages required.
#
# Language  : base R only (no packages; Bioconductor is not needed).
#
# OUTPUT STYLE
#   - All console output via cat() (base R data.frames only).
#   - check(tag, cond, fail_msg): stops the run on failure.
#   - note(tag, msg): informational; never stops.
#
# MISMATCH SEMANTICS
#   - Every IUPAC letter stands for a SET of DNA bases:
#       A={A}  C={C}  G={G}  T={T}
#       R={A,G}  Y={C,T}  S={C,G}  W={A,T}  K={G,T}  M={A,C}
#       B={C,G,T}  D={A,G,T}  H={A,C,T}  V={A,C,G}  N={A,C,G,T}
#   - A position is a mismatch ONLY when the primer base set and the subject
#     base set are DISJOINT (share no base).
#   - Consequence: a degenerate base is NEVER counted as a mismatch when it
#     overlaps. Primer 'R' vs a column carrying 'A', 'G', 'R', 'M', 'W', 'D',
#     'V' or 'N' is a MATCH at that position. Primer 'R' vs a column 'Y' or
#     'C' is a mismatch (no shared base).
#   - A subject gap '-' is the empty set, so it is always a mismatch.
#   - Primers must be gap-free design sequences; a gap in a primer is an error.
#
# CONSENSUS RULES (population-level)
#   Two engines, selected per template via CONFIG:
#
#   "gaps_never_win" (base_consensus, DEFAULT):
#     Per column the most frequent NON-GAP symbol wins. A tie among non-gap
#     symbols -> IUPAC code of their union. Gaps never win: a gap-majority
#     column still returns a base. Zero-coverage column -> '-' with occ 0.
#
#   "gaps_can_win" (gap_inclusive_consensus):
#     Per column the most frequent symbol over bases, IUPAC codes AND the gap
#     wins. A tie that includes the gap -> '-'. A base-only tie -> IUPAC code
#     of their union. Zero-coverage column -> '-'. Same return shape as
#     base_consensus: list(cons, occ).
#
#   Both require uppercase IUPAC + gaps input; anything else is a hard error.
#
# MATCH CHOICE
#   - primer_match() returns ALL windows with mismatches <= max_mismatches.
#   - best_match() picks the window with fewest mismatches; ties broken by
#     most concrete (single-base) support, then leftmost start.
#   - primer_min_mismatches() reports the closest an unmapped primer came.
#
# HOW THE STEPS USE THIS FILE
#   step 0: read_fasta, fasta_dedup, check, note, CONFIG, consensus_rule_for,
#           UTILS_VERSION
#   step 1: consensus, consensus_rule_for, gap_positions, best_match,
#           primer_min_mismatches, strip_gaps, window_mismatch_counts,
#           pair_span, span_keep, match_support, gated_best_match,
#           CONFIG settings, UTILS_VERSION
#   step 2: to_wsl_path, reverse_complement, extract_window, strip_gaps,
#           base_set, check, note, CONFIG settings, UTILS_VERSION
#   step 3: extract_window, strip_gaps
#   step 4: check, note
#   step 5: check, note, to_wsl_path, log_row_make, CONFIG settings,
#           UTILS_VERSION
# ============================================================================

# ----------------------------------------------------------------------------
# CONFIG — the single config block for the whole pipeline.
# Steps must NEVER hardcode thresholds; they reference CONFIG$... only.
# A user with different primers or data edits this block (or reassigns
# CONFIG$... after sourcing, in a fresh session — everything downstream
# reads these values, so one change propagates cleanly).
# ----------------------------------------------------------------------------
CONFIG <- list(
  ## Step 2 suitability rule. PrimerMiner returns a per-primer penalty for
  ## each primer x template combination. For an ASSAY x species combination:
  ##     total = score_F + score_R        (the SUM of both primers)
  ##     total >  SCORE_THRESHOLD  ->  unsuitable for detection
  ##     total <= SCORE_THRESHOLD  ->  suitable, proceed to barcoding
  ## Which primer carries the penalty is irrelevant — only the sum matters.
  ## A total of EXACTLY SCORE_THRESHOLD is suitable: the rule is "over 100
  ## is unsuitable". Step 2 prints how many combinations land exactly on the
  ## boundary so this choice is visible in every run.
  SCORE_THRESHOLD = 100,
  
  ## Step 1 mapping tolerance: maximum number of mismatching positions
  ## (see MISMATCH SEMANTICS — degenerate overlaps are never mismatches)
  ## allowed when matching a primer to the gap-free consensus of a taxon MSA.
  ##     0  = strict IUPAC overlap
  ##     n  = also report primers binding with up to n non-overlapping
  ##          positions (e.g. 3 or 6 to explore tolerant binding)
  ## Every run reports min_mismatches for unmapped primers either way, so
  ## raising the tolerance later shows exactly what you gained. Indel-rich
  ## regions are the typical reason to raise this value.
  MAX_MAP_MISMATCHES = 7L,
  
  ## Step 1 occupancy reporting threshold. For every mapped primer window,
  ## step 1 reports win_min_occ (minimum per-column base count across the
  ## window) and win_n_low_occ (number of window columns with occ < MIN_OCC).
  ## NOTE-level only — it never excludes a column or a mapping result.
  ## 1 = flag nothing; 2 = flag columns supported by a single sequence only
  ## (the default, suitable for partially covered MSAs).
  MIN_OCC = 2L,
  
  ## Consensus rule selection. Which consensus engine step 1 computes for
  ## each (Taxa|Gene) template:
  ##     "gaps_never_win"  = base_consensus() (DEFAULT): a gap-majority
  ##                         column still yields a base.
  ##     "gaps_can_win"    = gap_inclusive_consensus(): a gap-majority
  ##                         column -> '-'. Use for indel-rich alignments
  ##                         where primer windows span variable regions.
  ## CONSENSUS_RULE applies to every template; CONSENSUS_RULE_OVERRIDES
  ## is a named list keyed "Taxa|Gene" that wins over the default (see the
  ## sample entries below). Step 0 publishes the effective rule per template
  ## into the manifest; step 1 selects from the manifest and checks it
  ## against the CONFIG-derived rule (drift = stop; re-run step 0).
  CONSENSUS_RULE = "gaps_never_win",
  ## Per-template overrides: a named list keyed "Taxa|Gene" that wins over
  ## the global default for that one template. Leave EMPTY (taxa-agnostic
  ## default); add entries only where a template genuinely needs
  ## gap-inclusive consensus, e.g.
  ##   CONSENSUS_RULE_OVERRIDES = list("Anurans|mt12S" = "gaps_can_win")
  CONSENSUS_RULE_OVERRIDES = list("Anurans|mt12S" = "gaps_can_win",
                                  "Anurans|mt16S" = "gaps_can_win"),
  
  ## Pair-aware span gate. For every complete pair, when one strand is
  ## anchored at 0 mm (exact = high-confidence), the partner primer's
  ## candidate windows must FACE that anchor, be non-overlapping with it,
  ## and sit at a DE-GAPPED (clean-consensus) span of at most MAX_AMP_SPAN
  ## bp:  span = R_end - F_start + 1.
  ## Candidate windows are re-ranked WITHIN the span under the existing
  ## tie-break (mm -> concrete support -> leftmost). If no in-span window
  ## passes MAX_MAP_MISMATCHES, the partner is a documented span-rejection
  ## (an unmapped row; the distant window is never scored). No gate is
  ## applied when the partner's anchor is tolerant (1+ mm) or unmapped.
  ## Set this to the maximum expected amplicon length in your dataset
  ## (de-gapped). Raise only for datasets with genuinely longer amplicons.
  MAX_AMP_SPAN = 700L,
  
  ## Step 2 only: number of adjacent bases folded into the PrimerMiner
  ## penalty. 2 = broader penalty context (default); 1 = minimal penalty
  ## semantics. Pick one and keep it consistent for the whole project.
  ## Step 2 prints its provenance in [2.7].
  PRIMERMINER_ADJACENT = 2L
)

# ----------------------------------------------------------------------------
# Core version string. Step 0 records this in the manifest; steps 1-5
# re-check it against the manifest value. Bump when behavior or the test
# battery changes; re-run step 0 to publish the new version.
# ----------------------------------------------------------------------------
UTILS_VERSION <- "2.0"



# ----------------------------------------------------------------------------
# audit helpers
# ----------------------------------------------------------------------------
## check: PASS prints and continues; FAIL prints and STOPS the run.
check <- function(tag, cond, fail_msg = "see message") {
  if (isTRUE(cond)) {
    cat(sprintf("  [CHECK] %s : PASS\n", tag))
    return(invisible(TRUE))
  }
  cat(sprintf("  [CHECK] %s : FAIL\n", tag))
  stop("CHECK FAILED — ", tag, ": ", fail_msg, call. = FALSE)
}

## note: informational; never stops the run.
note <- function(tag, msg = "") {
  cat(sprintf("  [NOTE ] %s : %s\n", tag, msg))
  invisible(msg)
}

## log_row_make: the ONLY way any step appends to run_audit_log.csv.
## Fixed 7-column schema for ALL steps, so the appended log can never
## drift out of column sync (mixed-schema rbind = hard crash). Step-
## specific counters go into `metrics` as a "k=v;k=v" string. `timestamp`
## defaults to completion time — the row records when the step ENDED.
log_row_make <- function(run_id, step, status, base,
                         parent_run = NA_character_, metrics = "",
                         timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")) {
  data.frame(run_id = run_id, parent_run_id = parent_run, step = step,
             timestamp = timestamp, status = status, base = base,
             metrics = metrics, stringsAsFactors = FALSE)
}


# ----------------------------------------------------------------------------
# IUPAC tables
# ----------------------------------------------------------------------------
BASE_SET <- list(
  A = c("A"), C = c("C"), G = c("G"), T = c("T"),
  R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"),
  B = c("C", "G", "T"), D = c("A", "G", "T"),
  H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T"))

## key for a base set (sorted, joined) — table lookup key
set_key <- function(bases) paste(sort(bases), collapse = "")

## reverse lookup: base-set key -> IUPAC code
IUPAC_BY_SET <- list()
for (nm in names(BASE_SET)) IUPAC_BY_SET[[set_key(BASE_SET[[nm]])]] <- nm

## Set of DNA bases for one IUPAC letter (case-insensitive).
## '-' (gap) -> character(0) = the empty set. Anything else -> hard error.
base_set <- function(ch) {
  ch <- toupper(ch)
  if (ch == "-") return(character(0))
  bs <- BASE_SET[[ch]]
  if (is.null(bs))
    stop("base_set: not a valid IUPAC DNA letter: '", ch, "'", call. = FALSE)
  bs
}

## The IUPAC code for a column of symbols (each symbol may itself be a code):
## union of the members' base sets, then decode. Always decodable — every
## subset of {A,C,G,T} has exactly one IUPAC code.
iupac_code <- function(symbols) {
  symbols <- unique(toupper(symbols))
  if (length(symbols) == 0L) stop("iupac_code: no symbols given", call. = FALSE)
  bad <- setdiff(symbols, names(BASE_SET))
  if (length(bad) > 0L)
    stop("iupac_code: not an IUPAC DNA symbol: ", paste(bad, collapse = ", "), call. = FALSE)
  bases <- unique(unlist(lapply(symbols, base_set)))
  code <- IUPAC_BY_SET[[set_key(bases)]]
  if (is.null(code))
    stop("iupac_code: base set {", paste(bases, collapse = ","),
         "} has no IUPAC code (should be impossible)", call. = FALSE)
  code
}

## Validate a character vector of IUPAC letters (+ optional gaps); names the
## offending characters in the error.
validate_iupac <- function(x, what, allow_gap = FALSE) {
  okchars <- c(names(BASE_SET), if (allow_gap) "-")
  chars <- unique(unlist(strsplit(x, "", fixed = TRUE)))
  bad <- setdiff(chars, okchars)
  if (length(bad) > 0L)
    stop(what, " contains non-IUPAC character(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# FASTA
# ----------------------------------------------------------------------------
## Robust text-level FASTA reader.
##   - one record per '>' header; wrapped or single-line bodies; CRLF and
##     blank lines tolerated; sequence output is ALWAYS uppercase.
##   - name = full header text after '>' (description included); downstream
##     steps parse taxonomy from it.
##   - duplicate names are KEPT here and FLAGGED (dup_names); the policy
##     (keep first, drop rest) lives in fasta_dedup() so the reader stays
##     pure and the decision stays visible where it happens.
##   - raggedness is REPORTED (ragged = TRUE); step 0 turns it into a
##     CHECK failure — it is never silently narrowed.
read_fasta <- function(path) {
  if (!file.exists(path)) stop("read_fasta: file not found: ", path, call. = FALSE)
  ln <- readLines(path, n = -1, warn = FALSE)
  ln <- sub("\r$", "", ln)                          # tolerate CRLF
  h <- which(startsWith(ln, ">"))
  if (length(h) == 0L) stop("read_fasta: no FASTA headers in ", path, call. = FALSE)
  rec_names <- sub("^>[ \t]*", "", ln[h])
  body <- which(!startsWith(ln, ">") & nzchar(trimws(ln)))
  own <- findInterval(body, h)                      # body line -> nearest header above
  if (any(own == 0L))
    stop("read_fasta: sequence lines appear before the first '>' header in ", path, call. = FALSE)
  seqs <- character(length(h))
  for (j in seq_along(h)) seqs[j] <- paste(ln[body[own == j]], collapse = "")
  if (any(nchar(seqs) == 0L))
    stop("read_fasta: empty sequence record(s) in ", path, call. = FALSE)
  seqs <- toupper(seqs)                             # deterministic output case
  validate_iupac(seqs, sprintf("FASTA sequence in %s ", path), allow_gap = TRUE)
  w <- nchar(seqs)
  list(names        = rec_names,
       seqs         = seqs,
       records      = length(h),
       unique_names = length(unique(rec_names)),
       width        = w,
       ragged       = (min(w) < max(w)),
       dup_names    = rec_names[duplicated(rec_names)])
}

## Duplicate-name policy: keep the FIRST record per name, drop the rest,
## print a NOTE naming what was dropped. Returns same list shape with
## cleared dup_names.
fasta_dedup <- function(f) {
  keep <- !duplicated(f$names)
  n_dropped <- sum(!keep)
  if (n_dropped > 0L)
    note("fasta-dedup", sprintf("%d duplicate record(s) dropped, first occurrence kept: %s",
                                n_dropped, paste(unique(f$names[duplicated(f$names)]), collapse = ", ")))
  else
    note("fasta-dedup", "no duplicate names")
  list(names        = f$names[keep],
       seqs         = f$seqs[keep],
       records      = sum(keep),
       unique_names = length(unique(f$names[keep])),
       width        = f$width[keep],
       ragged       = f$ragged,
       dup_names    = character(0))
}

# ----------------------------------------------------------------------------
# degenerate-aware matching
# ----------------------------------------------------------------------------
## Mismatch count for EVERY window of `primer` in `subject` (see the
## MISMATCH SEMANTICS block at the top). Returns data.frame(start, n_mismatches)
## over all 1-based window starts (0 rows if the primer cannot fit).
## This is the engine; primer_match / best_match / primer_min_mismatches
## are thin, named views of it.
window_mismatch_counts <- function(primer, subject) {
  primer  <- toupper(primer)
  subject <- toupper(subject)
  if (grepl("-", primer, fixed = TRUE))
    stop("primer '", primer, "' contains a gap; primers must be gap-free design sequences", call. = FALSE)
  validate_iupac(primer, "primer", allow_gap = FALSE)
  validate_iupac(subject, "subject", allow_gap = TRUE)
  Lp <- nchar(primer)
  Lt <- nchar(subject)
  if (Lp == 0L || Lt < Lp)
    return(data.frame(start = integer(0), n_mismatches = integer(0), stringsAsFactors = FALSE))
  nwin <- Lt - Lp + 1L
  schars <- strsplit(subject, "", fixed = TRUE)[[1L]]
  ssets  <- lapply(schars, base_set)                # subject gap -> character(0)
  pchars <- strsplit(primer, "", fixed = TRUE)[[1L]]
  ## per distinct primer letter: TRUE over the whole subject where that
  ## letter's base set overlaps the column's base set (degenerate overlap)
  cache <- list()
  for (ch in unique(pchars)) {
    ps <- base_set(ch)
    cache[[ch]] <- vapply(ssets, function(ss) length(intersect(ps, ss)) > 0L, logical(1))
  }
  mm <- integer(nwin)
  for (j in seq_len(Lp))
    mm <- mm + as.integer(!cache[[pchars[j]]][j:(j + nwin - 1L)])
  data.frame(start = seq_len(nwin), n_mismatches = mm, stringsAsFactors = FALSE)
}

## All windows with mismatches <= max_mismatches.
primer_match <- function(primer, subject, max_mismatches = 4L) {
  w <- window_mismatch_counts(primer, subject)
  w[w$n_mismatches <= max_mismatches, , drop = FALSE]
}

## Concrete (single-base) support per candidate start: positions where the
## subject column is a single base (A/C/G/T) that also overlaps the primer.
## Used as a tie-break by best_match and gated_best_match.
match_support <- function(subject, primer, starts) {
  Lp <- nchar(primer)
  pc <- strsplit(toupper(primer), "", fixed = TRUE)[[1L]]
  vapply(starts, function(s) {
    a <- strsplit(substr(subject, s, s + Lp - 1L), "", fixed = TRUE)[[1L]]
    sum(vapply(seq_len(Lp), function(j)
      a[j] %in% c("A", "C", "G", "T") &&
        length(intersect(base_set(a[j]), base_set(pc[j]))) > 0L,
      logical(1)))
  }, integer(1))
}

## The single window to anchor on: fewest mismatches among passing windows,
## then most concrete (single-base) support, then leftmost start. Returns
## 1 row, or 0 rows (unmapped).
## WHY the support tie-break: at tolerance 0, a window over a fully
## degenerate consensus block (N-run) is a vacuous 0-mm match. Without the
## support tie-break, an early N-block would steal the anchor from a later
## real binding site.
best_match <- function(primer, subject, max_mismatches = 0L) {
  w <- primer_match(primer, subject, max_mismatches)
  if (nrow(w) == 0L) return(w)
  if (nrow(w) == 1L) return(w)
  w[order(w$n_mismatches, -match_support(subject, primer, w$start),
          w$start)[1L], , drop = FALSE]
}

## How close the primer came anywhere in the subject (NA if it cannot fit).
## Step 1 prints this for every unmapped primer — even at the strict default.
primer_min_mismatches <- function(primer, subject) {
  w <- window_mismatch_counts(primer, subject)
  if (nrow(w) == 0L) NA_integer_ else min(w$n_mismatches)
}

# ----------------------------------------------------------------------------
# Pair-aware span gate — clean-coordinate helpers + gated matcher.
# All coordinates are 1-based CLEAN (gap-free) consensus positions. The gate
# is orthogonal to the consensus rule and never loosens MAX_MAP_MISMATCHES.
# Step 1 applies the gate per complete pair; utils provides the primitives.
# ----------------------------------------------------------------------------

## Pair-level signed facing span (bp) of the two FINAL windows:
## positive = facing (F 5' of R); negative = anti-facing (F > R).
pair_span <- function(f_start, f_end, r_start, r_end) r_end - f_start + 1L

## Span-keep mask for candidate clean starts `s` (1-based), candidate primer
## length Lp, anchor window [a_start, a_end] (clean), side = "F" (candidate
## is the F strand; anchor is the partner R window) or "R" (candidate is the
## R strand; anchor is the partner F window).
##
##   side "F": facing span(s) = a_end - s + 1L
##   side "R": facing span(s) = (s + Lp - 1L) - a_start + 1L
##
## keep = non-overlapping AND within the size gate:
##   side "F": (s + Lp - 1L) < a_start  AND  a_end - s + 1L <= max_span
##   side "R": s > a_end               AND  (s + Lp - 1L) - a_start + 1L <= max_span
##
## Non-overlap makes the implied product >= Lp + anchor_Lp (a candidate that
## nests/overlaps the exact partner site is never a valid PCR partner) and
## therefore the kept span range is [LpF + LpR, max_span].
span_keep <- function(s, Lp, a_start, a_end, side, max_span) {
  if (side == "F") {
    (s + Lp - 1L) < a_start & (a_end - s + 1L) <= max_span
  } else {
    s > a_end & ((s + Lp - 1L) - a_start + 1L) <= max_span
  }
}

## best_match with the span constraint.
##   cands        = primer_match() frame for ONE orientation frame
##                  (columns: start, n_mismatches) — the candidate strand's
##                  passing windows in clean coordinates
##   primer       = the probe sequence that produced cands (support tie-break)
##   subject      = the clean consensus
##   a_start/a_end= the partner's exact anchor window (clean coords)
##   side         = "F" or "R" (candidate side; see span_keep)
##   max_span     = CONFIG$MAX_AMP_SPAN
## Returns list(
##   best         = 0-or-1-row frame (gated winner; same columns as best_match),
##   ungated_best = 0-or-1-row frame (the winner without the gate — for
##                  span_gate_mm reporting; 0 rows if nothing passed),
##   n_excluded   = integer, passing windows removed by the gate)
## Behavior: filter cands by span_keep; empty filter -> best = 0 rows
## (step 1 turns that into a rejection); otherwise re-rank the survivors
## with the SAME tie-break as best_match (fewest mm -> most concrete support
## -> leftmost). The gate never loosens max_mismatches: the candidate set is
## whatever primer_match already passed.
gated_best_match <- function(cands, primer, subject, a_start, a_end,
                             side, max_span) {
  if (nrow(cands) == 0L)
    return(list(best = cands, ungated_best = cands, n_excluded = 0L))
  ord <- order(cands$n_mismatches,
               -match_support(subject, primer, cands$start),
               cands$start)
  ungated_best <- cands[ord[1L], , drop = FALSE]
  keep <- span_keep(cands$start, nchar(primer), a_start, a_end, side, max_span)
  n_excluded <- as.integer(sum(!keep))
  if (n_excluded == nrow(cands))
    return(list(best = cands[0L, , drop = FALSE],
                ungated_best = ungated_best,
                n_excluded = n_excluded))
  sub <- cands[keep, , drop = FALSE]
  o2 <- order(sub$n_mismatches, -match_support(subject, primer, sub$start),
              sub$start)
  list(best = sub[o2[1L], , drop = FALSE],
       ungated_best = ungated_best,
       n_excluded = n_excluded)
}

# ----------------------------------------------------------------------------
# consensus + coordinate helpers
# ----------------------------------------------------------------------------
## BASE CONSENSUS (gaps never win; pure base R).
## seqs: character vector of equal-width, uppercase, gapped IUPAC sequences.
## Per column the most frequent NON-GAP SYMBOL wins. A tie among non-gap
## symbols -> the IUPAC code of their union. Gaps never win: a gap-majority
## column still returns a base (or tied IUPAC code), which prevents indel-
## rich alignments from opening holes in the consensus inside primer windows.
## Zero-coverage column (all gaps) -> '-' with occ 0.
## Returns list(cons, occ) where occ = per-column base count (occupancy).
base_consensus <- function(seqs) {
  if (length(seqs) == 0L) stop("base_consensus: no sequences", call. = FALSE)
  w <- vapply(seqs, nchar, integer(1L))
  if (any(w != w[1L])) stop("base_consensus: ragged input (uniform width required)", call. = FALSE)
  validate_iupac(seqs, "base_consensus input", allow_gap = TRUE)
  m   <- matrix(unlist(strsplit(seqs, "", fixed = TRUE), use.names = FALSE),
                nrow = length(seqs), byrow = TRUE)
  occ <- integer(ncol(m))
  out <- character(ncol(m))
  for (j in seq_len(ncol(m))) {
    col   <- m[, j]
    bases <- col[col != "-"]            # NON-GAP symbols only -> gaps never win
    occ[j] <- length(bases)             # occupancy = base count
    if (length(bases) == 0L) {
      out[j] <- "-"                     # zero-coverage column
    } else {
      tb   <- table(bases)
      tied <- names(tb)[tb == max(tb)]
      out[j] <- if (length(tied) == 1L) tied else iupac_code(tied)
    }
  }
  list(cons = paste(out, collapse = ""), occ = occ)
}


## GAP-INCLUSIVE CONSENSUS (gaps CAN win; pure base R).
## Per column the MOST FREQUENT SYMBOL over bases, IUPAC codes AND the gap
## wins. A tie that includes the gap -> '-'; a base-only tie -> the IUPAC
## code of their union. Zero-coverage column -> '-' with occ 0. Same return
## shape as base_consensus(): list(cons, occ) with occ = per-column NON-GAP
## base count. A gap column under this rule may carry occ > 0 (gap-majority
## with a few bases); those columns are "holes" that step 1 strips from the
## clean consensus. A mapped window's MSA span may include them, so
## win_min_occ can be 0 for rows mapped under this rule (NOTE-level only).
gap_inclusive_consensus <- function(seqs) {
  if (length(seqs) == 0L) stop("gap_inclusive_consensus: no sequences", call. = FALSE)
  w <- vapply(seqs, nchar, integer(1L))
  if (any(w != w[1L])) stop("gap_inclusive_consensus: ragged input (uniform width required)", call. = FALSE)
  validate_iupac(seqs, "gap_inclusive_consensus input", allow_gap = TRUE)
  m   <- matrix(unlist(strsplit(seqs, "", fixed = TRUE), use.names = FALSE),
                nrow = length(seqs), byrow = TRUE)
  occ <- integer(ncol(m))
  out <- character(ncol(m))
  for (j in seq_len(ncol(m))) {
    col  <- m[, j]
    occ[j] <- sum(col != "-")           # occupancy = NON-GAP base count
    tb   <- table(col)                  # gap is a first-class symbol here
    tied <- names(tb)[tb == max(tb)]
    if ("-" %in% tied) {
      out[j] <- "-"                    # gap in the tie -> gap wins
    } else {
      out[j] <- if (length(tied) == 1L) tied else iupac_code(tied)
    }
  }
  list(cons = paste(out, collapse = ""), occ = occ)
}

## Consensus rule vocabulary.
CONSENSUS_RULES_VALID <- c("gaps_never_win", "gaps_can_win")

## Rule selector. key = "Taxa|Gene"; overrides = named list keyed
## "Taxa|Gene" (e.g. CONFIG$CONSENSUS_RULE_OVERRIDES); default = fallback
## rule (e.g. CONFIG$CONSENSUS_RULE). Returns the effective rule, validated
## against CONSENSUS_RULES_VALID (a bad value is a hard error, named).
consensus_rule_for <- function(key, overrides = NULL, default = "gaps_never_win") {
  if (!default %in% CONSENSUS_RULES_VALID)
    stop("consensus_rule_for: invalid default rule '", default,
         "' (valid: ", paste(CONSENSUS_RULES_VALID, collapse = ", "), ")", call. = FALSE)
  r <- if (!is.null(overrides) && key %in% names(overrides)) overrides[[key]] else default
  if (!r %in% CONSENSUS_RULES_VALID)
    stop("consensus_rule_for: invalid consensus rule '", r, "' for template '", key,
         "' (valid: ", paste(CONSENSUS_RULES_VALID, collapse = ", "), ")", call. = FALSE)
  r
}

## Consensus dispatcher. rule = "gaps_never_win" (default; calls
## base_consensus) or "gaps_can_win" (calls gap_inclusive_consensus).
## Same return shape for both: list(cons, occ).
consensus <- function(seqs, rule = "gaps_never_win") {
  if (!rule %in% CONSENSUS_RULES_VALID)
    stop("consensus: unknown rule '", rule, "' (valid: ",
         paste(CONSENSUS_RULES_VALID, collapse = ", "), ")", call. = FALSE)
  if (rule == "gaps_never_win") base_consensus(seqs) else gap_inclusive_consensus(seqs)
}

## Compatibility wrapper: returns only the consensus string (gaps_never_win
## rule) without the occupancy vector. New code should call consensus() and
## use $occ.
## x: character vector of equal-width gapped IUPAC strings, OR a FASTA file
## path (read via read_fasta + fasta_dedup).
consensus_of_alignment <- function(x) {
  if (length(x) == 1L && file.exists(x)) {
    f <- read_fasta(x)
    x <- fasta_dedup(f)$seqs
  }
  base_consensus(x)$cons
}

## 1-based gapped-column positions carrying a base (gap-free projection map).
## Accepts EITHER a single gapped string ("AC-GT") OR a per-column symbol
## vector (c("A","C","-","G","T")); both normalize to per-column symbols.
gap_positions <- function(consensus) {
  if (length(consensus) == 1L && nchar(consensus) > 1L)
    consensus <- strsplit(consensus, "", fixed = TRUE)[[1L]]
  which(consensus != "-")
}

## Remove gaps (vectorized).
strip_gaps <- function(x) gsub("-", "", x, fixed = TRUE)

## Slice gapped windows [start, end] (inclusive, 1-based) out of aligned
## sequences (vectorized). Used by steps 2-3 to pull the primer-binding
## region / amplicon region out of the gapped MSA.
extract_window <- function(seqs, start, end) {
  if (length(start) != 1L || length(end) != 1L)
    stop("extract_window: start/end must be single integers", call. = FALSE)
  if (end < start)
    stop("extract_window: end must be >= start", call. = FALSE)
  if (end > max(nchar(seqs)))
    stop(sprintf("extract_window: end %d exceeds sequence width %d", end, max(nchar(seqs))), call. = FALSE)
  substring(seqs, start, end)
}

# ----------------------------------------------------------------------------
# reverse complement (IUPAC-aware, vectorized, output uppercase)
#   A<->T  C<->G  R<->Y  S<->S  W<->W  K<->M  B<->V  D<->H  H<->D  V<->B  N<->N
# ----------------------------------------------------------------------------
COMPLEMENT <- c(A = "T", C = "G", G = "C", T = "A",
                R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
                B = "V", V = "B", D = "H", H = "D", N = "N", "-" = "-")

reverse_complement <- function(x) {
  x <- toupper(x)
  validate_iupac(x, "reverse_complement input", allow_gap = TRUE)
  vapply(strsplit(x, "", fixed = TRUE),
         function(s) paste(rev(COMPLEMENT[s]), collapse = ""), character(1))
}

# ----------------------------------------------------------------------------
# path helpers
# ----------------------------------------------------------------------------
## Convert a Windows path (C:/a/b, backslashes tolerated) to its WSL /mnt
## form; POSIX paths pass through unchanged. Pure string operation, no file
## IO. Used by step 2 to derive the WSL spelling of the base directory.
to_wsl_path <- function(p) {
  p <- gsub("\\", "/", p, fixed = TRUE)
  m <- regmatches(p, regexec("^([A-Za-z]):(.*)$", p))[[1L]]
  if (length(m) == 3L) paste0("/mnt/", tolower(m[2L]), m[3L]) else p
}

# ============================================================================
# SELF-TEST BATTERY — runs on source; hard-stops on any failure.
# 115 algorithm tests with hand-computed expected values (no project data).
# If any test fails, the run stops before any pipeline step can execute.
# Pure base R; no external packages required.
# ============================================================================
utils_self_test <- function() {
  cat("\n===================== utils.R SELF-TEST (hard gate) =====================\n")
  n_ok <- 0L; n_fail <- 0L; failed <- character(0)
  t <- function(label, ok) {
    if (isTRUE(ok)) {
      n_ok <<- n_ok + 1L
      cat(sprintf("  [TEST %-52s] ok\n", label))
    } else {
      n_fail <<- n_fail + 1L
      failed <<- c(failed, label)
      cat(sprintf("  [TEST %-52s] ** FAIL **\n", label))
    }
    invisible(NULL)
  }
  errs <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
  
  ## ---- base_set / iupac_code ----
  t("base_set: A",                          identical(base_set("A"), "A"))
  t("base_set: R = {A,G}",                  identical(base_set("R"), c("A","G")))
  t("base_set: H = {A,C,T}",                identical(base_set("H"), c("A","C","T")))
  t("base_set: V = {A,C,G} (distinct from H)", identical(base_set("V"), c("A","C","G")))
  t("base_set: N = {A,C,G,T}",              identical(base_set("N"), c("A","C","G","T")))
  t("base_set: case-insensitive",           identical(base_set("r"), c("A","G")))
  t("base_set: gap -> empty set",           identical(base_set("-"), character(0)))
  t("base_set: invalid character errors",   errs(base_set("X")))
  t("iupac_code: C+T -> Y",                 identical(iupac_code(c("C","T")), "Y"))
  t("iupac_code: C+G+T -> B",               identical(iupac_code(c("C","G","T")), "B"))
  t("iupac_code: A+C+G+T -> N",             identical(iupac_code(c("A","C","G","T")), "N"))
  t("iupac_code: R alone -> R (passthrough)", identical(iupac_code("R"), "R"))
  
  ## ---- strict matching (max_mismatches = 0) ----
  m <- function(p, s) primer_match(p, s, 0L)$start
  t("match: ACGT in AAACGTTTACGT -> 3,9",   identical(m("ACGT", "AAACGTTTACGT"), c(3L, 9L)))
  t("match: absent pattern -> none",        identical(m("ACGTACGTACGT", "AAACGTTTACGT"), integer(0)))
  t("match: single A -> 1,2,3,9",           identical(m("A", "AAACGTTTACGT"), c(1L,2L,3L,9L)))
  t("match: A vs R (overlap: match)",       identical(m("A", "R"), 1L))
  t("match: A vs Y (disjoint: no match)",   identical(m("A", "Y"), integer(0)))
  t("match: A vs N (wildcard: match)",      identical(m("A", "N"), 1L))
  t("match: R in AG -> 1,2",                identical(m("R", "AG"), c(1L, 2L)))
  t("match: R in C -> none",                identical(m("R", "C"), integer(0)))
  t("match: N in ACGT -> 1,2,3,4",          identical(m("N", "ACGT"), c(1L,2L,3L,4L)))
  t("match: RAC in RACYG -> 1",             identical(m("RAC", "RACYG"), 1L))
  t("match: ACY in RACYG -> 2",             identical(m("ACY", "RACYG"), 2L))
  t("match: subject gap blocks the window", identical(m("ACGT", "AC-GTACGT"), 6L))
  t("match: degenerate subject: RY vs GT",  identical(m("RY", "GT"), 1L))
  t("match: degenerate subject: RY vs AC",  identical(m("RY", "AC"), 1L))
  t("match: H vs G (H={A,C,T}: no)",        identical(m("H", "G"), integer(0)))
  t("match: V vs G (V={A,C,G}: yes)",       identical(m("V", "G"), 1L))
  t("match: H vs T (yes)",                  identical(m("H", "T"), 1L))
  t("match: lowercase pattern ok",          identical(m("acgt", "aaacgtttacgt"), c(3L, 9L)))
  t("match: primer with a gap errors",      errs(primer_match("AC-GT", "ACGTA")))
  
  ## ---- over-match canaries (must be 0 hits at tolerance 0) ----
  t("canary: (ACGT)x4 in T20 -> 0",         identical(m(paste(rep("ACGT", 4L), collapse = ""), strrep("T", 20L)), integer(0)))
  t("canary: A16 in (AC)x8 -> 0",           identical(m(strrep("A", 16L), "ACACACACACACACACAC"), integer(0)))
  t("canary: pattern longer than subject",  identical(m("ACGTACGTACGT", "ACGT"), integer(0)))
  
  ## ---- vacuous N-block canaries (real site must beat an early N-block) ----
  bm_nb <- best_match("ACGT", paste0(strrep("N", 10L), "ACGTACGT"), 0L)
  t("best_match: real site beats early N-block",
    nrow(bm_nb) == 1L && bm_nb$start[1] == 11L)
  t("best_match: all-N subject still maps (leftmost)",
    identical(best_match("ACGT", strrep("N", 12L), 0L)$start, 1L))
  
  ## ---- mismatch tolerance ----
  t("tol=0: ACGT in TAAGTTACG -> none",     nrow(primer_match("ACGT", "TAAGTTACG", 0L)) == 0L)
  r1 <- primer_match("ACGT", "TAAGTTACG", 1L)
  t("tol=1: ACGT in TAAGTTACG -> start 2, 1 mm",
    nrow(r1) == 1L && r1$start[1] == 2L && r1$n_mismatches[1] == 1L)
  t("min_mm: ACGT in TAAGTTACG == 1",       identical(primer_min_mismatches("ACGT", "TAAGTTACG"), 1L))
  t("min_mm: ACGT in T20 == 3",             identical(primer_min_mismatches("ACGT", strrep("T", 20L)), 3L))
  t("min_mm: pattern cannot fit -> NA",     is.na(primer_min_mismatches("ACGT", "ACG")))
  bm <- best_match("ACGT", "TTACGTACGT", 0L)
  t("best_match: tie -> leftmost (3 of 3,7)", nrow(bm) == 1L && bm$start[1] == 3L)
  
  ## ---- gregexpr(fixed) equivalence (pure-DNA, non-self-overlapping) ----
  ## For a pattern with no border (no proper prefix that is also a suffix),
  ## occurrences cannot overlap, so gregexpr's non-overlapping scan finds
  ## exactly the same windows as the all-window scan. Patterns WITH borders
  ## (e.g. "CGAATTCG" border "CG", "ACGTACGT" border "ACGT") are excluded
  ## deliberately: gregexpr skips a second overlapping hit by design, which
  ## the window scan correctly does not.
  set.seed(42L)
  txt <- paste(sample(c("A","C","G","T"), 200L, replace = TRUE), collapse = "")
  ok_g <- TRUE
  for (pt in c("ACGT", "A", "GTCGAC", "TACCGA")) {
    g <- gregexpr(pt, txt, fixed = TRUE)[[1L]]
    g <- if (g[1L] == -1L) integer(0L) else as.integer(g)
    ok_g <- ok_g && identical(as.integer(m(pt, txt)), g)
  }
  t("gregexpr(fixed) equivalence on random 200-mer", ok_g)
  
  ## ---- base_consensus (gaps never win) ----
  ## Contract: most frequent NON-GAP symbol per column; tie -> IUPAC union;
  ## gaps never win; zero-coverage column -> '-' with occ 0.
  ## Returns list(cons, occ). All comparisons use whole-string matching on
  ## $cons.
  b1 <- base_consensus(c("ACGT", "ACGT"))
  t("base_consensus: uniform symbols + occ",
    identical(b1$cons, "ACGT") && identical(b1$occ, c(2L, 2L, 2L, 2L)))
  b2 <- base_consensus(c("AACT", "ACCT", "ACGT"))
  t("base_consensus: majority, not union (A,C,C -> C not M)",
    identical(b2$cons, "ACCT") && identical(b2$occ, c(3L, 3L, 3L, 3L)))
  b3 <- base_consensus(c("AC--", "CA--"))
  t("base_consensus: A/C tie in col 1 -> M",
    identical(b3$cons, "MM--") && identical(b3$occ, c(2L, 2L, 0L, 0L)))
  b4 <- base_consensus(c("ACGT", "TGCA", "GTAC", "CATG"))
  t("base_consensus: full 4-way tie -> N",
    identical(b4$cons, "NNNN") && identical(b4$occ, c(4L, 4L, 4L, 4L)))
  b5 <- base_consensus(c("RN--", "YN--"))
  t("base_consensus: IUPAC tie R+Y -> N (union of base sets)",
    identical(b5$cons, "NN--") && identical(b5$occ, c(2L, 2L, 0L, 0L)))
  b6 <- base_consensus(c("A-----", "C-----"))
  t("base_consensus: gap-dominated col still yields a base (A/C tie -> M)",
    identical(b6$cons, "M-----") && identical(b6$occ, c(2L, 0L, 0L, 0L, 0L, 0L)))
  b7 <- base_consensus(c("A-CG", "A-CG"))
  t("base_consensus: all-gap column -> '-' with occ 0",
    identical(b7$cons, "A-CG") && identical(b7$occ, c(2L, 0L, 2L, 2L)))
  b8 <- base_consensus(c("ACGT", "A-GT"))
  ## hand-verified: col 1 A x2 -> A(2); col 2 C x1 vs gap x1 -> base wins
  ## (gaps never win) -> C, occ 1; col 3 G x2 -> G(2); col 4 T x2 -> T(2)
  t("base_consensus: base/gap tie -> base wins (C beats gap)",
    identical(b8$cons, "ACGT") && identical(b8$occ, c(2L, 1L, 2L, 2L)))
  b9 <- base_consensus(c("RAC", "RAC"))
  t("base_consensus: IUPAC passthrough (RAC)",
    identical(b9$cons, "RAC") && identical(b9$occ, c(2L, 2L, 2L)))
  b10 <- base_consensus("A-CT")
  t("base_consensus: single-seq passthrough, gap kept",
    identical(b10$cons, "A-CT") && identical(b10$occ, c(1L, 0L, 1L, 1L)))
  b11 <- base_consensus(c("A--A", "A-C-", "A-T-", "A-G-"))
  ## hand-verified (gaps never win): col 1 A x4 -> A(4);
  ## col 2 all gaps -> '-'(0);
  ## col 3 C/T/G 3-way tie, gap ignored -> B (union {C,G,T}), occ 3;
  ## col 4 A x1 beats gap x3 -> A, occ 1
  t("base_consensus: mixed cols -> 'A-BA' (col 3 tie -> B; col 4 lone A wins)",
    identical(b11$cons, "A-BA") && identical(b11$occ, c(4L, 0L, 3L, 1L)))
  
  t("base_consensus: ragged input errors",   errs(base_consensus(c("ACGT", "ACG"))))
  t("base_consensus: no sequences errors",   errs(base_consensus(character(0))))
  t("base_consensus: invalid character errors",
    errs(base_consensus(c("ACGX", "ACGT"))))
  t("consensus_of_alignment: wrapper returns the cons string",
    identical(consensus_of_alignment(c("ACGT", "A-GT")), "ACGT"))
  
  ## ---- gap_inclusive_consensus (gaps can win) ----
  ## Contract: most frequent symbol over bases + IUPAC codes + the gap;
  ## tie including the gap -> '-'; base-only tie -> IUPAC union;
  ## zero-coverage column -> '-' with occ 0. Returns list(cons, occ).
  g1 <- gap_inclusive_consensus(c("ACGT", "ACGT"))
  t("gap_incl: uniform symbols + occ",
    identical(g1$cons, "ACGT") && identical(g1$occ, c(2L, 2L, 2L, 2L)))
  g2seqs <- c("A-----", "C-----", "------", "------", "------", "------")
  g2 <- gap_inclusive_consensus(g2seqs)
  ## hand-verified: col 1 gap x4 beats A x1, C x1 -> '-' (gap wins), occ 2;
  ## cols 2-6 all gaps -> '-' occ 0
  t("gap_incl: gap-dominated col -> gap wins ('-'), occ kept",
    identical(g2$cons, "------") && identical(g2$occ, c(2L, 0L, 0L, 0L, 0L, 0L)))
  g3 <- gap_inclusive_consensus(c("AC", "A-"))
  ## hand-verified: col 1 A x2 -> A(2); col 2 C x1 vs gap x1 TIE including
  ## the gap -> '-' (gap wins), occ 1. (Same input, base_consensus gives "AC".)
  t("gap_incl: base/gap tie -> gap wins",
    identical(g3$cons, "A-") && identical(g3$occ, c(2L, 1L)))
  g4 <- gap_inclusive_consensus(c("AC--", "CA--"))
  t("gap_incl: base tie without gap -> IUPAC union (M)",
    identical(g4$cons, "MM--") && identical(g4$occ, c(2L, 2L, 0L, 0L)))
  g5 <- gap_inclusive_consensus(c("ACGT", "ACGT", "----"))
  ## hand-verified: each col base x2 beats gap x1 -> base, occ 2
  t("gap_incl: 2 bases vs 1 gap -> base wins",
    identical(g5$cons, "ACGT") && identical(g5$occ, c(2L, 2L, 2L, 2L)))
  g6 <- gap_inclusive_consensus(c("RN--", "YN--"))
  t("gap_incl: IUPAC tie R+Y -> N (union of base sets)",
    identical(g6$cons, "NN--") && identical(g6$occ, c(2L, 2L, 0L, 0L)))
  t("gap_incl: ragged input errors",         errs(gap_inclusive_consensus(c("ACGT", "ACG"))))
  t("gap_incl: no sequences errors",         errs(gap_inclusive_consensus(character(0))))
  t("gap_incl: invalid character errors",    errs(gap_inclusive_consensus(c("ACGX", "ACGT"))))
  
  ## ---- consensus dispatcher + rule selector ----
  t("consensus: gaps_never_win routes to base_consensus",
    identical(consensus(c("AC--", "CA--"), "gaps_never_win"), base_consensus(c("AC--", "CA--"))))
  t("consensus: gaps_can_win routes to gap_inclusive_consensus",
    identical(consensus(g2seqs, "gaps_can_win"), gap_inclusive_consensus(g2seqs)))
  t("consensus: default rule is gaps_never_win",
    identical(consensus(g2seqs), base_consensus(g2seqs)))
  t("consensus: invalid rule errors",        errs(consensus("ACGT", "bogus")))
  t("rule_for: default when no override",
    identical(consensus_rule_for("X|12S", list("TaxonA|12S" = "gaps_can_win"), "gaps_never_win"),
              "gaps_never_win"))
  t("rule_for: override wins",
    identical(consensus_rule_for("TaxonA|12S", list("TaxonA|12S" = "gaps_can_win"), "gaps_never_win"),
              "gaps_can_win"))
  t("rule_for: invalid override errors",
    errs(consensus_rule_for("X|12S", list("X|12S" = "bogus"), "gaps_never_win")))
  
  
  ## ---- reverse complement ----
  rc <- reverse_complement
  t("rc: ACGT -> ACGT (palindrome)",        identical(rc("ACGT"), "ACGT"))
  t("rc: AAAA -> TTTT",                     identical(rc("AAAA"), "TTTT"))
  t("rc: RYKW -> WMRY",                     identical(rc("RYKW"), "WMRY"))
  t("rc: BHDV -> BHDV (palindrome)",        identical(rc("BHDV"), "BHDV"))
  t("rc: lowercase input -> uppercase",     identical(rc("acgt"), "ACGT"))
  t("rc: idempotent rc(rc(x)) = x",         identical(rc(rc("RYKMBDHVNACGTSW")), "RYKMBDHVNACGTSW"))
  t("rc: gap preserved",                    identical(rc("AC-GT"), "AC-GT"))
  t("rc: invalid character errors",         errs(rc("ACGX")))
  
  ## ---- small helpers ----
  t("gap_positions",                        identical(gap_positions("AC-GT"), c(1L,2L,4L,5L)))
  t("strip_gaps",                           identical(strip_gaps("AC-GT"), "ACGT"))
  t("extract_window",                       identical(extract_window(c("ACGTACGT","TTTTTTTT"), 3L, 6L), c("GTAC","TTTT")))
  t("extract_window: out-of-range errors",  errs(extract_window("ACGT", 1L, 9L)))
  
  ## ---- to_wsl_path (WSL path derivation) ----
  t("to_wsl_path: C:/Users/x/Code/ -> /mnt/c/Users/x/Code/",
    identical(to_wsl_path("C:/Users/x/Code/"), "/mnt/c/Users/x/Code/"))
  t("to_wsl_path: backslashes normalized",
    identical(to_wsl_path("C:\\Users\\x"), "/mnt/c/Users/x"))
  t("to_wsl_path: POSIX passthrough",       identical(to_wsl_path("/home/u/x"), "/home/u/x"))
  
  ## ---- FASTA reader + dedup ----
  tmp <- tempfile(fileext = ".fasta")
  con <- file(tmp, "w")
  for (l in c(">seq1 haplotype", "acgt", "acgt", "", ">seq2", "TTTTTTTT",
              ">seq1 haplotype", "aAaAaAaA"))
    writeLines(paste0(l, "\r"), con)        # CRLF endings, a blank line, lowercase
  close(con)
  f <- read_fasta(tmp)
  t("fasta: 3 records / 2 unique names",    f$records == 3L && f$unique_names == 2L)
  t("fasta: names + wrapped/CRLF/lowercase",
    identical(f$names, c("seq1 haplotype","seq2","seq1 haplotype")) &&
      identical(f$seqs, c("ACGTACGT","TTTTTTTT","AAAAAAAA")))
  t("fasta: uniform width, not ragged",     all(f$width == 8L) && !isTRUE(f$ragged))
  t("fasta: duplicate name flagged",        identical(f$dup_names, "seq1 haplotype"))
  f2 <- fasta_dedup(f)
  t("fasta_dedup: keep first, drop rest",
    f2$records == 2L && identical(f2$names, c("seq1 haplotype","seq2")) &&
      identical(f2$seqs, c("ACGTACGT","TTTTTTTT")))
  unlink(tmp)
  tmp2 <- tempfile(fileext = ".fasta"); writeLines("ACGT", tmp2)
  t("fasta: no header errors",              errs(read_fasta(tmp2))); unlink(tmp2)
  tmp3 <- tempfile(fileext = ".fasta"); writeLines(c("ACGT", ">seq1", "TTTT"), tmp3)
  t("fasta: line before header errors",     errs(read_fasta(tmp3))); unlink(tmp3)
  
  ## ---- audit helpers ----
  t("check: PASS path returns TRUE",        isTRUE(tryCatch(check("selftest", TRUE), error = function(e) NULL)))
  t("check: FAIL path stops the run",       errs(check("selftest", FALSE, "deliberate")))
  t("note: returns its message",            identical(tryCatch(note("selftest", "hello"), error = function(e) NULL), "hello"))
  
  lr0 <- log_row_make("r0", "0", "PASS", "C:/x", metrics = "n=1;m=2")
  t("log_row_make: fixed 7-column schema",
    nrow(lr0) == 1L &&
      identical(names(lr0), c("run_id", "parent_run_id", "step", "timestamp",
                              "status", "base", "metrics")) &&
      is.na(lr0$parent_run_id[1L]) && lr0$metrics[1L] == "n=1;m=2")
  lr1 <- rbind(lr0, log_row_make("r1", "1", "PASS", "C:/x", parent_run = "r0"))
  t("log_row_make: different steps rbind cleanly",
    nrow(lr1) == 2L && ncol(lr1) == 7L)
  
  ## ---- Span gate tests ----
  ## (1) pair_span: signed facing span of the two FINAL windows.
  t("span: pair_span facing 68 / anti-facing -1032",
    identical(pair_span(900L, 907L, 950L, 967L), 68L) &&
      identical(pair_span(2000L, 2017L, 950L, 967L), -1032L))
  
  ## (2) span_keep boundaries (hand-computed).
  ## Anchor R [950,967], F cand Lp 18, max 500: kept span range [36, 500].
  ## s = 468 kept (span 500); s = 467 dropped (span 501); s = 933 dropped
  ## (overlap); s = 932 kept (min span 36).
  ## Anchor F [100,117], R cand Lp 18: s = 582 kept (span 500); s = 583
  ## dropped (span 501); s = 117 dropped (overlap); s = 118 kept (min span
  ## 36).
  t("span: span_keep boundaries (F vs R anchor; R vs F anchor)",
    isTRUE(span_keep(468L, 18L, 950L, 967L, "F", 500L)) &&
      !isTRUE(span_keep(467L, 18L, 950L, 967L, "F", 500L)) &&
      !isTRUE(span_keep(933L, 18L, 950L, 967L, "F", 500L)) &&
      isTRUE(span_keep(932L, 18L, 950L, 967L, "F", 500L)) &&
      isTRUE(span_keep(582L, 18L, 100L, 117L, "R", 500L)) &&
      !isTRUE(span_keep(583L, 18L, 100L, 117L, "R", 500L)) &&
      !isTRUE(span_keep(117L, 18L, 100L, 117L, "R", 500L)) &&
      isTRUE(span_keep(118L, 18L, 100L, 117L, "R", 500L)))
  
  ## Gate fixtures (hand-computed): d8_P has no T, so every background T
  ## mismatches every primer position (18-mm background windows), and every
  ## flanking window of a designed region carries >= 7 mm (checked against
  ## the period-3 primer: shifts not divisible by 3 turn the matched
  ## positions into mismatches too), so primer_match(., 6L) sees ONLY the
  ## designed regions.
  d8_P      <- "ACGACGACGACGACGACG"
  d8_bg3000 <- strrep("T", 3000L)
  d8_bg500  <- strrep("T", 500L)
  d8_bg300  <- strrep("T", 300L)
  ## reg_high: 12 concrete matches + 6 T mismatches (6 mm, support 12)
  d8_pc <- strsplit(d8_P, "", fixed = TRUE)[[1L]]
  d8_rh <- d8_pc
  d8_rh[c(3L, 5L, 7L, 9L, 11L, 13L)] <- "T"
  d8_regHigh <- paste(d8_rh, collapse = "")
  ## reg_low: 6 concrete matches + 6 N (vacuous) + 6 T mismatches
  ## (6 mm, support 6)
  d8_rl <- character(18L)
  d8_rl[c(1L, 4L, 7L, 10L, 13L, 16L)] <- "N"
  d8_rl[c(2L, 5L, 8L, 11L, 14L, 17L)] <- d8_pc[c(2L, 5L, 8L, 11L, 14L, 17L)]
  d8_rl[c(3L, 6L, 9L, 12L, 15L, 18L)] <- "T"
  d8_regLow <- paste(d8_rl, collapse = "")
  
  ## (3) Distant-window exclusion: 3000-mer subject with two 6-mm candidate
  ## windows. The ungated best (position 597, support 12, facing span 2359)
  ## is outside the 500-bp span gate; the in-span window (position 2700,
  ## support 6, facing span 256) becomes the gated winner. One window excluded.
  d8_s3000 <- paste0(substr(d8_bg3000, 1L, 596L), d8_regHigh,
                     substr(d8_bg3000, 615L, 2699L), d8_regLow,
                     substr(d8_bg3000, 2718L, 3000L))
  d8_c3 <- primer_match(d8_P, d8_s3000, 6L)
  d8_g3 <- gated_best_match(d8_c3, d8_P, d8_s3000, 2938L, 2955L, "F", 500L)
  t("gate: distant window excluded by span (597 -> 2700, excluded 1)",
    nrow(d8_c3) == 2L && identical(d8_c3$start, c(597L, 2700L)) &&
      identical(d8_c3$n_mismatches, c(6L, 6L)) &&
      d8_g3$ungated_best$start[1L] == 597L &&
      d8_g3$best$start[1L] == 2700L &&
      identical(d8_g3$n_excluded, 1L))
  
  ## (4) No in-span candidate -> best = 0 rows (step 1 rejection); the
  ## ungated winner is still returned for span_gate_mm reporting.
  d8_g4 <- gated_best_match(d8_c3, d8_P, d8_s3000, 100L, 117L, "F", 500L)
  t("gate: no in-span candidate -> empty best (rejection)",
    nrow(d8_g4$best) == 0L &&
      nrow(d8_g4$ungated_best) == 1L &&
      d8_g4$ungated_best$start[1L] == 597L &&
      identical(d8_g4$n_excluded, 2L))
  
  ## (5) No-op case: ungated best already in span -> identical frames,
  ## n_excluded 0. 300-mer, one 6-mm window at 100, anchor R [200,217]
  ## (facing span 118 <= 500).
  d8_s300 <- paste0(substr(d8_bg300, 1L, 99L), d8_regHigh,
                    substr(d8_bg300, 118L, 300L))
  d8_c5 <- primer_match(d8_P, d8_s300, 6L)
  d8_g5 <- gated_best_match(d8_c5, d8_P, d8_s300, 200L, 217L, "F", 500L)
  t("gate: no-op (ungated best in span, excluded 0)",
    nrow(d8_c5) == 1L && d8_c5$start[1L] == 100L &&
      identical(d8_g5$best, d8_g5$ungated_best) &&
      identical(d8_g5$n_excluded, 0L))
  
  ## (6) Symmetric R-candidate vs F-anchor (mirror of 3): candidates 150
  ## (6 mm, support 6, facing span 68 — in span) and 600 (6 mm, support 12,
  ## facing span 518 — out of span); anchor F [100,117]. Ungated winner 600
  ## (support 12 > 6) -> gated winner 150.
  d8_s3000b <- paste0(substr(d8_bg3000, 1L, 149L), d8_regLow,
                      substr(d8_bg3000, 168L, 599L), d8_regHigh,
                      substr(d8_bg3000, 618L, 3000L))
  d8_c6 <- primer_match(d8_P, d8_s3000b, 6L)
  d8_g6 <- gated_best_match(d8_c6, d8_P, d8_s3000b, 100L, 117L, "R", 500L)
  t("gate: R-candidate vs F-anchor (600 -> 150, excluded 1)",
    nrow(d8_c6) == 2L && identical(d8_c6$start, c(150L, 600L)) &&
      d8_g6$ungated_best$start[1L] == 600L &&
      d8_g6$best$start[1L] == 150L &&
      identical(d8_g6$n_excluded, 1L))
  
  ## (7) Per-orientation gating: two frames gated against the same anchor.
  ## The out-of-span frame (597) empties; the in-span frame (2700) survives.
  ## Orientation selection between gated winners is handled in step 1.
  d8_g7a <- gated_best_match(d8_c3[d8_c3$start == 597L, , drop = FALSE],
                             d8_P, d8_s3000, 2938L, 2955L, "F", 500L)
  d8_g7b <- gated_best_match(d8_c3[d8_c3$start == 2700L, , drop = FALSE],
                             d8_P, d8_s3000, 2938L, 2955L, "F", 500L)
  t("gate: per-orientation frames gated independently",
    nrow(d8_g7a$best) == 0L &&
      nrow(d8_g7b$best) == 1L && d8_g7b$best$start[1L] == 2700L &&
      identical(d8_g7b$n_excluded, 0L))
  
  ## (8) Tie-break still applies WITHIN the gated set: two in-span, equal
  ## mm candidates -> higher concrete support wins (300 beats the leftmost
  ## 150); two equal-support candidates -> leftmost (150).
  d8_s500  <- paste0(substr(d8_bg500, 1L, 149L), d8_regLow,
                     substr(d8_bg500, 168L, 299L), d8_regHigh,
                     substr(d8_bg500, 318L, 500L))
  d8_c8 <- primer_match(d8_P, d8_s500, 6L)
  d8_g8 <- gated_best_match(d8_c8, d8_P, d8_s500, 400L, 417L, "F", 500L)
  d8_s500b <- paste0(substr(d8_bg500, 1L, 149L), d8_regHigh,
                     substr(d8_bg500, 168L, 299L), d8_regHigh,
                     substr(d8_bg500, 318L, 500L))
  d8_c8b <- primer_match(d8_P, d8_s500b, 6L)
  d8_g8b <- gated_best_match(d8_c8b, d8_P, d8_s500b, 400L, 417L, "F", 500L)
  t("gate: tie-break within span (support, then leftmost)",
    nrow(d8_c8) == 2L && d8_g8$best$start[1L] == 300L &&
      nrow(d8_c8b) == 2L && d8_g8b$best$start[1L] == 150L &&
      identical(d8_g8b$n_excluded, 0L))
  
  ## ---- CONFIG shape ----
  t("CONFIG: SCORE_THRESHOLD numeric",      is.numeric(CONFIG$SCORE_THRESHOLD))
  t("CONFIG: MAX_MAP_MISMATCHES integer >= 0",
    is.integer(CONFIG$MAX_MAP_MISMATCHES) && CONFIG$MAX_MAP_MISMATCHES >= 0L)
  t("CONFIG: MIN_OCC integer >= 1",
    is.integer(CONFIG$MIN_OCC) && CONFIG$MIN_OCC >= 1L)
  t("CONFIG: CONSENSUS_RULE valid",
    CONFIG$CONSENSUS_RULE %in% CONSENSUS_RULES_VALID)
  t("CONFIG: CONSENSUS_RULE_OVERRIDES values valid",
    all(unlist(CONFIG$CONSENSUS_RULE_OVERRIDES) %in% CONSENSUS_RULES_VALID))
  t("CONFIG: MAX_AMP_SPAN integer >= 1",
    is.integer(CONFIG$MAX_AMP_SPAN) && CONFIG$MAX_AMP_SPAN >= 1L)
  
  cat("===========================================================================\n")
  if (n_fail > 0L)
    stop("utils self-test FAILED on ", n_fail, " case(s): ", paste(failed, collapse = " ; "),
         "\n      The core utilities are untrustworthy — do not run any pipeline step.", call. = FALSE)
  cat(sprintf("utils self-test: %d / %d passed -> core verified, safe to run steps\n\n",
              n_ok, n_ok + n_fail))
  invisible(TRUE)
}

utils_self_test()

# ============================================================================
# END utils.R
# ============================================================================