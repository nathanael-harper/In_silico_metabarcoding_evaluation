# ============================================================================
# step0_setup_audit.r — STEP 0: INPUT SETUP + AUDIT  (Windows R)
# ----------------------------------------------------------------------------
# Primer mapping / DNA metabarcoding pipeline.
# This is the setup and input-validation step, preceding five analysis
# steps (steps 1-5).
#
# WHAT THIS STEP DOES (and only this)
#   Validates the three inputs, builds the assay pair table, inventories
#   every MSA (each FASTA record = one species record for that taxon x gene),
#   and hard-stops with named offenders if anything does not match. When it
#   passes, steps 1-5 may run — with any primer set, any genes, any MSA
#   dimensions. No expected counts are hardcoded anywhere in this file.
#
#   Consensus rule contract: for every (Taxa, Gene) template, the
#   population-consensus rule is selected from CONFIG$CONSENSUS_RULE (global
#   default) with per-template CONFIG$CONSENSUS_RULE_OVERRIDES (keys
#   'Taxa|Gene'), and PUBLISHED in the manifest as section 'consensus_rule'
#   (one row per template). Step 1 reads the rules back from the manifest,
#   CHECKs they equal the CONFIG-derived selection, and dispatches consensus()
#   accordingly — so the engine choice is part of the continuity contract.
#
#   Span gate parameter: CONFIG$MAX_AMP_SPAN is published in the manifest as
#   config|max_pair_span so step 1 can CHECK it against its own CONFIG.
#   [0.3] NOTEs any complete pair whose primer lengths sum to more than
#   MAX_AMP_SPAN: the span gate can never be satisfied for such a pair, so it
#   can only map with both strands tolerant (ungated) or not at all. NOTE
#   only — never stops the run. No mapping behavior changes in step 0; the
#   gate itself lives in step 1.
#
# INPUTS (exactly three files, all under BASE):
#   1. Primer file (primers_input2.csv)  columns: name, seq, gene
#                                    one row per primer strand; name ends in
#                                    _F or _R; seq = gap-free IUPAC design
#                                    sequence
#   2. Taxon directory (taxon_template_directory_2.csv)
#                                    columns: Taxa, Gene, filepath
#                                    filepath relative to BASE; each row =
#                                    one MSA; each MSA record = one species
#   3. The MSA FASTA files named in (2)
#
# OUTPUTS (BASE/results/, overwritten on each run; log is appended):
#   00_run_manifest.csv        run id, timestamp, config (consensus rule
#                              default + overrides, max_pair_span), input
#                              fingerprints, consensus_rule section (one row
#                              per template) — the continuity contract that
#                              steps 1-5 gate on
#   01_pair_table.csv          one row per complete assay pair (F + R);
#                              columns: assay, gene, f_name, f_seq, r_name,
#                              r_seq
#   02_template_inventory.csv  one row per MSA (records, width,
#                              consensus_rule, ...)
#   run_audit_log.csv          one row per completed step, keyed by run_id
#
# AUDITS (numbered; CHECK stops the run on failure, NOTE never stops):
#   [0.1]  BASE is a directory
#   [0.2]  primer file: exists, schema, rows, empty cells, unique names,
#          _F/_R suffix, gap-free IUPAC
#   [0.3]  assay pairs: complete vs incomplete (named), gene conflicts;
#          long-primer NOTE (Lp_F + Lp_R > MAX_AMP_SPAN)
#   [0.4]  gene targets + complete pairs per gene
#   [0.5]  taxon directory: exists, schema, duplicates, files exist
#   [0.6]  MSA parse per file: records, unique names, width (ragged = stop);
#          single-record MSA = NOTE; consensus rule validated and selected
#   [0.7]  duplicate FASTA names: keep first, drop rest (NOTE)
#   [0.8]  orphan templates: (Taxa, Gene) with no primers = stop
#   [0.9]  coverage grid: every (Taxa x gene) has an MSA
#   [0.10] outputs written + re-read (round-trip check); manifest
#          consensus_rule section validated
#
# RUN:  Rscript step0_setup_audit.r     (or paste the whole file; every
#       conditional is braced, so a paste cannot orphan an else)
#       Steps 1, 3, 4, 5: Windows R, same or fresh session.
#       Step 2: Ubuntu/WSL R (PrimerMiner dependency).
# ============================================================================

## ---- BASE — the ONE line to edit (steps 1-5 carry the same line and
## ---- CHECK it against the manifest they read) ----
BASE <- "C:/Users/Nathanael/Desktop/Vernal_Pool_Parameters1/Final_run"

source(file.path(BASE, "utils.R"))    # runs the self-test battery; hard-stops on any failure

## ---- Consensus rule config, normalized ONCE (the banner, the [0.6]
## ---- selection and the manifest packing all read these) ----
OV      <- as.list(CONFIG$CONSENSUS_RULE_OVERRIDES)   # named list: "Taxa|Gene" -> rule
ov_keys <- sort(names(OV))
ov_str  <- if (length(OV) == 0L) {
  "(none)"
} else {
  paste(sprintf("%s=%s", ov_keys, as.character(unlist(OV[ov_keys]))), collapse = ";")
}

run_id    <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat("============================================================================\n")
cat("STEP 0 - INPUT SETUP + AUDIT\n")
cat("  run_id :", run_id, "\n")
cat("  base   :", BASE, "\n")
## The span-gate parameter is printed with the other config values so
## every run banner carries the full config contract.
cat("  config : score_threshold =", CONFIG$SCORE_THRESHOLD,
    " | max_map_mismatches =", CONFIG$MAX_MAP_MISMATCHES,
    " | max_amp_span =", CONFIG$MAX_AMP_SPAN, "\n")
cat("           consensus_rule =", CONFIG$CONSENSUS_RULE,
    " | overrides =", ov_str, "\n")
cat("============================================================================\n")

## ---- [0.1] BASE -----------------------------------------------------------
cat("\n[0.1] base directory\n")
check("[0.1] BASE is a directory", dir.exists(BASE),
      paste0("BASE does not exist or is not a directory: ", BASE))
OUT <- file.path(BASE, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- [0.2] primer file ----------------------------------------------------
cat("\n[0.2] primer file (primers_input2.csv)\n")
PRIM_CSV <- file.path(BASE, "primers_input2.csv")
check("[0.2] file exists", file.exists(PRIM_CSV),
      paste0("not found: ", PRIM_CSV))
prim <- read.csv(PRIM_CSV, stringsAsFactors = FALSE, check.names = FALSE)
prim[] <- lapply(prim, trimws)
check("[0.2] schema is exactly name, seq, gene",
      identical(names(prim), c("name", "seq", "gene")),
      paste0("columns found: ", paste(names(prim), collapse = ", ")))
check("[0.2] at least one row", nrow(prim) >= 1L, "primer file is empty")
n_empty <- sum(prim$name == "" | prim$seq == "" | prim$gene == "")
check("[0.2] no empty name/seq/gene cells", n_empty == 0L,
      paste0(n_empty, " row(s) with an empty cell"))
dups <- prim$name[duplicated(prim$name)]
check("[0.2] primer names unique", length(dups) == 0L,
      paste0("duplicate name(s): ", paste(unique(dups), collapse = ", ")))
bad_sfx <- prim$name[!grepl("_(F|R)$", prim$name)]
check("[0.2] every name ends in _F or _R", length(bad_sfx) == 0L,
      paste0("bad name(s): ", paste(bad_sfx, collapse = ", ")))
bad_iupac <- character(0)
for (i in seq_len(nrow(prim))) {
  err <- tryCatch({ validate_iupac(prim$seq[i], "primer", allow_gap = FALSE); NULL },
                  error = function(e) conditionMessage(e))
  if (!is.null(err)) bad_iupac <- c(bad_iupac, prim$name[i])
}
check("[0.2] every seq is gap-free IUPAC", length(bad_iupac) == 0L,
      paste0("primer(s): ", paste(bad_iupac, collapse = ", ")))
prim$strand <- substr(prim$name, nchar(prim$name), nchar(prim$name))
prim$pair   <- sub("_(F|R)$", "", prim$name)
cat(sprintf("    %d primer rows loaded\n", nrow(prim)))

## ---- [0.3] assay pairs ----------------------------------------------------
cat("\n[0.3] assay pairs\n")
pair_keys <- unique(prim$pair)                       # first-appearance order
nF <- vapply(pair_keys, function(p) sum(prim$strand[prim$pair == p] == "F"), integer(1))
nR <- vapply(pair_keys, function(p) sum(prim$strand[prim$pair == p] == "R"), integer(1))
pairs <- data.frame(pair = pair_keys, n_F = nF, n_R = nR, stringsAsFactors = FALSE)
pairs$gene <- vapply(pair_keys,
                     function(p) prim$gene[prim$pair == p][1L], character(1))
pairs$complete <- pairs$n_F == 1L & pairs$n_R == 1L
gene_conflict <- pair_keys[vapply(pair_keys,
                                  function(p) length(unique(prim$gene[prim$pair == p])) > 1L,
                                  logical(1))]
check("[0.3] no pair spans two gene targets", length(gene_conflict) == 0L,
      paste0("pair(s): ", paste(gene_conflict, collapse = ", ")))
check("[0.3] at least one complete pair", sum(pairs$complete) >= 1L,
      "no pair has both an _F and an _R primer")
if (any(!pairs$complete))
  for (i in which(!pairs$complete))
    note("[0.3]", sprintf("incomplete pair %s (%s): F=%d R=%d -> excluded",
                          pairs$pair[i], pairs$gene[i], pairs$n_F[i], pairs$n_R[i]))
cat(sprintf("    %d complete | %d incomplete (of %d pairs)\n",
            sum(pairs$complete), sum(!pairs$complete), nrow(pairs)))

## Long-primer guard, once the complete pairs are built.
## span_keep requires the implied product to be >= Lp_F + Lp_R, so a
## complete pair whose primer lengths sum to more than MAX_AMP_SPAN can
## NEVER satisfy the span gate: it can only map with both strands tolerant
## (ungated) or not at all. NOTE only (never stops).
for (i in which(pairs$complete)) {
  p    <- pairs$pair[i]
  fseq <- prim$seq[prim$pair == p & prim$strand == "F"]
  rseq <- prim$seq[prim$pair == p & prim$strand == "R"]
  lp_sum <- nchar(fseq[1L]) + nchar(rseq[1L])
  if (lp_sum > CONFIG$MAX_AMP_SPAN)
    note("[0.3]", sprintf("pair %s (%s): primer lengths sum to %d bp > MAX_AMP_SPAN=%d — the span gate can never be satisfied for this pair; it can only map with both strands tolerant (ungated) or not at all",
                          p, pairs$gene[i], lp_sum, CONFIG$MAX_AMP_SPAN))
}

## ---- [0.4] gene targets ---------------------------------------------------
cat("\n[0.4] gene targets\n")
genes <- sort(unique(prim$gene))
for (g in genes)
  cat(sprintf("    %-8s %d complete pair(s)\n", g,
              sum(pairs$gene == g & pairs$complete)))

## ---- [0.5] taxon template directory ---------------------------------------
cat("\n[0.5] taxon template directory (taxon_template_directory_2.csv)\n")
TAXON_CSV <- file.path(BASE, "taxon_template_directory_2.csv")
check("[0.5] file exists", file.exists(TAXON_CSV),
      paste0("not found: ", TAXON_CSV))
taxon <- read.csv(TAXON_CSV, stringsAsFactors = FALSE, check.names = FALSE)
taxon[] <- lapply(taxon, trimws)
check("[0.5] schema is exactly Taxa, Gene, filepath",
      identical(names(taxon), c("Taxa", "Gene", "filepath")),
      paste0("columns found: ", paste(names(taxon), collapse = ", ")))
check("[0.5] at least one row", nrow(taxon) >= 1L, "directory is empty")
check("[0.5] no empty Taxa/Gene/filepath cells",
      !any(taxon$Taxa == "" | taxon$Gene == "" | taxon$filepath == ""),
      "empty cell(s) in directory")
dupkey <- taxon$Taxa[duplicated(paste(taxon$Taxa, taxon$Gene))]
check("[0.5] no duplicate (Taxa, Gene) rows", length(dupkey) == 0L,
      paste0("duplicate(s): ", paste(unique(dupkey), collapse = ", ")))
taxon$full_path <- file.path(BASE, taxon$filepath)
missing_f <- taxon$full_path[!file.exists(taxon$full_path)]
check("[0.5] every MSA file exists under BASE", length(missing_f) == 0L,
      paste0("missing: ", paste(missing_f, collapse = ", ")))

## ---- [0.6]/[0.7] MSA parse + duplicate-name policy + consensus rule
cat("\n[0.6] MSA parse (each FASTA record = one species record) + consensus rule selection\n")

## Consensus rule config: the default and every override must be in the
## rule vocabulary, and override keys must be unique 'Taxa|Gene' strings.
## (utils' self-test already validated the CONFIG shape at source time; this
## re-asserts the vocabulary so the manifest is only ever published with a
## rule the matcher can dispatch on.)
check("[0.6] CONFIG$CONSENSUS_RULE in rule vocabulary",
      CONFIG$CONSENSUS_RULE %in% CONSENSUS_RULES_VALID,
      paste0("CONFIG$CONSENSUS_RULE is '", CONFIG$CONSENSUS_RULE,
             "', expected one of: ", paste(CONSENSUS_RULES_VALID, collapse = " / ")))
check("[0.6] CONFIG$CONSENSUS_RULE_OVERRIDES valid (named 'Taxa|Gene' keys, unique, valid rules)",
      length(OV) == 0L ||
        (all(!is.na(names(OV))) && !any(names(OV) == "") &&
           !any(duplicated(names(OV))) && all(grepl("\\|", names(OV))) &&
           all(vapply(OV, function(r) is.character(r) && length(r) == 1L &&
                        r %in% CONSENSUS_RULES_VALID, logical(1)))),
      paste0("each override must be a unique 'Taxa|Gene' key mapped to exactly one of: ",
             paste(CONSENSUS_RULES_VALID, collapse = " / "),
             " | current: ", ov_str))

TEMPLATES <- vector("list", nrow(taxon))
inv <- data.frame(
  Taxa = taxon$Taxa, Gene = taxon$Gene, filepath = taxon$filepath,
  records_raw = integer(nrow(taxon)), records = integer(nrow(taxon)),
  width = integer(nrow(taxon)), dups_dropped = integer(nrow(taxon)),
  first_name = character(nrow(taxon)), last_name = character(nrow(taxon)),
  consensus_rule = character(nrow(taxon)),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(taxon))) {
  f <- read_fasta(taxon$full_path[i])
  check(sprintf("[0.6] %s %s uniform width", taxon$Taxa[i], taxon$Gene[i]),
        !f$ragged,
        sprintf("%s has ragged widths %d..%d — an MSA must be equal-width",
                taxon$filepath[i], min(f$width), max(f$width)))
  f2 <- fasta_dedup(f)          # prints a [NOTE] naming any dropped duplicates
  TEMPLATES[[i]] <- f2
  inv$records_raw[i]  <- f$records
  inv$records[i]      <- f2$records
  inv$width[i]        <- f2$width[1L]
  inv$dups_dropped[i] <- f$records - f2$records
  inv$first_name[i]   <- f2$names[1L]
  inv$last_name[i]    <- f2$names[f2$records]
  ## Select this template's consensus rule (utils consensus_rule_for:
  ## the override wins when 'Taxa|Gene' matches, else the global default)
  tpl_key <- paste0(taxon$Taxa[i], "|", taxon$Gene[i])
  rule_i  <- tryCatch(
    consensus_rule_for(tpl_key, CONFIG$CONSENSUS_RULE_OVERRIDES, CONFIG$CONSENSUS_RULE),
    error = function(e) stop("[0.6] consensus_rule_for failed for '", tpl_key, "': ",
                             conditionMessage(e),
                             " — check the utils.R rule-selection API",
                             call. = FALSE))
  check(sprintf("[0.6] %s %s consensus rule valid", taxon$Taxa[i], taxon$Gene[i]),
        rule_i %in% CONSENSUS_RULES_VALID,
        paste0("consensus_rule_for returned '", rule_i,
               "', expected one of: ", paste(CONSENSUS_RULES_VALID, collapse = " / ")))
  inv$consensus_rule[i] <- rule_i
  cat(sprintf("    %-14s %-4s records=%-5d width=%-5d rule=%-14s first=%s | last=%s\n",
              taxon$Taxa[i], taxon$Gene[i], f2$records, f2$width[1L], rule_i,
              f2$names[1L], f2$names[f2$records]))
  if (f2$records == 1L)
    note("[0.6]", sprintf("%s %s has a single record — if this is a consensus file, point the directory at the real MSA",
                          taxon$Taxa[i], taxon$Gene[i]))
}

## ---- [0.8] orphan templates -------------------------------------------------
cat("\n[0.8] orphan templates (Taxa x Gene with no primers)\n")
orphan <- taxon[!taxon$Gene %in% genes, ]
check("[0.8] none", nrow(orphan) == 0L,
      paste0("template(s) for gene(s) not in the primer file: ",
             paste(paste(orphan$Taxa, orphan$Gene, sep = " x "), collapse = ", ")))

## ---- [0.9] coverage grid ----------------------------------------------------
cat("\n[0.9] coverage grid (every taxon x gene has an MSA)\n")
taxa <- sort(unique(taxon$Taxa))
cells <- expand.grid(Taxa = taxa, Gene = genes, stringsAsFactors = FALSE)
cells$have <- vapply(seq_len(nrow(cells)),
                     function(i) any(taxon$Taxa == cells$Taxa[i] &
                                       taxon$Gene == cells$Gene[i]),
                     logical(1))
for (i in seq_len(nrow(cells)))
  cat(sprintf("    %-14s x %-4s %s\n", cells$Taxa[i], cells$Gene[i],
              if (cells$have[i]) "OK" else "** MISSING **"))
missing_cells <- cells[!cells$have, ]
check("[0.9] grid complete", nrow(missing_cells) == 0L,
      paste0("no MSA for: ",
             paste(paste(missing_cells$Taxa, missing_cells$Gene, sep = " x "),
                   collapse = ", ")))

## ---- [0.10] outputs ----------------------------------------------------------
cat("\n[0.10] outputs (written to results/)\n")

## 01_pair_table.csv — one row per complete assay pair
pair_table <- data.frame(assay = character(0), gene = character(0),
                         f_name = character(0), f_seq = character(0),
                         r_name = character(0), r_seq = character(0),
                         stringsAsFactors = FALSE)
for (i in which(pairs$complete)) {
  frow <- prim[prim$pair == pairs$pair[i] & prim$strand == "F", ]
  rrow <- prim[prim$pair == pairs$pair[i] & prim$strand == "R", ]
  pair_table <- rbind(pair_table, data.frame(
    assay = pairs$pair[i], gene = pairs$gene[i],
    f_name = frow$name, f_seq = frow$seq,
    r_name = rrow$name, r_seq = rrow$seq, stringsAsFactors = FALSE))
}
write.csv(pair_table, file.path(OUT, "01_pair_table.csv"), row.names = FALSE)

## 02_template_inventory.csv (includes per-template consensus rule)
write.csv(inv, file.path(OUT, "02_template_inventory.csv"), row.names = FALSE)

## 00_run_manifest.csv — the continuity contract steps 1-5 gate on
mani <- data.frame(section = character(0), item = character(0),
                   value = character(0), stringsAsFactors = FALSE)
addm <- function(s, i, v) mani <<- rbind(mani, data.frame(
  section = s, item = i, value = as.character(v), stringsAsFactors = FALSE))
addm("run", "run_id", run_id)
addm("run", "timestamp", timestamp)
addm("run", "base", BASE)
addm("run", "utils_version", UTILS_VERSION)
addm("config", "score_threshold", CONFIG$SCORE_THRESHOLD)
addm("config", "max_map_mismatches", CONFIG$MAX_MAP_MISMATCHES)
## Publish the span-gate parameter next to the other config rows;
## step 1 [1.1] CHECKs it against its own CONFIG$MAX_AMP_SPAN.
addm("config", "max_pair_span", CONFIG$MAX_AMP_SPAN)
addm("config", "consensus_rule", CONFIG$CONSENSUS_RULE)
addm("config", "consensus_rule_overrides", ov_str)
addm("inputs", "primer_file", "primers_input2.csv")
addm("inputs", "primer_rows", nrow(prim))
addm("inputs", "primer_genes", paste(genes, collapse = ";"))
addm("inputs", "pairs_complete", sum(pairs$complete))
addm("inputs", "pairs_incomplete", sum(!pairs$complete))
addm("inputs", "taxon_file", "taxon_template_directory_2.csv")
addm("inputs", "taxon_rows", nrow(taxon))
for (i in seq_len(nrow(taxon))) {
  key <- paste0(taxon$Taxa[i], "|", taxon$Gene[i])
  f2 <- TEMPLATES[[i]]
  first_seq <- f2$seqs[1L]
  last_seq  <- f2$seqs[f2$records]
  addm("templates", key,
       paste0("records=", f2$records, ";width=", f2$width[1L],
              ";first_name=", f2$names[1L], ";last_name=", f2$names[f2$records],
              ";first_head40=", substr(first_seq, 1L, 40L),
              ";last_tail40=",
              substr(last_seq, pmax(1L, nchar(last_seq) - 39L), nchar(last_seq))))
  ## Per-template consensus rule (step 1 reads these back and CHECKs
  ## them against its own CONFIG-derived selection)
  addm("consensus_rule", key, inv$consensus_rule[i])
}
addm("coverage", "grid",
     sprintf("%d taxa x %d genes = %d cells, all covered",
             length(taxa), length(genes), nrow(cells)))
addm("status", "step0", "PASS")
write.csv(mani, file.path(OUT, "00_run_manifest.csv"), row.names = FALSE)

## run_audit_log.csv — appended (history across runs, keyed by run_id)
LOG <- file.path(OUT, "run_audit_log.csv")
log_row <- log_row_make(run_id, "0", "PASS", BASE,
                        metrics = paste0("n_primers=", nrow(prim),
                                         ";n_pairs_complete=", sum(pairs$complete),
                                         ";n_genes=", length(genes),
                                         ";n_templates=", nrow(taxon),
                                         ";n_msa_records=", sum(inv$records),
                                         ";n_templates_gaps_can_win=",
                                         sum(inv$consensus_rule == "gaps_can_win"),
                                         ";n_templates_gaps_never_win=",
                                         sum(inv$consensus_rule == "gaps_never_win")))
if (file.exists(LOG))
  log_row <- rbind(read.csv(LOG, stringsAsFactors = FALSE), log_row)
write.csv(log_row, LOG, row.names = FALSE)

## Round-trip: re-read everything we just wrote
pt_chk   <- read.csv(file.path(OUT, "01_pair_table.csv"), stringsAsFactors = FALSE)
inv_chk  <- read.csv(file.path(OUT, "02_template_inventory.csv"), stringsAsFactors = FALSE)
mani_chk <- read.csv(file.path(OUT, "00_run_manifest.csv"), stringsAsFactors = FALSE)
check("[0.10] pair table round-trip",
      nrow(pt_chk) == sum(pairs$complete) &&
        identical(names(pt_chk), c("assay", "gene", "f_name", "f_seq", "r_name", "r_seq")),
      paste0("expected ", sum(pairs$complete), " rows, found ", nrow(pt_chk)))
check("[0.10] inventory round-trip",
      nrow(inv_chk) == nrow(taxon) && all(inv_chk$records == inv$records),
      "inventory differs after re-read")
check("[0.10] manifest round-trip",
      any(mani_chk$section == "status" & mani_chk$item == "step0" &
            mani_chk$value == "PASS"),
      "manifest status row missing")

## max_pair_span config row round-trip (the span-gate parameter must
## survive the write/re-read exactly; step 1 [1.1] gates on this same row)
mps <- mani_chk$value[mani_chk$section == "config" & mani_chk$item == "max_pair_span"]
check("[0.10] manifest config|max_pair_span round-trip",
      length(mps) == 1L && identical(mps[1L], as.character(CONFIG$MAX_AMP_SPAN)),
      paste0("expected config|max_pair_span = ", CONFIG$MAX_AMP_SPAN,
             " in the re-read manifest, found: ",
             if (length(mps) == 1L) mps[1L] else "(missing)"))

## consensus_rule section round-trip (one row per template, valid
## vocabulary, values exactly match this session's selection)
cr_keys  <- paste0(taxon$Taxa, "|", taxon$Gene)
cr       <- mani_chk[mani_chk$section == "consensus_rule", , drop = FALSE]
check("[0.10] manifest consensus_rule: one row per template, valid vocabulary",
      nrow(cr) == nrow(taxon) && all(cr_keys %in% cr$item) &&
        all(cr$value %in% CONSENSUS_RULES_VALID),
      sprintf("expected %d rows (one per template) with values in {%s}; found %d row(s)",
              nrow(taxon), paste(CONSENSUS_RULES_VALID, collapse = ", "), nrow(cr)))
cr_match <- TRUE
for (i in seq_len(nrow(taxon))) {
  v <- cr$value[cr$item == cr_keys[i]]
  if (length(v) != 1L || v != inv$consensus_rule[i]) cr_match <- FALSE
}
check("[0.10] manifest consensus_rule values match this session", cr_match,
      "a per-template rule in the manifest differs from the session selection — re-run step 0")

cat("\n============================================================================\n")
cat(sprintf("STEP 0 PASSED | run_id %s | %d primers | %d pairs | %d genes | %d MSAs | %d species records\n",
            run_id, nrow(prim), sum(pairs$complete), length(genes),
            nrow(taxon), sum(inv$records)))
cat(sprintf("  consensus rules (published in manifest): %d gaps_can_win | %d gaps_never_win\n",
            sum(inv$consensus_rule == "gaps_can_win"),
            sum(inv$consensus_rule == "gaps_never_win")))
cat(sprintf("  span gate (published in manifest): MAX_AMP_SPAN = %d bp\n",
            CONFIG$MAX_AMP_SPAN))
cat("Outputs: results/00_run_manifest.csv | 01_pair_table.csv | 02_template_inventory.csv | run_audit_log.csv\n")
cat("Next: step1_primer_mapping.r (Windows R — same or fresh session)\n")
cat("      step2_primerminer.r runs on Ubuntu/WSL R (PrimerMiner dependency)\n")
cat("============================================================================\n")
