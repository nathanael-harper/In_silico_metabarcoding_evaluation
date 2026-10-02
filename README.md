# DNA Metabarcoding Assay Evaluation Pipeline

This is supplementary material for: 



Assay Selection, Sequencing Depth, and PCR Stochasticity Affect the Detection of Amphibian Communities using Environmental DNA Metabarcoding 

Nathanael B. J. Harper1, Michael D. J. Lynch1, Cailyn M. Zamora1, Yuwei Xie2, Philip J. Ankley2, John P. Giesy2, Andrew C. Doxey1, Mark R. Servos1, Paul M. Craig1 and Barbara A. Katzenback1* 
1 Department of Biology, University of Waterloo, Waterloo, Ontario, Canada, N2L 3G1 
2 Department of Veterinary Biomedical Sciences, University of Saskatchewan, Saskatoon, Saskatchewan, Canada, S7N 5A2.  

This in silico pipeline evaluates the performance of DNA metabarcoding assays by mapping primers to reference databases, scoring primer-template mismatches, and analyzing barcode resolution (unique vs shared barcode). 
It produces both template-level diagnostics and species-level resolution summaries.

## 1. Prerequisites
*   **R Version:** Latest stable version.
*   **Operating System:** 
    *   **Windows:** Required for Steps 0, 1, 3, 4, 5, and all Supplemental Figures.
    *   **Ubuntu/WSL:** **Required** for Step 2 (`step2_primerminer.r`) due to `PrimerMiner` dependencies.
*   **R Packages:**
    *   **Core:** `ggplot2`, `dplyr`, `tidyr`, `ggnewscale` (for Supplemental Figures).
    *   **Step 2:** `PrimerMiner` (Bioconductor package, installed via WSL).
    *   *Note: The core pipeline logic relies on a local `utils.R` file (not included here) which must be present in the `BASE` directory.*

## 2. Configuration
Before running, you must edit the `BASE` path in each script to point to your working directory.
*   **Windows Path Example:** `BASE <- "C:/Users/Name/Project/Final_run"`
*   **WSL Path Example:** `BASE <- "/mnt/c/Users/Name/Project/Final_run"`

## 3. Input Requirements
Place the following files in your `BASE` directory (or referenced paths):

1.  **`primers_input2.csv`**: Primer sequences (columns: `name`, `seq`, `gene`).
2.  **`taxon_template_directory_2.csv`**: A manifest of reference sequences (columns: `Taxa`, `Gene`, `filepath`).
3.  **MSA Files**: The FASTA files referenced in the taxon directory.
4.  **`GR_spp.csv`**: (Required for Supplemental Figures 3 & 4) A list of Scientific Names used for fuzzy matching species data.

## 4. Execution Workflow
Run the steps in order. Each step creates a `results/` folder with intermediate CSVs.

### Phase 1: Core Analysis (Template Level)
1.  **Step 0 (`step0.R`):** Validates inputs and builds the pair table.
2.  **Step 1 (`step1.R`):** Maps primers to consensus sequences.
3.  **Step 2 (`step2.R`):** **[WSL ONLY]** Scores primers using PrimerMiner.
4.  **Step 3 (`step3.R`):** Extracts the actual amplicon sequences.
5.  **Step 4 (`step4.R`):** Analyzes barcode resolution and collisions.
6.  **Step 5 (`step5.R`):** Prepares data for the Hybrid Plot.

### Phase 2: Visualization (Supplemental Figures)
7.  **Supplemental 1 (`supplemental_figure_1.R`):** Generates the **Hybrid Jitter Plot** (Template-level score distribution).
8.  **Supplemental 2 (`supplemental_figure_2.R`):** Generates the **Stacked Barplot** (Template-level fate distribution).
9.  **Supplemental 3 (`supplemental_figure_3.r`):** Generates the **Species-Averaged Score Plot** (Weighted by haplotype count).
10. **Supplemental 4 (`Supplemental_figure_4.r`):** Generates the **Species-Level Resolution Plot** (Averages fates per species).

## 5. Script Descriptions  Position

### Core Pipeline (Steps 0-5)

#### **Step 0: Input Setup & Audit (`step0.R`)**
*   **Function:** Validates primers, taxon directory, and MSA files. Builds the assay pair table and audits the "consensus rule" (how to handle gaps).

#### **Step 1: Primer Mapping (`step1.R`)**
*   **Function:** Computes population-level consensus sequences. Maps primers against these sequences and applies a "Span Gate" to enforce maximum amplicon size limits.

#### **Step 2: PrimerMiner Scoring (`step2.R`)**
*   **Function:** Slices amplicon windows and scores them using the PrimerMiner algorithm (penalizing mismatches at the 3' end). Classifies every template into five categories: **Suitable, Unsuitable, No Coverage, Incomplete Coverage, or Unmapped**.

#### **Step 3: Amplicon Extraction (`step3.R`)**
*   **Function:** Cuts the specific sequence out of the MSA for every record, de-gaps the sequences, and joins this data with the status from Step 2.

#### **Step 4: Barcode Resolution (`step4.R`)**
*   **Function:** Analyzes extracted amplicons to identify **Collisions** (where two different species have the exact same sequence) and calculates **Resolution** (percentage of species uniquely identifiable).

#### **Step 5: Hybrid-Plot Labels (`step5.R`)**
*   **Function:** Aggregates the "Fate" of templates (No Coverage vs. Good vs. Unmapped) per assay and taxon group, generating formatted label strings for the X-axis.

### Supplemental Figures (Species Level)

#### **Supplemental Figure 1 (`supplemental_figure_1.r`)**
*   **Function:** Plots the PrimerMiner score for every assay-template combination, assuming each template is a unique taxa. It parses sequence IDs to extract haplotype counts, ensuring that common haplotypes weigh more heavily in the average score.

#### **Supplemental Figure 2 (`Supplemental_figure_2.r`)**
*   **Function:** Determines the "Fate" of every species for every assay-template combination (e.g., Unique Barcode, Shared Barcode), assuming each template is a unique taxa. It uses a priority hierarchy (Suitable > Unsuitable > Incomplete) to classify species, then plots the proportion of species falling into each category.

#### **Supplemental Figure 3 (`supplemental_figure_3.r`)**
*   **Function:** Calculates a **weighted average PrimerMiner score** for every assay-species combination, assuming multiple templates may map to a single species. It parses sequence IDs to extract haplotype counts, ensuring that common haplotypes weigh more heavily in the average score.

#### **Supplemental Figure 4 (`Supplemental_figure_4.r`)**
*   **Function:** Determines the "Fate" of every species for every assay-species combination (e.g., Unique Barcode, Shared Barcode), assuming multiple templates may map to a single species. It uses a priority hierarchy (Suitable > Unsuitable > Incomplete) to classify species, then plots the proportion of species falling into each category.

## 6. Key Outputs
All outputs are saved to `results/`.

| File ID | Description |
| :--- | :--- |
| `00_run_manifest.csv` | Audit log of the run (inputs, config, and versions). |
| `03_primer_mapping.csv` | Coordinates where primers bind. |
| `06_assay_template_status.csv` | The "Master Table" classifying every template (Suitable/Unsuitable/etc.). |
| `07_amplicons.csv` | The extracted DNA sequences. |
| `08_barcode_resolution_summary.csv` | Summary of how many species are uniquely identifiable. |
| `12_stacked_plot_data.csv` | Data for the stacked barplot figure. |
| `14_assay_axis_label_strings.csv` | Formatted labels for the hybrid plot. |
| `SuppFig3_primer_species_summary.csv` | Intermediate data for Species-Averaged Score Plot. |
| `species_fate_audit.csv` | Intermediate data for Species-Level Resolution Plot. |
| `figures/SuppFig3_...png` | The Species-Averaged Score Plot. |
| `figures/SuppFig4_...png` | The Species-Level Resolution Plot. |
| `figures/Primer_Full_Distribution...png` | The Hybrid Jitter Plot. |

## 7. Troubleshooting
*   **Run Stops:** Check the console output for the `[CHECK]` or `[1.1]` error message. The pipeline is strict; if a file changes or config drifts, it will stop to ensure data integrity.
*   **PrimerMiner Error:** Ensure you are running Step 2 in a WSL/Ubuntu session and that the package is installed.
*   **Missing Results:** Ensure `step0` has passed before attempting `step1`. The manifest created in Step 0 is required for all subsequent steps.
*   **Fuzzy Matching Errors:** Ensure `GR_spp.csv` is accessible for Supplemental Figures 3 and 4.
