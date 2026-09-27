# DSA4262-Ranier-Goh-Workflow-Code

Code for the long-read RNA-seq workflow in Task 5 of DSA4262 Assignment 1.

The workflow is adapted from the [Introduction to Genomics 3 Nextflow pipeline](https://github.com/GoekeLab/sg-nex-data/blob/master/docs/colab/workflow_longReadRNASeq.nf). An R script launches Nextflow and summarises QC results; Nextflow manages the analysis steps.

## Files

- `task5_from_R.R` — R controller for running both scenarios and summarising QC.
- `workflow_task5.nf` — Nextflow pipeline for alignment, BAM conversion, QC, and transcript discovery and quantification.

## Workflow

1. Align reads to GRCh38 using Minimap2, with different parameters for direct RNA and cDNA reads.
2. Convert alignments to sorted, indexed BAM files using Samtools.
3. Check primary read counts and mapping percentages.
4. Analyse all four samples jointly using Bambu, producing transcript counts, gene counts, and a GTF file.

Two scenarios are supported:

- **Scenario 1:** Bambu with GRCh38.91 reference annotations.
- **Scenario 2:** Bambu without reference annotations, using Nextflow’s `-resume` option to reuse unchanged tasks.

## Requirements

The analysis environment used:

- Nextflow 26.04.6
- Minimap2 2.26-r1175
- Samtools 1.19.2
- R with Bambu 3.13.1
- Bash on Ubuntu

All command-line tools must be available on `PATH`.

## Input files

The scripts expect the following directory structure on the analysis machine:

```text
~/workshop/
├── task5_from_R.R
├── nextflow/
│   └── workflow_task5.nf
├── fastq/
│   ├── SGNex_A549_directRNA_replicate1_run1.fastq.gz
│   ├── SGNex_Hct116_cDNA_replicate1_run6.fastq.gz
│   ├── SGNex_Hct116_directRNA_replicate6_run1.fastq.gz
│   └── SGNex_K562_cDNA_replicate1_run3.fastq.gz
└── reference/
    ├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
    ├── Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai
    └── Homo_sapiens.GRCh38.91.gtf
```

If both scripts are stored at the repository root, place the `.nf` file in `~/workshop/nextflow/` and the `.R` file in `~/workshop/` before running. Input datasets must be supplied separately.

## Running the workflow

Start R on the Ubuntu analysis machine:

```r
setwd("~/workshop")
source("task5_from_R.R")

# Scenario 1: with annotations
run_task5(1)

# Run only after Scenario 1 succeeds
run_task5(2)

# Summarise QC after both runs succeed
summarise_task5()
```

The R controller prints the Nextflow commands and saves execution reports and logs. Keep the `work/`, `.nextflow/`, and `task5_R_state.rds` files intact so completed tasks can be reused.

## Outputs

Results are saved in timestamped directories:

```text
results_task5_R_<timestamp>/
├── bam/
├── qc/
├── bambu_with/
└── bambu_without/

task5_reports_R_<timestamp>/
```

Outputs include:

- Sorted BAM files and indexes
- Per-sample QC results and a combined QC table
- Transcript and gene read counts
- Bambu output annotations in GTF format
- HTML execution reports, task traces, and console logs
- An R-generated QC summary CSV and mapping-percentage PDF

## QC criteria and limitations

Samples are flagged if they have fewer than **10,000 primary reads** or a **primary mapping percentage below 70%**. Secondary and supplementary alignments are excluded from these calculations. These are exploratory thresholds, and flagged samples remain included in Bambu analysis.

The workflow retains the workshop’s permissive `NDR=1` setting. Bambu recommended 0.128 during the annotated analysis, so discovery thresholds should be reviewed before interpreting novel transcripts.

Full-dataset analysis can take considerable time, particularly in Bambu. Resumed-run elapsed times exclude earlier computation and should not be presented as uninterrupted end-to-end runtimes.

## Acknowledgements

The original Nextflow workflow was developed by Jonathan Göke and provided through the GoekeLab genomics workshop. Adaptations include four-sample processing, protocol-specific alignment, BAM sorting and indexing, QC, an annotation-free mode, and an R controller.

Sequencing data are from the [Singapore Nanopore Expression Project (SG-NEx)](https://github.com/GoekeLab/sg-nex-data).
