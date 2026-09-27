#!/usr/bin/env nextflow
// Adapted from Jonathan Goke's GoekeLab course workflow:
// https://github.com/GoekeLab/sg-nex-data/blob/master/docs/colab/workflow_longReadRNASeq.nf
// ADDITIONS: four explicit samples, protocol-aware alignment, sorting/indexing,
// mapping QC, joint Bambu, annotation switch, resource limits, separate outputs.
nextflow.enable.dsl=2

params.fastq_dir = "${launchDir}/fastq"
params.refFa = "${launchDir}/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.refGtf = "${launchDir}/reference/Homo_sapiens.GRCh38.91.gtf"
params.outdir = "${launchDir}/results_task5"
params.annotation_mode = 'with'
// Preserve the course NDR setting in BOTH scenarios; this is permissive.
params.ndr = 1
// Exploratory assignment QC criteria, NOT validated universal cutoffs.
params.min_primary_reads = 10000
params.min_mapping_pct = 70

process MINIMAP2_ALIGN {
    tag sample
    cpus 8
    memory '32 GB'
    maxForks 2
    input:
    tuple val(sample), val(protocol), path(reads)
    path refFa
    output:
    tuple val(sample), path("${sample}.sam")
    script:
    def extra = protocol == 'directRNA' ? '-uf -k14' : ''
    """
    minimap2 -t ${task.cpus} -ax splice ${extra} "${refFa}" "${reads}" > ${sample}.sam
    """
}

process SAM_TO_BAM {
    tag sample
    cpus 4
    memory '8 GB'
    maxForks 2
    publishDir "${params.outdir}/bam", mode: 'copy'
    input:
    tuple val(sample), path(reads_sam)
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
    """
    set -euo pipefail
    samtools view -u "${reads_sam}" | samtools sort -@ 2 -m 1G -o ${sample}.bam -
    samtools index ${sample}.bam
    samtools quickcheck -v ${sample}.bam
    """
}

process QC {
    tag sample
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
    tuple val(sample), path(bam), path(bai)
    output:
    path "${sample}.qc.tsv", emit: summary
    path "${sample}.flagstat.txt", emit: flagstat
    script:
    """
    set -euo pipefail
    samtools quickcheck -v "${bam}"
    samtools flagstat "${bam}" > ${sample}.flagstat.txt
    # Exclude secondary (256) and supplementary (2048) records.
    total=\$(samtools view -c -F 2304 "${bam}")
    # Also exclude unmapped (4) records from the numerator.
    mapped=\$(samtools view -c -F 2308 "${bam}")
    awk -v sample="${sample}" -v total="\$total" -v mapped="\$mapped" \\
        -v minreads=${params.min_primary_reads} -v minpct=${params.min_mapping_pct} '
        BEGIN {
            pct = total > 0 ? 100 * mapped / total : 0;
            depth = total >= minreads ? "PASS" : "FAIL";
            mapping = total > 0 && pct >= minpct ? "PASS" : "FAIL";
            overall = depth == "PASS" && mapping == "PASS" ? "PASS" : "FAIL";
            printf "sample\\tprimary_reads\\tprimary_mapped\\tmapping_pct\\tdepth_qc\\tmapping_qc\\toverall\\tmin_primary_reads\\tmin_mapping_pct\\n";
            printf "%s\\t%d\\t%d\\t%.3f\\t%s\\t%s\\t%s\\t%d\\t%.3f\\n", sample, total, mapped, pct, depth, mapping, overall, minreads, minpct;
        }' > ${sample}.qc.tsv
    """
}

process QC_SUMMARY {
    cpus 1
    memory '1 GB'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
    path summaries
    output:
    path 'qc_summary.tsv'
    script:
    """
    awk 'FNR == 1 && NR != 1 {next} {print}' ${summaries} > qc_summary.tsv
    """
}

process BAMBU {
    tag annotation_mode
    cpus 4
    memory '96 GB'
    publishDir "${params.outdir}/bambu_${annotation_mode}", mode: 'copy'
    input:
    path bams
    path indexes
    path refFa
    path refFai
    path annotation_files
    val annotation_mode
    output:
    path 'counts_transcript.txt'
    path 'counts_gene.txt'
    path 'extended_annotations.gtf'
    path 'run_settings.txt'
    path 'sessionInfo.txt'
    script:
    def annotationCode = annotation_mode == 'with' ?
        "annotations <- prepareAnnotations(\"${annotation_files}\")" :
        'annotations <- NULL'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(bambu))
    ${annotationCode}
    bam_files <- sort(list.files(pattern = "[.]bam\$", full.names = TRUE))
    stopifnot(length(bam_files) == 4L)
    writeLines(c("annotation_mode=${annotation_mode}", "NDR=${params.ndr}",
                 "ncore=${task.cpus}", bam_files), "run_settings.txt")
    se <- bambu(reads = bam_files, annotations = annotations,
                genome = "${refFa}", NDR = ${params.ndr},
                discovery = TRUE, quant = TRUE, ncore = ${task.cpus})
    writeBambuOutput(se, path = "./")
    capture.output(sessionInfo(), file = "sessionInfo.txt")
    """
}

workflow {
    if (!(params.annotation_mode in ['with', 'without'])) {
        error "--annotation_mode must be with or without"
    }
    def ref_fa = file(params.refFa, checkIfExists: true)
    def ref_fai = file("${params.refFa}.fai", checkIfExists: true)
    def annotation = params.annotation_mode == 'with' ?
        file(params.refGtf, checkIfExists: true) : []
    def samples = [
        ['SGNex_A549_directRNA_replicate1_run1', 'directRNA'],
        ['SGNex_Hct116_cDNA_replicate1_run6', 'cDNA'],
        ['SGNex_Hct116_directRNA_replicate6_run1', 'directRNA'],
        ['SGNex_K562_cDNA_replicate1_run3', 'cDNA']
    ]
    reads_ch = Channel.fromList(samples).map { sample, protocol ->
        tuple(sample, protocol, file("${params.fastq_dir}/${sample}.fastq.gz", checkIfExists: true))
    }
    MINIMAP2_ALIGN(reads_ch, ref_fa)
    SAM_TO_BAM(MINIMAP2_ALIGN.out)
    QC(SAM_TO_BAM.out)
    QC_SUMMARY(QC.out.summary.collect(sort: true))
    // Stable ordering prevents completion order from changing Bambu inputs.
    bams_ch = SAM_TO_BAM.out.map { sample, bam, bai -> bam }.collect(sort: true)
    indexes_ch = SAM_TO_BAM.out.map { sample, bam, bai -> bai }.collect(sort: true)
    BAMBU(bams_ch, indexes_ch, ref_fa, ref_fai, annotation, params.annotation_mode)
}
