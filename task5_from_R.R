# Run in R on the Ubuntu VM, not in Windows R.
# This R controller launches the course-derived Nextflow pipeline.
# The Nextflow pipeline and this controller must both accompany the submission.
setwd(path.expand('~/workshop'))
nf <- 'nextflow/workflow_task5.nf'
stopifnot(file.exists(nf))

# Fix deferred evaluation of the Bambu task-specific publication directory.
code <- readLines(nf, warn = FALSE)
old <- 'publishDir "${params.outdir}/bambu_${annotation_mode}", mode: \'copy\''
new <- 'publishDir { "${params.outdir}/bambu_${annotation_mode}" }, mode: \'copy\''
if (any(grepl(old, code, fixed = TRUE))) {
    code <- vapply(code, function(x) {
        if (grepl(old, x, fixed = TRUE)) sub(old, new, x, fixed = TRUE) else x
    }, character(1))
    writeLines(code, nf)
}

state_file <- 'task5_R_state.rds'
if (file.exists(state_file)) {
    task5 <- readRDS(state_file)
} else {
    stamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
    task5 <- list(run_name = paste0('task5_R_', stamp),
                  reports = file.path(getwd(), paste0('task5_reports_R_', stamp)),
                  outputs = file.path(getwd(), paste0('results_task5_R_', stamp)),
                  scenario1_complete = FALSE, scenario2_complete = FALSE)
    saveRDS(task5, state_file)
}
dir.create(task5$reports, recursive = TRUE, showWarnings = FALSE)

run_task5 <- function(scenario = 1L) {
    stopifnot(scenario %in% c(1L, 2L))
    state <- readRDS(state_file)
    if (scenario == 2L && !isTRUE(state$scenario1_complete)) {
        stop('Complete Scenario 1 successfully before starting Scenario 2.')
    }
    label <- paste0('scenario', scenario)
    html <- file.path(state$reports, paste0('report_', label, '.html'))
    trace <- file.path(state$reports, paste0('trace_', label, '.tsv'))
    # Keep previous attempt reports instead of overwriting them.
    for (f in c(html, trace)) {
        if (file.exists(f)) {
            ok <- file.rename(f, paste0(f, '.', format(Sys.time(), '%Y%m%d_%H%M%S'), '.previous'))
            if (!ok) stop('Could not preserve previous report: ', f)
        }
    }
    attempt_marker <- file.path(state$reports, 'scenario1_attempted')
    resume_args <- if (scenario == 2L || file.exists(attempt_marker)) {
        c('-resume', state$run_name)
    } else c('-name', state$run_name)
    args <- c('run', nf, resume_args,
              '--annotation_mode', if (scenario == 1L) 'with' else 'without',
              '--outdir', state$outputs,
              '-with-report', html, '-with-trace', trace, '-ansi-log', 'false')
    command <- paste('nextflow', paste(shQuote(args), collapse = ' '))
    writeLines(command, file.path(state$reports, paste0(label, '.command.txt')))
    cat('\n', command, '\n\n', sep = '')
    if (scenario == 1L) file.create(attempt_marker)
    log <- file.path(state$reports, paste0(label, '.console.log'))
    status <- system2('bash', c('-o', 'pipefail', '-c',
                      shQuote(paste(command, '2>&1 | tee', shQuote(log)))))
    if (status != 0L) stop('Nextflow failed. Keep its error output; do not start the next scenario.')
    state[[paste0(label, '_complete')]] <- TRUE
    saveRDS(state, state_file)
    cat('\nSUCCESS: ', label, '\nReport: ', html, '\n', sep = '')
    invisible(html)
}

summarise_task5 <- function() {
    state <- readRDS(state_file)
    qc <- read.delim(file.path(state$outputs, 'qc', 'qc_summary.tsv'))
    print(qc, row.names = FALSE)
    write.csv(qc, file.path(state$reports, 'qc_summary.csv'), row.names = FALSE)
    pdf(file.path(state$reports, 'qc_mapping.pdf'), width = 11, height = 7)
    par(mar = c(14, 5, 2, 1))
    barplot(qc$mapping_pct, names.arg = qc$sample, las = 2,
            ylim = c(0, 100), ylab = 'Primary reads mapped (%)',
            col = ifelse(qc$mapping_qc == 'PASS', 'steelblue', 'tomato'))
    abline(h = qc$min_mapping_pct[1], lty = 2)
    dev.off()
    second_trace <- file.path(state$reports, 'trace_scenario2.tsv')
    if (file.exists(second_trace)) {
        trace <- read.delim(second_trace, check.names = FALSE)
        print(trace[, intersect(c('name', 'status', 'duration', 'realtime'), names(trace))], row.names = FALSE)
    }
    cat('\nOpen both report_scenario*.html files in: ', state$reports,
        '\nUse their overall elapsed times for the runtime comparison.\n', sep = '')
    invisible(qc)
}

cat('Ready. Run run_task5(1), then after success run_task5(2).\n')
cat('After both succeed, run summarise_task5().\n')
