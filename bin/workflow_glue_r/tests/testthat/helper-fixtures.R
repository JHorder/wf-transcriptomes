# Helpers

#' Get workflow_glue_r repository root directory from environment variable.
#' @return Normalized path to workflow_glue_r repository root directory.
#' @export
workflow_glue_r_repo_root <- function() {
    normalizePath(Sys.getenv("WORKFLOW_GLUE_R_REPO_ROOT"))
}


#' Get path to a test fixture file within the workflow_glue_r test data directory.
#' @param ... Path components within the test data directory.
#' @return Normalized path to the specified test fixture file.
#' @export
workflow_glue_r_fixture <- function(...) {
    file.path(Sys.getenv("TEST_DATA"), "workflow_glue_r", ...)
}

#' Run a system command and capture its output and status.
#' @param command Command to run.
#' @param args Arguments to the command.
#' @param stdout Whether to capture standard output (TRUE or file path).
#' @param stderr Whether to capture standard error (TRUE or file path).
#' @return A list with 'status' (exit status) and 'output' (captured output).
#' @export
run_system2 <- function(command, args, stdout = TRUE, stderr = TRUE) {
    output <- system2(command, args, stdout = stdout, stderr = stderr)
    if (is.character(output)) {
        status <- attr(output, "status")
        if (is.null(status)) {
            status <- 0L
        }
        return(list(status = status, output = output))
    }

    captured <- character(0)
    if (is.character(stderr) && file.exists(stderr)) {
        captured <- readLines(stderr, warn = FALSE)
    }

    list(status = as.integer(output), output = captured)
}
#' Run an R script using Rscript.
#' @param script_name Name of the R script to run.
#' @param args Arguments to pass to the R script.
#' @return A list with 'status' (exit status) and 'output' (captured output).
#' @export
run_rscript <- function(script_name, args) {
    command <- file.path(R.home("bin"), "Rscript")
    script <- file.path(workflow_glue_r_repo_root(), "bin", script_name)
    run_system2(command, c(script, args), stdout = TRUE, stderr = TRUE)
}

#' Create a test transcript-level SummarizedExperiment with synthetic count data.
#'
#' SummarizedExperiment is a Bioconductor container that holds:
#' - assays: matrices of counts/CPM (rows=transcripts, cols=samples)
#' - rowRanges: genomic coordinates and metadata for each transcript
#' - colData: sample metadata (not used here)
#'
#' This fixture creates realistic differential expression patterns:
#' - tx1: high in control, low in treated (downregulated)
#' - tx2: low in control, high in treated (upregulated)
#' - tx3: similar across conditions (not DE)
#' - tx4: low counts in both (filtered in real DE analysis)
#'
#' @param include_geneid Whether to include GENEID metadata column in rowRanges.
#' @param sample_names Optional vector of sample names to use as column names.
#' @return A SummarizedExperiment object with synthetic transcript-level data.
#' @export
make_test_tx_se <- function(include_geneid = TRUE, sample_names = NULL) {
    if (is.null(sample_names)) {
        sample_names <- c(
            "control_rep1",
            "control_rep2",
            "control_rep3",
            "treated_rep1",
            "treated_rep2",
            "treated_rep3"
        )
    }

    counts <- vapply(sample_names, function(sample_name) {
        condition <- sub("_rep[0-9]+$", "", sample_name)
        replicate_id <- suppressWarnings(as.integer(sub("^.*_rep", "", sample_name)))
        if (is.na(replicate_id)) {
            replicate_id <- match(sample_name, sample_names)
        }
        offset <- ((replicate_id - 1L) %% 3L) - 1L

        base_counts <- switch(
            condition,
            control = c(120, 18, 80, 15),
            treated = c(25, 95, 76, 14),
            treated2 = c(35, 85, 60, 40),
            baseline = c(120, 18, 80, 15),
            c(80, 40, 70, 20)
        )
        pmax(base_counts + c(offset, -offset, offset, 0), 1)
    }, numeric(4))
    dimnames(counts) <- list(c("tx1", "tx2", "tx3", "tx4"), sample_names)
    cpm <- t(t(counts) / colSums(counts)) * 1e6

    row_ranges <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(1, 101, 201, 301), width = 50),
        strand = rep("+", 4),
        TXNAME = rownames(counts),
        eqClassById = IRanges::CharacterList(list(c("1", "2"), "3", "4", "5"))
    )
    if (include_geneid) {
        S4Vectors::mcols(row_ranges)$GENEID <- c("gene1", "gene1", "gene2", "gene2")
    }

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts, CPM = cpm),
        rowRanges = row_ranges
    )
}

#' Create bambu-style transcript row ranges for output-writing tests.
#'
#' `bambu::writeToGTF()` expects a `GRangesList` shaped like the object
#' returned by `bambu::prepareAnnotations()`. This helper writes a minimal GTF
#' and returns that object with metadata columns used by test assertions.
#'
#' @param out_dir Directory where temporary annotation fixture is written.
#' @return A GRangesList with tx1-tx4 transcript entries and metadata.
#' @export
make_test_bambu_row_ranges <- function(out_dir) {
    gtf <- file.path(out_dir, "annotation.gtf")
    writeLines(
        c(
            paste(
                "chr1", "test", "transcript", "1", "50", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "1", "50", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx1"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "101", "150", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx2";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "101", "150", ".", "+", ".",
                'gene_id "gene1"; transcript_id "tx2"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "201", "250", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx3";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "201", "250", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx3"; exon_number "1";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "transcript", "301", "350", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx4";',
                sep = "\t"
            ),
            paste(
                "chr1", "test", "exon", "301", "350", ".", "+", ".",
                'gene_id "gene2"; transcript_id "tx4"; exon_number "1";',
                sep = "\t"
            )
        ),
        gtf
    )

    row_ranges <- bambu::prepareAnnotations(gtf)[c("tx1", "tx2", "tx3", "tx4")]
    S4Vectors::mcols(row_ranges)$TXNAME <- names(row_ranges)
    S4Vectors::mcols(row_ranges)$GENEID <- c("gene1", "gene1", "gene2", "gene2")
    S4Vectors::mcols(row_ranges)$eqClassById <- IRanges::CharacterList(list(c("1", "2"), "3", "4", "5"))

    row_ranges
}

#' Create a test gene-level SummarizedExperiment by aggregating a transcript-level SummarizedExperiment.
#' @param sample_names Optional vector of sample names to use as column names.
#' @return A SummarizedExperiment object with synthetic gene-level data.
#' @export
make_test_gene_se <- function(sample_names = NULL) {
    if (is.null(sample_names)) {
        sample_names <- c(
            "control_rep1",
            "control_rep2",
            "control_rep3",
            "treated_rep1",
            "treated_rep2",
            "treated_rep3"
        )
    }

    tx_se <- make_test_tx_se(sample_names = sample_names)
    tx_counts <- SummarizedExperiment::assays(tx_se)$counts
    counts <- rbind(
        gene1 = tx_counts["tx1", ] + tx_counts["tx2", ],
        gene2 = tx_counts["tx3", ] + tx_counts["tx4", ]
    )
    cpm <- t(t(counts) / colSums(counts)) * 1e6

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts, CPM = cpm),
        rowData = S4Vectors::DataFrame(GENEID = rownames(counts))
    )
}

#' Create a test DE analysis input bundle with synthetic transcript and gene
#' SummarizedExperiments and a sample sheet.
#' @param out_dir Directory to write the output files.
#' @param levels Optional vector of condition levels to use in the sample sheet.
#' @return A list with paths to the generated transcript RDS, gene RDS, and sample sheet CSV.
#' @export
write_de_fixture_bundle <- function(out_dir, levels = c("control", "treated")) {
    sample_names <- unlist(lapply(levels, function(level) {
        sprintf("%s_rep%d", level, seq_len(3))
    }))

    tx_se <- make_test_tx_se(sample_names = sample_names)
    gene_se <- make_test_gene_se(sample_names = sample_names)

    tx_path <- file.path(out_dir, "transcripts.rds")
    gene_path <- file.path(out_dir, "genes.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(out_dir, "sample_sheet.csv")
    sample_rows <- do.call(
        rbind,
        lapply(levels, function(level) {
            data.frame(
                alias = sprintf("%s_rep%d", level, seq_len(3)),
                condition = level,
                batch = c("b1", "b2", "b1"),
                stringsAsFactors = FALSE
            )
        })
    )
    utils::write.csv(sample_rows, sample_sheet, row.names = FALSE, quote = FALSE)

    list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet
    )
}

# Expecters

#' Expect a system command to succeed.
#' @param command Command to run.
#' @param args Arguments to the command.
#' @param stdout Whether to capture standard output (TRUE or file path).
#' @param stderr Whether to capture standard error (TRUE or file path).
#' @return Invisibly returns the result of the command.
#' @export
expect_command_success <- function(command, args, stdout = TRUE, stderr = TRUE) {
    result <- run_system2(command, args, stdout = stdout, stderr = stderr)
    testthat::expect_equal(
        result$status,
        0L,
        info = paste(c(command, args, result$output), collapse = "\n")
    )
    invisible(result)
}

#' Expect a BAM fixture can be built with minimap2 and samtools.
#' @param reference Path to the reference FASTA file.
#' @param reads Path to the FASTQ file containing reads to align.
#' @param out_dir Directory to write the output BAM file and intermediate files.
#' @param alias Sample alias to use in output file names.
#' @return Path to the generated BAM file.
#' @export
expect_bam_fixture_built <- function(reference, reads, out_dir, alias = "sampleA") {
    sam_path <- file.path(out_dir, sprintf("%s.sam", alias))
    bam_path <- file.path(out_dir, sprintf("%s.aligned.sorted.bam", alias))
    minimap2_stderr <- file.path(out_dir, sprintf("%s.minimap2.stderr.txt", alias))
    sort_stderr <- file.path(out_dir, sprintf("%s.samtools-sort.stderr.txt", alias))
    index_stderr <- file.path(out_dir, sprintf("%s.samtools-index.stderr.txt", alias))

    expect_command_success(
        "minimap2",
        c("-ax", "splice", "-uf", reference, reads),
        stdout = sam_path,
        stderr = minimap2_stderr
    )
    testthat::expect_true(file.exists(sam_path))
    testthat::expect_gt(file.info(sam_path)$size, 0)

    expect_command_success(
        "samtools",
        c("sort", "-o", bam_path, sam_path),
        stderr = sort_stderr
    )
    expect_command_success(
        "samtools",
        c("index", bam_path),
        stderr = index_stderr
    )

    bam_path
}

#' Create a degenerate DE fixture where both condition groups have identical
#' zero-heavy counts, simulating the control-vs-control case from the bug reports.
#'
#' Properties:
#' - Every transcript has at least one zero (forces DEXSeq poscounts path)
#' - Every gene has at least one zero (forces DESeq2 poscounts path via
#'   de_choose_size_factor_type)
#' - Both condition groups have identical counts (no differential signal,
#'   stresses dispersion estimation)
#'
#' @param out_dir Directory to write the output files.
#' @return A list with paths to transcript RDS, gene RDS, and sample sheet CSV.
#' @export
write_degenerate_de_fixture_bundle <- function(out_dir) {
    sample_names <- c(
        "ctrl_rep1", "ctrl_rep2", "ctrl_rep3",
        "case_rep1", "case_rep2", "case_rep3"
    )

    # Identical counts across both groups; zeros in every transcript.
    tx_counts <- matrix(
        c(
            0, 0, 4, 0, 0, 4,  # tx1: zero in rep1 and rep2 of both groups
            3, 4, 3, 3, 4, 3,  # tx2: all nonzero but no group difference
            2, 2, 0, 2, 2, 0,  # tx3: zero in rep3
            0, 1, 1, 0, 1, 1   # tx4: zero in rep1
        ),
        nrow = 4, byrow = TRUE,
        dimnames = list(c("tx1", "tx2", "tx3", "tx4"), sample_names)
    )
    tx_cpm <- t(t(tx_counts) / pmax(colSums(tx_counts), 1)) * 1e6

    row_ranges <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 4),
        ranges = IRanges::IRanges(start = c(1, 101, 201, 301), width = 50),
        strand = rep("+", 4),
        TXNAME = c("tx1", "tx2", "tx3", "tx4"),
        GENEID = c("gene1", "gene1", "gene2", "gene2"),
        eqClassById = IRanges::CharacterList(list(c("1", "2"), "3", "4", "5"))
    )

    tx_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = tx_counts, CPM = tx_cpm),
        rowRanges = row_ranges
    )

    # Gene counts: gene1 = tx1+tx2, gene2 = tx3+tx4.
    # Each gene has at least one zero per group (triggers poscounts for DGE).
    gene_counts <- rbind(
        gene1 = tx_counts["tx1", ] + tx_counts["tx2", ],
        gene2 = tx_counts["tx3", ] + tx_counts["tx4", ]
    )
    gene_cpm <- t(t(gene_counts) / pmax(colSums(gene_counts), 1)) * 1e6

    gene_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = gene_counts, CPM = gene_cpm),
        rowData = S4Vectors::DataFrame(GENEID = rownames(gene_counts))
    )

    tx_path <- file.path(out_dir, "transcripts.rds")
    gene_path <- file.path(out_dir, "genes.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(out_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = sample_names,
            condition = rep(c("control", "case"), each = 3),
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet
    )
}

#' Per-condition replicate counts for genes designed around the DGE
#' prefilter boundaries with --de_min_count 10 --de_min_samples 3.
#'
#' Expected outcome in a control vs treated contrast:
#' - retained: target_only, reference_only, exact_boundary, rounds_up,
#'   rounds_half_even, mixed_conditions
#' - excluded: singleton, below_count, below_samples, rounds_down,
#'   third_only, third_tips_balance, all_zero
#' In a control vs treated2 contrast, third_only and third_tips_balance are
#' retained because their support comes from treated2 samples.
gene_filter_designed_counts <- function() {
    zeros <- c(0, 0, 0)
    list(
        # Mirrors the motivating case: nonzero in one sample, zero in five.
        singleton = list(control = zeros, treated = c(0, 50, 0), treated2 = zeros),
        target_only = list(control = zeros, treated = c(25, 30, 20), treated2 = zeros),
        reference_only = list(control = c(15, 12, 18), treated = zeros, treated2 = zeros),
        exact_boundary = list(control = c(10, 10, 10), treated = c(9, 9, 9), treated2 = c(9, 9, 9)),
        below_count = list(control = c(9, 9, 9), treated = c(9, 9, 9), treated2 = c(9, 9, 9)),
        below_samples = list(control = c(40, 40, 9), treated = c(9, 0, 0), treated2 = zeros),
        # Unrounded 9.6 is below 10; DESeq2 sees round(9.6) = 10.
        rounds_up = list(control = c(9.6, 9.6, 9.6), treated = zeros, treated2 = zeros),
        # R rounds half to even, so 9.5 becomes 10, matching DESeq2's input.
        rounds_half_even = list(control = zeros, treated = c(9.5, 9.5, 9.5), treated2 = zeros),
        rounds_down = list(control = c(9.4, 9.4, 9.4), treated = c(9.4, 9.4, 9.4), treated2 = c(9.4, 9.4, 9.4)),
        mixed_conditions = list(control = c(12, 0, 0), treated = c(11, 0, 13), treated2 = zeros),
        third_only = list(control = zeros, treated = zeros, treated2 = c(30, 35, 40)),
        # Only two qualifying samples unless treated2 samples leak into the
        # control vs treated contrast.
        third_tips_balance = list(control = c(10, 10, 0), treated = zeros, treated2 = c(10, 10, 10)),
        all_zero = list(control = zeros, treated = zeros, treated2 = zeros)
    )
}

#' Build gene and transcript counts for DGE prefilter tests.
#'
#' Thirty deterministic, well-expressed background genes (five of them
#' 4-fold up in `treated`) give DESeq2 enough data to fit, alongside the
#' designed boundary genes from `gene_filter_designed_counts()`. Each gene is
#' split into two transcripts so DEXSeq can run on the same samples.
#'
#' @param levels Condition levels, three replicates each.
#' @return A list with `gene_counts`, `tx_counts`, `tx_meta`,
#'   `sample_names` and `conditions`.
make_gene_filter_counts <- function(levels = c("control", "treated")) {
    conditions <- rep(levels, each = 3)
    sample_names <- sprintf("%s_rep%d", conditions, rep(seq_len(3), length(levels)))
    size_factors <- rep(c(1, 1.15, 0.85), length(levels)) *
        rep(seq(0.9, 1.1, length.out = length(levels)), each = 3)

    n_background <- 30
    background <- t(vapply(seq_len(n_background), function(g) {
        base_mean <- 50 * 1.13^g
        fold_change <- ifelse(g <= 5 & conditions == "treated", 4, 1)
        noise <- 1 + 0.2 * sin(g * 1.7 + seq_along(conditions) * 2.3)
        round(base_mean * size_factors * fold_change * noise)
    }, numeric(length(conditions))))
    rownames(background) <- sprintf("bg%02d", seq_len(n_background))

    designed_spec <- gene_filter_designed_counts()
    designed <- t(vapply(designed_spec, function(spec) {
        unlist(lapply(levels, function(level) spec[[level]]), use.names = FALSE)
    }, numeric(length(conditions))))
    rownames(designed) <- names(designed_spec)

    gene_counts <- rbind(background, designed)
    colnames(gene_counts) <- sample_names

    # Two transcripts per gene; background genes 1-4 switch isoform usage in
    # the second condition level to give DEXSeq some signal.
    tx_prop <- t(vapply(seq_len(nrow(gene_counts)), function(g) {
        base_prop <- 0.3 + 0.2 * (g %% 3)
        ifelse(g <= 4 & conditions == levels[[2]], 1 - base_prop, base_prop)
    }, numeric(length(conditions))))
    tx_counts <- rbind(gene_counts * tx_prop, gene_counts * (1 - tx_prop))
    gene_ids <- rep(rownames(gene_counts), 2)
    tx_names <- paste0(gene_ids, rep(c("_tx1", "_tx2"), each = nrow(gene_counts)))
    tx_order <- order(gene_ids, tx_names)
    tx_counts <- tx_counts[tx_order, , drop = FALSE]
    dimnames(tx_counts) <- list(tx_names[tx_order], sample_names)

    list(
        gene_counts = gene_counts,
        tx_counts = tx_counts,
        tx_meta = data.frame(
            TXNAME = tx_names[tx_order],
            GENEID = gene_ids[tx_order],
            stringsAsFactors = FALSE
        ),
        sample_names = sample_names,
        conditions = conditions
    )
}

#' Write a DE input bundle for DGE prefilter tests.
#'
#' @param out_dir Directory to write the output files.
#' @param levels Condition levels, three replicates each.
#' @return A list with paths to the transcript RDS, gene RDS, and sample sheet CSV.
write_gene_filter_fixture_bundle <- function(out_dir, levels = c("control", "treated")) {
    fixture <- make_gene_filter_counts(levels)

    tx_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = fixture$tx_counts),
        rowData = S4Vectors::DataFrame(fixture$tx_meta)
    )
    gene_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = fixture$gene_counts),
        rowData = S4Vectors::DataFrame(
            GENEID = rownames(fixture$gene_counts),
            gene_name = toupper(rownames(fixture$gene_counts))
        )
    )

    tx_path <- file.path(out_dir, "transcripts.rds")
    gene_path <- file.path(out_dir, "genes.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(out_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = fixture$sample_names,
            condition = fixture$conditions,
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet
    )
}
