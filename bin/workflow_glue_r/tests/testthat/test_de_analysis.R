###
# Design parsing and validation
#
# Tests covering the experiment-design checks owned by supeRglue de_analysis:
# covariate parsing, required sample-sheet columns, alias matching, reference
# level selection, and the transcript metadata needed for DEXSeq grouping.

# Covariates CLI arg is comma-separated string with possible whitespace/empty values.
testthat::test_that("covariates parsed and trimmed", {
    testthat::expect_equal(
        workflow_glue_r_parse_csv_list(" batch, sex ,, site "),
        c("batch", "sex", "site")
    )
})


testthat::test_that("DTU transcript output renames contrast-specific log2fold column", {
    contrast_name <- "treated_control"
    dex_df <- data.frame(
        featureID = "tx1",
        groupID = "gene1",
        log2fold_treated_control = 1.5,
        pvalue = 0.01,
        padj = 0.05,
        exonBaseMean = 100,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    dtu_tx <- de_extract_dtu_transcript_table(dex_df, contrast_name)

    testthat::expect_equal(
        names(dtu_tx),
        c(
            "featureID",
            "groupID",
            "log2FoldChange",
            "pvalue",
            "padj",
            "exonBaseMean"
        )
    )
})

# Sample aliases CAN have hyphens (they're not used in formulas, just for matching).
# Sample sheet rows can be in different order than SE columns - should reorder automatically.
testthat::test_that("non-syntactic aliases allowed, sheets reordered", {
    sample_names <- c(
        "1control",
        "control-2",
        "control _rep 3",
        "04",
        "treated-5",
        "treated-rep.6"
    )
    tx_se <- make_test_tx_se(sample_names = sample_names)
    gene_se <- make_test_gene_se(sample_names = sample_names)
    sample_df <- data.frame(
        barcode = sprintf("barcode%02d", seq_along(sample_names)),
        sample = paste0("sample_", seq_along(sample_names)),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        alias = sample_names,
        stringsAsFactors = FALSE
    )
    sample_df <- sample_df[c(6, 3, 5, 1, 4, 2), , drop = FALSE]
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control"
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    validated <- de_validate_inputs(tx_se, gene_se, sample_df, argv)

    testthat::expect_equal(validated$sample_df$alias, sample_names)
    testthat::expect_equal(
        as.character(validated$sample_df$condition),
        rep(c("control", "treated"), each = 3)
    )
})

# If --reference_level not specified, default to "control" if it exists in condition column.
# Sample aliases must exactly match SE column names (no missing, no extra).
testthat::test_that("reference level defaults to control", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = NULL
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )

    validated <- de_validate_inputs(tx_se, gene_se, sample_df, argv)
    testthat::expect_equal(validated$reference_level, "control")

    bad_alias_df <- sample_df
    bad_alias_df$alias[1] <- "missing_alias"
    testthat::expect_error(
        de_validate_inputs(tx_se, gene_se, bad_alias_df, argv),
        "Sample sheet aliases do not match"
    )
})


# Count matrices must be numeric, non-NA, non-zero totals per sample.
# Transcript and gene SE sample names must match exactly (same order, same values).
testthat::test_that("unusable count data rejected", {
    tx_se <- make_test_tx_se()
    gene_se <- make_test_gene_se()
    argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv",
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control"
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())
    sample_df <- data.frame(
        alias = colnames(tx_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )

    duplicated_tx_se <- tx_se
    colnames(duplicated_tx_se)[2] <- colnames(duplicated_tx_se)[1]
    testthat::expect_error(
        de_validate_inputs(duplicated_tx_se, gene_se, sample_df, argv),
        "Transcript RDS sample names must be unique"
    )

    mismatched_gene_se <- gene_se
    colnames(mismatched_gene_se)[1] <- "missing_from_transcripts"
    testthat::expect_error(
        de_validate_inputs(tx_se, mismatched_gene_se, sample_df, argv),
        "Transcript and gene RDS sample names must match"
    )

    zero_count_tx_se <- tx_se
    SummarizedExperiment::assay(zero_count_tx_se, "counts")[, "control_rep1"] <- 0
    testthat::expect_error(
        de_validate_inputs(zero_count_tx_se, gene_se, sample_df, argv),
        "zero total counts"
    )

    character_count_tx_se <- tx_se
    tx_counts_char <- SummarizedExperiment::assay(character_count_tx_se, "counts")
    storage.mode(tx_counts_char) <- "character"
    SummarizedExperiment::assay(character_count_tx_se, "counts") <- tx_counts_char
    testthat::expect_error(
        de_validate_inputs(character_count_tx_se, gene_se, sample_df, argv),
        "Count matrices must be numeric"
    )

    na_count_gene_se <- gene_se
    SummarizedExperiment::assay(na_count_gene_se, "counts")[1, 1] <- NA
    testthat::expect_error(
        de_validate_inputs(tx_se, na_count_gene_se, sample_df, argv),
        "must not contain NA values"
    )
})

# DRIMSeq (DTU analysis) requires GENEID metadata column to group transcripts by gene.
# Fail fast if bambu output is missing this column.
testthat::test_that("transcript SE without GENEID rejected", {
    fixture_dir <- tempfile("de-no-geneid-")
    dir.create(fixture_dir)

    tx_se <- make_test_tx_se(include_geneid = FALSE)
    gene_se <- make_test_gene_se()
    tx_path <- file.path(fixture_dir, "tx.rds")
    gene_path <- file.path(fixture_dir, "gene.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = colnames(tx_se),
            condition = rep(c("control", "treated"), each = 3),
            batch = rep(c("b1", "b2", "b1"), 2),
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    argv <- list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet,
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control",
        out_dir = file.path(fixture_dir, "out")
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    testthat::expect_error(
        main_run_de_analysis(argv),
        "Transcript rowData must contain GENEID"
    )
})

# DE analysis requires ≥2 replicates per condition level and ≥2 condition levels.
# RDS files must be valid R objects (not corrupt/malformed).
testthat::test_that("underspecified designs rejected", {
    fixture_dir <- tempfile("de-bad-data-")
    dir.create(fixture_dir)

    underspecified_samples <- c("control_rep1", "treated_rep1", "treated_rep2")
    tx_se <- make_test_tx_se(sample_names = underspecified_samples)
    gene_se <- make_test_gene_se(sample_names = underspecified_samples)
    tx_path <- file.path(fixture_dir, "tx.rds")
    gene_path <- file.path(fixture_dir, "gene.rds")
    saveRDS(tx_se, tx_path)
    saveRDS(gene_se, gene_path)

    sample_sheet <- file.path(fixture_dir, "sample_sheet.csv")
    utils::write.csv(
        data.frame(
            alias = underspecified_samples,
            condition = c("control", "treated", "treated"),
            batch = c("b1", "b1", "b2"),
            stringsAsFactors = FALSE
        ),
        sample_sheet,
        row.names = FALSE,
        quote = FALSE
    )

    argv <- list(
        transcript_rds = tx_path,
        gene_rds = gene_path,
        sample_sheet = sample_sheet,
        condition_column = "condition",
        covariates = "batch",
        reference_level = "control",
        out_dir = file.path(fixture_dir, "out")
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    testthat::expect_error(
        suppressWarnings(main_run_de_analysis(argv)),
        "fewer than 2 replicates"
    )

    single_condition_sheet <- file.path(fixture_dir, "single_condition.csv")
    utils::write.csv(
        data.frame(
            alias = colnames(make_test_tx_se()),
            condition = rep("control", 6),
            batch = rep(c("b1", "b2", "b1"), 2),
            stringsAsFactors = FALSE
        ),
        single_condition_sheet,
        row.names = FALSE,
        quote = FALSE
    )
    valid_tx_path <- file.path(fixture_dir, "valid_tx.rds")
    valid_gene_path <- file.path(fixture_dir, "valid_gene.rds")
    saveRDS(make_test_tx_se(), valid_tx_path)
    saveRDS(make_test_gene_se(), valid_gene_path)
    argv$transcript_rds <- valid_tx_path
    argv$gene_rds <- valid_gene_path
    argv$sample_sheet <- single_condition_sheet
    argv$out_dir <- file.path(fixture_dir, "single-out")
    testthat::expect_error(
        main_run_de_analysis(argv),
        "requires at least two condition levels"
    )

    malformed_rds <- file.path(fixture_dir, "not-an-rds.rds")
    writeLines("not an RDS file", malformed_rds)
    argv$transcript_rds <- malformed_rds
    argv$sample_sheet <- sample_sheet
    argv$out_dir <- file.path(fixture_dir, "malformed-out")
    testthat::expect_error(
        main_run_de_analysis(argv)
    )
})

###
# Fallback helpers and metadata wiring
#
# Fixture-driven tests for DESeq2/DEXSeq helper behaviour and metadata.
# These avoid dependency injection and exercise the real package code paths.

testthat::test_that("de_choose_size_factor_type uses poscounts for zero-heavy matrices", {
    count_mat <- matrix(
        c(
            0, 4, 5,
            3, 0, 2,
            7, 1, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    testthat::expect_warning(
        sf_type <- de_choose_size_factor_type(count_mat, "test matrix"),
        "using DESeq2 size-factor estimation with sfType='poscounts'"
    )
    testthat::expect_equal(sf_type, "poscounts")
})

testthat::test_that("de_run_deseq_with_fallback returns structured metadata", {
    testthat::skip_if_not_installed("DESeq2")

    gene_se <- make_test_gene_se()
    sample_df <- data.frame(
        alias = colnames(gene_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    rownames(sample_df) <- sample_df$alias

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(gene_se, "counts"),
        colData = sample_df,
        design = ~ batch + condition
    )
    out_dir <- tempfile("deseq-fallback-")
    dir.create(out_dir)

    result <- suppressWarnings(de_run_deseq_with_fallback(
        dds = dds,
        contrast_name = "condition_treated_vs_control",
        out_dir = out_dir
    ))

    testthat::expect_true(!is.null(result$dds))
    testthat::expect_true(result$deseq2_dispersion_fallback$method_used %in% c(
        "parametric",
        "gene-wise"
    ))
    testthat::expect_true(is.logical(result$deseq2_dispersion_fallback$applied))
    if (isTRUE(result$deseq2_dispersion_fallback$applied)) {
        testthat::expect_true(file.exists(file.path(
            out_dir,
            result$deseq2_dispersion_fallback$diagnostic_file
        )))
    }
})

testthat::test_that("de_run_deseq_with_fallback rethrows non-recoverable errors", {
    testthat::skip_if_not_installed("DESeq2")

    out_dir <- tempfile("deseq-fallback-error-")
    dir.create(out_dir)

    testthat::expect_error(
        de_run_deseq_with_fallback(
            dds = list(not = "a DESeqDataSet"),
            contrast_name = "condition_treated_vs_control",
            out_dir = out_dir
        )
    )
})

testthat::test_that("de_estimate_dispersions_with_fallback reports method metadata", {
    testthat::skip_if_not_installed("DESeq2")

    gene_se <- make_test_gene_se()
    sample_df <- data.frame(
        alias = colnames(gene_se),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("b1", "b2", "b1"), 2),
        stringsAsFactors = FALSE
    )
    rownames(sample_df) <- sample_df$alias
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = SummarizedExperiment::assay(gene_se, "counts"),
        colData = sample_df,
        design = ~ batch + condition
    )
    dds <- DESeq2::estimateSizeFactors(dds)

    result <- suppressWarnings(suppressMessages(
        de_estimate_dispersions_with_fallback(dds, "DESeq2")
    ))

    testthat::expect_true(result$method_used %in% c(
        "parametric",
        "local",
        "mean",
        "gene-wise"
    ))
    testthat::expect_true(is.logical(result$fallback_applied))
    testthat::expect_true(!is.null(result$object))
})

testthat::test_that("de_is_recoverable_dexseq_error recognises expected messages", {
    testthat::expect_true(de_is_recoverable_dexseq_error(
        "all gene-wise dispersion estimates are within 2 orders of magnitude"
    ))
    testthat::expect_true(de_is_recoverable_dexseq_error("model matrix is not full rank"))
    testthat::expect_true(de_is_recoverable_dexseq_error("replacement has 1 row, data has 0"))
    # Internal object-shape mismatches should still be reported as failures,
    # but they are not known statistical fallback patterns.
    testthat::expect_false(de_is_recoverable_dexseq_error(
        "nb of cols in 'assay' (7) must equal nb of rows in 'colData' (12)"
    ))
    testthat::expect_true(de_is_recoverable_dexseq_error(
        "newsplit: out of vertex space"
    ))
    testthat::expect_true(de_is_recoverable_dexseq_error(
        "every gene contains at least one zero, cannot compute log geometric means"
    ))
    testthat::expect_false(de_is_recoverable_dexseq_error("random unrelated failure"))
})

testthat::test_that("de_dexseq_failure_hint identifies degenerate-data shape mismatch", {
    hint <- de_dexseq_failure_hint(
        "nb of cols in 'assay' (7) must equal nb of rows in 'colData' (12)"
    )
    testthat::expect_match(hint, "degenerate or near-identical sample groups")
    row_sums_hint <- de_dexseq_failure_hint(
        "Error in rowSums(numerator) : 'x' must be an array of at least two dimensions"
    )
    testthat::expect_match(row_sums_hint, "gene-level DTU summarisation")
    testthat::expect_null(de_dexseq_failure_hint("random unrelated failure"))
})

testthat::test_that("de_postprocess_dexseq_result degrades gracefully on gene-level DTU failure", {
    contrast_dir <- tempfile("dexseq-postprocess-")
    dir.create(contrast_dir)

    dex_res <- list(
        dxr = data.frame(
            featureID = c("tx1", "tx2"),
            groupID = c("gene1", "gene1"),
            log2fold_case_control = c(1, -1),
            pvalue = c(0.01, 0.2),
            padj = c(0.02, 0.3),
            exonBaseMean = c(10, 20),
            check.names = FALSE
        ),
        dxd = structure(list(), class = "mock_dxd")
    )
    contrast_samples <- data.frame(
        alias = c("s1", "s2"),
        condition = c("case", "control"),
        stringsAsFactors = FALSE
    )
    tx_counts <- matrix(
        c(10, 0, 4, 8),
        nrow = 2,
        dimnames = list(c("tx1", "tx2"), c("s1", "s2"))
    )

    result <- testthat::expect_warning(
        de_postprocess_dexseq_result(
            dex_res = dex_res,
            contrast_dir = contrast_dir,
            target_level = "case",
            reference_level = "control",
            contrast_samples = contrast_samples,
            condition_column = "condition",
            tx_counts = tx_counts,
            contrast_name = "case_control",
            per_gene_qvalue_fn = function(...) {
                stop("'x' must be an array of at least two dimensions", call. = FALSE)
            },
            plot_writer = function(...) NULL
        ),
        "DEXSeq post-processing failed"
    )

    testthat::expect_false(result$ok)
    testthat::expect_match(
        result$failure_hint,
        "gene-level DTU summarisation"
    )
    testthat::expect_true(file.exists(file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt")))
    failure_text <- readLines(file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt"))
    testthat::expect_true(any(grepl("DTU Analysis Failed", failure_text, fixed = TRUE)))
    testthat::expect_true(any(grepl("at least two dimensions", failure_text, fixed = TRUE)))
})

testthat::test_that("de_run_dexseq_result records rank-deficiency covariate drops", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    tx_se <- make_test_tx_se(
        sample_names = c(
            "control_rep1", "control_rep2", "control_rep3",
            "treated_rep1", "treated_rep2", "treated_rep3"
        )
    )
    tx_counts <- SummarizedExperiment::assay(tx_se, "counts")
    tx_meta <- as.data.frame(SummarizedExperiment::rowData(tx_se))
    coldata <- data.frame(
        alias = colnames(tx_counts),
        condition = rep(c("control", "treated"), each = 3),
        batch = rep(c("control", "treated"), each = 3),
        stringsAsFactors = FALSE
    )

    seen_messages <- character(0)
    result <- withCallingHandlers(
        tryCatch(
            de_run_dexseq_result(
                tx_counts,
                tx_meta,
                coldata,
                condition_column = "condition",
                covariates = c("batch")
            ),
            error = function(err) err
        ),
        message = function(m) {
            seen_messages <<- c(seen_messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )

    testthat::expect_true(any(grepl(
        "retrying without it",
        seen_messages,
        fixed = TRUE
    )))

    if (inherits(result, "error")) {
        testthat::expect_match(
            conditionMessage(result),
            "model matrix is not full rank",
            fixed = TRUE
        )
    } else {
        testthat::expect_equal(result$dexseq_covariates_dropped, "batch")
        testthat::expect_true(result$dexseq_dispersion_method %in% c(
            "parametric",
            "local",
            "mean",
            "gene-wise"
        ))
        testthat::expect_true(nrow(as.data.frame(result$dxr)) > 0)
    }
})

###
# Contrast planning and output writing
#
# Use small synthetic fixtures with the real DESeq2/DEXSeq path to verify
# workflow-level planning and output writing:
# - One contrast created per non-reference condition level (treated vs control, treated2 vs control)
# - Each contrast subsets to only reference + target samples (not all samples)
# - Output directories created with correct naming

# Multi-level design expands to multiple pairwise contrasts (all vs reference).
# Each contrast subsets samples to just reference + target level.
testthat::test_that("contrasts expanded and samples subsetted", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    fixture_dir <- tempfile("de-multi-")
    dir.create(fixture_dir)

    levels <- c("control", "treated", "treated2")
    bundle <- write_de_fixture_bundle(fixture_dir, levels = levels)

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    suppressWarnings(suppressMessages(main_run_de_analysis(argv)))

    treated_dir <- file.path(argv$out_dir, "condition_treated_vs_control")
    treated2_dir <- file.path(argv$out_dir, "condition_treated2_vs_control")

    treated_samples <- utils::read.delim(
        file.path(treated_dir, "samples_used.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    treated2_samples <- utils::read.delim(
        file.path(treated2_dir, "samples_used.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    testthat::expect_true(dir.exists(treated_dir))
    testthat::expect_true(dir.exists(treated2_dir))
    testthat::expect_equal(
        sort(treated_samples$alias),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated_rep1", "treated_rep2", "treated_rep3"))
    )
    testthat::expect_equal(
        sort(treated2_samples$alias),
        sort(c("control_rep1", "control_rep2", "control_rep3", "treated2_rep1", "treated2_rep2", "treated2_rep3"))
    )
})

# Transcript IDs may contain '|' (e.g., Ensembl IDs like ENST00000123.4|ENSG00000456.7).
# TSV reading defaults to using '|' as separator - verify workflow preserves these IDs.
testthat::test_that("pipe characters in transcript IDs preserved", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    fixture_dir <- tempfile("de-pipe-ids-")
    dir.create(fixture_dir)
    bundle <- write_de_fixture_bundle(fixture_dir, levels = c("control", "treated"))

    tx_se <- readRDS(bundle$transcript_rds)
    pipe_ids <- c(
        "ENSMUST000001|ENSMUSG000001|chr1:1-50",
        "ENSMUST000002|ENSMUSG000001|chr1:101-150",
        "ENSMUST000003|ENSMUSG000002|chr1:201-250",
        "ENSMUST000004|ENSMUSG000002|chr1:301-350"
    )
    rownames(tx_se) <- pipe_ids
    S4Vectors::mcols(SummarizedExperiment::rowRanges(tx_se))$TXNAME <- pipe_ids
    saveRDS(tx_se, bundle$transcript_rds)

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = "batch",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    suppressWarnings(suppressMessages(main_run_de_analysis(argv)))

    contrast_dir <- file.path(argv$out_dir, "condition_treated_vs_control")
    dtu_tx <- utils::read.delim(
        file.path(contrast_dir, "results_dtu_transcript.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    dexseq <- utils::read.delim(
        file.path(contrast_dir, "results_dexseq.tsv"),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    de_qc <- jsonlite::read_json(
        file.path(argv$out_dir, "de_qc_stats.json"),
        simplifyVector = TRUE
    )
    contrast_qc <- de_qc$contrasts[["condition_treated_vs_control"]]

    if (identical(contrast_qc$dtu_status, "SUCCESS")) {
        testthat::expect_equal(dtu_tx$featureID, pipe_ids)
        testthat::expect_equal(dexseq$featureID, pipe_ids)
    } else {
        testthat::expect_true(file.exists(file.path(
            contrast_dir,
            "DTU_ANALYSIS_FAILED.txt"
        )))
        testthat::expect_true(nrow(dtu_tx) == 0)
        testthat::expect_true(nrow(dexseq) == 0)
    }
})

###
# CLI integration tests
#
# End-to-end test with real DESeq2/DRIMSeq (not mocked):
# - Run `supeRglue de_analysis` on synthetic count data
# - Verify output directories and files exist for downstream report generation

testthat::test_that("CLI integration produces expected outputs", {
    fixture_dir <- tempfile("de-cli-")
    dir.create(fixture_dir)
    bundle <- write_de_fixture_bundle(fixture_dir, levels = c("control", "treated"))
    tx_se <- readRDS(bundle$transcript_rds)
    tx_row_data <- as.data.frame(SummarizedExperiment::rowData(tx_se))
    tx_row_data$gene_name <- c("GENEA", "GENEA", "GENEB", "GENEB")
    tx_row_data$transcript_name <- c("TXA", "TXB", "TXC", "TXD")
    SummarizedExperiment::rowData(tx_se) <- S4Vectors::DataFrame(tx_row_data)
    saveRDS(tx_se, bundle$transcript_rds)

    gene_se <- readRDS(bundle$gene_rds)
    gene_row_data <- as.data.frame(SummarizedExperiment::rowData(gene_se))
    gene_row_data$gene_name <- c("GENEA", "GENEB")
    SummarizedExperiment::rowData(gene_se) <- S4Vectors::DataFrame(gene_row_data)
    saveRDS(gene_se, bundle$gene_rds)

    out_dir <- file.path(fixture_dir, "out")

    result <- run_rscript(
        "supeRglue",
        c(
            "de_analysis",
            "--transcript_rds", bundle$transcript_rds,
            "--gene_rds", bundle$gene_rds,
            "--sample_sheet", bundle$sample_sheet,
            "--covariates", "batch",
            "--reference_level", "control",
            "--out_dir", out_dir
        )
    )

    testthat::expect_equal(
        result$status,
        0L,
        info = paste(result$output, collapse = "\n")
    )

    contrast_dir <- file.path(out_dir, "condition_treated_vs_control")
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dge.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_transcript.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_gene.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dexseq.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "samples_used.tsv")))

    dge <- utils::read.delim(file.path(contrast_dir, "results_dge.tsv"), check.names = FALSE)
    dtu_tx <- utils::read.delim(file.path(contrast_dir, "results_dtu_transcript.tsv"), check.names = FALSE)
    dexseq <- utils::read.delim(file.path(contrast_dir, "results_dexseq.tsv"), check.names = FALSE)
    de_qc <- jsonlite::read_json(
        file.path(out_dir, "de_qc_stats.json"),
        simplifyVector = TRUE
    )

    testthat::expect_gt(nrow(dge), 0)
    testthat::expect_true(all(c("GENEID", "gene_name", "log2FoldChange", "padj") %in% names(dge)))
    testthat::expect_true(all(c("featureID", "groupID", "gene_name", "transcript_name", "padj") %in% names(dtu_tx)))
    testthat::expect_true(all(stats::na.omit(dge$gene_name) %in% c("GENEA", "GENEB")))
    testthat::expect_true(all(stats::na.omit(dtu_tx$gene_name) %in% c("GENEA", "GENEB")))
    testthat::expect_true(all(stats::na.omit(dtu_tx$transcript_name) %in% c("TXA", "TXB", "TXC", "TXD")))
    testthat::expect_true(all(c("featureID", "groupID", "padj") %in% names(dexseq)))
    testthat::expect_true("analysis_fallbacks" %in% names(de_qc))
    testthat::expect_true("contrasts" %in% names(de_qc))
    contrast_qc <- de_qc$contrasts[["condition_treated_vs_control"]]
    testthat::expect_true("deseq2_size_factor_method" %in% names(contrast_qc))
    testthat::expect_true("deseq2_dispersion_fallback" %in% names(contrast_qc))
    testthat::expect_true("dge_status" %in% names(contrast_qc))
    testthat::expect_true("dexseq_size_factor_method" %in% names(contrast_qc))
    testthat::expect_true("dexseq_dispersion_method" %in% names(contrast_qc))
    testthat::expect_true("dexseq_covariates_dropped" %in% names(contrast_qc))
})

###
# Degenerate input (control-vs-control / zero-heavy counts)
#
# Exercises the graceful-degradation paths added for three cases:
# - newsplit/lfproc OOM in DESeq2 local regression
# - DEXSeq colData dimension mismatch on degenerate input
# - size-factor estimation failure when every gene has a zero
#
# The process must run to completion and produce placeholder outputs + a
# machine-readable status in de_qc_stats.json rather than crashing.

testthat::test_that("degenerate control-vs-control data completes gracefully", {
    testthat::skip_if_not_installed("DESeq2")
    testthat::skip_if_not_installed("DEXSeq")

    fixture_dir <- tempfile("de-degenerate-")
    dir.create(fixture_dir)
    bundle <- write_degenerate_de_fixture_bundle(fixture_dir)

    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            covariates = NULL,
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())

    # Must complete without throwing — any DESeq2/DEXSeq failure should be
    # caught and degraded to placeholder outputs rather than crashing the process.
    result <- tryCatch(
        suppressWarnings(suppressMessages(main_run_de_analysis(argv))),
        error = function(err) err
    )
    testthat::expect_false(
        inherits(result, "error"),
        info = if (inherits(result, "error")) conditionMessage(result) else ""
    )

    # QC JSON must exist with per-contrast status fields.
    qc_path <- file.path(argv$out_dir, "de_qc_stats.json")
    testthat::expect_true(file.exists(qc_path))
    de_qc <- jsonlite::read_json(qc_path, simplifyVector = TRUE)

    contrast_name <- "condition_case_vs_control"
    contrast_qc <- de_qc$contrasts[[contrast_name]]
    testthat::expect_false(is.null(contrast_qc))
    testthat::expect_true(contrast_qc$dge_status %in% c("SUCCESS", "FAILED"))
    testthat::expect_true(contrast_qc$dtu_status %in% c("SUCCESS", "FAILED"))

    # DGE results file must always exist (empty placeholder on failure).
    contrast_dir <- file.path(argv$out_dir, contrast_name)
    dge_tsv <- file.path(contrast_dir, "results_dge.tsv")
    testthat::expect_true(file.exists(dge_tsv))

    if (identical(contrast_qc$dge_status, "FAILED")) {
        # Failure log and placeholder PDF written.
        testthat::expect_true(file.exists(file.path(contrast_dir, "DGE_ANALYSIS_FAILED.txt")))
        testthat::expect_true(file.exists(file.path(contrast_dir, "results_dge.pdf")))
        # Empty TSV must have the expected column headers.
        dge_df <- utils::read.delim(dge_tsv, check.names = FALSE)
        testthat::expect_equal(nrow(dge_df), 0)
        testthat::expect_true(all(c("GENEID", "log2FoldChange", "padj") %in% names(dge_df)))
    }

    # DTU outputs must always exist.
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_transcript.tsv")))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_gene.tsv")))

    if (identical(contrast_qc$dtu_status, "FAILED")) {
        testthat::expect_true(file.exists(file.path(contrast_dir, "DTU_ANALYSIS_FAILED.txt")))
    }

    # analysis_fallbacks must record failed counts.
    fallbacks <- de_qc$analysis_fallbacks
    testthat::expect_true(!is.null(fallbacks$failed_dge_contrasts))
    testthat::expect_true(!is.null(fallbacks$failed_dtu_contrasts))
})

###
# DGE gene expression prefilter
#
# Optional --de_min_count / --de_min_samples filter applied per contrast,
# after selecting the contrast's samples and before DESeq2 fitting. Fixtures
# from make_gene_filter_counts() place genes on either side of the 10 / 3
# boundaries; see gene_filter_designed_counts() for the expected outcomes.

gene_filter_10_3 <- list(enabled = TRUE, min_count = 10L, min_samples = 3L)

gene_filter_two_level_counts <- function() {
    make_gene_filter_counts(c("control", "treated"))$gene_counts
}

read_de_tsv <- function(path) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

run_gene_filter_de <- function(
    levels = c("control", "treated"),
    de_min_count = NULL,
    de_min_samples = NULL
) {
    fixture_dir <- tempfile("de-gene-filter-")
    dir.create(fixture_dir)
    bundle <- write_gene_filter_fixture_bundle(fixture_dir, levels = levels)
    argv <- c(
        bundle,
        list(
            condition_column = "condition",
            reference_level = "control",
            out_dir = file.path(fixture_dir, "out")
        )
    )
    argv$de_min_count <- de_min_count
    argv$de_min_samples <- de_min_samples
    argv <- workflow_glue_r_normalise_args(argv, de_analysis_arg_spec())
    suppressWarnings(suppressMessages(main_run_de_analysis(argv)))
    list(
        out_dir = argv$out_dir,
        qc = jsonlite::read_json(
            file.path(argv$out_dir, "de_qc_stats.json"),
            simplifyVector = TRUE
        )
    )
}

# DESeq2 results for a count matrix via the unchanged v2.0.3 fitting path,
# keyed by GENEID for comparison with results_dge.tsv.
direct_deseq2_results <- function(count_mat, target_level = "treated") {
    conditions <- sub("_rep[0-9]+$", "", colnames(count_mat))
    coldata <- data.frame(
        alias = colnames(count_mat),
        condition = stats::relevel(factor(conditions), ref = "control"),
        stringsAsFactors = FALSE
    )
    out_dir <- tempfile("direct-deseq2-")
    dir.create(out_dir)
    run <- suppressWarnings(suppressMessages(de_run_deseq2_result(
        count_mat,
        coldata,
        target_level,
        "control",
        "condition",
        character(0),
        out_dir,
        "direct"
    )))
    result <- as.data.frame(run$result)
    result$GENEID <- rownames(result)
    result
}

expect_dge_matches <- function(observed, expected) {
    observed <- observed[match(expected$GENEID, observed$GENEID), , drop = FALSE]
    testthat::expect_equal(observed$GENEID, expected$GENEID)
    for (column in c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) {
        testthat::expect_equal(
            observed[[column]],
            expected[[column]],
            tolerance = 1e-8,
            info = column
        )
    }
}

testthat::test_that("gene prefilter settings: 0/0 disables, pairs enable, partial or invalid rejected", {
    disabled <- de_gene_filter_settings()
    testthat::expect_false(disabled$enabled)
    testthat::expect_identical(de_gene_filter_settings(0L, 0L), disabled)
    testthat::expect_identical(de_gene_filter_settings(NULL, NULL), disabled)

    enabled <- de_gene_filter_settings(10L, 3L)
    testthat::expect_true(enabled$enabled)
    testthat::expect_identical(enabled$min_count, 10L)
    testthat::expect_identical(enabled$min_samples, 3L)

    testthat::expect_error(de_gene_filter_settings(10L, 0L), "must either both be 0")
    testthat::expect_error(de_gene_filter_settings(0L, 3L), "must either both be 0")
    testthat::expect_error(de_gene_filter_settings(-1, 3L), "--de_min_count must be a non-negative integer")
    testthat::expect_error(de_gene_filter_settings(10L, 2.5), "--de_min_samples must be a non-negative integer")
    testthat::expect_error(de_gene_filter_settings("10", 3L), "--de_min_count must be a non-negative integer")
})

testthat::test_that("gene prefilter arguments default to disabled and reject bad values", {
    base_argv <- list(
        transcript_rds = "transcripts.rds",
        gene_rds = "genes.rds",
        sample_sheet = "sample_sheet.csv"
    )
    defaults <- workflow_glue_r_normalise_args(base_argv, de_analysis_arg_spec())
    testthat::expect_identical(defaults$de_min_count, 0L)
    testthat::expect_identical(defaults$de_min_samples, 0L)

    parsed <- workflow_glue_r_normalise_args(
        argparser::parse_args(
            de_analysis_arg_parser(),
            argv = c(
                "--transcript_rds", "t.rds", "--gene_rds", "g.rds", "--sample_sheet", "s.csv",
                "--de_min_count", "10", "--de_min_samples", "3"
            )
        ),
        de_analysis_arg_spec()
    )
    testthat::expect_identical(parsed$de_min_count, 10L)
    testthat::expect_identical(parsed$de_min_samples, 3L)

    testthat::expect_error(
        workflow_glue_r_normalise_args(c(base_argv, list(de_min_count = -1)), de_analysis_arg_spec()),
        "--de_min_count must be a non-negative integer"
    )
    testthat::expect_error(
        workflow_glue_r_normalise_args(c(base_argv, list(de_min_samples = 2.5)), de_analysis_arg_spec()),
        "--de_min_samples must be a non-negative integer"
    )
    # Settings are checked before any input is read.
    testthat::expect_error(
        run_de_analysis_cli(c(
            "--transcript_rds", "missing.rds", "--gene_rds", "missing.rds",
            "--sample_sheet", "missing.csv", "--de_min_count", "10"
        )),
        "must either both be 0"
    )
})

testthat::test_that("disabled gene prefilter returns the count matrix untouched", {
    counts <- gene_filter_two_level_counts()
    result <- de_apply_gene_expression_filter(counts, de_gene_filter_settings())

    # Unrounded values and every row are passed on exactly as before.
    testthat::expect_identical(result$counts, counts)
    testthat::expect_null(result$keep)
    testthat::expect_equal(result$n_input, nrow(counts))
    testthat::expect_equal(result$n_retained, nrow(counts))
    testthat::expect_equal(result$n_removed, 0)
})

testthat::test_that("gene prefilter thresholds are inclusive at both boundaries", {
    result <- de_apply_gene_expression_filter(gene_filter_two_level_counts(), gene_filter_10_3)

    # Exactly 10 counts in exactly 3 samples passes.
    testthat::expect_true(result$keep[["exact_boundary"]])
    # 9 counts everywhere fails the count threshold.
    testthat::expect_false(result$keep[["below_count"]])
    # 10+ counts in only 2 samples fails the sample threshold.
    testthat::expect_false(result$keep[["below_samples"]])
    testthat::expect_equal(sum(result$meets_min_count["below_samples", ]), 2)
})

testthat::test_that("gene prefilter rounds counts before applying thresholds", {
    counts <- gene_filter_two_level_counts()
    result <- de_apply_gene_expression_filter(counts, gene_filter_10_3)

    # 9.6 rounds to 10; unrounded it would fail the threshold.
    testthat::expect_true(all(counts["rounds_up", 1:3] < 10))
    testthat::expect_true(result$keep[["rounds_up"]])
    # R rounds half to even, matching the round() applied for DESeq2.
    testthat::expect_true(result$keep[["rounds_half_even"]])
    testthat::expect_false(result$keep[["rounds_down"]])

    # Retained counts are the rounded values DESeq2 receives.
    testthat::expect_identical(result$counts, round(counts)[result$keep, , drop = FALSE])
    testthat::expect_equal(unname(result$counts["rounds_up", 1:3]), c(10, 10, 10))
})

testthat::test_that("10/3 prefilter removes singleton genes and keeps one-condition genes", {
    result <- de_apply_gene_expression_filter(gene_filter_two_level_counts(), gene_filter_10_3)

    testthat::expect_false(result$keep[["singleton"]])
    testthat::expect_false(result$keep[["all_zero"]])
    # Expression in all three samples of either condition is enough.
    testthat::expect_true(result$keep[["target_only"]])
    testthat::expect_true(result$keep[["reference_only"]])
    # Qualifying samples may come from both conditions.
    testthat::expect_true(result$keep[["mixed_conditions"]])
    testthat::expect_true(all(result$keep[grepl("^bg", names(result$keep))]))
    testthat::expect_equal(result$n_retained, 36)
    testthat::expect_equal(result$n_removed, 7)
})

testthat::test_that("gene prefilter rejects de_min_samples larger than a contrast", {
    sizes <- c(condition_treated_vs_control = 6L, condition_treated2_vs_control = 5L)

    testthat::expect_silent(de_validate_gene_filter_contrasts(gene_filter_10_3, sizes))
    testthat::expect_silent(
        de_validate_gene_filter_contrasts(de_gene_filter_settings(10L, 5L), sizes)
    )
    testthat::expect_error(
        de_validate_gene_filter_contrasts(de_gene_filter_settings(10L, 6L), sizes),
        "condition_treated2_vs_control \\(5 samples\\).*at most 5"
    )
    # Disabled settings are never checked against contrast sizes.
    testthat::expect_silent(de_validate_gene_filter_contrasts(de_gene_filter_settings(), c(x = 0L)))

    fixture_dir <- tempfile("de-gene-filter-impossible-")
    dir.create(fixture_dir)
    bundle <- write_gene_filter_fixture_bundle(fixture_dir)
    argv <- workflow_glue_r_normalise_args(
        c(bundle, list(
            reference_level = "control",
            de_min_count = 10L,
            de_min_samples = 7L,
            out_dir = file.path(fixture_dir, "out")
        )),
        de_analysis_arg_spec()
    )
    testthat::expect_error(
        suppressWarnings(main_run_de_analysis(argv)),
        "--de_min_samples \\(7\\) exceeds the number of samples in contrast: condition_treated_vs_control \\(6 samples\\)"
    )
    # Rejected before any contrast is fitted.
    testthat::expect_false(dir.exists(file.path(argv$out_dir, "condition_treated_vs_control")))
})

testthat::test_that("disabled gene prefilter preserves the unfiltered DESeq2 analysis", {
    run <- run_gene_filter_de()
    contrast <- "condition_treated_vs_control"
    contrast_dir <- file.path(run$out_dir, contrast)
    counts <- gene_filter_two_level_counts()

    expect_dge_matches(
        read_de_tsv(file.path(contrast_dir, "results_dge.tsv")),
        direct_deseq2_results(counts)
    )
    contrast_qc <- run$qc$contrasts[[contrast]]
    testthat::expect_false(run$qc$gene_expression_filter$enabled)
    testthat::expect_false(contrast_qc$gene_filter$enabled)
    testthat::expect_equal(contrast_qc$genes_tested, nrow(counts))
    testthat::expect_equal(contrast_qc$dge_total_genes, nrow(counts))
    testthat::expect_equal(contrast_qc$gene_filter$genes_removed, 0)
    testthat::expect_null(contrast_qc$gene_filter$audit_file)
    testthat::expect_false(file.exists(file.path(contrast_dir, de_gene_filter_audit_filename)))

    # Explicit 0/0 is the same disabled state as the defaults.
    explicit <- run_gene_filter_de(de_min_count = 0L, de_min_samples = 0L)
    testthat::expect_equal(
        read_de_tsv(file.path(explicit$out_dir, contrast, "results_dge.tsv")),
        read_de_tsv(file.path(contrast_dir, "results_dge.tsv")),
        tolerance = 1e-8
    )
})

testthat::test_that("enabled gene prefilter fits DESeq2 on retained genes only and leaves DTU unchanged", {
    unfiltered <- run_gene_filter_de()
    filtered <- run_gene_filter_de(de_min_count = 10L, de_min_samples = 3L)
    contrast <- "condition_treated_vs_control"
    contrast_dir <- file.path(filtered$out_dir, contrast)

    audit <- read_de_tsv(file.path(contrast_dir, de_gene_filter_audit_filename))
    testthat::expect_equal(
        names(audit),
        c(
            "GENEID", "gene_name", "n_samples_meeting_min_count",
            "n_target_samples_meeting_min_count", "n_reference_samples_meeting_min_count",
            "filter_status"
        )
    )
    status <- stats::setNames(audit$filter_status, audit$GENEID)
    testthat::expect_equal(
        sort(names(status)[status == "excluded"]),
        sort(c(
            "singleton", "below_count", "below_samples", "rounds_down",
            "third_only", "third_tips_balance", "all_zero"
        ))
    )
    singleton <- audit[audit$GENEID == "singleton", ]
    testthat::expect_equal(singleton$n_samples_meeting_min_count, 1)
    testthat::expect_equal(singleton$n_target_samples_meeting_min_count, 1)
    target_only <- audit[audit$GENEID == "target_only", ]
    testthat::expect_equal(target_only$n_target_samples_meeting_min_count, 3)
    testthat::expect_equal(target_only$n_reference_samples_meeting_min_count, 0)

    # Only retained genes are tested, and the fit matches DESeq2 run on the
    # filtered matrix, so size factors and dispersions used filtered data.
    dge <- read_de_tsv(file.path(contrast_dir, "results_dge.tsv"))
    retained <- audit$GENEID[audit$filter_status == "retained"]
    testthat::expect_setequal(dge$GENEID, retained)
    counts <- gene_filter_two_level_counts()
    expect_dge_matches(dge, direct_deseq2_results(round(counts)[retained, , drop = FALSE]))
    # Not a post-hoc subset of the unfiltered fit.
    unfiltered_dge <- read_de_tsv(file.path(unfiltered$out_dir, contrast, "results_dge.tsv"))
    background <- sprintf("bg%02d", 1:30)
    testthat::expect_false(isTRUE(all.equal(
        dge$pvalue[match(background, dge$GENEID)],
        unfiltered_dge$pvalue[match(background, unfiltered_dge$GENEID)]
    )))

    contrast_qc <- filtered$qc$contrasts[[contrast]]
    testthat::expect_true(filtered$qc$gene_expression_filter$enabled)
    testthat::expect_equal(filtered$qc$gene_expression_filter$min_count, 10)
    testthat::expect_equal(filtered$qc$gene_expression_filter$min_samples, 3)
    testthat::expect_equal(contrast_qc$gene_filter$genes_input, 43)
    testthat::expect_equal(contrast_qc$gene_filter$genes_retained, 36)
    testthat::expect_equal(contrast_qc$gene_filter$genes_removed, 7)
    testthat::expect_equal(contrast_qc$genes_tested, 36)
    testthat::expect_equal(contrast_qc$dge_total_genes, 36)
    testthat::expect_equal(contrast_qc$dge_status, "SUCCESS")
    testthat::expect_equal(
        contrast_qc$dge_significant_fdr05,
        sum(dge$padj < 0.05, na.rm = TRUE)
    )
    testthat::expect_equal(
        contrast_qc$gene_filter$audit_file,
        file.path(contrast, de_gene_filter_audit_filename)
    )
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dge.pdf")))

    # DTU sees every transcript whether or not the gene prefilter is enabled.
    unfiltered_qc <- unfiltered$qc$contrasts[[contrast]]
    testthat::expect_equal(contrast_qc$transcripts_tested, unfiltered_qc$transcripts_tested)
    testthat::expect_equal(contrast_qc$transcripts_tested, 86)
    testthat::expect_equal(contrast_qc$dtu_status, unfiltered_qc$dtu_status)
    testthat::expect_equal(
        read_de_tsv(file.path(contrast_dir, "results_dexseq.tsv")),
        read_de_tsv(file.path(unfiltered$out_dir, contrast, "results_dexseq.tsv")),
        tolerance = 1e-8
    )
})

testthat::test_that("samples from an excluded third condition do not affect the prefilter", {
    run <- run_gene_filter_de(
        levels = c("control", "treated", "treated2"),
        de_min_count = 10L,
        de_min_samples = 3L
    )
    audit_for <- function(contrast) {
        audit <- read_de_tsv(file.path(run$out_dir, contrast, de_gene_filter_audit_filename))
        audit[match(c("third_only", "third_tips_balance", "target_only"), audit$GENEID), ]
    }

    # control vs treated: treated2 support must not count.
    treated <- audit_for("condition_treated_vs_control")
    testthat::expect_equal(treated$n_samples_meeting_min_count, c(0, 2, 3))
    testthat::expect_equal(treated$filter_status, c("excluded", "excluded", "retained"))

    # control vs treated2: the same genes pass on treated2 samples.
    treated2 <- audit_for("condition_treated2_vs_control")
    testthat::expect_equal(treated2$n_samples_meeting_min_count, c(3, 5, 0))
    testthat::expect_equal(treated2$filter_status, c("retained", "retained", "excluded"))

    dge_treated <- read_de_tsv(file.path(run$out_dir, "condition_treated_vs_control", "results_dge.tsv"))
    testthat::expect_false(any(c("third_only", "third_tips_balance") %in% dge_treated$GENEID))
    testthat::expect_equal(run$qc$contrasts$condition_treated_vs_control$gene_filter$genes_retained, 36)
    testthat::expect_equal(run$qc$contrasts$condition_treated2_vs_control$gene_filter$genes_retained, 35)
})

testthat::test_that("no genes passing the prefilter is a reported DGE failure, not zero DEGs", {
    run <- run_gene_filter_de(de_min_count = 100000L, de_min_samples = 3L)
    contrast <- "condition_treated_vs_control"
    contrast_dir <- file.path(run$out_dir, contrast)
    contrast_qc <- run$qc$contrasts[[contrast]]

    testthat::expect_equal(contrast_qc$dge_status, "FAILED")
    testthat::expect_equal(contrast_qc$dge_failure_reason, "no_genes_passed_expression_filter")
    testthat::expect_match(contrast_qc$dge_failure_message, "No genes passed the gene expression prefilter")
    testthat::expect_equal(contrast_qc$genes_tested, 0)
    testthat::expect_equal(contrast_qc$gene_filter$genes_retained, 0)
    testthat::expect_equal(contrast_qc$gene_filter$genes_removed, 43)
    testthat::expect_equal(run$qc$analysis_fallbacks$failed_dge_contrasts, 1)
    testthat::expect_equal(run$qc$analysis_fallbacks$dge_no_genes_passed_filter_contrasts, 1)
    testthat::expect_equal(
        run$qc$analysis_fallbacks$dge_no_genes_passed_filter_contrast_names,
        contrast
    )

    failure <- readLines(file.path(contrast_dir, "DGE_ANALYSIS_FAILED.txt"))
    testthat::expect_true(any(grepl("No genes passed the gene expression prefilter", failure, fixed = TRUE)))
    testthat::expect_true(any(grepl(
        "rounded count >= 100000 in >= 3 of the 6 contrast samples", failure, fixed = TRUE
    )))
    testthat::expect_true(any(grepl("Genes before prefilter: 43", failure, fixed = TRUE)))
    testthat::expect_true(any(grepl("DGE RESULTS CANNOT BE INTERPRETED", failure, fixed = TRUE)))

    dge <- read_de_tsv(file.path(contrast_dir, "results_dge.tsv"))
    testthat::expect_equal(nrow(dge), 0)
    testthat::expect_true(all(c("GENEID", "log2FoldChange", "padj") %in% names(dge)))
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dge.pdf")))
    audit <- read_de_tsv(file.path(contrast_dir, de_gene_filter_audit_filename))
    testthat::expect_true(all(audit$filter_status == "excluded"))

    summary <- readLines(file.path(contrast_dir, "contrast_qc_summary.txt"))
    testthat::expect_true(any(grepl("Significant (FDR < 0.05): N/A", summary, fixed = TRUE)))
    testthat::expect_true(any(grepl("no genes passed the gene expression prefilter", summary, fixed = TRUE)))
    overall <- readLines(file.path(run$out_dir, "de_overall_summary.txt"))
    testthat::expect_true(any(grepl("no genes passed the prefilter; no test performed", overall, fixed = TRUE)))

    # DTU still runs on all transcripts.
    testthat::expect_equal(contrast_qc$transcripts_tested, 86)
    testthat::expect_true(file.exists(file.path(contrast_dir, "results_dtu_transcript.tsv")))
})

testthat::test_that("CLI propagates prefilter settings to QC JSON, summaries and audit table", {
    fixture_dir <- tempfile("de-gene-filter-cli-")
    dir.create(fixture_dir)
    bundle <- write_gene_filter_fixture_bundle(fixture_dir)
    out_dir <- file.path(fixture_dir, "out")
    base_args <- c(
        "de_analysis",
        "--transcript_rds", bundle$transcript_rds,
        "--gene_rds", bundle$gene_rds,
        "--sample_sheet", bundle$sample_sheet,
        "--reference_level", "control"
    )

    result <- run_rscript(
        "supeRglue",
        c(base_args, "--de_min_count", "10", "--de_min_samples", "3", "--out_dir", out_dir)
    )
    testthat::expect_equal(result$status, 0L, info = paste(result$output, collapse = "\n"))
    testthat::expect_true(any(grepl(
        "condition_treated_vs_control: gene expression prefilter retained 36 of 43 genes (7 removed).",
        result$output,
        fixed = TRUE
    )))

    de_qc <- jsonlite::read_json(file.path(out_dir, "de_qc_stats.json"), simplifyVector = TRUE)
    testthat::expect_equal(
        de_qc$gene_expression_filter[c("enabled", "min_count", "min_samples")],
        list(enabled = TRUE, min_count = 10L, min_samples = 3L)
    )
    gene_filter <- de_qc$contrasts$condition_treated_vs_control$gene_filter
    testthat::expect_equal(
        gene_filter[c("genes_input", "genes_retained", "genes_removed")],
        list(genes_input = 43L, genes_retained = 36L, genes_removed = 7L)
    )
    testthat::expect_true(file.exists(file.path(out_dir, gene_filter$audit_file)))

    summary <- readLines(file.path(out_dir, "condition_treated_vs_control", "contrast_qc_summary.txt"))
    testthat::expect_true(any(grepl(
        "Rule: keep genes with rounded count >= 10 in >= 3 of the 6 contrast samples",
        summary,
        fixed = TRUE
    )))
    testthat::expect_true(any(grepl("Genes tested: 36", summary, fixed = TRUE)))
    overall <- readLines(file.path(out_dir, "de_overall_summary.txt"))
    testthat::expect_true(any(grepl(
        "Gene expression prefilter (DGE only): keep genes with rounded count >= 10 in >= 3",
        overall,
        fixed = TRUE
    )))
    testthat::expect_true(any(grepl("DGE genes tested: 36 of 43 (7 removed by prefilter)", overall, fixed = TRUE)))

    # Expected failures: system2() warns about the non-zero exit status.
    partial <- suppressWarnings(run_rscript(
        "supeRglue",
        c(base_args, "--de_min_samples", "3", "--out_dir", file.path(fixture_dir, "partial"))
    ))
    testthat::expect_false(identical(partial$status, 0L))
    testthat::expect_true(any(grepl("must either both be 0", partial$output, fixed = TRUE)))

    negative <- suppressWarnings(run_rscript(
        "supeRglue",
        c(base_args, "--de_min_count", "-1", "--de_min_samples", "3", "--out_dir", file.path(fixture_dir, "neg"))
    ))
    testthat::expect_false(identical(negative$status, 0L))
    testthat::expect_true(any(grepl("--de_min_count must be a non-negative integer", negative$output, fixed = TRUE)))
})
