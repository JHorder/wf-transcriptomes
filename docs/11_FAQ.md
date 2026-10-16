### Does the workflow support both cDNA and direct RNA?

Yes. Use `--direct_rna` for direct RNA data. cDNA is the default mode.

### Do I need both `--ref_genome` and `--ref_annotation` in fixed-annotation mode?

Yes. `bambu` still uses the genome together with the imported annotation.

### Does the workflow create one shared transcriptome or one per sample?

Both. The joint cohort model is the primary result for reporting and DE/DTU, and
the workflow also emits independent per-sample transcriptomes.

### Can I run DE/DTU without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--de_analysis`.

### How do I remove low-expression genes before differential gene expression?

Set `--de_min_count` and `--de_min_samples` together, for example
`--de_min_count 10 --de_min_samples 3`. Within each contrast, a gene is passed to
`DESeq2` only if at least 3 of that contrast's samples have a rounded count of at
least 10. See the pipeline overview for details. The prefilter is disabled by
default and does not affect DTU.

### How is the gene expression prefilter different from DESeq2's independent filtering?

The prefilter decides which genes `DESeq2` models at all, and runs before size
factors and dispersions are estimated. `DESeq2`'s independent filtering is
unchanged: it runs after testing, uses the mean of normalised counts to choose
which tested genes receive an adjusted p-value, and sets `padj` to `NA` for the
rest.

### Can I enable the prefilter on a completed analysis without re-running everything?

Yes. Launch the workflow again from the same launch directory (the one
containing the `.nextflow/` folder of the original run), with the same work
directory, inputs, profile and upstream options, add the prefilter options and
`-resume`, and write to a new output directory:

```
nextflow run <this workflow> \
    <original input, reference, sample sheet and DE options> \
    --de_min_count 10 --de_min_samples 3 \
    --out_dir output_prefilter_10_3 \
    -profile <original profile> \
    -resume
```

Both options are part of the differential analysis task's command, so changing
either one makes Nextflow re-run that task and the report, while alignment,
`bambu`, `SQANTI3` and other upstream tasks are reused from the cache. Cached
tasks still publish their outputs, so the new output directory receives a full
set of results, including copies of large files such as BAMs.

+ `-resume` continues the most recent run in the launch directory. If other runs
  have happened since, find the original run name with `nextflow log` and use
  `-resume <run name>`.
+ Keep the same `-w` work directory if one was set, and do not run
  `nextflow clean` on it beforehand.
+ Changing the input files, the profile or container settings, or upstream
  options can invalidate cached tasks.

### Which copy of the differential analysis R code does the workflow run?

The differential analysis task runs `supeRglue de_analysis`. `bin/supeRglue`
loads the R files in `bin/workflow_glue_r/R/` from the copy of the workflow
being run; the container only provides R and the R packages (`DESeq2`,
`DEXSeq`, ...). Nextflow makes the workflow's `bin/` directory available inside
the container, so running the workflow from a checkout that contains a change is
enough; no package installation or container rebuild is needed.

Nextflow's task cache does not track these R files. If you change only files
in `bin/workflow_glue_r/R/` and re-run with `-resume`, the cached differential
analysis task is reused. To force that task to re-run while reusing the rest,
pass a config file with `-c` containing:

```
process {
    withName: runDifferentialAnalysis { cache = false }
}
```

### What changed from the previous workflow version?

The current workflow uses `bambu`, `SQANTI3`, `DESeq2`, and `DEXSeq`.
The main transcriptome result is now one shared `bambu` model built from all
samples together, with separate per-sample `bambu` outputs published alongside
it. The most important differences are summarised in this FAQ section.

#### Why do the results now mention `bambu` and `SQANTI3` instead of StringTie or GffCompare?

The workflow now uses a different set of transcript analysis tools. It builds
its transcript models with `bambu`, classifies them with `SQANTI3`,
and performs DE/DTU from the cohort `bambu` outputs.

#### Why is `--ref_annotation` now required even in fixed-annotation mode?

`bambu` still requires the annotation together with the genome in both
`discover` and `fixed_annotation` modes. Fixed-annotation mode means
annotation-driven quantification, not annotation-only execution.

#### Can I still use `--ref_transcriptome`?

No. `--ref_transcriptome` has been removed from the workflow interface.

Short old-to-new example:

```text
Previous workflow version:
  --transcriptome_source precomputed --ref_transcriptome transcripts.fa

Current workflow:
  --transcriptome_mode fixed_annotation \
  --ref_genome genome.fa \
  --ref_annotation annotation.gtf
```

#### Why do I now get both cohort and per-sample transcriptomes?

The workflow now treats both as important outputs. The shared cohort model is
the main transcriptome used for reporting and optional DE/DTU, while the
per-sample transcriptomes are provided for looking at each sample separately.

#### Why are DE results under `de_analysis/<contrast>/`?

The workflow now writes one subdirectory per comparison instead of publishing
one flat DE result set. This makes the output folder clearer when more than one
comparison is present.

Short old-to-new example:

```text
Previous workflow version:
  de_analysis/results_dge.tsv

Current workflow:
  de_analysis/<contrast>/results_dge.tsv
```

#### Can I still run fixed-annotation quantification without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--ref_genome`,
`--ref_annotation`, and any optional DE/DTU settings.

#### What should I expect to differ in report contents and output filenames?

Expect the report and output folder to emphasise:

+ the joint cohort `bambu` transcriptome under `cohort/`
+ the per-sample `bambu` transcriptomes under `samples/<alias>/`
+ `SQANTI3` results under cohort and per-sample directories
+ contrast-specific DE/DTU outputs under `de_analysis/<contrast>/`

If your question is not answered here, please start a discussion on the
[community](https://community.nanoporetech.com/).
