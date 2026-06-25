rule compare_methods:
    input:
        mixscape_stats=str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_stats_merged.tsv.gz"),
        mixscale_cell=str(RESULTS_DIR / "{dataset}" / "mixscale" / "cell_scores.tsv.gz"),
        mixscale_de=str(RESULTS_DIR / "{dataset}" / "mixscale" / "perturbation_de.tsv.gz"),
        ps_cell=str(RESULTS_DIR / "{dataset}" / "ps" / "cell_scores.tsv.gz"),
        ps_summary=str(RESULTS_DIR / "{dataset}" / "ps" / "perturbation_summary.tsv.gz")
    output:
        long=str(RESULTS_DIR / "{dataset}" / "comparison" / "comparison_long.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "comparison" / "comparison_summary.tsv"),
        meta=str(RESULTS_DIR / "{dataset}" / "comparison" / "method_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "comparison")
    conda:
        CFG["conda_env"]
    shell:
        (
            "python {BASE_DIR}/scripts/10_compare_methods.py "
            "--dataset {wildcards.dataset} "
            "--mixscape-stats {input.mixscape_stats} "
            "--mixscale-cell-scores {input.mixscale_cell} "
            "--mixscale-de {input.mixscale_de} "
            "--ps-cell-scores {input.ps_cell} "
            "--ps-summary {input.ps_summary} "
            "--outdir {params.outdir}"
        )
