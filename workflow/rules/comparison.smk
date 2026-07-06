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
        dataset=lambda wc: wc.dataset,
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "comparison")
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "10_compare_methods.py")
