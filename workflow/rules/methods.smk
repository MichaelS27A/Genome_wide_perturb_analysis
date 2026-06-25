rule install_r_method_packages:
    output:
        marker=r_bootstrap_marker()
    params:
        rscript=CFG.get("rscript_bin", "Rscript")
    conda:
        CFG["r_env"]
    shell:
        (
            "{params.rscript} {BASE_DIR}/scripts/11_install_r_method_packages.R "
            "--marker {output.marker}"
        )


rule run_mixscale_method:
    input:
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        audit_ready=maybe_audit_ready_input,
        r_pkgs=(lambda wc: r_bootstrap_marker() if r_bootstrap_enabled() else [])
    output:
        done=str(RESULTS_DIR / "{dataset}" / "mixscale" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "mixscale" / "cell_scores.tsv.gz"),
        de=str(RESULTS_DIR / "{dataset}" / "mixscale" / "perturbation_de.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "mixscale" / "method_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "mixscale"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        min_cells=CFG.get("methods", {}).get("mixscale", {}).get("min_cells_per_perturbation", 30),
        max_pert=CFG.get("methods", {}).get("mixscale", {}).get("max_perturbations", 100),
        max_cells=CFG.get("methods", {}).get("mixscale", {}).get("max_cells", 0),
        seed=CFG.get("methods", {}).get("mixscale", {}).get("random_seed", 0),
        rscript=CFG.get("rscript_bin", "Rscript")
    resources:
        mem_mb=180000,
        runtime=2400
    conda:
        CFG["r_env"]
    shell:
        (
            "{params.rscript} {BASE_DIR}/scripts/08_run_mixscale.R "
            "--h5ad {input.h5ad} "
            "--outdir {params.outdir} "
            "--pert-col {params.pert_col} "
            "--control-label {params.control} "
            "--min-cells-per-perturbation {params.min_cells} "
            "--max-perturbations {params.max_pert} "
            "--max-cells {params.max_cells} "
            "--random-seed {params.seed}"
        )


rule run_ps_method:
    input:
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        audit_ready=maybe_audit_ready_input,
        r_pkgs=(lambda wc: r_bootstrap_marker() if r_bootstrap_enabled() else [])
    output:
        done=str(RESULTS_DIR / "{dataset}" / "ps" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "ps" / "cell_scores.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "ps" / "perturbation_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "ps" / "method_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        min_cells=CFG.get("methods", {}).get("ps", {}).get("min_cells_per_perturbation", 30),
        max_pert=CFG.get("methods", {}).get("ps", {}).get("max_perturbations", 100),
        max_cells=CFG.get("methods", {}).get("ps", {}).get("max_cells", 0),
        seed=CFG.get("methods", {}).get("ps", {}).get("random_seed", 0),
        rscript=CFG.get("rscript_bin", "Rscript")
    resources:
        mem_mb=180000,
        runtime=2400
    conda:
        CFG["r_env"]
    shell:
        (
            "{params.rscript} {BASE_DIR}/scripts/09_run_ps_score.R "
            "--h5ad {input.h5ad} "
            "--outdir {params.outdir} "
            "--pert-col {params.pert_col} "
            "--control-label {params.control} "
            "--min-cells-per-perturbation {params.min_cells} "
            "--max-perturbations {params.max_pert} "
            "--max-cells {params.max_cells} "
            "--random-seed {params.seed}"
        )
