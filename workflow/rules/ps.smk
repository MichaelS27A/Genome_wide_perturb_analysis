def ps_chunk_cell_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "cell_scores.tsv.gz"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


def ps_chunk_summary_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "perturbation_summary.tsv.gz"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


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


rule run_chunk_ps_method:
    input:
        manifest=lambda wc: checkpoints.build_chunk_manifests.get(dataset=wc.dataset).output.manifest,
        chunk_cells=chunk_cells_path,
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        audit_ready=maybe_audit_ready_input,
        r_pkgs=(lambda wc: r_bootstrap_marker() if r_bootstrap_enabled() else [])
    output:
        done=str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "cell_scores.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "perturbation_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "ps" / "chunk_runs" / "chunk_{chunk}" / "method_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps" / "chunk_runs" / f"chunk_{wc.chunk}"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        min_cells=CFG.get("methods", {}).get("ps", {}).get("min_cells_per_perturbation", 30),
        max_pert=CFG.get("methods", {}).get("ps", {}).get("max_perturbations", 0),
        max_cells=CFG.get("methods", {}).get("ps", {}).get("max_cells", 0),
        seed=CFG.get("methods", {}).get("ps", {}).get("random_seed", 0),
        rscript=CFG.get("rscript_bin", "Rscript")
    resources:
        mem_mb=180000,
        runtime=660
    conda:
        CFG["r_env"]
    shell:
        (
            "{params.rscript} {BASE_DIR}/scripts/09_run_ps_score.R "
            "--h5ad {input.h5ad} "
            "--chunk-cells {input.chunk_cells} "
            "--chunk-id {wildcards.chunk} "
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
        cell_scores=ps_chunk_cell_inputs,
        summary=ps_chunk_summary_inputs
    output:
        done=str(RESULTS_DIR / "{dataset}" / "ps" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "ps" / "cell_scores.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "ps" / "perturbation_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "ps" / "method_meta.json")
    params:
        dataset=lambda wc: wc.dataset,
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps")
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "13_merge_ps_chunk_results.py")


rule run_ps_postprocess_de:
    input:
        cell_scores=str(RESULTS_DIR / "{dataset}" / "ps" / "cell_scores.tsv.gz"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "done.txt"),
        de=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "perturbation_differential_genes.tsv.gz"),
        selected=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "selected_ps_cells.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "perturbation_selection_counts.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "de_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps" / "de"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        score_col=CFG.get("ps_de", {}).get("score_column", "ps_score"),
        score_mode=CFG.get("ps_de", {}).get("score_mode", "top_positive"),
        score_quantile=CFG.get("ps_de", {}).get("score_quantile", 0.90),
        min_selected=CFG.get("ps_de", {}).get("min_selected_cells", 20),
        max_control_cells=CFG.get("ps_de", {}).get("max_control_cells", 50000),
        n_top_de_genes=CFG.get("ps_de", {}).get("n_top_de_genes", 200),
        random_seed=CFG.get("ps_de", {}).get("random_seed", 0)
    resources:
        mem_mb=180000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "17_run_ps_postprocess_de.py")
