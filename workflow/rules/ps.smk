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
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps")
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    shell:
        (
            "python {BASE_DIR}/scripts/13_merge_ps_chunk_results.py "
            "--dataset {wildcards.dataset} "
            "--chunk-cell-scores {input.cell_scores} "
            "--chunk-summary {input.summary} "
            "--outdir {params.outdir}"
        )
