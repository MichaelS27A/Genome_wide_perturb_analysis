def mixscale_chunk_cell_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "cell_scores.tsv.gz"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


def mixscale_chunk_de_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "perturbation_de.tsv.gz"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


rule run_chunk_mixscale_method:
    input:
        manifest=lambda wc: checkpoints.build_chunk_manifests.get(dataset=wc.dataset).output.manifest,
        chunk_cells=chunk_cells_path,
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        audit_ready=maybe_audit_ready_input
    output:
        done=str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "cell_scores.tsv.gz"),
        de=str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "perturbation_de.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "mixscale" / "chunk_runs" / "chunk_{chunk}" / "method_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "mixscale" / "chunk_runs" / f"chunk_{wc.chunk}"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        min_cells=CFG.get("methods", {}).get("mixscale", {}).get("min_cells_per_perturbation", 30),
        max_pert=CFG.get("methods", {}).get("mixscale", {}).get("max_perturbations", 0),
        max_cells=CFG.get("methods", {}).get("mixscale", {}).get("max_cells", 0),
        seed=CFG.get("methods", {}).get("mixscale", {}).get("random_seed", 0),
        normalize_target_sum=CFG.get("mixscale", {}).get("normalize_target_sum", 10000),
        logfc_threshold=CFG.get("mixscale", {}).get("logfc_threshold", 0.10),
        pval_cutoff=CFG.get("mixscale", {}).get("pval_cutoff", 0.05),
        min_de_genes=CFG.get("mixscale", {}).get("min_de_genes", 5),
        max_de_genes=CFG.get("mixscale", {}).get("max_de_genes", 100),
        batch_size=CFG.get("mixscale", {}).get("batch_size", 0),
        auto_batch_max_elements=CFG.get("mixscale", {}).get("auto_batch_max_elements", 800000000),
        auto_batch_size=CFG.get("mixscale", {}).get("auto_batch_size", 2000),
        csc_max_genes=lambda wc: CFG["datasets"][wc.dataset].get(
            "csc_max_genes", CFG.get("mixscale", {}).get("csc_max_genes", 0)
        ),
        csc_max_total_nnz=lambda wc: CFG["datasets"][wc.dataset].get(
            "csc_max_total_nnz", CFG.get("mixscale", {}).get("csc_max_total_nnz", 0)
        ),
        use_hvg_for_pca=CFG.get("mixscale", {}).get("use_hvg_for_pca", False)
    resources:
        mem_mb=200000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "08_run_mixscale.py")


rule run_mixscale_method:
    input:
        cell_scores=mixscale_chunk_cell_inputs,
        de=mixscale_chunk_de_inputs
    output:
        done=str(RESULTS_DIR / "{dataset}" / "mixscale" / "done.txt"),
        cell_scores=str(RESULTS_DIR / "{dataset}" / "mixscale" / "cell_scores.tsv.gz"),
        de=str(RESULTS_DIR / "{dataset}" / "mixscale" / "perturbation_de.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "mixscale" / "method_meta.json")
    params:
        dataset=lambda wc: wc.dataset,
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "mixscale")
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "12_merge_mixscale_chunk_results.py")
