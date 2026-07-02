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
        pca_gene_flag="--use-hvg-for-pca" if CFG.get("mixscale", {}).get("use_hvg_for_pca", False) else ""
    resources:
        mem_mb=180000,
        runtime=660
    conda:
        CFG["conda_env"]
    shell:
        (
            "export OMP_NUM_THREADS=1; "
            "export OPENBLAS_NUM_THREADS=1; "
            "export MKL_NUM_THREADS=1; "
            "export VECLIB_MAXIMUM_THREADS=1; "
            "export NUMEXPR_NUM_THREADS=1; "
            "python {BASE_DIR}/scripts/08_run_mixscale.py "
            "--h5ad {input.h5ad} "
            "--chunk-cells {input.chunk_cells} "
            "--chunk-id {wildcards.chunk} "
            "--outdir {params.outdir} "
            "--pert-col {params.pert_col} "
            "--control-label {params.control} "
            "--min-cells-per-perturbation {params.min_cells} "
            "--max-perturbations {params.max_pert} "
            "--max-cells {params.max_cells} "
            "--random-seed {params.seed} "
            "--normalize-target-sum {params.normalize_target_sum} "
            "--mixscale-logfc-threshold {params.logfc_threshold} "
            "--mixscale-pval-cutoff {params.pval_cutoff} "
            "--mixscale-min-de-genes {params.min_de_genes} "
            "--mixscale-max-de-genes {params.max_de_genes} "
            "--batch-size {params.batch_size} "
            "{params.pca_gene_flag}"
        )


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
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "mixscale")
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    shell:
        (
            "python {BASE_DIR}/scripts/12_merge_mixscale_chunk_results.py "
            "--dataset {wildcards.dataset} "
            "--chunk-cell-scores {input.cell_scores} "
            "--chunk-de {input.de} "
            "--outdir {params.outdir}"
        )
