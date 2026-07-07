checkpoint build_chunk_manifests:
    input:
        source=lambda wc: CFG["datasets"][wc.dataset].get(
            "guide_calls_csv", CFG["datasets"][wc.dataset]["h5ad"]
        ),
        audit_ready=maybe_audit_ready_input
    output:
        manifest=str(RESULTS_DIR / "{dataset}" / "manifest.json"),
        summary=str(RESULTS_DIR / "{dataset}" / "chunk_summary.tsv"),
        chunkdir=directory(str(RESULTS_DIR / "{dataset}" / "chunks"))
    params:
        source_type=lambda wc: "csv" if CFG["datasets"][wc.dataset].get("guide_calls_csv") else "h5ad",
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        ppc=CFG["chunking"]["perturbations_per_chunk"],
        min_cells=CFG["chunking"]["min_cells_per_perturbation"],
        max_controls=CFG["chunking"]["max_controls_per_chunk"],
        control_seed=CFG["chunking"].get("control_sample_seed", 0),
        read_chunk=CFG["chunking"]["csv_read_chunk_size"]
    resources:
        mem_mb=350000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "04_build_mixscape_chunk_manifests.py")


rule run_chunk_mixscape:
    input:
        manifest=lambda wc: checkpoints.build_chunk_manifests.get(dataset=wc.dataset).output.manifest,
        chunk_cells=chunk_cells_path,
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "done.txt"),
        stats=str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "perturbation_stats.tsv"),
        effects=str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "perturbation_effects_pca.tsv"),
        labels=str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "cell_mixscape_labels.tsv.gz"),
        selected_cells=str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "selected_perturbed_cells.tsv.gz")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "chunk_runs" / f"chunk_{wc.chunk}"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        pca_dims=CFG["mixscape"]["pca_dims"],
        batch_size=CFG["mixscape"].get("batch_size", 0),
        auto_batch_max_elements=CFG["mixscape"].get("auto_batch_max_elements", 800000000),
        auto_batch_size=CFG["mixscape"].get("auto_batch_size", 2000),
        csc_max_genes=CFG["mixscape"].get("csc_max_genes", 1000),
        csc_max_total_nnz=CFG["mixscape"].get("csc_max_total_nnz", 120000000),
        normalize_target_sum=CFG["mixscape"].get("normalize_target_sum", 10000),
        logfc_threshold=CFG["mixscape"].get("logfc_threshold", 0.10),
        pval_cutoff=CFG["mixscape"].get("pval_cutoff", 0.05),
        write_subset=CFG["mixscape"]["write_subset_h5ad"],
        use_hvg_for_pca=CFG["mixscape"].get("use_hvg_for_pca", False)
    resources:
        mem_mb=200000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "05_run_mixscape_chunk.py")


rule merge_dataset_results:
    input:
        stats=chunk_stats_inputs,
        effects=chunk_effects_inputs,
        selected_cells=chunk_selected_inputs
    output:
        stats_merged=str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_stats_merged.tsv.gz"),
        effects_merged=str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_effects_pca_merged.tsv.gz"),
        clusters=str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_clusters.tsv"),
        selected_cells_merged=str(RESULTS_DIR / "{dataset}" / "merged" / "selected_perturbed_cells_merged.tsv.gz"),
        selection_summary=str(RESULTS_DIR / "{dataset}" / "merged" / "mixscape_selection_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "merged" / "merge_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "merged"),
        n_clusters=CFG["clustering"]["n_clusters"]
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "06_merge_mixscape_chunk_results.py")


rule postprocess_dataset:
    input:
        selected_cells=lambda wc: str(RESULTS_DIR / wc.dataset / "merged" / "selected_perturbed_cells_merged.tsv.gz"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        meta=str(RESULTS_DIR / "{dataset}" / "postprocess" / "postprocess_meta.json"),
        pseudobulk=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_pseudobulk.h5ad"),
        umap_leiden=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_umap_leiden.tsv.gz"),
        de=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_differential_genes.tsv.gz"),
        long_table=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_long_table.tsv.gz"),
        gene_long_table=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_gene_long_table.tsv.gz")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "postprocess"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        min_selected=CFG["postprocess"]["min_selected_cells"],
        max_controls=CFG["postprocess"]["max_control_cells_for_de"],
        n_top=CFG["postprocess"]["n_top_de_genes"],
        seed=CFG["postprocess"]["random_seed"]
    resources:
        mem_mb=180000,
        runtime=660
    conda:
        CFG["postprocess_conda_env"]
    script:
        str(BASE_DIR / "scripts" / "07_perturbation_embedding_and_de.py")
