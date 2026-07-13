def _target_cluster_de_path(dataset, stream):
    if stream == "global_de":
        return str(RESULTS_DIR / dataset / "global_de" / "perturbation_differential_genes.tsv.gz")
    if stream == "mixscape":
        return str(RESULTS_DIR / dataset / "postprocess" / "perturbation_differential_genes.tsv.gz")
    if stream == "mixscale":
        return str(RESULTS_DIR / dataset / "mixscale" / "perturbation_de.tsv.gz")
    if stream == "ps":
        return str(RESULTS_DIR / dataset / "ps" / "de" / "perturbation_differential_genes.tsv.gz")
    raise ValueError(f"Unsupported target-random clustering DEG stream: {stream}")


def _target_cluster_de_inputs(wc):
    return [_target_cluster_de_path(wc.dataset, s) for s in TARGET_CLUSTER_DE_STREAMS]


def _target_cluster_ora_path(dataset, stream):
    if stream == "global_de":
        return str(RESULTS_DIR / dataset / "global_de" / "ora" / "ora_terms.tsv.gz")
    if stream == "mixscape":
        return str(RESULTS_DIR / dataset / "postprocess" / "ora" / "ora_terms.tsv.gz")
    if stream == "ps":
        return str(RESULTS_DIR / dataset / "ps" / "de" / "ora" / "ora_terms.tsv.gz")
    raise ValueError(f"Unsupported ORA stream for target-random clustering: {stream}")


def _target_cluster_ora_inputs(wc):
    if not ORA_ENABLED:
        return []
    streams = [s for s in TARGET_CLUSTER_DE_STREAMS if s in ORA_STREAMS and s in {"global_de", "mixscape", "ps"}]
    return [_target_cluster_ora_path(wc.dataset, s) for s in streams]


rule run_target_random_clustering:
    input:
        de_tables=_target_cluster_de_inputs,
        ora_tables=_target_cluster_ora_inputs,
        mixscape_selected=str(RESULTS_DIR / "{dataset}" / "merged" / "selected_perturbed_cells_merged.tsv.gz"),
        mixscale_cells=str(RESULTS_DIR / "{dataset}" / "mixscale" / "cell_scores.tsv.gz"),
        chunk_summary=lambda wc: checkpoints.build_chunk_manifests.get(dataset=wc.dataset).output.summary,
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        target_list=lambda wc: CFG["datasets"][wc.dataset]["target_perturbation_list"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "done.txt"),
        selected=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "selected_perturbations.tsv"),
        chunk_summary_copy=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "chunk_summary.tsv"),
        expression_matrix=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "pseudobulk_expression_matrix.tsv.gz"),
        deg_matrix=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "deg_profile_matrix.tsv.gz"),
        clusters_expression=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "perturbation_clusters_expression.tsv"),
        clusters_deg=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "perturbation_clusters.tsv"),
        heatmap_expression=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "heatmap_expression.png"),
        heatmap_deg=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "heatmap_deg_profile.png"),
        volcano_manifest=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "volcano_plots.tsv"),
        ora_terms_targets=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "ora_terms_selected_targets.tsv.gz"),
        ora_summary_targets=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "ora_summary_selected_targets.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "clustering_meta.json"),
        methods_summary=str(RESULTS_DIR / "{dataset}" / "target_random_clustering" / "method_clustering_summary.tsv")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "target_random_clustering"),
        dataset=lambda wc: wc.dataset,
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        deg_streams=",".join(TARGET_CLUSTER_DE_STREAMS),
        ora_streams=",".join([s for s in TARGET_CLUSTER_DE_STREAMS if s in ORA_STREAMS and s in {"global_de", "mixscape", "ps"}]) if ORA_ENABLED else "",
        random_match_target_count=CFG.get("target_random_clustering", {}).get("random_controls_match_target_count", True),
        n_random=CFG.get("target_random_clustering", {}).get("n_random_controls", 0),
        min_abs_effect=CFG.get("target_random_clustering", {}).get("min_abs_effect", 0.25),
        max_pval_adj=CFG.get("target_random_clustering", {}).get("max_pval_adj", 0.05),
        rank_significant_max=CFG.get("target_random_clustering", {}).get("rank_significant_max", 100),
        top_genes_per_pert=CFG.get("target_random_clustering", {}).get("top_genes_per_perturbation", 0),
        min_feature_genes=CFG.get("target_random_clustering", {}).get("min_feature_genes", 20),
        assignment_logic=CFG.get("target_random_clustering", {}).get("assignment_logic", "intersection"),
        mixscale_score_threshold=CFG.get("target_random_clustering", {}).get("mixscale_score_threshold", 0.0),
        max_control_cells=CFG.get("target_random_clustering", {}).get("max_control_cells", 50000),
        expression_normalize_target_sum=CFG.get("target_random_clustering", {}).get("expression_normalize_target_sum", 10000),
        expression_log1p=CFG.get("target_random_clustering", {}).get("expression_log1p", True),
        heatmap_max_genes=CFG.get("target_random_clustering", {}).get("heatmap_max_genes", 200),
        heatmap_gene_order=CFG.get("target_random_clustering", {}).get("heatmap_gene_order", "variance"),
        volcano_max_targets=CFG.get("target_random_clustering", {}).get("volcano_max_targets", 20),
        recurrent_gene_fraction=CFG.get("target_random_clustering", {}).get("recurrent_gene_fraction", 0.25),
        volcano_stream_priority=",".join(CFG.get("target_random_clustering", {}).get("volcano_stream_priority", ["global_de", "mixscape", "ps"])),
        n_clusters=CFG.get("clustering", {}).get("n_clusters", 20),
        linkage=CFG.get("target_random_clustering", {}).get("linkage", "average"),
        distance=CFG.get("target_random_clustering", {}).get("distance", "euclidean"),
        seed=CFG.get("target_random_clustering", {}).get("random_seed", 0),
    resources:
        mem_mb=32000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "18_cluster_target_vs_random_de.py")
