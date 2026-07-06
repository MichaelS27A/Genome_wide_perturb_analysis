rule run_global_de_all_cells:
    input:
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        audit_ready=maybe_audit_ready_input
    output:
        done=str(RESULTS_DIR / "{dataset}" / "global_de" / "done.txt"),
        de=str(RESULTS_DIR / "{dataset}" / "global_de" / "perturbation_differential_genes.tsv.gz"),
        group_counts=str(RESULTS_DIR / "{dataset}" / "global_de" / "group_cell_counts.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "global_de" / "de_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "global_de"),
        pert_col=lambda wc: CFG["datasets"][wc.dataset].get("pert_col", "gene_target"),
        control=lambda wc: CFG["datasets"][wc.dataset].get("control_label", "Non-Targeting"),
        method=CFG.get("global_de", {}).get("method", "wilcoxon"),
        normalize_target_sum=CFG.get("global_de", {}).get("normalize_target_sum", 10000),
        log1p=CFG.get("global_de", {}).get("log1p", True),
        n_top_de_genes=CFG.get("global_de", {}).get("n_top_de_genes", 2000),
        min_cells=CFG.get("global_de", {}).get("min_cells_per_perturbation", 30),
        max_control_cells=CFG.get("global_de", {}).get("max_control_cells", 0),
        max_cells_per_pert=CFG.get("global_de", {}).get("max_cells_per_perturbation", 0),
        seed=CFG.get("global_de", {}).get("random_seed", 0),
    resources:
        mem_mb=200000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "15_run_global_de_all_cells.py")


rule run_global_de_ora:
    input:
        de=str(RESULTS_DIR / "{dataset}" / "global_de" / "perturbation_differential_genes.tsv.gz"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "global_de" / "ora" / "done.txt"),
        terms=str(RESULTS_DIR / "{dataset}" / "global_de" / "ora" / "ora_terms.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "global_de" / "ora" / "ora_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "global_de" / "ora" / "ora_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "global_de" / "ora"),
        gmt_files=lambda wc: [str(resolve_from_root(p)) for p in CFG.get("ora", {}).get("gmt_files", [])],
        fdr_alpha=CFG.get("ora", {}).get("fdr_alpha", 0.05),
        min_abs_logfc=CFG.get("ora", {}).get("min_abs_logfc", 0.25),
        min_deg_genes=CFG.get("ora", {}).get("min_deg_genes", 10),
        max_terms_per_direction=CFG.get("ora", {}).get("max_terms_per_direction", 50),
        min_term_size=CFG.get("ora", {}).get("min_term_size", 5),
        max_term_size=CFG.get("ora", {}).get("max_term_size", 5000),
    resources:
        mem_mb=64000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "16_run_ora_from_global_de.py")


rule run_mixscape_postprocess_ora:
    input:
        de=str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_differential_genes.tsv.gz"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "postprocess" / "ora" / "done.txt"),
        terms=str(RESULTS_DIR / "{dataset}" / "postprocess" / "ora" / "ora_terms.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "postprocess" / "ora" / "ora_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "postprocess" / "ora" / "ora_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "postprocess" / "ora"),
        gmt_files=lambda wc: [str(resolve_from_root(p)) for p in CFG.get("ora", {}).get("gmt_files", [])],
        fdr_alpha=CFG.get("ora", {}).get("fdr_alpha", 0.05),
        min_abs_logfc=CFG.get("ora", {}).get("min_abs_logfc", 0.25),
        min_deg_genes=CFG.get("ora", {}).get("min_deg_genes", 10),
        max_terms_per_direction=CFG.get("ora", {}).get("max_terms_per_direction", 50),
        min_term_size=CFG.get("ora", {}).get("min_term_size", 5),
        max_term_size=CFG.get("ora", {}).get("max_term_size", 5000),
    resources:
        mem_mb=64000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "16_run_ora_from_global_de.py")


rule run_ps_de_ora:
    input:
        de=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "perturbation_differential_genes.tsv.gz"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        done=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "ora" / "done.txt"),
        terms=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "ora" / "ora_terms.tsv.gz"),
        summary=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "ora" / "ora_summary.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "ps" / "de" / "ora" / "ora_meta.json")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "ps" / "de" / "ora"),
        gmt_files=lambda wc: [str(resolve_from_root(p)) for p in CFG.get("ora", {}).get("gmt_files", [])],
        fdr_alpha=CFG.get("ora", {}).get("fdr_alpha", 0.05),
        min_abs_logfc=CFG.get("ora", {}).get("min_abs_logfc", 0.25),
        min_deg_genes=CFG.get("ora", {}).get("min_deg_genes", 10),
        max_terms_per_direction=CFG.get("ora", {}).get("max_terms_per_direction", 50),
        min_term_size=CFG.get("ora", {}).get("min_term_size", 5),
        max_term_size=CFG.get("ora", {}).get("max_term_size", 5000),
    resources:
        mem_mb=64000,
        runtime=660
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "16_run_ora_from_global_de.py")
