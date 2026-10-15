rule audit_dataset_structure:
    input:
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        json=str(RESULTS_DIR / "{dataset}" / "inspection" / "{dataset}.json"),
        summary=str(RESULTS_DIR / "{dataset}" / "inspection" / "summary.tsv")
    params:
        dataset=lambda wc: wc.dataset,
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "inspection")
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "00b_batch_inspect_dataset_structures.py")


rule validate_dataset_audit:
    input:
        json=str(RESULTS_DIR / "{dataset}" / "inspection" / "{dataset}.json"),
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"],
        calls=lambda wc: [CFG["datasets"][wc.dataset]["guide_calls_csv"]] if CFG["datasets"][wc.dataset].get("guide_calls_csv") else []
    output:
        ready=str(RESULTS_DIR / "{dataset}" / "inspection" / "ready.ok")
    params:
        dataset=lambda wc: wc.dataset,
        require_ready=CFG.get("audit", {}).get("require_ready", True),
        require_raw=CFG.get("audit", {}).get("require_raw_source", True),
        pert_col=lambda wc: CFG["datasets"][wc.dataset]["pert_col"],
        control=lambda wc: CFG["datasets"][wc.dataset]["control_label"]
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "00c_validate_audit_ready.py")
