rule audit_dataset_structure:
    input:
        h5ad=lambda wc: CFG["datasets"][wc.dataset]["h5ad"]
    output:
        json=str(RESULTS_DIR / "{dataset}" / "inspection" / "{dataset}.json"),
        summary=str(RESULTS_DIR / "{dataset}" / "inspection" / "summary.tsv")
    params:
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "inspection")
    conda:
        CFG["conda_env"]
    shell:
        (
            "python {BASE_DIR}/scripts/00b_batch_inspect_dataset_structures.py "
            "--dataset {wildcards.dataset}={input.h5ad} "
            "--outdir {params.outdir}"
        )


rule validate_dataset_audit:
    input:
        json=str(RESULTS_DIR / "{dataset}" / "inspection" / "{dataset}.json")
    output:
        ready=str(RESULTS_DIR / "{dataset}" / "inspection" / "ready.ok")
    params:
        require_ready="--require-ready" if CFG.get("audit", {}).get("require_ready", True) else "",
        require_raw="--require-raw-source" if CFG.get("audit", {}).get("require_raw_source", True) else ""
    conda:
        CFG["conda_env"]
    shell:
        (
            "python {BASE_DIR}/scripts/00c_validate_audit_ready.py "
            "--audit-json {input.json} "
            "--dataset {wildcards.dataset} "
            "{params.require_ready} {params.require_raw} && "
            "echo ok > {output.ready}"
        )
