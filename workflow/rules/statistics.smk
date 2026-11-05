"""Optional raw-count DESeq2 and distribution tests, independent of Mixscape."""


def statistics_inputs(wc):
    dataset = CFG["datasets"][wc.dataset]
    return [dataset["h5ad"]] + ([dataset["guide_calls_csv"]] if dataset.get("guide_calls_csv") else [])


rule build_pseudobulk_counts:
    input:
        data=statistics_inputs,
        support=[str(BASE_DIR / "scripts" / name) for name in ("statistics_support.py", "workflow_support.py")],
        audit_ready=maybe_audit_ready_input
    output:
        counts=str(RESULTS_DIR / "{dataset}" / "pseudobulk_counts" / "counts.h5ad"),
        samples=str(RESULTS_DIR / "{dataset}" / "pseudobulk_counts" / "samples.tsv.gz"),
        qc=str(RESULTS_DIR / "{dataset}" / "pseudobulk_counts" / "pool_qc.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "pseudobulk_counts" / "aggregation_meta.json")
    params:
        dataset=lambda wc: json.dumps(CFG["datasets"][wc.dataset]),
        settings=json.dumps(CFG["pseudobulk_de"]),
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "pseudobulk_counts")
    resources:
        mem_mb=16000,
        runtime=660
    conda:
        CFG["statistics_conda_env"]
    script:
        str(BASE_DIR / "scripts" / "build_pseudobulk_counts.py")


rule run_pseudobulk_de:
    input:
        counts=str(RESULTS_DIR / "{dataset}" / "pseudobulk_counts" / "counts.h5ad"),
        support=str(BASE_DIR / "scripts" / "statistics_support.py")
    output:
        de=str(RESULTS_DIR / "{dataset}" / "pseudobulk_de" / "differential_genes.tsv.gz"),
        qc=str(RESULTS_DIR / "{dataset}" / "pseudobulk_de" / "contrast_qc.tsv.gz"),
        agreement=str(RESULTS_DIR / "{dataset}" / "pseudobulk_de" / "guide_agreement.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "pseudobulk_de" / "de_meta.json")
    params:
        dataset=lambda wc: json.dumps(CFG["datasets"][wc.dataset]),
        settings=json.dumps(CFG["pseudobulk_de"]),
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "pseudobulk_de")
    resources:
        mem_mb=16000,
        runtime=660
    conda:
        CFG["statistics_conda_env"]
    script:
        str(BASE_DIR / "scripts" / "run_pseudobulk_de.py")


rule run_distribution_tests:
    input:
        data=statistics_inputs,
        support=[str(BASE_DIR / "scripts" / name) for name in ("statistics_support.py", "workflow_support.py")],
        audit_ready=maybe_audit_ready_input
    output:
        energy=str(RESULTS_DIR / "{dataset}" / "distribution_tests" / "energy_distance.tsv.gz"),
        ad=str(RESULTS_DIR / "{dataset}" / "distribution_tests" / "anderson_darling_permutation.tsv.gz"),
        qc=str(RESULTS_DIR / "{dataset}" / "distribution_tests" / "contrast_qc.tsv.gz"),
        meta=str(RESULTS_DIR / "{dataset}" / "distribution_tests" / "distribution_meta.json"),
        basis=str(RESULTS_DIR / "{dataset}" / "distribution_tests" / "basis.npz")
    params:
        dataset=lambda wc: json.dumps(CFG["datasets"][wc.dataset]),
        settings=json.dumps(CFG["distribution_tests"]),
        energy_enabled=ENERGY_DISTANCE_ENABLED,
        ad_enabled=AD_PERMUTATION_ENABLED,
        outdir=lambda wc: str(RESULTS_DIR / wc.dataset / "distribution_tests")
    resources:
        mem_mb=16000,
        runtime=660
    conda:
        CFG["statistics_conda_env"]
    script:
        str(BASE_DIR / "scripts" / "run_distribution_tests.py")
