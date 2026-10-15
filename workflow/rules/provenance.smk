rule record_run:
    input:
        sources=WORKFLOW_SOURCES,
        data_files=sorted({
            dataset[key]
            for name, dataset in CFG["datasets"].items()
            if name in DATASETS
            for key in ("h5ad", "guide_calls_csv", "target_perturbation_list")
            if dataset.get(key)
        })
    output:
        str(RESULTS_DIR / "run_provenance.json")
    params:
        resolved_config=json.dumps(CFG, sort_keys=True),
        repository=str(ROOT_DIR)
    conda:
        CFG["conda_env"]
    script:
        str(BASE_DIR / "scripts" / "record_run.py")
