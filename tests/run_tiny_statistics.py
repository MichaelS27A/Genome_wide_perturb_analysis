"""Execute the optional statistics branches and check their planted mean effects."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from make_tiny_statistics import make_example
from build_pseudobulk_counts import run as build
from run_pseudobulk_de import run as de
from run_distribution_tests import run as distributions


def execute(config, workflow):
    cfg = yaml.safe_load(config.read_text())
    if workflow:
        subprocess.run([sys.executable, "-m", "snakemake", "-s", str(ROOT / "workflow/Snakefile"),
                        "--directory", str(ROOT), "--configfile", str(config),
                        "--cores", "1", "--rerun-incomplete"], check=True)
    else:
        dataset = cfg["datasets"]["TINY_STATS"]
        base = Path(cfg["results_dir"]) / "TINY_STATS"
        build(dataset, cfg["pseudobulk_de"], base / "pseudobulk_counts")
        de(base / "pseudobulk_counts/counts.h5ad", dataset, cfg["pseudobulk_de"], base / "pseudobulk_de")
        distributions(dataset, cfg["distribution_tests"], True, True, base / "distribution_tests")
    return Path(cfg["results_dir"]) / "TINY_STATS"


def verify(base):
    results = pd.read_csv(base / "pseudobulk_de/differential_genes.tsv.gz", sep="\t")
    primary = results.loc[results["analysis_level"].eq("target")]
    self_effect = primary.loc[primary["perturbation"].eq(primary["gene"])]
    assert len(self_effect) >= 2 and (self_effect["log2FoldChange"] < -1).all()
    induced = primary.loc[primary["perturbation"].eq("TARGET_A") & primary["gene"].isin([f"G{i}" for i in range(2, 10)])]
    assert (induced["log2FoldChange"] > 1).all() and (induced["padj"] < .05).all()
    qc = pd.read_csv(base / "pseudobulk_de/contrast_qc.tsv.gz", sep="\t")
    assert qc.loc[qc["analysis_level"].isin(["target", "ntc_null"]), "status"].eq("ok").all()
    assert qc["analysis_level"].eq("ntc_null").any()
    assert qc.loc[qc["status"].eq("ok"), "dispersion_fit"].eq("mean").all()
    energy = pd.read_csv(base / "distribution_tests/energy_distance.tsv.gz", sep="\t")
    assert (energy["pvalue"] >= .01).all() and len(energy) == 12
    ad = pd.read_csv(base / "distribution_tests/anderson_darling_permutation.tsv.gz", sep="\t")
    assert len(ad) == 12 * 16 and ad["padj"].between(0, 1).all()
    meta = json.loads((base / "distribution_tests/distribution_meta.json").read_text())
    assert meta["paper_reproduction"] is False and len(meta["basis"]) == 4


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mode", choices=["biological_sample", "guide_batch", "kolf_resample"],
                        default="biological_sample")
    parser.add_argument("--workflow", action="store_true")
    parser.add_argument("--repeat", action="store_true")
    args = parser.parse_args()
    config = make_example(args.outdir, args.mode)
    first = execute(config, args.workflow)
    verify(first)
    if args.repeat:
        cfg = yaml.safe_load(config.read_text())
        cfg["results_dir"] = str(Path(cfg["results_dir"]).parent / "repeat_results")
        config.write_text(yaml.safe_dump(cfg, sort_keys=False))
        second = execute(config, args.workflow)
        verify(second)
        for relative in ("pseudobulk_de/differential_genes.tsv.gz", "pseudobulk_de/guide_agreement.tsv.gz",
                         "distribution_tests/energy_distance.tsv.gz",
                         "distribution_tests/anderson_darling_permutation.tsv.gz"):
            pd.testing.assert_frame_equal(pd.read_csv(first / relative, sep="\t"),
                                          pd.read_csv(second / relative, sep="\t"),
                                          rtol=1e-6, atol=1e-6)
    print(f"PASS: optional statistics ({args.mode}); repeat={args.repeat}")
