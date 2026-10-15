"""Record resolved settings, code hashes, input identity and installed versions."""

import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path


def record(snk):
    root = Path(snk.params.repository)
    cfg = json.loads(snk.params.resolved_config)
    files = {}
    for path in snk.input.sources:
        p = Path(path)
        files[str(p.relative_to(root))] = hashlib.sha256(p.read_bytes()).hexdigest()
    inputs = {}
    for name, dataset in cfg["datasets"].items():
        if not dataset.get("enabled", True):
            continue
        inputs[name] = {}
        for key in ("h5ad", "guide_calls_csv", "target_perturbation_list"):
            if dataset.get(key):
                p = Path(dataset[key])
                st = p.stat()
                inputs[name][key] = {"path": str(p), "bytes": st.st_size, "mtime_ns": st.st_mtime_ns}
        inputs[name]["data_id"] = dataset.get("data_id")
        inputs[name]["declared_h5ad_sha256"] = dataset.get("h5ad_sha256")
    versions = {}
    for package in ("snakemake", "pertpy", "anndata", "scanpy", "numpy", "pandas",
                    "scipy", "scikit-learn", "h5py", "umap-learn", "igraph"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    out = Path(snk.output[0])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "resolved_config": cfg, "source_sha256": files, "inputs": inputs,
        "python": platform.python_version(), "package_versions": versions,
        "input_checksum_note": "Declared checksums are recorded, not verified. Sizes and mtimes are not content hashes.",
    }, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    record(snakemake)
