"""
uniPort integration pipeline (from Uniport_integration.py).
Reads 10x MTX datasets from config, runs uniPort, saves h5ad with X_uniport for scib.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import scanpy as sc
import uniport as up
import yaml

def read_10x_folder(path: str, prefix: str):
    ad = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True)
    ad.obs_names = [f"{prefix}_{bc}" for bc in ad.obs_names]
    ad.obs["source"] = prefix
    return ad

def prep_common_hvg(adata_cm, n_top_genes: int):
    sc.pp.normalize_total(adata_cm)
    sc.pp.log1p(adata_cm)
    sc.pp.highly_variable_genes(adata_cm, n_top_genes=n_top_genes, subset=True)
    up.batch_scale(adata_cm)
    return adata_cm

def prep_specific_hvg(ad, n_top_genes: int):
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes, subset=True)
    up.batch_scale(ad)
    return ad

def main():
    out_h5ad = snakemake.output.h5ad
    out_dir = Path(snakemake.params.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    figdir = out_dir / "figs"
    figdir.mkdir(parents=True, exist_ok=True)

    with open(snakemake.input.config) as f:
        cfg = yaml.safe_load(f)

    data_cfg = cfg["data"]
    datasets = [(d["path"], d["prefix"]) for d in data_cfg["datasets"]]
    u_cfg = cfg.get("uniport", {})
    n_hvg_common = u_cfg.get("n_hvg_common", 2000)
    n_hvg_specific = u_cfg.get("n_hvg_specific", 2000)
    mode = u_cfg.get("mode", "h")
    lambda_s = u_cfg.get("lambda_s", 1.0)
    source_name = u_cfg.get("source_name", "source")
    use_rep_name = u_cfg.get("use_rep", "X")

    umap_cfg = cfg.get("umap", {})
    umap_min_dist = float(umap_cfg.get("min_dist", 0.3))
    leiden_cfg = cfg.get("leiden", {})
    leiden_res = float(leiden_cfg.get("resolution", 0.02))

    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2

    print("\n[1] Reading raw 10x datasets ...")
    adatas_raw = []
    for path, prefix in datasets:
        ad = read_10x_folder(path, prefix)
        adatas_raw.append(ad)
        print(prefix, ad.shape)

    print("\n[2] Assigning domain_id ...")
    for i, ad in enumerate(adatas_raw):
        ad.obs["domain_id"] = i
        ad.obs["domain_id"] = ad.obs["domain_id"].astype("category")

    print("\n[3] Building adata_cm (common genes, common HVG) ...")
    adata_cm = sc.concat(adatas_raw, join="inner", axis=0, index_unique=None)
    print("adata_cm (before HVG):", adata_cm.shape)
    adata_cm = prep_common_hvg(adata_cm, n_hvg_common)
    print("adata_cm (after common HVG):", adata_cm.shape)

    print("\n[4] Preparing specific HVG for each dataset ...")
    adatas = [prep_specific_hvg(ad.copy(), n_hvg_specific) for ad in adatas_raw]
    print("Specific HVG shapes:", [a.shape for a in adatas])

    print("\n[5] Running uniPort ...")
    n = len(adatas)
    adata_int = up.Run(
        adatas=adatas,
        adata_cm=adata_cm,
        mode=mode,
        lambda_s=lambda_s,
        use_rep=[use_rep_name] * n,
        source_name=source_name,
        outdir=str(out_dir),
    )
    print("adata_int:", adata_int)

    print("\n[6] UMAP / Leiden on latent ...")
    sc.pp.neighbors(adata_int, use_rep="latent")
    sc.tl.umap(adata_int, min_dist=umap_min_dist)
    sc.tl.leiden(adata_int, resolution=leiden_res)
    sc.pl.umap(
        adata_int,
        color=["source", "domain_id", "leiden"],
        wspace=0.4,
        save="_uniport_latent_overview.pdf",
    )

    # Store latent as X_uniport for scib evaluation
    adata_int.obsm["X_uniport"] = np.asarray(adata_int.obsm["latent"])

    adata_int.write_h5ad(out_h5ad)
    print("Saved:", out_h5ad)

if __name__ == "__main__":
    main()
