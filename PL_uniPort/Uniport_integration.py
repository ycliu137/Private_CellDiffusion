
import os
import sys
import scanpy as sc
import numpy as np
import uniport as up

DATASETS = [
    ("/home/ziyinghuang/project/survey/dataset/NM/63", "63"),
    ("/home/ziyinghuang/project/survey/dataset/NM/19", "19"),
    ("/home/ziyinghuang/project/survey/dataset/NM/39", "39"),
    ("/home/ziyinghuang/project/survey/dataset/NM/29", "29"),
]

OUTDIR = "/home/ziyinghuang/project/diffusion/uniport_nm4_out"
FIGDIR = os.path.join(OUTDIR, "figs")

N_HVG_COMMON = 2000
N_HVG_SPECIFIC = 2000

LEIDEN_RES = 0.02
UMAP_MIN_DIST = 0.3

def read_10x_folder(path: str, prefix: str):

    ad = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        make_unique=True
    )

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

#################################################################################

sc.settings.figdir = FIGDIR
sc.settings.verbosity = 2

print("\n[1] Reading raw 10x datasets ...")
adatas_raw = []
for path, prefix in DATASETS:
    ad = read_10x_folder(path, prefix)
    adatas_raw.append(ad)
    print(prefix, ad.shape)


print("\n[2] Assigning domain_id ...")
for i, ad in enumerate(adatas_raw):
    ad.obs["domain_id"] = i
    ad.obs["domain_id"] = ad.obs["domain_id"].astype("category")

print("\n[3] Building adata_cm (common genes, common HVG) ...")
adata_cm = sc.concat(
    adatas_raw,
    join="inner",
    axis=0,
    index_unique=None  # 不改写 obs_names
)

print("adata_cm (before HVG):", adata_cm.shape)
adata_cm = prep_common_hvg(adata_cm, N_HVG_COMMON)
print("adata_cm (after common HVG):", adata_cm.shape)


print("\n[4] Preparing specific HVG for each dataset ...")
adatas = []
for ad in adatas_raw:
    ad_hvg = prep_specific_hvg(ad.copy(), N_HVG_SPECIFIC)
    adatas.append(ad_hvg)
print("Specific HVG shapes:", [a.shape for a in adatas])

print("\n[5] Running uniPort (mode='h') ...")
n = len(adatas)
adata_int = up.Run(
    adatas=adatas,
    adata_cm=adata_cm,
    mode="h",
    lambda_s=1.0,
    use_rep=["X"] * n,
    source_name="source",
    outdir=OUTDIR
)

print("adata_int:", adata_int)

print("\n[6] UMAP/Leiden on latent ...")
sc.pp.neighbors(adata_int, use_rep="latent")
sc.tl.umap(adata_int, min_dist=UMAP_MIN_DIST)
sc.tl.leiden(adata_int, resolution=LEIDEN_RES)
sc.pl.umap(adata_int, color=["source", "domain_id", "leiden"], wspace=0.4, save="_uniport_latent_overview.pdf")

adata_int_path = os.path.join(OUTDIR, "uniport_integrated_hvg.h5ad")
adata_int.write_h5ad(adata_int_path)
print("Saved:", adata_int_path)