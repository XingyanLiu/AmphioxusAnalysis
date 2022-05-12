#!/usr/bin/env python
# coding: utf-8
"""GRN inference for amphioxus and zebrafish"""
import os
import sys
from pathlib import Path
import time
import logging
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
ROOT = Path('../')
sys.path.append(str(ROOT))

from src import preproc as pp
# from src import CONSTANTS as C
from src import geneinfo as G

import scanpy as sc
sc.set_figure_params(fontsize=13)


def save_pickle(obj, fpath):
    """ save the object into a .pickle file
    """
    with open(fpath, 'wb') as f:
        pickle.dump(obj, f)
    print('object saved into:\n\t', fpath)


def load_pickle(fp):
    """ load the object from a .pickle file

    Examples
    --------
    >>> load_pickle('dpair.pickle')
    """
    with open(fp, 'rb') as f:
        res = pickle.load(f)
    return res


# Load TFs
tfs_zb = G.load_tfs('zebrafish', onlylist=True)

df_tf = G.load_amzb_gmap_tfs()
tfs_am = set(df_tf['QueryID'])
# tfs_zb = set(df_tf['Symbol'])

print('amphTFs: {}\namph2zebrTFs: {}\nzebrTFs: {}'.format(len(tfs_am), len(set(df_tf['Symbol'])), len(tfs_zb)))


resdir = Path('./grn_res')
figdir = resdir / 'figs'
os.makedirs(figdir, exist_ok=True)


# # Main for SCENIC
# 
# ## Load adata

"""amphioxus"""
tag_sp = 'am'
stg = ["G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0"][4]
fn = f'/Users/xingyan/Downloads/BoiData/003-amph/amph-emb/raw-{stg}.h5ad'
tfset = set(tfs_am)
tf_names = sorted(tfset)

"""zebrafish"""
# tag_sp = 'zb'
# stg = ['4hpf', '6hpf', '8hpf', '10hpf', '14hpf', '18hpf', '24hpf'][1]
# fn = f'/Users/xingyan/Downloads/BoiData/003-amph/zebrafish/raw-{stg}.h5ad'
# tfset = set(tfs_zb)
# tf_names = sorted(tfset)

"""load data"""
adata_raw = sc.read(fn)
adata_raw


adata = adata_raw.copy()

sc.pp.filter_genes(adata, min_cells=1)
print(adata.shape)

pp.normalize_default(adata, target_sum=None, copy=False)
sc.tl.rank_genes_groups(adata, groupby='lineage', pts=True)


deg_info = pp.get_marker_info_table(adata, cut_padj=0.05, cut_pts=0.01)
deg_uniq = deg_info['names'].unique()

print(deg_info['group'].value_counts())
print('total unique DEGs:', len(deg_uniq))

deg_info.to_csv(resdir / f'DEGinfo-{tag_sp}-{stg}.csv', index=False)

genes_used = sorted(tfset.union(deg_uniq).intersection(adata.var_names))
tf_used = sorted(tfset.intersection(deg_uniq))
# tf_used = sorted(tfset.intersection(adata.var_names))

print(f'unsed genes: {len(genes_used)}')
print(f'unsed TFs: {len(tf_used)}')


""" 1: Inference of co-expression modules
"""
subadt = adata_raw[:, genes_used]

ex_matrix = subadt.to_df()
gene_name = subadt.var_names

# 两条命令解决
t0 = time.time()
adjacencies = grnboost2(ex_matrix, gene_names=gene_name, 
                        tf_names=tf_used, verbose=True) #耗时
t1 = time.time() 
print(t1 - t0)

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
print(time.time() - t1)

adjacencies.to_csv(resdir / f'SCENIC.edge-{tag_sp}-{stg}.csv', index=False)
print(adjacencies.shape)
adjacencies.head()


save_pickle(modules, resdir / f'SCENIC.module-{tag_sp}-{stg}.pkl')


# =============================================================================
"""zebrafish"""
tag_sp = 'zb'
stg = ['4hpf', '6hpf', '8hpf', '10hpf', '14hpf', '18hpf', '24hpf'][1]
fn = f'/Users/xingyan/Downloads/BoiData/003-amph/zebrafish/raw-{stg}.h5ad'
tfset = set(tfs_zb)
tf_names = sorted(tfset)


"""load data"""
adata_raw = sc.read(fn)


def sep_zb_notochord(df, inplace=True):
    """separate Notochord cluster"""
    if inplace:
        print('changed inplace!')
    else:
        df = df.copy()
    df['lineage'] = df['TissueName'].astype(str)
    is_notochord = df['ClusterName'].apply(lambda x: 'noto' in x)
    df.loc[is_notochord, 'lineage'] = 'Notochord'
    return None if inplace else df


if tag_sp == 'zb':
    sep_zb_notochord(adata_raw.obs)
    print(adata_raw.obs['lineage'].value_counts())
    
""" process """
adata = adata_raw.copy()

sc.pp.filter_genes(adata, min_cells=1)
print(adata.shape)

pp.normalize_default(adata, target_sum=None, copy=False)
sc.tl.rank_genes_groups(adata, groupby='lineage', pts=True)

"""DEGs"""

deg_info = pp.get_marker_info_table(adata, cut_padj=0.05, cut_pts=0.005)
deg_uniq = deg_info['names'].unique()

print(deg_info['group'].value_counts())
print('total unique DEGs:', len(deg_uniq))

deg_info.to_csv(resdir / f'DEGinfo-{tag_sp}-{stg}.csv', index=False)

genes_used = sorted(tfset.union(deg_uniq).intersection(adata.var_names))
tf_used = sorted(tfset.intersection(deg_uniq))
# tf_used = sorted(tfset.intersection(adata.var_names))

print(f'unsed genes: {len(genes_used)}')
print(f'unsed TFs: {len(tf_used)}')


""" 1: Inference of co-expression modules
"""
subadt = adata_raw[:, genes_used]

ex_matrix = subadt.to_df()
gene_name = subadt.var_names

t0 = time.time()
adjacencies = grnboost2(ex_matrix, gene_names=gene_name, 
                        tf_names=tf_used, verbose=True) #耗时
t1 = time.time() 
print(t1 - t0)

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
print(time.time() - t1)

""" Save results """
adjacencies.to_csv(resdir / f'SCENIC.edge-{tag_sp}-{stg}.csv', index=False)
print("adjacencies.shape:", adjacencies.shape)

save_pickle(modules, resdir / f'SCENIC.module-{tag_sp}-{stg}.pkl')





