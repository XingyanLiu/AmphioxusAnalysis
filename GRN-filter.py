#!/usr/bin/env python
# coding: utf-8
"""Filter TF-target pairs based on their DE lineages"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('../')
sys.path.append(str(ROOT))
from src import geneinfo as G

df_gmap = G.load_amzb_gmap_combined()
df_tfmap = G.load_amzb_gmap_tfs()

# result directory
resdir = Path('./grn_res')
figdir = resdir / 'figs'
os.makedirs(figdir, exist_ok=True)

stg_am = 'G3'
stg_zb = ['6hpf', '8hpf'][1]

tag_grn = ['SCENIC', 'SCENIC.module',
           'SCENIC.wt75top5tgt', 'SCENIC.wt90top5tgt',
           'SCENIC.top50top5tgt'][2]

common_edges = pd.read_csv(resdir / f'{tag_grn}.common-{stg_am}-{stg_zb}.csv')
print(common_edges.shape)
common_edges.head()

# DEGs
deginfo_am = pd.read_csv(resdir / f'DEGinfo-am-{stg_am}.csv')
deginfo_zb = pd.read_csv(resdir / f'DEGinfo-zb-{stg_zb}.csv')

deginfo_zb.head()


def make_name2groups(deginfo, sep=' & '):
    return deginfo.groupby('names')['group'].apply(
        lambda x: sep.join(map(str, x))).to_dict()


def fill_de_lineage(common_edges, deginfo1, deginfo2, prefix='DE_group_'):
    """ Supplement the source of differential expression of each gene """
    source_lins_am = make_name2groups(deginfo1)
    source_lins_zb = make_name2groups(deginfo2)

    src_lin1 = common_edges[['TF1', 'target1']].applymap(
        lambda x: source_lins_am.get(x, ''))
    src_lin1.columns = prefix + src_lin1.columns
    src_lin2 = common_edges[['TF2', 'target2']].applymap(
        lambda x: source_lins_zb.get(x, ''))
    src_lin2.columns = prefix + src_lin2.columns

    return pd.concat([common_edges, src_lin1, src_lin2], axis=1)


def common_source_group(s1, s2, sep=' & '):
    set1 = s1.split(sep)
    set2 = s2.split(sep)
    inter = set(set1).intersection(set2)
    return sep.join(inter)


df_all = fill_de_lineage(common_edges, deginfo_am, deginfo_zb, )
df_all.head()

# save results
_fp = resdir / f'detail-{tag_grn}.common-{stg_am}-{stg_zb}.xlsx'
_fp_csv = str(_fp).replace('.xlsx', '.csv')
print(_fp)
print(_fp_csv)

df_all.set_index(['TF1', 'TF2']).to_excel(_fp)
df_all.set_index(['TF1', 'TF2']).to_csv(_fp_csv)

# TF and target genes should come from the same lineage
cmm1 = df_all[['DE_group_TF1', 'DE_group_target1']].apply(
    lambda x: common_source_group(*x), axis=1)
cmm2 = df_all[['DE_group_TF2', 'DE_group_target2']].apply(
    lambda x: common_source_group(*x), axis=1)

df_flt = df_all[(cmm1 != '') & (cmm2 != '')]

# save results
_fp = resdir / f'flt-{tag_grn}.common-{stg_am}-{stg_zb}.xlsx'
_fp_csv = str(_fp).replace('.xlsx', '.csv')
print(_fp)
print(_fp_csv)

df_flt.set_index(['TF1', 'TF2']).to_excel(_fp)
df_flt.set_index(['TF1', 'TF2']).to_csv(_fp_csv)
