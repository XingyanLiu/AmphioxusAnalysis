#!/usr/bin/env python
# coding: utf-8

import os
import sys
from pathlib import Path
import logging
import pandas as pd
import numpy as np
# import scanpy as sc

ROOT = Path('../')
sys.path.append(str(ROOT))

# from src import CONSTANTS as C
from src import geneinfo as G

# result directory 
resdir = Path('./grn_res')
figdir = resdir / 'figs'
os.makedirs(figdir, exist_ok=True)


df_gmap = G.load_amzb_gmap_combined()
df_tfmap = G.load_amzb_gmap_tfs()


# ## Option 1. DataFrame based

""" (DataFrame based)
load tf-target edges as a DataFrame
"""
stg_am = 'G3'
stg_zb = '8hpf' # '6hpf'

tag_grn = ['corrSp', 'corrPs', 'SCENIC', ][-1]

edges_am0 = pd.read_csv(resdir / f'{tag_grn}.edge-am-{stg_am}.csv', )#.iloc[:, -3:]
edges_zb0 = pd.read_csv(resdir / f'{tag_grn}.edge-zb-{stg_zb}.csv', )

from src.utils_grn import species_common_edges

common_edges = species_common_edges(edges_am0, edges_zb0, df_gmap, df_tfmap)

common_edges.to_excel(resdir / f'{tag_grn}.common-{stg_am}-{stg_zb}.xlsx', index_label=('TF1', 'TF2'))
common_edges.to_csv(resdir / f'{tag_grn}.common-{stg_am}-{stg_zb}.csv', index_label=('TF1', 'TF2'))
common_edges


# ## Option 2. Module based

"""load tf modules and filter to edge-dfs """

stg_am = 'G3'
stg_zb = ['6hpf', '8hpf'][-1]

tag_grn =  'SCENIC.module'

from src.base import load_pickle
modules_am = load_pickle(resdir / f'{tag_grn}-am-{stg_am}.pkl')
modules_zb = load_pickle(resdir / f'{tag_grn}-zb-{stg_zb}.pkl')

print(len(modules_am), len(modules_zb))


from src.utils_grn import edge_from_modules

contexts = {'top50',
            'top5perTarget',
           }
tag_flt = 'top50top5tgt'
# =====================================
contexts = {
            'weight>90.0%',
            'top5perTarget',
           }
tag_flt = 'wt90top5tgt'
# =====================================
contexts = {
            'weight>75.0%',
            'weight>90.0%',
            'top5perTarget',
           }
tag_flt = 'wt75top5tgt'
# =====================================


edges_am0 = edge_from_modules(modules_am, contexts=contexts, union=True)
edges_zb0 = edge_from_modules(modules_zb, contexts=contexts, union=True)

print(edges_am0.isna().any())
print(edges_zb0.isna().any())

from src.utils_grn import species_common_edges

common_edges = species_common_edges(edges_am0, edges_zb0, df_gmap, df_tfmap)

print(common_edges.isna().any())

# In[13]:

tag_grn_flt = tag_grn.replace('module', tag_flt)
common_edges.to_excel(resdir / f'{tag_grn_flt}.common-{stg_am}-{stg_zb}.xlsx', index_label=('TF1', 'TF2'))
common_edges.to_csv(resdir / f'{tag_grn_flt}.common-{stg_am}-{stg_zb}.csv', index_label=('TF1', 'TF2'))




