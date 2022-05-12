# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 16:01:40 2020

@author: xyliu
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc

import matplotlib as mpl

WORKDIR = Path(r'D:\Users\xyliu\003')
os.chdir(WORKDIR)

from src.pipe_handler import PipeHandler

DATADIR = Path(r'E:\Users\xyliu\data003\amph')
datadir = DATADIR / 'afterQC_formal'

NPC_LIST = [5, 8, 10, 10, 10,
            15, 15, 20, 25, 25,
            25, 25, 30, 30, 30]

Params = dict(qc=True, min_genes=100, max_genes=3000,
              rmv_mito=False,  # mito_perc=0.02,
              counts_level=None,
              plot_hvgs=True, min_mean=0.04, min_disp=0.25,
              do_regress=False, batch_key=None,
              n_comps=50,
              metric='cosine',
              nneigh=20, n_pcs=30,
              min_dist=0.3,
              de=True, plot_de=True,
              cluster=True,
              save_middle=True,
              save_by_default=True
              )
# In[]
''' Stagewise -2 (formal)
'''

qc = False
regress_cnt = False
norm_rev = True
batch_key = 'primer'
nneigh_set = [30] * 6 + [20] * 4 + [10] * 5

resdir = DATADIR / 'res-scanpy' / 'stgwise'

StageNamesAll = ("2cell", "4cell", "8cell", "32cell", "256cell",
                 "B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0" )
StageNames = ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0" )
stages = StageNames

for i, sname in enumerate(stages):

    params = Params.copy()
    params['qc'] = qc
    params['log1first'] = norm_rev
    params['counts_level'] = 300
    params['batch_key'] = batch_key
    params['nneigh'] = nneigh_set[i]
    params['n_pcs'] = NPC_LIST[i]
    if regress_cnt:
        params['do_regress'] = 'n_counts'
    print(params)

    adt = sc.read_h5ad(datadir / f'{sname}_afterQC.h5ad')

    ph = PipeHandler(adt, name=sname, resdir=resdir)
    ph.AutoPipe(**params)

    tt = f'{ph.name} ({ph.adata.shape[0]} cells)'
    ph.vis(color='primer', title=tt, save=f'_primer_{ph.name}.png')
    ph.vis(color='leiden', title=tt, save=f'_leiden_{ph.name}.png')
    ph.save_dims('umap', ftype='csv')
