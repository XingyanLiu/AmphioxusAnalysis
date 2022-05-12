# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2021-12-10
@File: CONSTANTS.py
@Project: AmphioxusDevTree
"""
import os
from pathlib import Path


# CMAP_CLASS = sc.pl.palettes.zeileis_26
CMAP_STAGE0 = 'Spectral'
CMAP_STAGE = 'plasma_r'  # 'Spectral'

# ================== Paths ====================
DATADIR = Path(r'E:/Users/xyliu/data003/amph')
DATADIR_RAW = DATADIR / 'RAW'
FORMALDIR = DATADIR / 'afterQC_formal'

FORMALDIR_NEW = DATADIR / '_formal_data'
FORMALDIR_NEW_ENB = FORMALDIR_NEW / 'embryos'
FORMALDIR_NEW_ADT = FORMALDIR_NEW / 'adults'

LINEDIR = DATADIR / 'lineage'
DATADIR_AD = DATADIR / 'Adult'
GENEDIE = DATADIR / 'genes'

# scATAC-seq related
DATADIR_GA = Path("E:/lxy_pro/003_data/amph/GA")


# ======= compared with zebrafish ======
lin_pairs_zbam = [
    ('Germline', 'Primordial germ cells'),
    ('Endoderm', 'Endoderm'),
    ('Mesoderm', 'Mesoderm'),
    ('Notochord', 'Notochord'),
    ('Epidermal', 'Epithelial ectoderm'),
    ('Hindbrain / Spinal Cord', 'Neural ectoderm'),
    ('Forebrain / Optic', 'Neural ectoderm'),
    ('Midbrain', 'Neural ectoderm'),
    ('Neural Crest', 'Neural ectoderm'),
]
LineageMapZA = dict(lin_pairs_zbam)

# ================== names ====================
# Stages0 = ('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9',
#            'E10', 'E11', 'E12', 'E13', 'E14', 'E15')
usePolyT = ('E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E10',)

Stages = ('E6', 'E7', 'E8', 'E9',
          'E10', 'E11', 'E12', 'E13', 'E14')
StageNames = ("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
StageNamesAll = ("2cell", "4cell", "8cell", "32cell", "256cell",
                 "B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")

Tissues = ('T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7',
           'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14')

Primers = ('random', 'polyT')

StageNameDict = dict(
    E6="B",
    E7="G3",
    E8="G4",
    E9="G5",
    E10="G6",
    E11="N0",
    E12="N1",
    E13="N3",
    E14="L0"
)
TissueNameDict = {
    'T1': 'Tentacles',
    'T2': 'Neural Tube',
    'T3': 'Male Gonad',
    'T4': 'Endostyle',
    'T5': 'Epidermis',
    'T6': 'Hepatic diverticulum',
    'T7': 'Fore End',
    'T8': 'Hind End',
    'T9': 'Gill Branchia',
    'T10': 'End Gut',
    'T11': 'Rostral Side',
    'T12': 'Female Gonad',
    'T13': 'Muscle',
    'T14': 'Notochord'
}
TissueNamesUse = [
    'Female Gonad',
    'Hind End',
    'End Gut',
    'Fore End',
    'Rostral Side',
    'Gill Branchia',
    'Male Gonad',
    'Notochord',
    'Epidermis',
    'Neural Tube',
    'Tentacles',
    'Endostyle',
    # 'Muscle',
    # 'Hepatic diverticulum'
]

LineageOrdAll = (
    "B_0",
    "B_1",
    "B_2",
    "Primordial germ cells",
    "Epithelial ectoderm",
    "Neural ectoderm",
    "Notochord",
    "Mesoderm",
    "Unassigned",
    "Endoderm",
)

LineageOrd = LineageOrdAll[-7:]
LineageAnno = dict(
    E7_0="Epithelial ectoderm",
    E7_1="Endoderm",
    E7_2="Mesoderm",
    E7_3="Notochord",
    E7_4="Neural ectoderm",
    E7_5="Unassigned",
    E7_6="Primordial germ cells"
)

LineageAnnoAll = dict(
    E6_0='B_0',
    E6_1='B_1',
    E6_2='B_2',
    E6_3="Primordial germ cells",
    E7_0="Epithelial ectoderm",
    E7_1="Endoderm",
    E7_2="Mesoderm",
    E7_3="Notochord",
    E7_4="Neural ectoderm",
    E7_5="Unassigned",
    E7_6="Primordial germ cells"
)

LineageColorMap = {
    "B_0": "#596e79",
    "B_1": "#b3b3b3",
    "B_2": "#c7b198",
    "B_3": "#40bad5",
    "Primordial germ cells": "#40bad5",
    "Epithelial ectoderm": "#984ea3",
    "Neural ectoderm": "#36622b",
    "Notochord": "#035aa6",
    "Mesoderm": "#fcbf1e",
    "Unassigned": "#af0404",
    "Endoderm": "#dd7631"
}

MarkerDictShow = {
    'bf_00013097': 'DNMT1',
    'bf_00019348': 'Hu_Elav',  # neural ectoderm
    'bf_00016583': 'Tbx2',
    'bf_00004829': 'Nanos',
    'bf_00008183': 'LRIG3',
    'bf_00002546': 'Sox17',  # endoderm
    'bf_00006571': 'Wnt8',  # Mesoderm
    'bf_00016125': 'FoxJ1',  # Epithelial ectoderm
    'bf_00023325': 'Keratin',  # Epithelial ectoderm
    'bf_00020985': 'CHRD',  # notochord
    'bf_00012984': 'Cdx',
}




