# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2022-02-09
@File: geneinfo.py
@Project: AmphioxusDevTree
"""
import json
import os
import sys
from pathlib import Path
from typing import Union, Optional, Sequence, Mapping
import time
import logging
import numpy as np
import pandas as pd
ROOT = Path(__file__).parents[1]
DIR_RESOURCE = ROOT / 'resources'
DIR_GMAP = DIR_RESOURCE / 'gene_mappings'


def load_amph_annos(**kwds):
    annos = pd.read_csv(DIR_RESOURCE / 'Annos1.csv', **kwds)
    annos.index = annos['ID'].values
    print(annos.columns)
    return annos


def load_amph_gmap(sp='zebrafish', with_id=False):
    """ 对于其他物种, 相当于直接读取表格
    >>> # amphioxus-zebrafish TF-homologies
    >>> df_tfmap = load_amph_gmap('zebrTF6k')
    """
    df_gmap = pd.read_csv(DIR_GMAP / f'bf__v__{sp}.csv')
    if sp == 'zebrafish':
        if not with_id:
            df_gmap = df_gmap.iloc[:, [0, 2]]
        externals = pd.DataFrame([
            ['bf_00002546', 'sox17'],
            ['bf_00001714', 'sox2'],
            ['bf_00001715', 'sox2'],
        ], columns=['amph_id', f'{sp}_name'])
        df_gmap = df_gmap.append(externals, ignore_index=True)
    return df_gmap


def load_amzb_gmap_tfs():
    """文昌鱼与斑马鱼之间潜在的同源TF映射表"""
    df_tfmap = load_amph_gmap('zebrTF6k')
    print('load candidate TF-homologies between amphioxus and zebrafish')
    print(f'\tShape: {df_tfmap.shape}')
    print(f'\tColumns: {df_tfmap.columns}')
    return df_tfmap


def load_amzb_gmap_combined():
    print('load combined gene-mappings from different sources')
    df_tfmap = load_amph_gmap('zebrTF6k')
    df_tfmap = df_tfmap.iloc[:, [0, 1]]
    df_gmap = pd.read_csv(DIR_GMAP / f'bf__v__zebrafish.csv').iloc[:, [0, 2]]
    # 合并两个来源的映射
    df_tfmap.columns = df_gmap.columns
    df_gmap = pd.concat([df_gmap, df_tfmap][::-1], ignore_index=True)
    df_gmap.drop_duplicates(inplace=True)
    print(f'\tShape: {df_gmap.shape}')
    print(f'\tColumns: {df_gmap.columns}')
    return df_gmap


def load_gmap_details(sp1, sp2, dir_gmap=None):
    """ mart_imputed """
    if dir_gmap is None:
        dir_gmap = Path('E:/lxy_pro/004/resources/mart_exports/imputed')
    df = pd.read_csv(dir_gmap / f'mart_imputed-{sp1}2{sp2}.csv')
    return df


def load_gmap(sp1, sp2, dir_gmap=None):
    if dir_gmap is None:
        dir_gmap = Path('E:/lxy_pro/004/resources/mart_exports/exported_gene_matches')
    df = pd.read_csv(dir_gmap / f'gene_matches_{sp1}2{sp2}.csv')
    return df


def load_gene_info(
        sp, dir_info=None,
        as_id2name=False,
        key_id='Gene stable ID', key_sym='Gene name',
):
    if dir_info is None:
        dir_info = Path('E:/lxy_pro/004/resources/mart_exports/gene_info')
    df = pd.read_csv(dir_info / f'gene_info-{sp}.csv')
    print(df.columns)
    if as_id2name:
        return dict(df[[key_id, key_sym]].values)
    return df


def load_tfs(species='human', onlylist=False):
    sp_dict = {
        'human': 'Homo_sapiens',
        'mouse': 'Mus_musculus',
        'chick': 'Gallus_gallus',
        'zebrafish': 'Danio_rerio',
        'pig': 'Sus_scrofa'
    }
    dir_tf = ROOT / 'resources/TF'
    sp_name = sp_dict.get(species, species)
    df_tf = pd.read_csv(dir_tf / f'{sp_name}_TF.txt', sep='\t')
    if onlylist:
        col = ['Ensembl', 'Symbol'][-1]
        return df_tf[col].tolist()
    return df_tf


def __test__():
    df = load_amph_gmap('mouse')
    df = load_gmap_details('chick', 'zebrafish')
    df = load_gmap('chick', 'zebrafish')
    df = load_tfs(onlylist=True)
    logging.info(df)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(filename)s-%(lineno)d-%(funcName)s(): '
               '%(levelname)s\n %(message)s')
    t = time.time()

    __test__()

    print('Done running file: {}\nTime: {}'.format(
        os.path.abspath(__file__), time.time() - t,
    ))
