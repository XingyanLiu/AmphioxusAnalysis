# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2022-02-15
@File: utils_grn.py
@Project: AmphioxusDevTree
"""
import os
from pathlib import Path
from typing import Union, Optional, Sequence, Mapping
import time
import logging
import numpy as np
import pandas as pd


# import preprocess as pp


def filter_edges_by_homo(
        edge_df1: pd.DataFrame, edge_df2: pd.DataFrame,
        gmap: Union[pd.DataFrame],
        tfmap: Union[pd.DataFrame] = None,
):
    if tfmap is None:
        tfmap = gmap

    def _filter(edges, candi_tfs, candi_genes):
        print(f'n_edges before: {edges.shape[0]}')
        edges_flt = subset_matches(edges, candi_tfs, candi_genes)
        print(f'n_edges after: {edges_flt.shape[0]}')
        return edges_flt

    edges1 = _filter(
        edge_df1,
        candi_tfs=tfmap.iloc[:, 0].tolist(),
        candi_genes=gmap.iloc[:, 0].tolist())
    edges2 = _filter(
        edge_df2,
        candi_tfs=tfmap.iloc[:, 1].tolist(),
        candi_genes=gmap.iloc[:, 1].tolist())

    tfmap_used = subset_matches(
        tfmap, edges1.iloc[:, 0].unique(),
        edges2.iloc[:, 0].unique(), )
    print(f'used tf-homo-pairs: {tfmap_used.shape[0]}')
    return edges1, edges2, tfmap_used


def species_common_edges(
        edge_df1: pd.DataFrame, edge_df2: pd.DataFrame,
        gmap: Union[pd.DataFrame],
        tfmap: Union[pd.DataFrame] = None,
):
    """

    :param edge_df1: pd.DataFrame with 3 columns ['source', 'target', 'weight']
    :param edge_df2: pd.DataFrame with 3 columns ['source', 'target', 'weight']
    :param gmap: pd.DataFrame with 2 columns, each row represents a homologous
        gene pair.
    :param tfmap: pd.DataFrame with 2 columns, each row represents a homologous
        gene pair (restricted to transcription factors).
    :return:
    DataFrame (storing intersected TF-target pairs),
        with multi-index ('source1', 'source2'),
        and columns ['target1', 'target2', 'weight1', 'weight2']
    """
    if tfmap is None:
        tfmap = gmap

    edges1, edges2, tfmap_used = filter_edges_by_homo(
        edge_df1, edge_df2, gmap, tfmap
    )

    def _to_multi_index_srs(df):
        df = df.set_index(list(df.columns[: 2])).iloc[:, 0]
        return df

    # for each pair of homo-TFs, check their targets
    edges1 = _to_multi_index_srs(edges1)  # 确保是 series 而不是 df
    edges2 = _to_multi_index_srs(edges2)

    tf_pairs = list(map(tuple, tfmap_used.iloc[:, [0, 1]].values))
    # gmap = df_gmap.iloc[:, :2]

    tgt_columns = ['target1', 'target2']
    col1, col2 = tgt_columns

    record = {}
    for tf1, tf2 in tf_pairs:
        weights1 = edges1.loc[tf1]  # series
        targets1 = [x for x in weights1.index if x != tf1]  # weights1.index

        weights2 = edges2.loc[tf2]  # series
        targets2 = [x for x in weights2.index if x != tf2]  # weights2.index

        gmap_tgt = subset_matches(
            gmap, targets1, targets2, union=False).copy()  # df
        gmap_tgt.columns = tgt_columns
        gmap_tgt['weight1'] = gmap_tgt[col1].apply(
            lambda x: weights1.get(x, np.nan))  # 按理说不会有 NaN
        gmap_tgt['weight2'] = gmap_tgt[col2].apply(
            lambda x: weights2.get(x, np.nan))  # 按理说不会有 NaN

        # 保留相关性方向一致的 target
        _same_sign = (gmap_tgt['weight1'] * gmap_tgt['weight2']) > 0.

        # print(_same_sign.sum())
        record[(tf1, tf2)] = gmap_tgt[_same_sign]

    record = pd.concat(record, )
    record.index = record.index.droplevel(2)
    return record


def edge_from_modules(
        modules,
        contexts=None,  # used to filter unused modules
        union=False,
):
    """

    :param modules: a list of ctxcore.genesig.Regulon (pySECNIC), output by
        grnboost2
    :param contexts: set
    :return: pd.DataFrame

    :examples
    >>> adjacencies = grnboost2(ex_matrix, gene_names=gene_name,
    ...                 tf_names=tf_used, verbose=True) #耗时
    >>> modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    >>> edge_from_modules(modules, contexts={'weight>75.0%', 'top10perTarget'})
    """
    if contexts is None:
        contexts = {'top50',
                    'weight>50.0%', 'weight>75.0%',  # 'weight>90.0%',
                    'top50perTarget', 'top10perTarget'}
    record = {}
    for m in modules:
        if len(m.context.intersection(contexts)) <= 0:
            continue
        _tf = m.transcription_factor
        if _tf is None:
            continue
        if _tf in record.keys():
            if union:
                record[_tf] = record[_tf].union(m)
            else:
                record[_tf] = record[_tf].intersection(m)
        else:
            record[_tf] = m

    for _tf in record.keys():
        record[_tf] = pd.Series(record[_tf].gene2weight).to_frame()
    df = pd.concat(record, ).reset_index()
    df.columns = ['source', 'target', 'weight']

    return df


def subset_matches(df_match: pd.DataFrame,
                   left: Sequence,
                   right: Sequence,
                   union: bool = False,
                   cols: Union[None, Sequence[str]] = None,
                   indicators=False) -> pd.DataFrame:
    """ Take a subset of token matches (e.g., gene homologies)
    (copies from 'preprocess.py')
    Parameters
    ----------
    df_match: pd.DataFrame
        a dataframe with at least 2 columns
    left:
        list-like, for extracting elements in the first column.
    right:
        list-like, for extracting elements in the second column.
    union:
        whether to take union of both sides of genes
    cols:
        if specified, should be 2 column names for left and right genes,
        respectively.
    indicators:
        if True, only return the indicators.
    """
    if cols is None:
        cols = df_match.columns[: 2]

    c1, c2 = cols
    tmp = df_match[c1].isin(left).to_frame(c1)
    tmp[c2] = df_match[c2].isin(right)

    keep = tmp.max(1) if union else tmp.min(1)

    return keep if indicators else df_match[keep]


