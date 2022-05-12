# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2022-05-12
@File: _zscore.py
@Project: AmphioxusDevTree
"""
import os
import sys
from pathlib import Path
from typing import Union, Optional, Sequence, Mapping
import time
import logging
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
from sklearn.preprocessing import StandardScaler


def zscore(X, with_mean=True, scale=True, ):
    """ For each column of X, do centering (z-scoring)
    """
    # code borrowed from `scanpy.pp._simple`
    scaler = StandardScaler(with_mean=with_mean, copy=True).partial_fit(X)
    if scale:
        # user R convention (unbiased estimator)
        e_adjust = np.sqrt(X.shape[0] / (X.shape[0] - 1))
        scaler.scale_ *= e_adjust
    else:
        scaler.scale_ = np.array([1] * X.shape[1])
    X_new = scaler.transform(X)
    if isinstance(X, pd.DataFrame):
        X_new = pd.DataFrame(X_new, index=X.index, columns=X.columns)
    return X_new


def group_zscore(X: Union[np.ndarray, pd.DataFrame],
                 labels: Union[Sequence, np.ndarray],
                 with_mean: bool = True,
                 scale: bool = True,
                 max_value: float = None):
    """
    For each column of X, do within-group centering (z-scoring)

    Parameters
    ----------
    X: np.ndarray or pd.DataFrame
        A matrix of shape (n_samples, n_features), each row of X is an
        observation, wile each column is a feature
    labels: np.ndarray
        the group labels, of shape (n_samples,)
    with_mean: boolean, True by default
        If True, center the data before scaling, and X shoud be a dense matrix.
    scale: bool
        whether to scale with standard deviation
    max_value: float
        if given, the absolute values of the result matrix will be
        clipped at this value.

    Returns
    -------
    the scaled data matrix
    """
    isdf = False
    if isinstance(X, pd.DataFrame):
        isdf = True
        index, columns, X = X.index, X.columns, X.values
    X = X.astype(np.float).copy()
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    for lb in unique_labels:
        ind = labels == lb
        if sum(ind) == 1:
            logging.warning(f'ignoring class {lb} with only one sample.')
            continue
        X[ind, :] = zscore(X[ind, :], with_mean=with_mean, scale=scale)

    if max_value is not None:
        X[X > max_value] = max_value
        logging.info(f'... clipping at max_value {max_value}')

    if isdf:
        X = pd.DataFrame(X, index=index, columns=columns)
    return X


def group_zscore_adata(adt: sc.AnnData,
                       groupby: str = 'batch',
                       key: str = 'counts',
                       key_new: str = None,
                       max_value: float = None,
                       with_mean: bool = True,
                       cover: bool = True,
                       **kwds):
    """Calculate z-scores for each group of observations in an ``AnnData`` object

    Parameters
    ----------
    adt: AnnData
    groupby: str
        A key from adt.obs, from which the labels are take
    key: str, {'X_pca', 'count'}
        can be a key from adt.obsm, e.g. `key='X_pca'`
        If key == 'counts', then do scaling on `adt.X`
        and cover the old count matrix, ignoring the `cover` parameter
    key_new: str
        used when ``key != 'counts'``
    with_mean: boolean, True by default
        If True, center the data before scaling, and X shoud be a dense matrix.
    max_value: float
        if given, the absolute values of the result matrix will be
        clipped at this value.
    cover: bool
        whether to cover the old X with the scored X
    """
    labels = adt.obs[groupby]
    if key == 'counts':
        logging.info(
            'Z-score scaling on count matrix, transformed into a dense array')
        if sparse.issparse(adt.X):
            X = adt.X.toarray()
        else:
            X = adt.X
        if not cover: adt = adt.copy()
        adt.X = group_zscore(
            X, labels, with_mean=with_mean,
            max_value=max_value, **kwds)
    else:
        if cover:
            key_new = key
        else:
            key_new = key + '_new' if key_new is None else key_new
        adt.obsm[key_new] = group_zscore(
            adt.obsm[key], labels,
            with_mean=with_mean, max_value=None, **kwds)

    return adt if not cover else None


def wrapper_scale(adata, zero_center=True, max_value=None,
                  groupby=None, copy=False, **kwds):
    """
    Wrapper function for centering and scaling data matrix `X` in sc.AnnData,
    extended for within-batch processing.

    Examples
    --------
    >>> wrapper_scale(adata, groupby='batch')
    """
    if groupby is not None:
        logging.info(f'doing within-group scaling, group by [ {groupby} ]')
        return group_zscore_adata(adata,
                                  max_value=max_value,
                                  groupby=groupby,
                                  key='counts',
                                  with_mean=zero_center,
                                  cover=not copy,
                                  **kwds)
    else:
        logging.info('using the build-in function `sc.pp.scale(..)`')
        return sc.pp.scale(adata, zero_center=zero_center,
                           max_value=max_value, copy=copy)


