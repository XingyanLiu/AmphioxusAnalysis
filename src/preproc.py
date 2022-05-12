# -*- coding: UTF-8 -*-
"""
@Author: Xingyan Liu
@CreateDate: 2022-05-12
@File: _utils.py
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
import scanpy as sc
from scipy import sparse, io
from ._zscore import *


def check_dirs(path):
    if os.path.exists(path):
        print('already exists:\n\t%s' % path)
    else:
        os.makedirs(path)
        print('a new directory made:\n\t%s' % path)


def change_names(names: Sequence,
                 foo_change,  #: Callable,
                 **kwargs):
    return list(map(foo_change, names, **kwargs))


def save_named_mtx(adata, dirname, field=None, raw=True, backup_npz=True,
                   **kwds):
    check_dirs(dirname)
    """ better set field='integer' if possible """
    if adata.raw is not None:
        adata = adata.raw

    mtx_file = '%s/matrix.mtx' % (dirname)
    bcd_file = '%s/barcodes.tsv' % (dirname)
    gene_file = '%s/genes.tsv' % (dirname)

    genes = adata.var_names.to_frame(index=False, name='genes')
    barcodes = adata.obs_names.to_frame(index=False, name='barcodes')
    logging.info(adata.X.data[:5])
    genes.to_csv(gene_file, index=False, header=False)
    barcodes.to_csv(bcd_file, index=False, header=False)
    if backup_npz:
        sparse.save_npz(f'{dirname}/matrix.npz', adata.X.T)
    io.mmwrite(mtx_file, adata.X.T, field=field, **kwds)

    print('Matrix saved to directory `%s`' % dirname)



def normalize_default(adata: sc.AnnData,
                      target_sum=None,
                      copy: bool = False,
                      log_only: bool = False,
                      force_return: bool = False, ):
    """ Normalizing datasets with default settings (total-counts normalization
    followed by log(x+1) transform).

    Parameters
    ----------
    adata
        ``AnnData`` object
    target_sum
        scale factor of total-count normalization
    copy
        whether to copy the dataset
    log_only
        whether to skip the "total-counts normalization" and only perform
        log(x+1) transform
    force_return
        whether to return the data, even if changes are made for the
        original object

    Returns
    -------
    ``AnnData`` or None

    """
    if copy:
        adata = adata.copy()
        logging.info('A copy of AnnData made!')
    else:
        logging.info('No copy was made, the input AnnData will be changed!')
    logging.info('normalizing datasets with default settings.')
    if not log_only:
        logging.info(f'performing total-sum normalization, target_sum={target_sum}...')
        sc.pp.normalize_total(adata, target_sum=target_sum)
    else:
        logging.info('skipping total-sum normalization')
    sc.pp.log1p(adata)
    return adata if copy or force_return else None


def normalize_log_then_total(
        adata, target_sum=None,
        copy=False, force_return=False):
    """ For SplitSeq data, performing log(x+1) BEFORE total-sum normalization
        will results a better UMAP visualization (e.g. clusters would be less
        confounded by different total-counts ).
    """
    if copy:
        adata = adata.copy()
    logging.info(
        'normalizing datasets with default settings (log1p first).'
        f'target_sum = {target_sum}')
    sc.pp.log1p(adata)  # x <- log(x + 1)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    return adata if copy or force_return else None


def top_markers_from_df(marker_df, n=5, groups=None, unique=True, ):
    """
    Parameters
    ----------
    marker_df:
        a data-frame with cluster names as columns, and genes as values
    n: int
    groups:
        a list of cluster names (column names)
    unique: bool
        whether to flatten the results into a unique gene-list, default is True.

    Returns
    -------
    A dataframe (``union=False``) or a flattened marker list (``union=True``)
    """
    groups = marker_df.columns if groups is None else groups
    top = marker_df[groups].iloc[: n]
    #    top = marker_df[groups].iloc[: n].values.T.flatten()
    if unique:
        top = top.values.T.flatten()
        print(f'{len(top)} genes before taking unique')
        top = pd.unique(top)
        print(f'taking total of {len(top)} unique differential expressed genes')
    return top


def top_markers_from_adata(adata: sc.AnnData,
                           n=5, groups=None,
                           unique=True,
                           cut_padj=0.05,
                           key='rank_genes_groups'):
    df_info = get_marker_info_table(adata, groups, key=key, cut_padj=cut_padj)

    if n is None:
        if unique:
            return df_info['names'].unique().tolist()
        return df_info['names'].tolist()
    else:
        names = df_info.groupby('group').apply(lambda x: x.head(n))['names']
        return names.unique().tolist() if unique else names.tolist()


def get_marker_info_table(
        adata, groups=None, key='rank_genes_groups',
        cut_padj: Optional[float] = 0.05,
        cut_logfc: Optional[float] = 0.25,
        cut_pts: Optional[float] = None,
):
    result = adata.uns[key]
    if groups is None:
        groups = result['names'].dtype.names

    dfs = []
    cols = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores', ]
    # cols = [c for c in cols if c in result.keys()]
    flag_pts = 'pts' in result.keys()

    for group in groups:
        _df = pd.DataFrame({
            key: result[key][group] for key in cols
        })
        _df['group'] = group

        if flag_pts:
            # expression proportions, avoid row mismatching
            _df['pts'] = _df['names'].map(result['pts'][group])
            _df['pts_rest'] = _df['names'].map(result['pts_rest'][group])
            if cut_pts is not None:
                _df = _df[_df['pts'] >= cut_pts]
        if cut_padj is not None:
            _df = _df[_df['pvals_adj'] <= cut_padj]
        if cut_logfc is not None:
            _df = _df[_df['logfoldchanges'] >= cut_logfc]
        dfs.append(_df.copy())  # [['group'] + cols])
    df = pd.concat(dfs, axis=0, keys=groups)
    if flag_pts:
        cols += ['pts', 'pts_rest']
    return df[['group'] + cols]


def get_marker_name_table(adata, key='rank_genes_groups'):
    return pd.DataFrame(adata.uns[key]['names'])


def take_group_labels(labels: Sequence, group_names: Sequence,
                      indicate=False, remove=False):
    """
    Parameters
    ----------
    labels: list-like
    group_names:
        names of groups that you want to take out
    indicate: bool
        if True, return a Series of bool indicators of the groups
        else, return the labels.
    remove: bool
        False by default, set as True if you want to keep groups that
        NOT in the `given group_names`
    """
    if isinstance(group_names, (str, int)):
        group_names = [group_names]
    if remove:
        indicators = change_names(labels, lambda x: x not in group_names)
    else:
        indicators = change_names(labels, lambda x: x in group_names)
    if indicate:
        return np.array(indicators)
    else:
        return np.array(labels)[indicators]


def take_adata_groups(adata: sc.AnnData,
                      key: str,
                      group_names: Sequence,
                      onlyx: bool = False,
                      copy: bool = False):
    """ Take given groups from an AnnData object """
    indicators = take_group_labels(adata.obs[key], group_names,
                                   indicate=True)
    if copy:
        return adata[indicators, :].X.copy() if onlyx else adata[indicators,
                                                           :].copy()
    else:
        return adata[indicators, :].X if onlyx else adata[indicators, :]


# ========================================================================
def ccaAdt(adata, groupby, key_add='X_cca', n_comps=50,
           method='plssvd', copy=False, **kwds):
    '''

    testing code:
    B.ccaAdt(ph.adata, groupby='RNAcap', n_comps=2, copy=True)
    '''
    labels = adata.obs[groupby]
    X = adata.X
    print('Computing cross-decomposition on 2 groups by `%s`' % groupby)
    X_cca, CCs = cca(X, labels, n_comps=n_comps, method=method,
                     return_info=False, **kwds)
    if copy:
        adata = adata.copy()
    adata.obsm[key_add] = X_cca
    for lb in CCs.keys():
        adata.varm['CCs_' + lb] = CCs[lb]

    return adata if copy else None


def cca(X, labels, n_comps=2, method='plssvd', return_model=False, scale=False,
        **kwds):
    '''
    inputs
    ======
    Suppose that n_cells = p + q.

    X : array-like, shape = [n_cells, n_genes] (will be transposed during computation)
    labels : array-like, shape = [n_cells, 1]
    scale : bool, False by default !

    return
    =======
    X_cca : array, [n_cells, n_components] (i.e. [p + q, n_components])
    x_weights_ : array, [p, n_components]
    y_weights_ : array, [q, n_components]
    x_scores_ : array, [n_genes, n_components]
    y_scores_ : array, [n_genes, n_components]

    if method is not 'plssvd':
        x_loadings_ : array, [p, n_components]
        y_loadings_ : array, [q, n_components]

    NOTE:
        * `n_features` here means `n_cells`
        * `n_samples here means `n_genes`
        * len(pd.unique(labels)) should be 2

    '''
    ## data
    n_cells = X.shape[0]
    X = X.T  # should be transposed
    labels = pd.Series(labels)
    unique_lbs = labels.unique()
    ind1 = (labels == unique_lbs[0])
    print(f'unique labels: {unique_lbs}')
    ## model construction
    from sklearn.cross_decomposition import PLSSVD, CCA, PLSCanonical
    print('Using method [ %s ], with [ %d ] components' % (method, n_comps))
    if method == 'plssvd':
        md = PLSSVD(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_weights_, md.y_weights_
    elif method == 'cca':  # NOT applicatable !!!
        md = CCA(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_loadings_, md.y_loadings_
    #        X1_, X2_ = md.x_rotations_,  md.y_rotations_
    elif method == 'plsc':
        md = PLSCanonical(n_components=n_comps, scale=scale, **kwds)
        md.fit(X[:, ind1], X[:, ~ind1])
        X1_, X2_ = md.x_weights_, md.y_weights_
    else:
        raise ValueError(
            'the argument `method` should be either `plssvd` or `cca`')

    X_cca = np.zeros((n_cells, n_comps))
    X_cca[ind1, :] = X1_
    X_cca[~ ind1, :] = X2_

    if return_model:
        return X_cca, md
    else:
        CCs = {str(unique_lbs[0]): md.x_scores_,
               str(unique_lbs[1]): md.y_scores_}
        return X_cca, CCs


def partialPCAdt(adata, based_on, key_add='X_pca', n_comps=50,
                 copy=False, **kwds):
    '''
    adata.X.shape: n_cells x n_genes
    based_on: a tuple, or a list of length 2, e.g. (key, class_name)
        based_on[0]: key for the class labels
        based_on[1]: class name of the subset of X to fit the PCA
    '''
    labels = adata.obs[based_on[0]]
    X = adata.X
    labels = pd.Series(labels)
    ind1 = (labels == based_on[1])
    print('Computing PCA on class `{}` grouped by `{}`'.format(*based_on))

    X1 = X[ind1, :]
    X2 = X[~ ind1, :]
    X1_pca, X2_pca, pca = partialPCA(X1, X2, n_comps=n_comps, **kwds)
    if copy:
        adata = adata.copy()

    X_pca = np.zeros((X.shape[0], n_comps))
    X_pca[ind1, :] = X1_pca
    X_pca[~ ind1, :] = X2_pca
    adata.obsm[key_add] = X_pca
    adata.uns['pca'] = {}
    adata.uns['pca']['variance'] = pca.noise_variance_
    adata.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
    adata.varm['PCs'] = pca.components_.T

    return adata if copy else None


def partialPCA(X1, X2, n_comps=2, **kwds):
    '''
    fit PCA on X, and transform both X and Y based on the fitted components
    X1: array-like, shape (n_samples1, n_features)
    X2: array-like, shape (n_samples2, n_features)

    return
    X_pca: array-like, shape (n_samples1 + n_samples2, n_features)
    estimator: `PCA` object

    '''
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_comps, **kwds)
    X1_pca = pca.fit_transform(X1)
    X2_pca = pca.transform(X2)

    return X1_pca, X2_pca, pca

