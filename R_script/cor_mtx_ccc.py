# The python scripts which apply a new correlation calculation method #
# Test on my own PC with the files with 1000 samples and all gene sets #
# A time cost of 20 hours using single core #

import numpy as np
import pandas as pd
from ccc.coef import ccc
import os
import math
import warnings
from tqdm import tqdm
np.seterr(divide = 'ignore')
warnings.simplefilter(action='ignore', category=FutureWarning)


os.chdir('/Users/kimhan/Desktop/CCI_DK/ccc')
# Import the expression matrix
gene_expression_df = pd.read_csv('./gene_expression_matrix.csv', index_col=0)
print(gene_expression_df.shape)
gene_set_expression_df = pd.read_csv('./gene_set_expression_matrix.csv', index_col=0)
print(gene_set_expression_df.shape)
icn_df = pd.read_csv('./icn.csv', index_col=0)

ligand = icn_df.loc[:, 'source_genesymbol'].drop_duplicates().values.tolist()
print(len(ligand))
receptor = icn_df.loc[:, 'target_genesymbol'].drop_duplicates().values.tolist()
print(len(receptor))
ligand_receptor = ligand + receptor


def gm_mean(x, zero_propagate=True):
    if (np.array(x) < 0).any():
        return None
    if zero_propagate is True:
        if 0 in x:
            return 0
        else:
            g_mean = math.exp(np.mean(np.log(x)))
            return g_mean
    else:
        g_mean = math.exp(sum(np.log([i for i in x if i > 0]))/len(x))
        return g_mean


def complex_expression(unit_list):
    exp = gene_expression_df.loc[unit_list, :]
    complex_exp = exp.apply(gm_mean,axis=0)
    return complex_exp


cor_df = pd.DataFrame(columns= ligand_receptor,index= gene_set_expression_df.index.values)
loop = 0
with tqdm(total=cor_df.shape[1]) as pbar:
    for lr_symbol in ligand_receptor:
        loop += 1
        print(lr_symbol)
        if '_' in lr_symbol:
            unit_list = lr_symbol.split(sep='_')
            if set(unit_list).issubset(set(gene_expression_df.index.values)):
                lr_score = complex_expression(unit_list)
            else:
                continue
        else:
            if lr_symbol in gene_expression_df.index:
                lr_score = gene_expression_df.loc[lr_symbol]
            else:
                continue
        for gene_set_symbol in gene_set_expression_df.index.values:
            gene_set_score = gene_set_expression_df.loc[gene_set_symbol]
            cor = ccc(np.log(lr_score+1),np.log(gene_set_score+1))  # logarithmic transformation and psedo-count 1
            cor_df.at[gene_set_symbol,lr_symbol] = cor
        pbar.update(loop)
    print('All is done!')

cor_df.to_csv('cor_matrix_ccc.csv')
print('The correlation matrix is saved in /home/projects/kvs_ccc/cor_matrix/cor_matrix_ccc.csv')

# Other methods to calculate the correlation in python #
# pearson and spearman
# cor_pearson = lr_score.corr(gene_set_score,method='pearson')
# cor_spearman = lr_score.corr(gene_set_score,method='spearman')
