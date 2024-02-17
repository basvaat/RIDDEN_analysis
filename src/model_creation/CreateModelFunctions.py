# iterating 
import pandas as pd
import numpy as np
import statsmodels.api as sm
from src.config import *

def fit_linear_model_and_get_coefficients_by_receptor(y, X):
    # y = signature
    # X = model_matrix
    # fit linear model, y:signature, X:model matrix
    coeff_m = pd.DataFrame(columns=y.columns, index=[*X.columns])
    for gene in list(y.columns):
        for receptor in list(X.columns):
            Xj = X.loc[y.index, receptor]
            Xj = sm.add_constant(Xj)

            yi = y.loc[:, gene]

            assert(all(Xj.index == yi.index))
            model = sm.OLS(yi, Xj) 
            results = model.fit()
            coeff_s = results.params
            assert len(coeff_s) == 2
            assert coeff_s.index[1]==receptor

            coeff_m.loc[receptor, gene] = coeff_s.loc[receptor]

    return coeff_m
