# -*- coding: utf-8 -*-
"""
Code to automatically run all clustering/regression permutations in sentinel-ssc.py.

Created on Wed Jul  6 11:33:31 2022

@author: valencig
"""

import os
from tqdm import tqdm
    
# Set working directory
wd = 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\'
print('Setting working directory as', 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python')
os.chdir(wd)

function_calls = [
    # Train Cluster
    'python sentinel-ssc_no_standardizing.py --task cluster --mode train --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\exports\\GEE_raw\\ssc_harmonized.csv" --cluster_type kMeans --reg_vars raw_bands',
    # Infer cluster labels
    'python sentinel-ssc_no_standardizing.py --task cluster --mode infer --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\exports\\GEE_raw\\ssc_harmonized.csv" --cluster_type kMeans --reg_vars raw_bands',
    # Train all 4 regressors
    'python sentinel-ssc_no_standardizing.py --task regression --mode train --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type linear --reg_vars full_bands --holdout 0.3',
    'python sentinel-ssc_no_standardizing.py --task regression --mode train --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type lasso --reg_vars full_bands --holdout 0.3',
    'python sentinel-ssc_no_standardizing.py --task regression --mode train --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type ridge --reg_vars full_bands --holdout 0.3',
    'python sentinel-ssc_no_standardizing.py --task regression --mode train --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type elasticNet --reg_vars full_bands --holdout 0.3',
    # Infer SSC from regressors
    'python sentinel-ssc_no_standardizing.py --task regression --mode infer --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type linear --reg_vars full_bands --holdout 0.',
    'python sentinel-ssc_no_standardizing.py --task regression --mode infer --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type lasso --reg_vars full_bands --holdout 0.',
    'python sentinel-ssc_no_standardizing.py --task regression --mode infer --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type ridge --reg_vars full_bands --holdout 0.',
    'python sentinel-ssc_no_standardizing.py --task regression --mode infer --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\clusters\\clustered_kMeans_raw_bands.csv" --reg_type elasticNet --reg_vars full_bands --holdout 0.',
    # Evaluate all regressions
    'python sentinel-ssc_no_standardizing.py --task evaluate --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\regression\\reg_linear.csv" --reg_vars full_bands',
    'python sentinel-ssc_no_standardizing.py --task evaluate --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\regression\\reg_lasso.csv" --reg_vars full_bands',
    'python sentinel-ssc_no_standardizing.py --task evaluate --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\regression\\reg_elasticNet.csv" --reg_vars full_bands',
    'python sentinel-ssc_no_standardizing.py --task evaluate --csv "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\no_standard\\regression\\reg_ridge.csv" --reg_vars full_bands'
    ]

for function_call in tqdm(function_calls):
    os.system(function_call)
