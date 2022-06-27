# -*- coding: utf-8 -*-
"""
Code to automatically run clustering and regression on harmonized Landsat/USGS SSC data from sentinel-ssc.R

Created on Fri Jun 24 13:23:38 2022

@author: valencig
"""
import argparse
import pandas as pd
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from sklearn import cluster as cl
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from yellowbrick.cluster import KElbowVisualizer

from pickle import dump, load
import json

from argparse import Namespace

### TO DO/ISSUES
# - Set working directory in R
# - Split training and test data (for regression)
# - Save clustering centers to file
# - Elbow graph to determine k - value: 
    # https://towardsdatascience.com/cheat-sheet-to-implementing-7-methods-for-selecting-optimal-number-of-clusters-in-python-898241e1d6ad
    # https://github.com/ivalencius/Machine_learning_INF264/blob/main/Homework/6/Hw6_INF264_template_part_1.ipynb
# - Add tag to specify regression variables
# - Cluster args hardcoded in
# - https://scikit-learn.org/stable/modules/clustering.html#
# - Make code pretty with docs and tqdm
# - Determine max value of sensor 

def get_regressors(name):
    if name == 'raw bands':
        return  ['B1','B2','B3','B4', 'B5', 'B6', 'B7', 'B2.B1','B3.B1','B4.B1',
                 'B5.B1', 'B7.B1', 'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2', 'B4.B3', 
                 'B5.B3','B7.B3','B5.B4', 'B7.B4', 'B7.B5', 'B1.2', 'B2.2', 'B3.2', 
                 'B4.2', 'B5.2', 'B7.2']
    elif name == 'with drainage':
        ...
    elif name == 'with width':
        ...
    elif name == 'drainage and width':
        ...

def standardize(csv_path, reg_vars, mode):
    csv_path = 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\tmp_vars\\ex_landsat.csv'
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        raise ValueError("Unknown CSV")
    # Extract data you wish to regress
    data = df[reg_vars].to_numpy()
    ssc = df['SSC_mgL'].to_numpy()
    if mode == 'train':
        # Standardize data
        scaler = preprocessing.RobustScaler() # JUSTIFY ROBUST SCALER
        data_std = scaler.fit_transform(data)
        dump(scaler, open('algs\\scaler.pkl', 'wb'))
    if mode == 'infer':
        scaler = load(open('algs\\scaler.pkl', 'rb'))
        data_std = scaler.transform(data)
    return data_std, ssc, scaler

def false_color_clust(centers, reg_vars, cluster_type):
    # Need to scale inputs into 0-255 range so need to know maximum sensor value
    max_val = 1000 # FIX, dependent on sensor
    # Get indexes of RGB in reg vars
    r = reg_vars.index('B3') 
    g = reg_vars.index('B2') 
    b = reg_vars.index('B1')
    # Get number of clusters
    num_clust = np.shape(centers)[0]
    # Make plots
    Cols = 3
    # Compute Rows required
    Rows = num_clust // Cols 
    Rows += num_clust % Cols
    # Create a Position index
    Position = range(1,num_clust + 1)
    plt.clf()
    fig = plt.figure(1)
    for i in range(num_clust):
    #for i in range(1):
        ax = fig.add_subplot(Rows,Cols,Position[i])
        rgb = (int(255 / max_val * centers[i, r]), int(255 / max_val * centers[i, g]), int(255 / max_val * centers[i, b]))
        ax.imshow([[rgb]], vmin=0, vmax=255)
        ax.grid(False)
        ax.axis('off')
        ax.title.set_text('Clust '+str(i))
        
    fig.suptitle('False Color RGB of Cluster Centers', y=0.80)
    plt.savefig('figures\\'+cluster_type+'_cluster_viz.pdf')
    
def get_cluster(cluster_type, num_clust = None):
    if num_clust == None:
        if cluster_type == 'kMeans':
            clust = cl.KMeans(n_init=20, verbose=0, random_state=0)
    elif num_clust != None:
        if cluster_type == 'kMeans':
            clust = cl.KMeans(n_clusters=num_clust, n_init=1, verbose=0, random_state=0)
    return clust

def load_cluster(cluster_type):
    if cluster_type == 'kMeans':
        clust = load(open('algs\\kmeans_model.pkl', 'wb'))
    return clust

def cluster(data, cluster_type, csv_path, mode, scaler):
    cluster_type = 'kMeans'
    csv_path = "D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\tmp_vars\\ex_landsat.csv"

    if mode == 'train':
        cluster = get_cluster(cluster_type)
        visualizer = KElbowVisualizer(cluster, k=(1,10))
        
        # Fit data and create elbow graph
        visualizer.fit(data)
        # Render elbow graph
        visualizer.show(outpath="figures\\"+cluster_type+'_elbow.pdf')
        # Select optimal k value
        k_optim = int(visualizer.elbow_value_)
        # Apply and train cluster on optimal k value
        cluster = get_cluster(cluster_type, num_clust = k_optim)
        cluster.fit(data)
        labels = cluster.predict(data)
        # Once cluster is trained: save it
        dump(cluster, open('algs\\kmeans_'+str(k_optim)+'.pkl', 'wb'))
        # Save cluster information
        clust_info = {
            'reg_vars' : reg_vars,
            'sum_sqr_dists' : cluster.inertia_,
            'centers' : cluster.cluster_centers_.tolist(),
            'centers_inv' : scaler.inverse_transform(cluster.cluster_centers_).tolist()
            }
        clust_info = {**cluster.get_params(), **clust_info}
        # Output clustering info to json file
        with open('clusters\\'+cluster_type+'_info.json', 'w') as file:
            file.write(json.dumps(clust_info))
        # Output pretty print data to txt for readability (not loading data)
        with open('clusters\\pprint_'+cluster_type+'_info.txt', 'w') as file:
                file.write(json.dumps(clust_info, indent=4, sort_keys=True))
        # Save data (append to non-standardized data)
        df = pd.read_csv(csv_path)
        df['cluster'] = labels
        cluster_file = 'clusters\\clustered_' +cluster_type+'_'+str(k_optim)+'.csv' 
        df.to_csv(cluster_file, index=False)
        # Generate false color plots of clusters
        false_color_clust(scaler.inverse_transform(cluster.cluster_centers_).astype(int), reg_vars, cluster_type)
        
    elif mode == 'infer':
        cluster = load_cluster(cluster_type)
        

if __name__ == "__main__":
    """
    Usage:
        python sentinel-ssc cluster --mode 'train' --csv 'PATH TO CSV' --cluster_type 'cluster' --reg_vars 'raw bands' # Run clustering algorithm
    """
    parser = argparse.ArgumentParser(description="Run a clustering regression pipeline on SSC data.")
    parser.add_argument(
        "task", metavar="task", default="cluster", choices=("cluster","regress"), type=str, help="task of workflow"
    )
    parser.add_argument("--mode", default="train", choices=("train","infer"), type=str, help="workflow mode")
    parser.add_argument("--csv", default="", type=str, help="path to CSV to import")
    parser.add_argument("--cluster_type", choices=("kMeans"), default="kMeans", type=str, help="clustering algorithm")
    parser.add_argument("--reg_vars", choices=('raw bands', 'with drainage', 'with width', 'drainage and width'), default='raw bands', type=str, help='variables to regress')
    args = parser.parse_args()
    
    # Set working directory
    wd = 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python'
    print('Setting working directory as', 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python')
    os.chdir(wd)
    
    args = Namespace(
        task = 'cluster',
        mode='train',
        csv="D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\tmp_vars\\ex_landsat.csv",
        cluster_type = 'kMeans',
        reg_vars = 'raw bands')

    
    # Extract regression variables
    reg_vars = get_regressors(args.reg_vars)
    # Extract and standardize data to numpy array
    data, ssc_vals, scaler = standardize(args.csv, reg_vars, args.mode)
        
    if args.task == "cluster":
        cluster(data=data, cluster_type=args.cluster_type, csv_path=args.csv, mode=args.mode, scaler=scaler)
    elif args.task == "regress":
        ... # neeed to use ssc_vals in regression and
    else:
        raise ValueError("Unknown mode.")