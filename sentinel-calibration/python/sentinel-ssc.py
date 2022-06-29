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
from pathlib import Path

import matplotlib.pyplot as plt
# from matplotlib.colors import Normalize
# from scipy.stats import gaussian_kde
# import mpl_scatter_density

# Model Pre-processing
from sklearn import preprocessing
# from sklearn.model_selection import train_test_split
from yellowbrick.cluster import KElbowVisualizer

# Import Models
from sklearn import cluster as cl
from sklearn.linear_model import LinearRegression, Lasso, Ridge, ElasticNet


from pickle import dump, load
import json
import glob

from argparse import Namespace

### TO DO/ISSUES
# - Split training and test data (for regression)
# - Elbow graph to determine k - value: 
    # https://towardsdatascience.com/cheat-sheet-to-implementing-7-methods-for-selecting-optimal-number-of-clusters-in-python-898241e1d6ad
    # https://github.com/ivalencius/Machine_learning_INF264/blob/main/Homework/6/Hw6_INF264_template_part_1.ipynb
# - Add tag to specify regression variables
# - Cluster args hardcoded in
# - Make code pretty with docs and tqdm
# - Determine max value of sensor 
# - Implement other clustering
# - Clustering and regression save to different folders

def get_regressors(name):
    if name == 'raw bands':
        return  ['B1','B2','B3','B4', 'B5', 'B6', 'B7', 'B2.B1','B3.B1','B4.B1',
                 'B5.B1', 'B7.B1', 'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2', 'B4.B3', 
                 'B5.B3','B7.B3','B5.B4', 'B7.B4', 'B7.B5', 'B1.2', 'B2.2', 'B3.2', 
                 'B4.2', 'B5.2', 'B7.2']
    elif name == 'station':
        return ['B1','B2','B3','B4', 'B5', 'B6', 'B7', 'B2.B1','B3.B1','B4.B1',
                 'B5.B1', 'B7.B1', 'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2', 'B4.B3', 
                 'B5.B3','B7.B3','B5.B4', 'B7.B4', 'B7.B5', 'B1.2', 'B2.2', 'B3.2', 
                 'B4.2', 'B5.2', 'B7.2', 'station_nm'] # Add station name as dummy regression variable

def standardize(csv_path, reg_vars, mode):
    try:
        df = pd.read_csv(csv_path)
    except Exception:
        raise ValueError("!!! Unknown CSV !!!")
    # Extract data you wish to regress
    data = df[reg_vars].to_numpy()
    # Check to see if data has SSC data
    try:
        ssc = df['SSC_mgL'].to_numpy()
    except:
        ssc = None
    # Check to see if data has pred_SSC data
    try:
        pred_ssc = df['pred_SSC_mgL'].to_numpy()
    except:
        pred_ssc = None
    # Check to see if data is clustered
    try:
        num_clust = max(df['cluster'])+1
    except:
        num_clust = None
    # Remove infinite ssc values (set to 0)
    ssc[np.isnan(ssc)] = 0
    data[np.isnan(ssc)] = 0
    # Remove nan ssc values (set to 0)
    ssc[~np.isfinite(ssc)] = 0
    data[~np.isfinite(ssc)] = 0
    if mode == 'train':
        # Standardize data
        scaler = preprocessing.RobustScaler() # JUSTIFY ROBUST SCALER
        data_std = scaler.fit_transform(data)
        dump(scaler, open('algs\\scaler.pkl', 'wb'))
    elif mode == 'evaluate' or mode == 'infer':
        scaler = load(open('algs\\scaler.pkl', 'rb'))
        data_std = scaler.transform(data)
    return data_std, ssc, pred_ssc, scaler, num_clust

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
            return cl.KMeans(n_init=20, verbose=0, random_state=0)
    elif num_clust != None:
        if cluster_type == 'kMeans':
            return cl.KMeans(n_clusters=num_clust, n_init=1, verbose=0, random_state=0)

def load_cluster(cluster_type):
    file = glob.glob('algs\\'+cluster_type+'*.pkl')[0]
    try:
        return load(open(file, 'rb'))
    except:
        print('!!! No Cluster Algorithm Found !!!')

def cluster(data, cluster_type, csv_path, mode, scaler):
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
        # Once cluster is trained: save it
        dump(cluster, open('algs\\kmeans_'+str(k_optim)+'.pkl', 'wb'))
        # Save cluster information
        clust_info = {
            'clusters' : k_optim,
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
        # Generate false color plots of clusters
        false_color_clust(scaler.inverse_transform(cluster.cluster_centers_).astype(int), reg_vars, cluster_type)
    elif mode == 'infer':
        # Need to implement
        cluster = load_cluster(cluster_type)
        # Save data (append to non-standardized data)
        df = pd.read_csv(csv_path)
        labels = cluster.predict(data)
        df['cluster'] = labels
        cluster_file = 'clusters\\clustered_' +cluster_type+'.csv' 
        df.to_csv(cluster_file, index=False)
        
def get_reg(reg_type):
    if reg_type == 'linear':
        return LinearRegression()
    if reg_type == 'lasso':
        return Lasso(max_iter=10**5)
    if reg_type == 'ridge':
        return Ridge(max_iter=10**5)
    if reg_type == 'elasticNet':
        return ElasticNet(max_iter=10**5, random_state=0)

def load_reg(reg_type):
    algs = []
    # Ensure first value in list is first cluster
    regression_algs = sorted(glob.glob('regression\\'+reg_type+'\\*.pkl'))
    for reg_alg in regression_algs:
        algs.append(load(open(reg_alg, 'rb')))
    return algs

def regress(data, ssc, num_clust, reg_type, csv_path, mode, scaler):
    # Assumes data already has clusters, implement automatically setting the
    if num_clust == None:
        print('NEED TO IMPLEMENT CLUSTERING FOR REGRESSION')
    # Create dictionary to store metrics
    reg_info = {'reg_vars' : reg_vars}
    reg_folder = 'regression\\'+reg_type+'\\'
    if not os.path.exists(reg_folder):
        os.makedirs(reg_folder)
    # Get cluster rows
    df = pd.read_csv(csv_path)
    clusters = df['cluster']
    if mode == 'train':
        for clust in range(num_clust):
            c = str(clust)
            # Extract data for clusters
            data_clust = data[clusters == clust]
            ssc_clust = ssc[clusters == clust]
            # Get regressor
            regressor = get_reg(reg_type)
            # Fit regressor
            reg = regressor.fit(data_clust, ssc_clust)
            # Save regressor
            dump(reg, open(reg_folder+c+'.pkl', 'wb'))
            # Map to dictionary
            reg_info[c+'_R2'] = reg.score(data_clust, ssc_clust)
            reg_info[c+'_coefficients'] = reg.coef_.tolist()
            reg_info[c+'_intercept'] = reg.intercept_
        # Output regression info to json file
        with open(reg_folder+'info.json', 'w') as file:
            file.write(json.dumps(reg_info))
        # Output pretty print data to txt for readability (not loading data)
        with open(reg_folder+'pprint_info.txt', 'w') as file:
                file.write(json.dumps(reg_info, indent=4, sort_keys=True))
        # Save data (append to non-standardized data)
        # df = pd.read_csv(csv_path)
        # df['pred_SSC_mgL'] = reg.predict(data)
        # reg_file = 'regression\\reg_' +reg_type+'.csv' 
        # df.to_csv(reg_file, index=False)
    elif mode == 'infer':
        # Get regression algorithms for each cluster
        regression_algs = load_reg(reg_type)
        # Create dummy variable to store predicted SSC
        pred_ssc = []
        # Get cluster rows
        df = pd.read_csv(csv_path)
        rows = np.shape(data)[0]
        clusters = df['cluster']
        for i in range(rows):
            # Extract regressor for this cluster
            reg = regression_algs[clusters[i]]
            # Extract data
            data_row = data[i, :]
            # Predict ssc
            ssc = reg.predict(data_row.reshape(1,-1))
            # Append predicted data
            pred_ssc.append(ssc)
        # Save data (append to non-standardized data)
        df['pred_SSC_mgL'] = np.array(pred_ssc)
        reg_file = 'regression\\reg_' +reg_type+'.csv' 
        df.to_csv(reg_file, index=False)
        

def evaluate(ssc, pred_ssc, num_clust, csv_path):
    assert np.shape(ssc) == np.shape(pred_ssc), "Number of predicted ssc measurements not equal to number of in-situ measurements"
    # Convert to int to make plots look better (no 10**-3 artifacts and such)
    ssc = ssc.astype(int)
    pred_ssc = pred_ssc.astype(int)
    # Extract csv basename
    csv_name = Path(csv_path).stem
    
    # # Set bounds to remove outliers
    # max_ssc = np.inf
    # # Filter by in-situ ssc
    # rows = ssc < max_ssc
    # ssc = ssc[rows]
    # pred_ssc = pred_ssc[rows]
    # # Filter by predicted ssc
    # rows = pred_ssc < max_ssc
    # ssc = ssc[rows]
    # pred_ssc = pred_ssc[rows]
    
    # Normal 1:1 plot
    fig, ax = plt.subplots()
    ax.scatter(ssc, pred_ssc, color='black', alpha = 0.05)
    plt.xscale('log')
    plt.yscale('log')
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    lims = [0, np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, '--', color='red')
    #ax.grid(False)
    #ax.axis('off')
    plt.xlabel('in situ SSC (mg/L)')
    plt.ylabel('Satellite-estimated SSC (mg/L)')
    plt.title(csv_name+' SSC Correlation')
    plt.savefig('figures\\'+csv_name+'_correlation.pdf')
    
    # Density plot -- NEED TO IMPLEMENT
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    # ax.scatter_density(ssc, pred_ssc, dpi=10, cmap='jet')
    # plt.show()
    

        

if __name__ == "__main__":
    """
    Usage:
        python sentinel-ssc cluster --mode 'train' --csv 'PATH TO CSV' --cluster_type 'cluster' --reg_vars 'raw bands' # Run clustering algorithm
    """
    parser = argparse.ArgumentParser(description="Run a clustering regression pipeline on SSC data.")
    parser.add_argument(
        "task", metavar="task", default="cluster", choices=("cluster","regression", "evaluate"), type=str, help="task of workflow")
    parser.add_argument("--mode", default="infer", choices=("train","infer"), type=str, help="workflow mode")
    parser.add_argument("--csv", default="", type=str, help="path to CSV to import")
    parser.add_argument("--cluster_type", choices=("kMeans"), default="kMeans", type=str, help="clustering algorithm")
    parser.add_argument("--reg_type", choices=("linear","lasso","ridge","elasticNet"), default="linear", type=str, help="regression algorithm")
    parser.add_argument("--reg_vars", choices=('raw bands', 'station'), default='raw bands', type=str, help='variables to regress')
    args = parser.parse_args()
    
    # Set working directory
    wd = 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python'
    print('Setting working directory as', 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python')
    os.chdir(wd)
    # Set up subfolders
    folders = ['algs\\','clusters\\','figures\\','regression\\']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    args = Namespace(
        task = 'evaluate',
        mode='infer',
        #csv="D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\tmp_vars\\ex_landsat.csv",
        #csv ="D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python\\clusters\\clustered_kMeans.csv",
        csv = "D:\\valencig\\Thesis\sentinel-ssc\\sentinel-calibration\\python\\regression\\reg_elasticNet.csv",
        #reg_type = "elasticNet",
        reg_vars = 'raw bands')
    
    # Extract regression variables
    reg_vars = get_regressors(args.reg_vars)
    # Extract and standardize data to numpy array
    data, ssc, pred_ssc, scaler, num_clust = standardize(args.csv, reg_vars, args.mode)
    
    if args.task == "cluster":
        cluster(data=data, cluster_type=args.cluster_type, csv_path=args.csv, mode=args.mode, scaler=scaler)
    elif args.task == "regression":
        regress(data=data, ssc=ssc, num_clust=num_clust, reg_type=args.reg_type, csv_path=args.csv, mode=args.mode, scaler=scaler)
    elif args.task == "evaluate":
        evaluate(ssc=ssc, pred_ssc=pred_ssc, num_clust=num_clust, csv_path=args.csv)
    else:
        raise ValueError("!!! Unknown mode or missing data !!!")