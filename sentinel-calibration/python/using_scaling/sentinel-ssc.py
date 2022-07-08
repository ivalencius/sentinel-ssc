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
import math

import matplotlib.pyplot as plt
from matplotlib import patches
# from matplotlib import colors
# from matplotlib.colors import Normalize
# from scipy.stats import gaussian_kde
# import mpl_scatter_density
# from sklearn.neighbors import KernelDensity

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
# - Log10 in relative errror plot
# - Relative errors just for holdout set??
# - Make code pretty with docs and tqdm
# - Determine max value of sensor 
# - Implement other clustering

def get_regressors(name):
    if name == 'raw_bands':
        return  ['B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11','B12']
    elif name == 'full_bands':
        return  ['B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B8A', 'B9', 'B11', 
                 'B12','B2.B1','B3.B1','B4.B1','B1.2','B2.2','B3.2','B4.2',
                 'B5.2','B6.2','B7.2','B8.2','B8A.2','B9.2','B11.2','B12.2',
                 'B2.B1','B3.B1','B4.B1','B5.B1','B6.B1','B7.B1','B8.B1',
                 'B8A.B1','B9.B1','B11.B1','B12.B1','B3.B2','B4.B2','B5.B2',
                 'B6.B2','B7.B2','B8.B2','B8A.B2','B9.B2','B11.B2','B12.B2',
                 'B4.B3','B5.B3','B6.B3','B7.B3','B8.B3','B8A.B3','B9.B3',
                 'B11.B3','B12.B3','B5.B4','B6.B4','B7.B4','B8.B4','B8A.B4',
                 'B9.B4','B11.B4','B12.B4','B6.B5','B7.B5','B8.B5','B8A.B5',
                 'B9.B5','B11.B5','B12.B5','B7.B6','B8.B6','B8A.B6','B9.B6',
                 'B11.B6','B12.B6','B8.B7','B8A.B7','B9.B7','B11.B7','B12.B7',	
                 'B8A.B8','B9.B8','B11.B8','B12.B8','B9.B8A','B11.B8A',
                 'B12.B8A','B11.B9','B12.B9','B12.B11']

    # elif name == 'station':
    #     return ['B1','B2','B3','B4', 'B5', 'B6', 'B7', 'B2.B1','B3.B1','B4.B1',
    #              'B5.B1', 'B7.B1', 'B3.B2', 'B4.B2', 'B5.B2', 'B7.B2', 'B4.B3', 
    #              'B5.B3','B7.B3','B5.B4', 'B7.B4', 'B7.B5', 'B1.2', 'B2.2', 'B3.2', 
    #              'B4.2', 'B5.2', 'B7.2', 'station_nm'] # Add station name as dummy regression variable
    
def data_extract(csv_path, reg_vars):
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
    # Set nan values to 0 and +inf and -inf to large and small #'s
    ssc = np.nan_to_num(ssc, nan=0.0)
    data = np.nan_to_num(data, nan=0.0, posinf=10**6)
    return data, ssc, pred_ssc, num_clust
    
def standardize(csv_path, reg_vars):
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
    # Set nan values to 0 and +inf and -inf to large and small #'s
    ssc = np.nan_to_num(ssc, nan=0.0)
    data = np.nan_to_num(data, nan=0.0, posinf=10**6)
    # Use whole dataset for cluster scaler
    if not os.path.exists('algs\\cluster_scaler.pkl'):
        # Standardize data
        scaler = preprocessing.RobustScaler() # JUSTIFY ROBUST SCALER
        data_std = scaler.fit_transform(data)
        dump(scaler, open('algs\\cluster_scaler.pkl', 'wb'))
    else:
        scaler = load(open('algs\\cluster_scaler.pkl', 'rb'))
        data_std = scaler.transform(data)
    return data_std, ssc, pred_ssc, scaler, num_clust

# def standardize_ssc(ssc, mode):
#     if mode == 'train':
#         # Standardize data
#         scaler = preprocessing.RobustScaler() # JUSTIFY ROBUST SCALER
#         ssc_std = scaler.fit_transform(ssc.reshape(-1,1))
#         dump(scaler, open('algs\\ssc_scaler.pkl', 'wb'))
#     elif mode == 'evaluate' or mode == 'infer':
#         scaler = load(open('algs\\ssc_scaler.pkl', 'rb'))
#         ssc_std = scaler.transform(ssc.reshape(-1,1))
#     return ssc_std

def false_color_clust(centers, reg_vars, cluster_type):
    # Need to scale inputs into 0-255 range so need to know maximum sensor value
    max_val = 2000 # FIX, dependent on sensor
    # Get indexes of RGB in reg vars
    r = reg_vars.index('B4') 
    g = reg_vars.index('B3') 
    b = reg_vars.index('B2')
    # Get number of clusters
    num_clust = np.shape(centers)[0]
    # Make plots
    Cols = 4
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
        # elif cluster_type == 'agglomerative':
        #     return cl.AgglomerativeClustering()
    elif num_clust != None:
        if cluster_type == 'kMeans':
            return cl.KMeans(n_clusters=num_clust, n_init=1, verbose=0, random_state=0)
        # elif cluster_type == 'agglomerative':
        #     return cl.AgglomerativeClustering(n_clusters=num_clust)

def load_cluster(cluster_type, reg_names):
    file = glob.glob('algs\\'+cluster_type+'_'+reg_names+'.pkl')[0]
    print(file)
    try:
        return load(open(file, 'rb'))
    except:
        print('!!! No Cluster Algorithm Found !!!')

def cluster(data, cluster_type, csv_path, mode, scaler, reg_vars, reg_names):
    if mode == 'train':
        cluster = get_cluster(cluster_type)
        visualizer = KElbowVisualizer(cluster, k=(2,10))
        
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
        dump(cluster, open('algs\\'+cluster_type+'_'+reg_names+'.pkl', 'wb'))
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
        with open('clusters\\'+cluster_type+'_'+reg_names+'_info.json', 'w') as file:
            file.write(json.dumps(clust_info))
        # Output pretty print data to txt for readability (not loading data)
        with open('clusters\\pprint_'+cluster_type+'_'+reg_names+'_info.txt', 'w') as file:
            file.write(json.dumps(clust_info, indent=4, sort_keys=True))
        # Generate false color plots of clusters
        false_color_clust(scaler.inverse_transform(cluster.cluster_centers_), reg_vars, cluster_type)
    elif mode == 'infer':
        # Need to implement
        cluster = load_cluster(cluster_type, reg_names)
        # Save data (append to non-standardized data)
        df = pd.read_csv(csv_path)
        labels = cluster.predict(data)
        # Create Histogram
        fig, ax = plt.subplots()
        ax.hist(labels, bins=max(labels)+1, density=True, histtype='bar', ec='black')
        ax.locator_params(axis='x', integer=True)
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Probability')
        plt.title('Cluster Distribution')
        plt.savefig('figures\\'+cluster_type+'_cluster_distribution.pdf')
        
        df['cluster'] = labels
        cluster_file = 'clusters\\clustered_' +cluster_type+'_'+reg_names+'.csv' 
        df.to_csv(cluster_file, index=False)
        
def get_reg(reg_type):
    if reg_type == 'linear':
        return LinearRegression()
    if reg_type == 'lasso':
        return Lasso(max_iter=10**5, alpha=14.0, selection='random')
    if reg_type == 'ridge':
        return Ridge(max_iter=10**5,)
    if reg_type == 'elasticNet':
        return ElasticNet(max_iter=10**5, random_state=0, alpha=14.0)

def load_reg(reg_type, reg_names):
    algs = []
    # Ensure first value in list is first cluster
    regression_algs = sorted(glob.glob('regression\\'+reg_type+'_'+reg_names+'\\*.pkl'))
    for reg_alg in regression_algs:
        algs.append(load(open(reg_alg, 'rb')))
    return algs

def regress(data, ssc, num_clust, reg_type, csv_path, mode, holdout, reg_vars, reg_names):
    # Assumes data already has clusters, implement automatically setting the
    if num_clust == None:
        print('NEED TO IMPLEMENT CLUSTERING FOR REGRESSION')
    # Get cluster rows
    df = pd.read_csv(csv_path)
    clusters = df['cluster'].to_numpy()
    if mode == 'train':
        # Create dictionary to store metrics
        reg_info = {'reg_vars' : reg_names,
                    'holdout_frac' : holdout}
        reg_folder = 'regression\\'+reg_type+'_'+reg_names+'\\'
        if not os.path.exists(reg_folder):
            os.makedirs(reg_folder)
            
        assert len(data) == len(ssc) and len(ssc) == len(clusters)
        num_samples = len(data) # rows x cols
        num_holdout = int(holdout * num_samples)
        np.random.seed(0)
        random_rows = np.random.choice(len(ssc), size=num_holdout, replace=False)
        data_holdout = np.copy(data[random_rows])
        ssc_holdout = np.copy(ssc[random_rows])
        clusters_holdout = np.copy(clusters[random_rows])
        # Set rows not selected to null
        data = np.delete(data, obj=random_rows, axis=0)
        ssc = np.delete(ssc, obj=random_rows, axis=0)
        clusters = np.delete(clusters, obj=random_rows, axis=0)
        # Standardize data
        scaler = preprocessing.RobustScaler() # JUSTIFY ROBUST SCALER
        # Only use training set to train scaler
        if os.path.exists('algs\\regression_scaler.pkl'):
            scaler = load(open('algs\\regression_scaler.pkl', 'rb'))
            data = scaler.transform(data)
            data_holdout = scaler.transform(data_holdout)
        else:
            data = scaler.fit_transform(data)
            data_holdout = scaler.transform(data_holdout)
            dump(scaler, open('algs\\regression_scaler.pkl', 'wb'))
        
        assert len(data_holdout) == len(ssc_holdout) and len(ssc_holdout) == len(clusters_holdout)
        
        R2_scores = []
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
            # Extract data from holdout set
            metric_data = data_holdout[clusters_holdout == clust]
            metric_ssc = ssc_holdout[clusters_holdout == clust]
            # Get metrics and export
            reg_info[c+'_train_samples'] = len(ssc_clust)
            reg_info[c+'_validation_samples'] = len(metric_ssc)
            reg_info[c+'_training_frac'] = len(ssc_clust) / (len(ssc_clust) + len(metric_ssc))
            reg_info[c+'_R2'] = reg.score(metric_data, metric_ssc)
            R2_scores.append(reg.score(metric_data, metric_ssc))
            # Undo normalization for coefficients
            reg_info[c+'_coefficients_unormalized'] = (reg.coef_ / scaler.scale_ ).tolist()
            reg_info[c+'_intercept'] = reg.intercept_
        # Get average R2 score
        reg_info['unweighted_R2'] = R2_scores
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
        # Standardize Data
        scaler = load(open('algs\\regression_scaler.pkl', 'rb'))
        data = scaler.transform(data)
        # Get regression algorithms for each cluster
        regression_algs = load_reg(reg_type, reg_names)
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
        df['pred_SSC_mgL'] = np.array(pred_ssc).astype(int)
        reg_file = 'regression\\reg_' +reg_type+'.csv' 
        df.to_csv(reg_file, index=False)
        
def relative_error_fxn(row):
    if row[0] != 0 and row[1] != 0:
        return math.log10(abs(row[1])/abs(row[0]))
    else:
        return 0
    
def false_color_reg(data, pred_ssc, num_clust, csv_path, reg_vars):
    df = pd.read_csv(csv_path)
    max_val = 2500
    # Get indexes of RGB in reg vars
    r = reg_vars.index('B4') 
    g = reg_vars.index('B3') 
    b = reg_vars.index('B2')
    # Set SSC levels
    ssc_levels = [50, 100, 250, 500, 750, 10**4]
    fig, ax = plt.subplots()
    for clust in range(num_clust):
        data_clust = df.loc[df['cluster'] == clust]
        for i in range(len(ssc_levels)):
            ssc_max = ssc_levels[i]
            if i == 0:
                ssc_min = 0
            else:
                ssc_min = ssc_levels[i-1]
            # Filter by max and min
            ssc_data = data_clust.loc[(data_clust['SSC_mgL'] <= ssc_max) & (data_clust['SSC_mgL'] >= ssc_min)]
            ssc_data = ssc_data[reg_vars].to_numpy()
            # Get average value of bands
            med = np.median(ssc_data, axis=0)
            rgb = ((1 / max_val * med[r]), (1 / max_val * med[g]), (1 / max_val * med[b]), 1)
            rectangle = patches.Rectangle((clust, i), 1, 1, facecolor=rgb)
            ax.add_patch(rectangle)
    # Make axes and ticks look categorical
    ax.set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
    ax.set_yticklabels(['0-50', '50-100', '100-250', '250-500', '500-750', '>750'])
    ax.set_ylabel('USGS Determined SSC (mg/L)')
    ax.set_xticks([0.5, 1.5, 2.5, 3.5])
    ax.set_xticklabels(['1', '2', '3', '4'])
    # Deal with removing grid and chaning outline to black
    plt.margins(x=0.04, y=0.04)
    plt.grid(False)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black') 
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.set_xlabel('River Grouping')
    plt.title('Median "True Color" Appearance')
    plt.savefig('figures\\ssc_level_viz.pdf')

def evaluate(data, ssc, pred_ssc, num_clust, csv_path, reg_vars):
    assert np.shape(ssc) == np.shape(pred_ssc), "Number of predicted ssc measurements not equal to number of in-situ measurements"
    # Ensure int ssc measurementes (no small decimal artifacts)
    ssc = ssc.astype(int)
    pred_ssc = pred_ssc.astype(int)
    # Extract csv basename
    csv_name = Path(csv_path).stem
    
    # Determin relative rror
    combined = np.concatenate([ssc.reshape(-1,1), pred_ssc.reshape(-1,1)], axis=1)
    log_arr = np.apply_along_axis(relative_error_fxn, -1, combined)
    relative_error = 10**(np.median(log_arr))-1
        
    # Normal 1:1 plot
    fig, ax = plt.subplots()
    ax.scatter(ssc, pred_ssc, color='black', alpha = 0.05)
    plt.xscale('log')
    plt.yscale('log')
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.grid(False)
    lims = [0, np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, '--', color='red')
    plt.xlabel('in situ SSC (mg/L)')
    plt.ylabel('Satellite-estimated SSC (mg/L)')
    plt.figtext(0.2, 0.9, 'Relative error = '+str(round(relative_error, 3)), ha='center', va='center',transform=ax.transAxes)
    plt.title(csv_name+' SSC Correlation')
    plt.savefig('figures\\'+csv_name+'_correlation.pdf')
    
    false_color_reg(data, pred_ssc, num_clust, csv_path, reg_vars)
    

if __name__ == "__main__":
    """
    Usage:
        python sentinel-ssc cluster --mode 'train' --csv 'PATH TO CSV' --cluster_type 'cluster' --reg_vars 'raw bands' # Run clustering algorithm
    """
    parser = argparse.ArgumentParser(description="Run a clustering regression pipeline on SSC data.")
    parser.add_argument("--task", choices=("cluster","regression", "evaluate"), type=str, help="task of workflow")
    parser.add_argument("--mode", default="infer", choices=("train","infer"), type=str, help="workflow mode")
    parser.add_argument("--csv", type=str, help="path to CSV to import")
    parser.add_argument("--cluster_type", choices=("kMeans", "NOT IMPLEMENTED"), type=str, help="clustering algorithm")
    parser.add_argument("--reg_type", choices=("linear","lasso","ridge","elasticNet"), type=str, help="regression algorithm")
    parser.add_argument("--reg_vars", choices=('raw_bands', 'full_bands'), type=str, help='variables to regress')
    parser.add_argument("--holdout", type=float, help='holdout percentage for validation set')
    args = parser.parse_args()
    
    # Check for correct args
    if args.task == 'cluster' and (args.mode is None or args.csv is None or args.cluster_type is None):
        parser.error("clustering requires --mode --csv -- cluster_type --reg_vars")
    elif args.task == 'regression' and args.mode == 'train' and (args.csv is None or args.reg_type is None or args.holdout is None):
        parser.error("training regression requires --csv -- reg_type --reg_vars --holdout")
    elif args.task == 'regression' and args.mode == 'infer' and (args.csv is None or args.reg_type is None):
        parser.error("inferring regression requires --csv -- reg_type")
    elif args.task == 'evaluate' and (args.csv is None):
        parser.error("evaluation requires --csv")
    
    # Set working directory
    # wd = 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python'
    # print('Setting working directory as', 'D:\\valencig\\Thesis\\sentinel-ssc\\sentinel-calibration\\python')
    # os.chdir(wd)
    # Set up subfolders
    folders = ['algs\\','clusters\\','figures\\','regression\\']
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
    
    # args = Namespace(
    #     task = 'regression',
    #     mode= 'infer',
    #     # csv = "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/exports/GEE_raw/transect/transect_harmonized.csv",
    #     csv = "D:/valencig/Thesis/sentinel-ssc/sentinel-calibration/python/clusters/clustered_kMeans_raw_bands.csv",
    #     # cluster_type = "kMeans",
    #     reg_type  = "elasticNet",
    #     reg_vars = 'full_bands',
    #     holdout = 0.
    #     )
    
    # Extract regression variables
    reg_vars = get_regressors(args.reg_vars)
    
    if args.task == "cluster":
        # Extract and standardize data to numpy array
        data, ssc, pred_ssc, scaler, num_clust = standardize(args.csv, reg_vars)
        cluster(data=data, cluster_type=args.cluster_type, csv_path=args.csv, mode=args.mode, scaler=scaler, reg_vars=reg_vars, reg_names=args.reg_vars)
    elif args.task == "regression":
        # Extract but don't standardize data
        data, ssc, pred_ssc, num_clust = data_extract(args.csv, reg_vars)
        regress(data=data, ssc=ssc, num_clust=num_clust, reg_type=args.reg_type, csv_path=args.csv, mode=args.mode, holdout=args.holdout, reg_vars=reg_vars, reg_names=args.reg_vars)
    elif args.task == "evaluate":
        # Extract but don't standardize data
        data, ssc, pred_ssc, num_clust = data_extract(args.csv, reg_vars)
        evaluate(data=data, ssc=ssc, pred_ssc=pred_ssc, num_clust=num_clust, csv_path=args.csv, reg_vars=reg_vars)
    else:
        raise ValueError("!!! Unknown mode !!!")