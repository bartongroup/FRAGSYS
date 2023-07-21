### MAIN PACKAGES ### 

import os
import math
import scipy
import pickle
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

### ML PACKAGES ###

import tensorflow as tf
from tensorflow import keras

from sklearn import metrics
from sklearn.utils import class_weight
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedStratifiedKFold

from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Sequential

from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

### PLOTTING ARGUMENTS

alpha = 0.75

PROPS_pred_KNN = {
    'boxprops':{'facecolor':'firebrick', 'edgecolor':'k', 'alpha' : alpha},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'},
    
}

PROPS_pred_ANN = {
    'boxprops':{'facecolor':'royalblue', 'edgecolor':'k', 'alpha' : alpha},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'},
    
}

PROPS_rand_KNN = {
    'boxprops':{'facecolor':'tomato', 'edgecolor':'k', 'alpha' : alpha},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'},
    
}

PROPS_rand_ANN = {
    'boxprops':{'facecolor':'cornflowerblue', 'edgecolor':'k', 'alpha' : alpha},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'},
    
}

PROPS_rand = {
    'boxprops':{'facecolor':'gold', 'edgecolor':'k', 'alpha' : alpha},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'},
    
}

flierprops_pred_KNN = dict(marker='o', markerfacecolor='firebrick', markersize=5,
                  linestyle='none', color="k", markeredgecolor='k')

flierprops_pred_ANN = dict(marker='o', markerfacecolor='royalblue', markersize=5,
                  linestyle='none', color="k", markeredgecolor='k')

flierprops_rand_KNN = dict(marker='o', markerfacecolor='tomato', markersize=5,
                  linestyle='none', color="k", markeredgecolor='k')

flierprops_rand_ANN = dict(marker='o', markerfacecolor='cornflowerblue', markersize=5,
                  linestyle='none', color="k", markeredgecolor='k')

flierprops_rand = dict(marker='o', markerfacecolor='gold', markersize=5,
                  linestyle='none', color="k", markeredgecolor='k')

### FUNCTIONS ###

def dump_pickle(data, f_out):
    """
    dumps pickle
    """
    with open(f_out, "wb") as f:
        pickle.dump(data, f)
        
def load_pickle(f_in):
    """
    loads pickle
    """
    with open(f_in, "rb") as f:
        data = pickle.load(f)
    return data

def plot_loss(hist_df, epoch_stop = None):
    """
    customised lineplot function to plot
    training loss and validation loss from    
    an artificial neural network training
    """
    n_epochs = len(hist_df)
    ax1 = sns.lineplot(x = list(range(n_epochs)), y = hist_df["loss"], c = "royalblue")
    sns.lineplot(x = list(range(n_epochs)), y = hist_df["val_loss"], c = "orange", ax = ax1)
    if epoch_stop != None:
        plt.axvline(x = epoch_stop, linestyle = "--", color = "red", linewidth = 1.5)

def plot_acc(hist_df, epoch_stop = None):
    """
    customised lineplot function to plot
    training accuracy and validation accuracy from    
    an artificial neural network training
    """
    n_epochs = len(hist_df)
    ax2 = sns.lineplot(x = list(range(n_epochs)), y = hist_df["accuracy"], c = "royalblue")
    sns.lineplot(x = list(range(n_epochs)), y = hist_df["val_accuracy"], c = "orange", ax = ax2)
    if epoch_stop != None:
        plt.axvline(x = epoch_stop, linestyle = "--", color = "red", linewidth = 1.5)

def randomise(df, seed1):
    """
    randomises rows and index of a dataframe
    using two different seeds so the table index 
    corresponds to random rows
    """
    seed2 = seed1 + 1
    df_rand = df.sample(frac = 1, random_state = seed1)
    random.seed(seed2)
    df_rand.index = random.sample(df_rand.index.tolist(), len(df))
    return df_rand

def get_confidences_df(preds_l, round_preds_l, vals_l):
    """
    preds_l: predictions probabilities
    round_preds_l: rounded predictions
    vals_l: real labels
    returns prediction confidence dataframe
    """
    confs = []
    for pred in preds_l:
        for row in pred:
            row = sorted(list(row), reverse = True) # sorts list of probabilities
            confidence = int(10*(row[0] - row[1])) # confidence score: goes from 0-10
            confs.append(confidence)
            
    vals = []
    preds = []
    ids = []
    for i in range(len(vals_l)):
        val = vals_l[i].label.tolist()
        idd = vals_l[i].index.tolist()
        pred = list(round_preds_l[i])
        vals.extend(val)
        preds.extend(pred)
        ids.extend(idd)
        
    conf_df = pd.DataFrame(list(zip(ids, confs, vals, preds)), columns = ["id", "conf", "val", "pred"]) # contains bs_id, real label, predicted label, and associated confidence
    
    un_confs = sorted(conf_df.conf.unique().tolist())
    covs = []
    accs = []
    for un_conf in un_confs:
        confi_df = conf_df[conf_df.conf >= un_conf]
        covs.append(len(confi_df)/len(conf_df))
        p_correct = len(confi_df[confi_df.val == confi_df.pred])/len(confi_df)
        accs.append(p_correct)
    
    conf_df_sum = pd.DataFrame(list(zip(un_confs, covs, accs)), columns = ["conf", "cov", "acc"]) # cumulative confidence dataframe with conf score, coverage, and average accuracy
    
    return conf_df, conf_df_sum

def plot_conf_acc_cov(conf_df_sum, max_xtick = 10, f_size = (5, 5), dpi = 100, out = None):
    """
    plots confidence vs coverage and accuracy
    for the cross-validation of a neural network
    """
    # create figure and axis objects with subplots()
    fig, ax = plt.subplots(figsize = f_size, dpi = dpi)
    # make a plot
    ax.plot(conf_df_sum["conf"],
            conf_df_sum["cov"],
            color="darkgreen", 
            marker="o",
           label = "Coverage", markeredgewidth = 1, markeredgecolor = "k")
    # set x-axis label
    ax.set_xlabel("Confidence score", fontsize = 14)
    # set y-axis label
    ax.set_ylabel("Coverage",
                  color="darkgreen",
                  fontsize=14)
    
    #ax.set_yticks(np.arange(0.3,1.1, 0.1))
    #ax.set_ylim(0.339, 1.01)
    #plt.legend()
    # twin object for two different y-axis on the sample plot
    ax2=ax.twinx()
    # make a plot with different y-axis using second axis object
    ax2.plot(conf_df_sum["conf"], conf_df_sum["acc"],color="purple",marker="D", label = "Average accuracy", markeredgewidth = 1, markeredgecolor = "k")
    ax2.set_ylabel("Average accuracy",color="purple",fontsize=14)
    ax2.set_xlabel("Average confidence",color="k",fontsize=14)
    #ax2.set_yticks(np.arange(0.90, 1.02, 0.02))
    #ax2.set_ylim(0.89, 1)
    #plt.legend()
    plt.xticks(range(0, max_xtick))
    #plt.xlim(-0.5, 9.5)
    if out != None:
        plt.savefig(out)
    plt.show()

def get_confidences_df_blind(preds, round_preds, vals):
    """
    plots confidence vs coverage and accuracy
    for the test set prediction of a neural network
    """
    confs = []
    for row in preds:
        row = sorted(list(row), reverse = True)
        confidence = int(10*(row[0] - row[1]))
        confs.append(confidence)
        
    conf_df = pd.DataFrame(list(zip(confs, list(vals), list(round_preds))), columns = ["conf", "val", "pred"])
    
    un_confs = sorted(conf_df.conf.unique().tolist())
    covs = []
    accs = []
    for un_conf in un_confs:
        confi_df = conf_df[conf_df.conf >= un_conf]
        covs.append(len(confi_df)/len(conf_df))
        p_correct = len(confi_df[confi_df.val == confi_df.pred])/len(confi_df)
        accs.append(p_correct)
    
    conf_df_sum = pd.DataFrame(list(zip(un_confs, covs, accs)), columns = ["conf", "cov", "acc"])
    
    return conf_df, conf_df_sum

### the end ###