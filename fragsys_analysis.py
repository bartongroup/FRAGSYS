### THE PACKAGES ###
import os
import copy
import math
import time
import scipy
import pickle
import random
import sklearn
import itertools
import matplotlib
import statistics
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from scipy import stats
from color_library import *
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from IPython.display import display
#from statannotations.Annotator import Annotator seaborn version inconsistency

### THE VARIABLES ###

sample_colors = list(itertools.islice(rgbs(), 200)) # list of colours

### THE FUNCTIONS ###

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

def plot_hist(data, feat_name, bns, col, show = True, out = None, dpi = 100, figsize = (5, 5), yticks = None):
    """
    customised histogram plotting function
    """
    if type(data) != list:
        data = data.tolist()

    st_median = round(statistics.median(data), 2)
    st_mean = round(statistics.mean(data), 2)
    st_max = round(max(data),2)
    st_min = round(min(data), 2)
    print("\tMIN = {}\tMEAN = {}\tMEDIAN = {}\tMAX = {}".format(st_min, st_mean, st_median, st_max))
    
    fig = plt.figure(figsize = figsize, dpi = dpi)
    sns.histplot(data, fill = True, color = col, stat = "count", element = "step", bins = bns)
    plt.ylabel("Count")
    plt.xlabel(feat_name)
    plt.axvline(x = st_median, linewidth = 1, color = "k", linestyle = "--")
    
    if yticks != None:
        plt.yticks(yticks)
    if out != None:
        plt.savefig(out)
    if show == False:
        plt.close()

def pearsonr_ci(x, y, alpha = 0.05):
    """
    calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    """

    r, p = stats.pearsonr(x,y)
    #r, p = stats.spearmanr(x, y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    
    print(
        "r = {}, p = {}, 95% CI = [{}, {}]".format(
            round(r, 2), round(p, 4), round(lo, 2), round(hi, 2)
        )
    )
    
    return r, p, lo, hi

def ezy_scatter(x, y, x_lab, y_lab, alpha = 0.05):
    """
    customisted scatter plotting function
    """
    m, b = np.polyfit(x, y, 1)
    r, p, lo, hi = pearsonr_ci(x, y, alpha)
    plt.scatter(x, y, linewidth = 0.5, edgecolor = "k")
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.plot(x, m*x+b, color = "orange")
    plt.show()
    
def bin_data_points(df, col_name, bin_list):
    """
    creates binned values list according to
    a set of values and predetermined bins
    """
    bins_list = []
    for i in range(len(df)):
        b = 0
        col_value = df.loc[i, col_name]
        for binn in bin_list:
            if col_value in range(binn[0], binn[1]):
                bins_list.append(binn)
                b = 1
                break
            else:
                continue
        if b == 0:
            print(i, bs_size)

    return bins_list

def regression_coordinates(df, xcol, ycol, alpha = 0.05):
    """
    obtains linear regression coordinats
    """
    x = df[xcol]
    y = df[ycol]
    m, b = np.polyfit(x, y, 1)
    return m*x+b

def plot_distrs(data, feat_name, labels, bns, cols, show = True, out = None, dpi = 100, figsize = (5, 5)):
    """
    customised histogram distribution plotting function
    """
    fig = plt.figure(figsize = figsize, dpi = dpi)
    sns.histplot(data[0], label = labels[0], fill = True, color = cols[0], stat = "proportion", element = "step", bins = bns)
    sns.histplot(data[1], label = labels[1], fill = True, color = cols[1], stat = "proportion", element = "step", bins = bns)
    plt.ylabel("p")
    plt.xlabel(feat_name)
    plt.legend(edgecolor = "k")
    if out != None:
        plt.savefig(out)
    if show == False:
        plt.close()
        
def run_tests(d):
    """
    runs a series of statistical tests to compare two distributions
    """
    a, b = d
    A, cvs, p = scipy.stats.anderson_ksamp([a, b])
    
    if A < cvs[0]:
        print("Anderson-Darling 2 sample test:\nA statistic = {}\tp > {}\n".format(round(A,2), round(p,3)))
    elif A > cvs[-1]:
        print("Anderson-Darling 2 sample test:\nA statistic = {}\tp < {}\n".format(round(A,2), round(p,3)))
    else:
        print("Anderson-Darling 2 sample test:\nA statistic = {}\tp = {}\n".format(round(A,2), round(p,3)))

    s, p = scipy.stats.ks_2samp(a, b, alternative="two-sided", mode = "exact")
    print("Kolmogorov-Smirnov two sample test result:\nKS statistic = {}\tp = {}\n".format(round(s,2), round(p,2)))

    s, p = scipy.stats.mannwhitneyu(a, b, alternative="two-sided")
    print("Mann-Whitney U rank test result:\nU = {}\tp = {}\n".format(round(s,2), round(p,2)))

    s, p = scipy.stats.ttest_ind(a, b)
    print("T-test result:\nT = {}\tp = {}".format(round(s,2), round(p,2)))

def get_rsa_profiles(bss_df, ress_df):
    """
    Returns a dictionary containint the binding site IDs as keys,
    and the site RSA profile as value. This is a list containing
    the RSA values of the residues forming the site.
    """
    prots = bss_df.protein.unique().tolist()
    rsa_profs = {}
    for prot in prots:
        prot_groups = sorted(bss_df.query('protein == @prot').group.unique().tolist())
        for group in prot_groups:
            prot_group_ress = ress_df.query('protein == @prot & group == @group')
            if len(prot_group_ress) == 0:
                print("Group {} of {} has 0 residues. Skipping!".format(group, prot))
                continue
            prot_bs_ids = sorted(bss_df.query('protein == @prot & group == @group').bs_id.unique().tolist())
            for prot_bs_id in prot_bs_ids:
                prot_bs_ress = prot_group_ress[prot_group_ress[prot_bs_id] == 1]
                if len(prot_bs_ress) == 0:
                    print("0 res at {} of group {} of {}".format(prot_bs_id, group, prot))
                    continue
                bs_rsas = prot_bs_ress.RSA.tolist()
                dk = "{}_{}_{}".format(prot, group, prot_bs_id)
                rsa_profs[dk] = sorted([round(el, 1) for el in bs_rsas])
    return rsa_profs

def get_U(a, b):
    """
    calculates U, Umax, Urel and associated p-value
    for a pair of RSA profiles
    """
    U, p = stats.mannwhitneyu(a, b, use_continuity = False)
    Umax = (len(a)*len(b))/2
    Urel = U/Umax
    return U, Umax, Urel, p

def get_UD_matrix(rsa_profs):
    """
    calculates the U distance matrix
    """
    ks = list(rsa_profs.keys())
    Urels = {i: {} for i in ks}
    for i in range(len(ks)):
        a = rsa_profs[ks[i]]
        U, Umax, Urel, p = get_U(a, a)
        Urels[ks[i]][ks[i]] = Urel
        for j in range(i+1, len(ks)):
            try:
                b = rsa_profs[ks[j]]
                U, Umax, Urel, p = get_U(a, b)
                Urels[ks[i]][ks[j]] = Urel
                Urels[ks[j]][ks[i]] = Urel
            except:
                print(a, b)
    Urels_df = pd.DataFrame(Urels)
    UD_df = abs(1 - Urels_df)
    return UD_df

def find_rsa_t_in_bs(rsa_prof, rsa_t = 25):
    """
    finds position where to draw the line
    on the vector representation of a binding
    site profile for a given RSA threshold
    """
    for pos, i in enumerate(rsa_prof):
        if i < rsa_t:
            continue
        else:
            return pos + 1
    return pos + 1

def plot_color_dd(dist_mat, clust_method, dist_t, sample_colors, out = None, show = True):
    """
    plots a dendrogram coloured according to the cluster labels
    generated by the linkage
    """
    z_c = scipy.cluster.hierarchy.linkage(
        scipy.spatial.distance.squareform(dist_mat),
        method = clust_method
    )

    ks = dist_mat.index.tolist()
    
    cutree = scipy.cluster.hierarchy.cut_tree(z_c, height = dist_t) # clusters coordinates
    cluster_ids = [int(cut) for cut in cutree]
    site_cluster_dict = {ks[i]: cluster_ids[i] for i in range(len(ks))}

    plt.rcParams['lines.linewidth'] = 3

    fig = plt.figure(figsize=(40, 15))
    plt.axhline(dist_t, linestyle = "--", linewidth = 3, color = "black")   
    dd = scipy.cluster.hierarchy.dendrogram(z_c, no_labels = True, labels = ks, leaf_font_size = 10, color_threshold = dist_t, above_threshold_color = "black")
    plt.xlabel("Binding site")
    plt.ylabel("UD")

    dd_labs = dd["ivl"]
    dd_c_labs = [site_cluster_dict[i] for i in dd_labs]
    ordered_labs = list(dict(zip(dd_c_labs,dd_c_labs)))
    n_clusters = len(ordered_labs)
    scipy.cluster.hierarchy.set_link_color_palette([matplotlib.colors.to_hex(sample_colors[i]) for i in ordered_labs])
    if out!= None:
        plt.savefig(out)
    if show == True:
        plt.show()
    return site_cluster_dict, cluster_ids, n_clusters

def process_heatmap_df(rsas_profs_hmap, rsas_profs_lens, cluster_dict, K):
    """
    does some data processing to prepare a simple dataframe
    to be plotted as we need to represent the ligand
    binding sites RSA profiles
    """
    heatmap_df = pd.DataFrame(rsas_profs_hmap).T
    sorted_dfs = []
    for c in range(0, K):
        c_labs = [k for k, v in cluster_dict.items() if v == c]
        sorted_c_keys = list({k: v for k, v in sorted({k2:v2 for k2, v2 in rsas_profs_lens.items() if cluster_dict[k2] == c}.items(), key=lambda item: item[1])}.keys())
        heatmap_df_c = heatmap_df.loc[c_labs,:]

        sorterIndex = dict(zip(sorted_c_keys, range(len(sorted_c_keys))))
        heatmap_df_c['Tm_Rank'] = heatmap_df_c.index.map(sorterIndex)
        heatmap_df_c.sort_values(['Tm_Rank'],
                ascending = [True], inplace = True)
        heatmap_df_c.drop(axis = 1, columns = 'Tm_Rank', inplace = True)
        sorted_dfs.append(heatmap_df_c)
    sorted_len_heatmap_df = pd.concat(sorted_dfs)
    return sorted_len_heatmap_df

def get_rsa_t_poss(sorted_len_heatmap_df, rsa_t = 33.33):
    """
    gets axis coordinates to draw RSA threshold line
    """
    rsa_t_poss = []
    for bs_id in sorted_len_heatmap_df.index.tolist():
        rsa_t_pos = find_rsa_t_in_bs(sorted_len_heatmap_df.loc[bs_id,:].dropna().tolist(), rsa_t)
        rsa_t_poss.append(rsa_t_pos)
        
    return rsa_t_poss

def plot_heatmaps(sorted_len_heatmap_df, rsa_t_poss, n_clusters, nks, fsize = (10, 30), dpi = 300, cmap = "cividis", out = None, show = False, rsa_t_lw = 5):
    """
    does the plotting itself of the ligand
    binding site RSA-based clusters
    """
    plt.figure(figsize = fsize, dpi = dpi)
    ax = sns.heatmap(sorted_len_heatmap_df, cmap = cmap, vmin = 0, vmax = 100, linewidth = 0, linecolor = "k", yticklabels = False)
    ax.set_xticks(np.arange(0, 45, 5))
    ax.set_xticklabels(np.arange(0, 45, 5), rotation = -90)
    ax.set_xlim(0, 45)
    ax.xaxis.set_ticks_position('top') # the rest is the same
    
    for i in range(n_clusters):
        if i == 0:
            start_index = 0
            end_index = nks[i]
        else:
            start_index = end_index
            end_index += nks[i]

        plt.plot(
            rsa_t_poss[start_index:end_index], np.arange(len(sorted_len_heatmap_df))[start_index:end_index] + 0.5,
            color = sample_colors[i], lw = rsa_t_lw
        )
        
    if out != None:
        plt.savefig(out)
        
    if show == False:
        plt.close()
        
def plot_clusters(rsas_profs_hmap, rsas_profs_lens, cluster_dict, nks, n_clusters, rsa_t = 33.33, fsize = (10, 30), dpi = 300, cmap = "cividis", out = None, show = False, rsa_t_lw = 5):
    """
    this functions does the data preprocessing
    necessary to plot the RSA-based ligand binding
    site clusters and does the plotting itself as well
    """
    
    sorted_len_heatmap_df = process_heatmap_df(rsas_profs_hmap, rsas_profs_lens, cluster_dict, n_clusters)
    
    rsa_t_poss = get_rsa_t_poss(sorted_len_heatmap_df, rsa_t)
    
    plot_heatmaps(sorted_len_heatmap_df, rsa_t_poss, n_clusters, nks, fsize, dpi, cmap, out, show, rsa_t_lw)

def plot_boxes(df, x, y, order, pairs, col_palette, flierprops, fsize = (5, 5), dpi = 100, annotate = False, ann_loc = "inside", out = None, show = True):
    """
    customised boxplot plotting function
    """
    plt.figure(figsize = fsize, dpi = dpi)
    
    ax = sns.boxplot(data = df, x = x, y = y, order = order, palette = col_palette, flierprops = flierprops)

    if annotate == True:
        annotator = Annotator(ax, pairs, data = df, x = x, y = y, order = order)
        annotator.configure(test = "Mann-Whitney", text_format = "star", loc = ann_loc)
        annotator.apply_and_annotate()
        
    if show == False:
        plt.close()

def plot_reg_boxes(x_reg, y_reg, x_box, y_box, binns, xlab, ylab, f_size = (7.5, 7.5), dpi = 100, out = None, show = True):
    """
    customised scatter + box + regression line plotting function
    """
    
    PROPS = {
    'boxprops':{'facecolor':'none', 'edgecolor':'k'},
    'medianprops':{'color':'k'},
    'whiskerprops':{'color':'k'},
    'capprops':{'color':'k'}
    }
    
    fig = plt.figure(figsize = f_size, dpi = dpi)
    sns.regplot(x = x_reg, y = y_reg, color = "orange", scatter = False)
    sns.swarmplot(x = x_box, y = y_box, order = binns, color = "royalblue", linewidth = 1, edgecolor = "k")
    ax = sns.boxplot(x = x_box, y = y_box, order = binns, palette = ["w"], **PROPS)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    
    if out != None:
        plt.savefig(out)

    if show == False:
        plt.close()

def get_OR(df, idx, t_row):
    """
    calculates enrichment in functional
    sites per RSA cluster
    """
    for i_row in idx:
        i_kf = df.loc[i_row, "kf"]
        i_tot = df.loc[i_row, "tot"]
        rest_kf = df.kf.sum() - i_kf
        rest_tot = df.tot.sum() - i_tot
        oddsr, pval = stats.fisher_exact([[i_kf, rest_kf], [i_tot, rest_tot]])
        vals = [i_kf, rest_kf, i_tot, rest_tot]
        print(i_row, vals)
        se_logor = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x), vals)))))
        logor = math.log(oddsr)
        df.loc[i_row, "oddsratio"] = round(oddsr, 2)
        df.loc[i_row, "log_oddsratio"] = round(logor, 2)
        df.loc[i_row, "pvalue"] = round(pval, 2)
        df.loc[i_row, "ci_dist"] = round(se_logor, 2)
    return df

### THE END ###