#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm, colors # colormap
import logging, ete3, argparse
import numpy as np, pandas as pd, seaborn as sns
from sklearn import manifold, metrics, cluster, neighbors, decomposition, preprocessing
#import skbio, parasail, dendropy, datetime, time, codecs, joypy
import sys, gzip, bz2, re, glob, pickle, collections, subprocess, os, errno, random, itertools, pathlib

logger = logging.getLogger(__name__)
#logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

def generate_time_heatmap (df0, date_col = None, group_col = None, use_max = True, label_interval = None):
    if date_col is None: date_col = "collection_date"
    if group_col is None: group_col = "lineage"
    if label_interval is None: label_interval = 7 #  once every week
    df = df0
    df["date"] = pd.to_datetime(df["collection_date"], infer_datetime_format=False, errors='coerce')
    ## see also df.pivot(a,b,c) which takes df with cols a,b,c and create axb matrix with c values
    df = df.groupby(["date","adm2"]).size().unstack() ## real one will use cluster_id, not adm 

    idx = pd.date_range(df.index.min(),df.index.max(), freq="1D") # creates uniform interval
    df.fillna(0,inplace=True) ## otherwise nan is treated differently
    df = df.reindex(idx, fill_value = 0) # does not accept inplace
    ## df.asfreq('D') # simpler alternative to idx[] above, to interpolate df with regular ("D"aily) intervals

    recent = {}
    if use_max: # sort by most recent case
        for column in df.columns:
            res = max (((df.index - df.index.values.min())/peroba.np.timedelta64(1, 'D')).astype(int) * (df[column].astype(peroba.np.int64())>0))
            recent[column] = res
    else: # sort by weighted average of dates (weights are number of cases)
        for column in df.columns:
            x = ((df.index - df.index.values.min())/peroba.np.timedelta64(1, 'D')).astype(int)
            res = sum(x * df[column].astype(peroba.np.int64()))/sum(df[column].astype(peroba.np.int64()))
            recent[column] = res

    reclist = [k for k, v in sorted(recent.items(), key=lambda item: item[1], reverse=True)]
    labs = ["" if i%label_interval else a for i,a in enumerate(df.index.strftime('%Y-%m-%d'))]
    return df[reclist], labs

def new_date_column_with_jitter (df0, original_date_col = None, new_date_col = None, label_interval = None):
    """ this function can be merged with above? or we can recalc here the column order, using most recent dates:
    we can sort the float dates and sum first 3 or 5...
    """ 

    if new_date_col is None: new_date_col = "float_date"
    if original_date_col is None: original_date_col = "collection_date"
    if label_interval is None: label_interval = 7 #  once every week
    tmp_col = new_date_col + "TMP"
    df[tmp_col] = pd.to_datetime(df0[original_date_col], infer_datetime_format=False, errors='coerce')
    df[new_date_col] = ((df[tmp_col] - df[tmp_col].min())/peroba.np.timedelta64(1, 'D')).astype(float) +  np.random.uniform(-0.35,0.35, len(df[tmp_col]))
    # below is WRONG since the x-axis does NOT follow table order!
    labs = ["" if i%label_interval else a for i,a in enumerate(df[tmp_col].dt.strftime('%Y-%m-%d'))] # notice the "dt"to convert
    return df[new_date_col], labs

def plotper_group_bubble():
    df=subm[6]
    df = df.groupby(["collection_datetime","adm2"]).size().to_frame(name="size")
    df = df.reset_index() # collection and adm2 are indices of series (that we transformed into frame)
    display(df)
    peroba.sns.set(font_scale=1); peroba.sns.set_palette("cubehelix", 8) 
    g = peroba.sns.scatterplot(y="adm2", x="collection_datetime", size="size", edgecolor="white", alpha=0.8, s=10, data=df)
    #g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='right') # nt working
     plt.xticks(rotation=45, horizontalalignment='right') # alternative to loop above (pyplot only)


## plot the heatmap with 
# df2 = df2.T ## transpose rows and cols
#g = sns.heatmap( df2, mask=df2.isnull(), square=False, cmap= 'Blues', linewidth=0.5
#    cbar_kws={'fraction' : 0.005, "ticks":[0,df.values.max()]})# fraction= shrink colorbar size, ticks=text,
## alternatives to labs[] above: to make seaborn decide number of tick labels
#x = g.set_xticklabels(df.index.strftime('%Y-%m-%d')) 
##  or, to have one tick per week: (strftime above can be removed but then seaborn may decide to drop some)
#x = g.get_xticklabels()
#for i, lab in enumerate(x):
#    if i%7:
#        lab.set_text("")
#    else:
#        lab.set_text(lab.get_text()[:10])
#x = g.set_xticklabels (x, rotation=30, horizontalalignment='right', fontweight='light')

## plot using the jitter (Although not reordered, can use order from heatmap...)
# g = peroba.sns.stripplot(y="lineage", x="date_float", edgecolor="black", jitter=0.3, linewidth=0.1, marker="o", alpha=0.8, s=2, data=df2)

def plot_genomes_sequenced_over_time (metadata, output_dir):
    df = metadata.copy() 
    logger.debug("counter: ", str(collections.Counter(df["submission_org_code"]).most_common(25)))
    ## *similar* (not equal) to df.pivot(a,b,c) that takes df with cols a,b,c and create axb matrix with c values
    df = df[df["submission_org_code"].str.contains("NORW", na=False)] 
    #df = df.groupby(["date","adm2"]).size().to_frame("size").reset_index()  #.unstack()
    df["region"] = df["adm2"].map({"NORFOLK":"Norfolk", "Norfolk":"Norfolk"}).fillna("others") # maps only Norfolk; other rows become NaN (which is filled by "others")

    plt.figure(figsize=(10,6)); sns.set()
    #rcParams['figure.figsize'] = 12,4
    sns.set_context("paper", rc={"figure.dpi":300, "font.size":8,
                                "axes.titlesize":8,"axes.labelsize":8, 
                                "xtick.labelsize":6, "ytick.labelsize":6})  
    sns.set_palette("cubehelix", 3)
    g = sns.countplot(x="collection_date", hue="region", data=df)
    g.set(title="Genomes sequenced at QIB over time", xlabel="", ylabel="Number of samples")
    leg = g.axes.legend()
    leg.set(title="sampling region")
    plt.setp(leg.get_texts(), fontweight='light', fontsize=8)
    plt.setp(leg.get_title(),  fontweight='normal', fontsize=8)
    x = g.get_xticklabels()
    for lab in x:
        lab.set_text(lab.get_text()[:10])
    x = g.set_xticklabels (x, rotation=30, horizontalalignment='right', fontweight='light')
    
    md_description = """
## Genomes sequenced at the QIB considered here PLEASE DO NOT USE 
These counts are **not** the total sums, they are based on the merged database (which remove data without sequence etc.)
"""
    md_description += "\n<br>![](genomes_over_time.pdf)\n<br>\n" # same directory as report.md, it doesn't need full path? 
    fname = os.path.join(output_dir,"genomes_over_time.pdf")
    g.figure.savefig(fname, format="pdf")  # or plt.savefig()
    g.figure.clf()
    return md_description 

def plot_jitter_lineages (metadata, output_dir=None):
    if output_dir is None: output_dir = cwd
    df = metadata.copy() 
    df["date_float"] = ((df["collection_date"] - df["collection_date"].min())/np.timedelta64(1, 'D')).astype(float) + np.random.normal(0,0.1, len(df["collection_date"]))
    df = df[df["adm2"].str.contains("folk",na=False)]
    print (df["date_float"])
    
    plt.figure(figsize=(10,6)); sns.set()
    sns.set_context("paper", rc={"figure.dpi":300, "font.size":8,
                                "axes.titlesize":8,"axes.labelsize":8, 
                                "xtick.labelsize":6, "ytick.labelsize":6})  
    #df2["date"].strftime('%Y-%m-%d')
    g = sns.stripplot (y="lineage", x="date_float", edgecolor="black", jitter=0.3, linewidth=0.1, marker="o", alpha=0.8, s=2, data=df)
    #x = g.set_xticklabels(df2["date"].dt.strftime('%Y-%m-%d'))  ## order in table is NOT PLOT ORDER!
    #x = g.get_xticklabels ()
    #for lab in x:
    #    print (lab)
    #    lab.set_text(float_to_date[lab.get_text()])
    #x = g.set_xticklabels ([1,2,3, 30])
    ## see also https://www.machinelearningplus.com/plots/top-50-matplotlib-visualizations-the-master-plots-python/
    md_description = """
## PLEASE DO NOT USE 
the days are a float number 
"""
    md_description += "\n<br>![](jitter_lineages.pdf)\n<br>\n" # same directory as report.md, it doesn't need full path? 
    fname = os.path.join(output_dir,"jitter_lineages.pdf")
    g.figure.savefig(fname, format="pdf")
    g.figure.clf()
    return md_description 

