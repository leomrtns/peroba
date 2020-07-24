#!/usr/bin/env python
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm, colors, patches
from mpl_toolkits.basemap import Basemap  # patches.Polygon, colors.Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection ## import collections may clash with collections.Counter one 
import warnings, shapely.geometry ## patchcollection doesnt understand geopandas geometry

import logging, ete3, argparse
import numpy as np, pandas as pd, seaborn as sns, geopandas as gpd
from sklearn import manifold, metrics, cluster, neighbors, decomposition, preprocessing
import sys, gzip, bz2, re, glob, pickle, collections, subprocess, os, errno, random, itertools, pathlib, datetime

logger = logging.getLogger(__name__)
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)
    
seaborn_rc = {"figure.dpi":300, "font.size":8, "axes.titlesize":8,"axes.labelsize":8, "xtick.labelsize":6, "ytick.labelsize":6}  
# basemap needs {dbf,fix,shp,prj, shx} files so suffix excluded, while geopandas needs just the shp file
uk_postcode_area_file = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/Districts") 
uk_gadm_nuts2_file = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/gadm36_GBR_2") 


def generate_time_heatmap (df0, date_col = None, group_col = None, use_max = True, label_interval = None):
    if date_col is None: date_col = "collection_date"
    if group_col is None: group_col = "peroba_uk_lineage"
    df = df0.copy()
    df["group"] = df[group_col].str.split('/').str[0]
#    df["date"] = pd.to_datetime(df[date_col], infer_datetime_format=True, yearfirst=True, errors='coerce')
    df["date"] = df[date_col]
    df.dropna(subset=["date"], inplace=True)
    ## see also df.pivot(a,b,c) which takes df with cols a,b,c and create axb matrix with c values
    df = df.groupby(["date","group"]).size().unstack() ## real one will use cluster_id, not adm 

    idx = pd.date_range(df.index.min(),df.index.max(), freq="1D") # creates uniform interval
    df.fillna(0,inplace=True) ## otherwise nan is treated differently
    df = df.reindex(idx, fill_value = 0) # does not accept inplace
    ## df.asfreq('D') # simpler alternative to idx[] above, to interpolate df with regular ("D"aily) intervals

    recent = {}
    if use_max: # sort by most recent case
        for column in df.columns:
            res = max (((df.index - df.index.values.min())/np.timedelta64(1, 'D')).astype(int) * (df[column].astype('int64')>0))
            recent[column] = res
    else: # sort by weighted average of dates (weights are number of cases)
        for column in df.columns:
            x = ((df.index - df.index.values.min())/np.timedelta64(1, 'D')).astype('int64')
            res = sum(x * df[column].astype(np.int64()))/sum(df[column].astype(np.int64()))
            recent[column] = res

    reclist = [k for k, v in sorted(recent.items(), key=lambda item: item[1], reverse=True)]
    if label_interval is None: label_interval = int (len(idx)/14)
    if label_interval < 1:     label_interval = 1
    labs = ["" if i%label_interval else a for i,a in enumerate(df.index.strftime('%Y-%m-%d'))]
    return df[reclist], labs

def plot_time_heatmap_seaborn (metadata, counter, output_dir):
    fname = f"heat{counter}.pdf"
    df2, labs = generate_time_heatmap (metadata)
    df2 = df2.T ## transpose rows and cols
    ratio = df2.shape[0]/df2.shape[1]
    if ratio > 4: ratio = 4
    if ratio < 0.25: ratio = 0.25
    #plt.figure(figsize=(dyn_width,dyn_height)); 
    plt.figure(figsize=(8, ratio * 8)); 
    sns.set(font_scale=1); sns.set_context("paper", rc=seaborn_rc);
    sns.set_palette("cubehelix", 8)
    g = sns.heatmap( df2, mask=df2.isnull(), square=False, cmap= 'Blues', linewidth=0.5, 
            cbar_kws={'fraction' : 0.005, "ticks":[0,df2.values.max()]})# fraction= shrink colorbar size, ticks=text,
    # alternatives to labs[] above: to make seaborn decide number of tick labels
    #x = g.set_xticklabels(df.index.strftime('%Y-%m-%d')) 
    #  or, to have one tick per week: (strftime above can be removed but then seaborn may decide to drop some)
    #x = g.get_xticklabels()
    #for i, lab in enumerate(x):
    #    if i%7:
    #        lab.set_text("")
    #    else:
    #        lab.set_text(lab.get_text()[:10])
    #x = g.get_yticklabels()
    #for i, lab in enumerate(x):
    #        lab.set_text(lab.get_text()[:16])
    x = g.set_xticklabels (labs, rotation=30, horizontalalignment='right', fontweight='light') ## x is to avoid printing
    g.figure.savefig(os.path.join(output_dir,fname), format="pdf")  # or plt.savefig()
    g.figure.clf()

    caption = f"Number of samples of each UK_lineage over time for cluster {counter}"
    md_description = f"\n![heat{counter}]({fname})*{caption}*\n<br>\n"

    return md_description

#def sum_heatmap_by_week (x, col_names):
    #x.reshape(-1, 4, 3).sum(axis=2)

def plot_time_heatmap (metadata, counter, output_dir, figdir, rlim = 80, clim = 120):
    df2, labs = generate_time_heatmap (metadata)
    df2 = df2.T ## transpose rows and cols

    row_names = df2.index
    #col_names = pd.to_datetime(df2.columns, infer_datetime_format=False, errors='coerce').strftime('%Y-%m-%d')
    col_names = labs
    df2 = df2.to_numpy()
    nrows = df2.shape[0];
    ncols = df2.shape[1];
    if nrows > rlim:
        df2 = df2[:rlim,:]  ## most recent at the top
        row_names = row_names[:rlim]
        nrows = rlim
    if ncols > clim:
        df2 = df2[:,-clim:]  ## most recent
        col_names = col_names[-clim:]
        ncols = clim

    ratio = nrows/ncols
    if ratio > 1.3: ratio = 1.3
    if ratio < 0.2: ratio = 0.2
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(12, 10 * ratio)) # I'm leaving ax in case we want several plots
    extent = (0, ncols, nrows, 0) # x0, x1, y0, y1
    im = ax.imshow(df2, aspect="auto", cmap="Blues", extent=extent, interpolation="nearest")
    cax = make_axes_locatable(ax).append_axes("right", size="2%", pad=0.03); cbar = plt.colorbar(im, cax=cax)
    cbar.ax.set_ylabel("counts", rotation=-90, va="bottom", fontsize=12)
    #ax.get_xticklabels().set_fontsize(14)

    ax.set_xticks(np.arange(0.5, ncols, 1)); # where labels go
    ax.set_yticks(np.arange(0.5, nrows, 1));
    ax.set_xticks(np.arange(ncols), minor=True); # just for break lines
    ax.set_yticks(np.arange(nrows), minor=True);
    ax.set_xticklabels(col_names)
    ax.set_yticklabels(row_names)
    ax.grid(which='minor', color='w', linestyle='-', linewidth=1)
    ax.grid(which='major', visible = False)
    if nrows < 5: yfsize = 18
    elif nrows < 12: yfsize = 12
    elif nrows < 24: yfsize = 8
    else: yfsize = 6
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=yfsize)

    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")
    #ax.set_title (f"UK lineages over time for Cluster {counter}", fontsize=20)
    fig.tight_layout()
    caption = f"Number of samples from each UK\_lineage over time for cluster {counter}. The colour intensity \
indicates the number of samples (from white to blue). Plot truncated to 80 most recent UK lineages and last 4 months"

    fname = f"{figdir}/heat{counter}.pdf"
    fig.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(100))  # or plt.savefig()
    pdf_desc = f"\n![{caption}]({fname})\n\n"
    fname = f"{figdir}/heat{counter}.png"
    fig.savefig(os.path.join(output_dir,fname), format="png", dpi=(100))  # or plt.savefig()
    html_desc = f"\n\n![{caption}]({fname})\n\n<br>"
    #   control size in pandocmd: "![](file.jpg){ width=50% }"

    return html_desc, pdf_desc 

def new_date_column_with_jitter (df0, original_date_col = None, new_date_col = None, label_interval = None):
    """ this function can be merged with above? or we can recalc here the column order, using most recent dates:
    we can sort the float dates and sum first 3 or 5...
    """ 

    if new_date_col is None: new_date_col = "float_date"
    if original_date_col is None: original_date_col = "collection_date"
    if label_interval is None: label_interval = 7 #  once every week
    tmp_col = new_date_col + "TMP"
    #df[tmp_col] = pd.to_datetime(df0[original_date_col], infer_datetime_format=True, yearfirst=True, errors='coerce')
    df[tmp_col] = df0[original_date_col]
    df[new_date_col] = ((df[tmp_col] - df[tmp_col].min())/np.timedelta64(1, 'D')).astype(float) +  np.random.uniform(-0.35,0.35, len(df[tmp_col]))
    # below is WRONG since the x-axis does NOT follow table order!
    labs = ["" if i%label_interval else a for i,a in enumerate(df[tmp_col].dt.strftime('%Y-%m-%d'))] # notice the "dt"to convert
    return df[new_date_col], labs

## plot using the jitter (Although not reordered, can use order from heatmap...)
# g = sns.stripplot(y="lineage", x="date_float", edgecolor="black", jitter=0.3, linewidth=0.1, marker="o", alpha=0.8, s=2, data=df2)

def plot_genomes_sequenced_over_time (metadata, output_dir, figdir):
    df = metadata.copy() 
    logger.debug("counter: ", str(collections.Counter(df["submission_org_code"]).most_common(25)))
    ## *similar* (not equal) to df.pivot(a,b,c) that takes df with cols a,b,c and create axb matrix with c values
    df = df[df["submission_org_code"].str.contains("NORW", na=False)] 
    #df = df.groupby(["date","adm2"]).size().to_frame("size").reset_index()  #.unstack()
    df["region"] = df["adm2"].map({"NORFOLK":"Norfolk", "Norfolk":"Norfolk"}).fillna("others") # maps only Norfolk; other rows become NaN (which is filled by "others")

    plt.figure(figsize=(12,8)); sns.set(); sns.set_context("paper", rc=seaborn_rc);
    #rcParams['figure.figsize'] = 12,4
    sns.set_palette("cubehelix", 3)
    g = sns.countplot(x="collection_date", hue="region", data=df)
    g.set(title="Genomes sequenced at QIB over time", xlabel="", ylabel="Number of samples")
    leg = g.axes.legend()
    leg.set(title="sampling region")
    plt.setp(leg.get_texts(), fontweight='light', fontsize=12)
    plt.setp(leg.get_title(),  fontweight='normal', fontsize=12)
    x = g.get_xticklabels()
    for i,lab in enumerate(x):
        if i%2:
            lab.set_text("")
        else:
            lab.set_text(lab.get_text()[:10])
    x = g.set_xticklabels (x, rotation=30, horizontalalignment='right', fontweight='light')
    
    desc = """
## Genomes sequenced at the QIB considered here 
These counts are **not** the total sums, since they are based on sequenced genomes of high phylogenetic quality\n<br> 
"""
    fname = f"{figdir}/genomes_over_time.pdf"
    g.figure.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(100))
    pdf_desc = f"{desc}![]({fname})\n\n"
    fname = f"{figdir}/genomes_over_time.png"
    g.figure.savefig(os.path.join(output_dir,fname), format="png", dpi=(100))
    html_desc = f"\n\n{desc}![]({fname})\n\n<br>"
    g.figure.clf()
    return html_desc, pdf_desc 

def plot_jitter_lineages_seaborn (metadata, output_dir=None):
    if output_dir is None: output_dir = cwd
    df = metadata.copy() 
    df["date_float"] = ((df["collection_date"] - df["collection_date"].min())/np.timedelta64(1, 'D')).astype(float) + np.random.normal(0,0.1, len(df["collection_date"]))
    df = df[df["adm2"].str.contains("folk",na=False)]
    plt.figure(figsize=(10,6)); sns.set(); sns.set_context("paper", rc=seaborn_rc);
    #df2["date"].strftime('%Y-%m-%d')
    g = sns.stripplot (y="peroba_lineage", x="date_float", edgecolor="black", jitter=0.3, linewidth=0.1, marker="o", alpha=0.8, s=2, data=df)
    #x = g.set_xticklabels(df2["date"].dt.strftime('%Y-%m-%d'))  ## order in table is NOT PLOT ORDER!
    #x = g.get_xticklabels ()
    #for lab in x:
    #    print (lab)
    #    lab.set_text(float_to_date[lab.get_text()])
    #x = g.set_xticklabels ([1,2,3, 30])
    ## see also https://www.machinelearningplus.com/plots/top-50-matplotlib-visualizations-the-master-plots-python/
    md_description = """
## Lineages over time
"""
    md_description += "\n<br>![](jitter_lineages.pdf)\n<br>\n" # same directory as report.md, it doesn't need full path? 
    fname = os.path.join(output_dir,"jitter_lineages.pdf")
    g.figure.savefig(fname, format="pdf")
    g.figure.clf()
    return md_description 

def plot_jitter_lineages (metadata, output_dir, figdir):
    df = metadata.copy() 
    df = df[df["submission_org_code"].str.contains("NORW", na=False)]
    df.dropna(subset=["collection_date"], inplace=True)
    #df["dtime"] = pd.to_datetime(df["collection_date"], infer_datetime_format=True, yearfirst=True, errors='coerce')
    df["dtime"] = df["collection_date"]
    df["peroba_lineage"] = df["peroba_lineage"].str.split('/').str[0] # only first inference
    df["peroba_lineage"].fillna("unclassified", inplace=True)

    sl_df = df.sort_values("dtime").groupby("peroba_lineage").tail(1)  # lineages, sorted by most recent dtime
    lineage_idx = {x:i for i,x in enumerate(sl_df["peroba_lineage"])} # map lineage to order in plot 
    y_label = [x for x in sl_df["peroba_lineage"]] 

    ## add Normal() jitter
    x = ((df["dtime"] - df["dtime"].min())/np.timedelta64(1, 'D')).astype(float) + np.random.normal(0,0.25, len(df["dtime"]))
    yfix = [lineage_idx[i] for i in df["peroba_lineage"]]
    y = yfix + np.random.normal(0,0.12, len(df["peroba_lineage"]))
    ## create x labels
    x_range = pd.date_range(df["dtime"].min(), df["dtime"].max(), freq="1D") # creates uniform interval
    label_interval = int (len(x_range)/14)
    x_label = ["" if i%label_interval else a for i,a in enumerate(x_range.strftime('%Y-%m-%d'))]

    my_cmap = colors.LinearSegmentedColormap.from_list("custom", [(x,cm.get_cmap("tab10")(x)) for x in np.linspace(0,1,10)])

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(12, 7)) # I'm leaving ax in case we want several plots
    ax.scatter (x,y, alpha=0.2, edgecolors=(0,0,0,0.1), cmap=my_cmap, c=yfix, s=30) 

    ax.set_xticks(np.arange(0, len(x_label), 1));
    ax.set_yticks(np.arange(0, len(y_label), 1));
    ax.set_xticklabels(x_label)
    ax.set_yticklabels(y_label)
    ax.grid(axis="y", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=14)

    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")
    ax.set_title (f"Lineages over time for all genomes sequenced at QIB", fontsize=16)
    [t.set_color(my_cmap(i/len(y_label))) for (i,t) in enumerate(ax.yaxis.get_ticklabels())]
    fig.tight_layout()
    caption= "Genomes sequenced at the QIB by lineage over time (sequences not included in phylogenetic pipeline\
are considered 'unclassified')"
    desc = "\n\n## Lineages over time\n\n"

    fname = f"{figdir}/jitter_lineages.pdf"
    fig.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(100))  # or plt.savefig()
    pdf_desc = f"{desc}![{caption}]({fname})\n\n" 
    fname = f"{figdir}/jitter_lineages.png"
    fig.savefig(os.path.join(output_dir,fname), format="png", dpi=(100))  # or plt.savefig()
    html_desc = f"{desc}![{caption}]({fname})\n\n<br>" 

    hdesc, pdesc = plot_jitter_uk_lineages (metadata, output_dir, figdir)
    html_desc += hdesc
    pdf_desc  += pdesc
    return html_desc, pdf_desc 

def plot_jitter_uk_lineages (metadata, output_dir, figdir):
    df = metadata.copy() 
    df = df[df["submission_org_code"].str.contains("NORW", na=False)]
    df.dropna(subset=["collection_date"], inplace=True)
    df["dtime"] = df["collection_date"]
    df["peroba_uk_lineage"] = df["peroba_uk_lineage"].str.split('/').str[0] # only first inference
    df["peroba_uk_lineage"].fillna("unclassified", inplace=True)

    sl_df = df.sort_values("dtime").groupby("peroba_uk_lineage").tail(1)  # lineages, sorted by most recent dtime
    lineage_idx = {x:i for i,x in enumerate(sl_df["peroba_uk_lineage"])} # map lineage to order in plot 
    y_label = [x for x in sl_df["peroba_uk_lineage"]] 

    ## add Normal() jitter
    x = ((df["dtime"] - df["dtime"].min())/np.timedelta64(1, 'D')).astype(float) + np.random.normal(0,0.2, len(df["dtime"]))
    yfix = [lineage_idx[i] for i in df["peroba_uk_lineage"]]
    y = yfix + np.random.normal(0,0.08, len(df["peroba_uk_lineage"]))
    ## create x labels
    x_range = pd.date_range(df["dtime"].min(), df["dtime"].max(), freq="1D") # creates uniform interval
    label_interval = int (len(x_range)/14)
    x_label = ["" if i%label_interval else a for i,a in enumerate(x_range.strftime('%Y-%m-%d'))]

    
    iteriter = 3 * [x for x in np.linspace(0,1,20)] # 5 * [] to iterate 5 times over colours, cyclically
    npoints = len(iteriter) -1  # i must be [0,1] (i.e. both inclusive)
    iteriter = [(i/npoints,cm.get_cmap("tab20b")(x)) for i,x in enumerate(iteriter)]
    my_cmap = colors.LinearSegmentedColormap.from_list("custom", iteriter)

    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(12, 14)) # I'm leaving ax in case we want several plots
    ax.scatter (x,y, alpha=0.25, edgecolors=(0,0,0,0.2), cmap=my_cmap, c=yfix, s=25) 

    ax.set_xticks(np.arange(0, len(x_label), 1));
    ax.set_yticks(np.arange(0, len(y_label), 1));
    ax.set_xticklabels(x_label)
    ax.set_yticklabels(y_label)
    ax.set_ylim(len(y_label)-75,len(y_label))
    ax.grid(axis="y", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(axis='x', which='major', labelsize=14)
    ax.tick_params(axis='y', which='major', labelsize=11)

    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")
    ax.set_title (f"UK lineages over time for all genomes sequenced at QIB", fontsize=16)
    [t.set_color(my_cmap(i/len(y_label))) for (i,t) in enumerate(ax.yaxis.get_ticklabels())]
    fig.tight_layout()
    caption= "Genomes sequenced at the QIB by UK_lineage over time (sequences not included in phylogenetic pipeline\
are considered 'unclassified')"

    fname = f"{figdir}/jitter_uk_lineages.pdf"
    fig.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(100))  # or plt.savefig()
    pdf_desc = f"![{caption}]({fname})\n\n" 
    fname = f"{figdir}/jitter_uk_lineages.png"
    fig.savefig(os.path.join(output_dir,fname), format="png", dpi=(100))  # or plt.savefig()
    html_desc = f"\n\n![{caption}]({fname})\n\n<br>" 
    return html_desc, pdf_desc 

def plot_bubble_per_cluster (metadata, counter, output_dir, figdir):
    df=metadata.copy()
    df.dropna(subset=["collection_date"], inplace=True)
    df = df.groupby(["collection_date","adm2"]).size().to_frame(name="size")
    df = df.reset_index() # collection and adm2 are indices of series (that we transformed into frame)
    if len(df["adm2"].unique()) == 0  or len(df["collection_date"].unique()) == 0:
         logger.warning("\ntoo few samples with information for bubble plot<br>\n") ## not all samples have collection_date, for instance
         return "\n", "\n"  # also happens for global (extended_mode) reports
    ratio = len(df["adm2"].unique())/len(df["collection_date"].unique())
    if ratio > 4: ratio = 4
    if ratio < 0.25: ratio = 0.25
    #plt.figure(figsize=(dyn_width,dyn_height)); 
    plt.figure(figsize=(12, ratio * 12)); 
    sns.set(font_scale=1.5); sns.set_context("paper", rc=seaborn_rc);
    sns.set_palette("cubehelix", 8)
    g = sns.scatterplot(y="adm2", x="collection_date", size="size", sizes=(30,150), edgecolor="white", alpha=0.7, data=df)
    #g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='right') # nt working
    plt.xticks(rotation=30, horizontalalignment='right') # alternative to loop above (pyplot only)
    plt.tight_layout()
    caption = f"Regions (*adm2*) where genomes from cluster {counter} were collected, over time"

    fname = f"{figdir}/bubble{counter}.pdf"
    g.figure.savefig(os.path.join(output_dir,fname), format="pdf")  # or plt.savefig()
    pdf_desc= f"\n![{caption}]({fname})\n\n"
    fname = f"{figdir}/bubble{counter}.png"
    g.figure.savefig(os.path.join(output_dir,fname), format="png")  # or plt.savefig()
    html_desc= f"\n\n![{caption}]({fname})\n\n<br>"
    g.figure.clf()
    return html_desc, pdf_desc 

def plot_postcode_map (metadata, counter, output_dir, figdir):
    df=metadata.copy()
    # counts = df.groupby(["adm2_private","peroba_lineage"]).size().unstack() # 2D: order is [rows, columns]
    ## casecounts is 1D: just counts per postcode (Series, not DataFrame) therefore we create a DF with column *name*
    casecounts = df.groupby(["adm2_private"]).size().to_frame(name="cnt") 
    casecounts.fillna(0, inplace=True)
    casecounts.reset_index(drop=False, inplace=True)
    casecounts.rename(columns={"adm2_private":"area"},inplace=True)
    logger.debug("cluster %s: nrows = %s, valid postcodes=%s",counter, df.shape, sum(casecounts["cnt"]))

    ## prepare geographical lines (land/sea etc.)
    fig, ax = plt.subplots(nrows=1, ncols=1,figsize=(6,8)) # I'm leaving ax in case we want several plots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    m = Basemap(resolution='h', projection='merc', lat_0=54.5, lon_0=-4.36, ax=ax,
             llcrnrlon=-0.92, llcrnrlat= 51.5, urcrnrlon=1.8, urcrnrlat=53.4)
    #m.drawmapboundary(fill_color='aliceblue') ## fillcolor is ocean; 'lake' below are... lakes!
    #m.fillcontinents(color='#ffffff',lake_color='aliceblue')
    #m.drawcoastlines()
    m.shadedrelief()

    m.readshapefile(uk_postcode_area_file, 'areas') # Areas.shp 
    poly = pd.DataFrame({
            'shapes': [patches.Polygon(np.array(shape), True) for shape in m.areas],
            'area': [area['name'] for area in m.areas_info]
        })
    poly = poly.merge(casecounts, on="area", how="left")

    colscale = [(0,(1,1,1))] + [((x+0.01)/1.01,cm.get_cmap("magma_r")(x)) for x in np.linspace(0, 1, 16)]
    cmap = colors.LinearSegmentedColormap.from_list("custom", colscale)
    # cmap = plt.get_cmap('Oranges')   
    pc = PatchCollection(poly.shapes, zorder=2)
    norm = colors.Normalize()

    pc.set_facecolor(cmap(norm(poly["cnt"].fillna(0).values)))
    pc.set_edgecolor("slategrey")
    pc.set_linewidth(0.2)
    ax.add_collection(pc)
    ax.set_title(f"Cluster {counter}")
    mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(poly["cnt"])
    plt.colorbar(mapper, shrink=0.9, ax=ax, orientation="horizontal", pad= 0.02)
    plt.gcf().set_rasterized(True)
    fig.tight_layout()
    caption = f"Number of samples per region (postal code) for cluster {counter}, for those samples with this information"

    fname = f"{figdir}/map{counter}.pdf"
    fig.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(400), pad_inches = 0)  # or plt.savefig()
    pdf_desc = f"\n![{caption}]({fname})\n\n"
    fname = f"{figdir}/map{counter}.png"
    fig.savefig(os.path.join(output_dir,fname), format="png", dpi=(400), pad_inches = 0)  # or plt.savefig()
    hdesc1 = f"![]({fname})"
    fig.clf()

    caption2, hdesc2, pdesc = plot_adm2_map (metadata, counter, output_dir, figdir)
    pdf_desc  += pdesc
    html_desc  = f"\n\n| {caption} | {caption2} |\n|-------|--------|\n|{hdesc1}|{hdesc2}|\n\n"
    return html_desc, pdf_desc 

def plot_adm2_map (metadata, counter, output_dir, figdir):
    df=metadata.copy()
    warnings.filterwarnings('ignore')
    nuts = gpd.read_file (f"{uk_gadm_nuts2_file}.shp")
    nuts.crs = "EPSG:4326"

    casecounts = df.groupby(["adm2"]).size().to_frame(name="cnt") 
    casecounts.fillna(0, inplace=True)
    casecounts.reset_index(drop=False, inplace=True)
    casecounts.rename(columns={"adm2":"NAME_2"},inplace=True)
    poly = nuts.merge(casecounts, on="NAME_2", how="left")
    poly["cnt"].fillna(0,inplace=True) # geopandas works on a dataframe

    fig, axs = plt.subplots(nrows=1, ncols=1,figsize=(5,6)) # I'm leaving ax in case we want several plots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    poly.plot(column="cnt", ax=axs, edgecolor='black', linewidth=0.1, #missing_kwds={'color': 'lightgrey'},
              legend=True, cmap='YlGn', scheme='natural_breaks')
    axs.set_xlim(-8.1, 1.9)
    axs.set_ylim(49.9,61)
    axs.set_aspect('equal') 
    #axs.set_axis_off()
    axs.set_title(f"Cluster {counter}")
    fig.tight_layout()
    axs.axis('off') # remove frame and optimise space
    caption = f"Number of samples per area (NUTS2) for cluster {counter}, for samples where this information is available" 

    plt.gcf().set_rasterized(True)
    fname = f"{figdir}/admap{counter}.pdf"
    fig.savefig(os.path.join(output_dir,fname), format="pdf", dpi=(400))  # or plt.savefig()
    pdf_desc = f"\n![{caption}]({fname})\n\n"
    fname = f"{figdir}/admap{counter}.png"
    fig.savefig(os.path.join(output_dir,fname), format="png", dpi=(400))  # or plt.savefig()
    html_desc = f" ![]({fname}) " # no nowlines, we'll have a table
    fig.clf()

    return caption, html_desc, pdf_desc 
