
import pandas as pd, seaborn as sns, numpy as np


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
    labs = ["" if i%label_interval else a for i,a in enumerate(df[tmp_col].dt.strftime('%Y-%m-%d'))] # notice the "dt"to convert
    return df[new_date_col], labs


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
