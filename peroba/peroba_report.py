
import logging, ete3, argparse
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams

from utils import *

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='perobaDB %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

cwd = os.getcwd()
template_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)), "../report")

#outdir = os.path.join(cwd, args.outdir)

def tex_formattted_string (string): # YAML header is translated into latex by pandoc; we have to tell it's already translated
    return "\"`" + string + "`{=latex}\""  # returns  "`/usr/file_name.txt`{=latex}"

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

def plot_jitter_lineages (metadata, output_dir):
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

def merge_metadata_with_csv (metadata0, csv0, tree, tree_leaves):
    """
    Local (NORW) sequences are named NORW-EXXXX, which is central_sample_id. We match those to COGUK when possible
    This adds all our local (stablished) sequences to the database, i.e. not only the ones to query
    The priority column will have a value from 0 (sequence can be safely excluded) to 255 (sequence must remain);
    therefore local sequences must start with at least 127
    The class column is to distinguish NORW_SEQ and NORW_MISSING (if sequence is present or not). It may also be
    NORW_QUERY, COGUK and GISAID
    New metadata has `date_sequenced` and `Postcode` that are useful
    """
    ## local csv may contain info on rejected samples (e.g. as 2020.05.12 there are 
    ## 452 rows but only 302 sequences, since 150 were not sequenced yet or QC rejected)
    metadata = metadata0.copy() ## to make sure we don't modify global metadata 
    csv = csv0.copy()
    matched = metadata[ metadata["central_sample_id"].isin(csv.index) ] # index of csv is "central_sample_id"
    csv["submission_org_code"] = "NORW"
    csv["submission_org"] = "Norwich"

    leaf_names = [x for x in tree_leaves.keys()] # some will be updated below to COGUK format, some are already in COGUK format
    seqs_not_in_tree = dict()  ## NORW sequences that should be in tree but are not
    # change leaf names whenever possible (i.e. they match to 'official' COGUK long names)
    for shortname, longname  in zip(matched["central_sample_id"], matched["sequence_name"]):
        if shortname in leaf_names:
            tree_leaves[str(shortname)].name = longname # leaf names E9999 --> England/E9999/2020
        else:
            seqs_not_in_tree[str(shortname)] = longname
    tree_leaves = {(leaf.name):leaf for leaf in tree.iter_leaves()} # dict {leaf_str_name: ete3_node}
    leaf_names = set (leaf_names + [x for x in tree_leaves.keys()]) # update with new leaf names (from COGUK)
#    norwseqnames = [x.id for x in seq_matrix]
    if len(seqs_not_in_tree):
        tbl_warning = [f"[WARNING]\t{x}\t{y}" for x,y in seqs_not_in_tree.items()]
        logger.warning("Samples from NORW not found on tree (excluded from phylo analysis due to low quality perhaps?):\n%s", "\n".join(tbl_warning))

    # remove csv rows not present in tree
    csv = csv[ csv.index.isin(leaf_names) ]
    logger.info("Number of NORW samples present in tree (according to CSV and global metadata): %s", len(csv))
    ## temporarily use central_sample_id as index, so that we can merge_by_index
    matched.reset_index(inplace=True) ## downgrades current index to a regular column
    matched.set_index ("central_sample_id", append=False, drop = True, inplace = True) # append creates a new col w/ index 
    # delete empty columns
    csv.dropna      (axis=1, how='all', inplace=True) # currently, useless columns
    matched.dropna  (axis=1, how='all', inplace=True) # with NA only 
    ## merge csv with corresponding elements from global metadata (note that these are just intersection with csv)
    csv = df_merge_metadata_by_index (csv, matched) 
    # replace receive leaf names in case it's NORW-E996C 
    csv["peroba_seq_uid"] = csv["peroba_seq_uid"].fillna(csv.index.to_series())
    csv["sequence_name"] = csv["sequence_name"].fillna(csv.index.to_series())
    #csv[metadata.index.names[0]] = csv[metadata.index.names[0]].fillna(csv.index.to_series()) # same as above

    ## revert index to same as global metadata ("peroba_seq_uid" usually)
    csv.reset_index (drop=False, inplace=True) ## drop=True means drop index completely, not even becomes a column
    csv.set_index (metadata.index.names, drop = True, inplace = True) # drop to avoid an extra 'peroba_seq_uid' column
    # merge with COGUK 
    csv = df_merge_metadata_by_index (csv, metadata) 
    csv['days_since_Dec19'] = csv['collection_date'].map(lambda a: get_days_since_2019(a, impute = True))
    csv["collection_date"] = pd.to_datetime(csv["collection_date"], infer_datetime_format=False, errors='coerce')

    # remove all rows not present in tree
    csv = csv[ csv.index.isin(leaf_names) ]
    # check if we have information about all leaves (o.w. we must prune the tree)
    len_csv = len(csv)
    len_tre = len(tree_leaves)
    logger.info("Number tree leaves mapped to the metadata: %s", len_csv)
    if len_csv < len_tre:
        logger.warning("Number of leaves (%s) is higher than matched metadata (%s)",len_tre, len_csv)
        id_meta = set(csv.index.array) # <PandasArray> object, works like an np.array
        id_leaves = set([x for x in tree_leaves.keys()]) 
        only_in_tree = list(id_leaves - id_meta)
        logger.warning("Unmapped leaves:\n%s", "\n".join(only_in_tree))
        for i in only_in_tree:
            del tree_leaves[i]
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates

    metadata0 = df_merge_metadata_by_index (csv, metadata0) 
    metadata0["collection_date"] = pd.to_datetime(metadata0["collection_date"], infer_datetime_format=False, errors='coerce')
    return metadata0, csv, tree, tree_leaves

def fix_csv_columns (csv, csv_cols=None):
    if csv_cols is None: 
        csv_cols = ["adm2", "Postcode", "submission_org_code", "date_sequenced",
                "source_age", "source_sex", "collecting_org", "ICU_admission", "PCR Ct value", "Repeat Sample ID"]
    csv_cols = [x for x in csv_cols if x in csv.columns]
    if csv_cols:
        if "submission_org_code" in csv_cols:
            csv.loc[~csv["submission_org_code"].str.contains("NORW", na=False), "date_sequenced"] = "nil" # only for NORW
        for col in csv_cols:
            csv[col].fillna("nil", inplace=True)
        if "ICU_admission" in csv_cols:
            csv["ICU_admission"] = csv["ICU_admission"].replace("Unknown", "nil", regex=True)
    logger.info ("Follow-up columns found in CSV: %s", " ".join(csv_cols))
    return csv, csv_cols

def ASR_subtrees (metadata0, tree0, csv_cols = None, reroot = True, method = None):
    """ ancestral state reconstruction assuming binary or not. Clusters are based on "locality": patient is from Norfolk
    *or* was sequenced here in NORW.
    One column at a time, so the n_threads doesn't help.
    """
    if method is None: method = ["DOWNPASS", "DOWNPASS"]
    if isinstance (method, str): method = [method, method]

    metadata = metadata0.copy()  # work with copies
    tree = tree0.copy()
    if reroot: ## this should be A,B 
        R = tree.get_midpoint_outgroup()
        tree.set_outgroup(R)

    md_description="""
    Sequence clusters are based on "locality": patients from Norfolk (field `adm2`) *or* patients that were sequenced here
    (submission org = `NORW`). <br>
    """
    ## subtrees will be defined by 'locality'
    yesno = (metadata["adm2"].str.contains("Norfolk", case=False)) | (metadata["submission_org_code"].str.contains("NORW", case=False))
    df = pd.DataFrame(yesno.astype(str), columns=["local"], index=metadata.index.copy())
    x = get_binary_trait_subtrees (tree, df, trait_column = "local", trait_value = "True", elements = 1, method=method[0])
    subtree, mono, result, trait_name = x  # most important is subtree
    print (mono)
    
    ## decorate the tree with ancestral states
    if (csv_cols):
        result = acr (tree, df, prediction_method = method[1], force_joint=False) ## annotates tree nodes with states (e.g. tre2.adm2)

    node2leaves = tree.get_cached_content(store_attr="name") # dict node:[leaves]
    submetadata = [metadata.loc[metadata.index.isin(node2leaves[x])] for x in subtree]
    # not currently used
    supermetadata = [metadata.loc[metadata.index.isin(node2leaves[x.up])] for x in subtree if x.up]

    return subtree, submetadata, md_description


def prepare_report_files (metadata, csv, tree, tree_leaves, input_dir, output_dir, title_date):
    mkd_file_name = os.path.join(output_dir,f"report_{title_date}.md")
    pdf_file_name = os.path.join(output_dir,f"report_{title_date}.pdf")
    titlepage_name = os.path.join(input_dir,"titlepage-fig.pdf")
    pagebackground_name = os.path.join(input_dir,"letterhead-fig.pdf")
    pandoc_template_name = os.path.join(input_dir,"eisvogel")
    
    md_description = """
---
title: "Phylogenomic Report on SARS-CoV2 in Norfolk area"
author: [Leonardo de Oliveira Martins]
date: "{title_date}"
keywords: [Markdown, SARS-CoV19]
titlepage: true,
titlepage-text-color: "FFFFFF"
titlepage-rule-color: "360049"
titlepage-rule-height: 0
titlepage-background: {titlepage}
page-background: {pagebackground}
page-background-opacity: "1"
fontsize: 11pt
papersize: a4
geometry:
- top=35mm
- bottom=30mm
- left=15mm
- right=15mm
- heightrounded
...
# Report for {title_date}
""".format(title_date = title_date, 
            titlepage=tex_formattted_string(titlepage_name), # protects from YAML substitution 
            pagebackground=tex_formattted_string(pagebackground_name))

    report_fw = open (mkd_file_name, "w")
    report_fw.write(md_description)

    ## prepare data (merging, renaming, main ASR)
    csv, csv_cols = fix_csv_columns (csv)
    metadata, csv, tree, tree_leaves = merge_metadata_with_csv (metadata, csv, tree, tree_leaves)
    sub_tree, sub_csv, description = ASR_subtrees (csv, tree, csv_cols)
    ## start plotting 
    md_description = plot_jitter_lineages (metadata, output_dir)
    report_fw.write(md_description)
    md_description = plot_genomes_sequenced_over_time (metadata, output_dir)
    report_fw.write(md_description)
    
    report_fw.close()

    runstr = f"pandoc {mkd_file_name} -o {pdf_file_name} --from markdown --template {pandoc_template_name} --listings"
    logger.info("running command:: %s", runstr)
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    logger.debug("output verbatim:\n%s",proc_run)
    
class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="peroba_report is the script that generates a PDF report given a tree and metadata for new genomes",
    usage='''peroba_report [options]''')

    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
            help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-m', '--metadata', metavar='csv.gz', required=True, help="metadata file formated by peroba")
    parser.add_argument('-c', '--csv', metavar='csv', required=True, help="csv table with metadata from NORW")
    parser.add_argument('-t', '--tree', metavar='treefile', required=True,  help="single treefile in newick format")
    parser.add_argument('-p', '--prefix', action="store", help="Date to be added to report")
    parser.add_argument('-i', '--input', action="store", help="Input directory with pandoc templates, if not using default")
    parser.add_argument('-o', '--output', action="store", help="Output database directory. Default: working directory")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: 
        output_d = os.path.join(cwd, args.output)
        pathlib.Path(output_d).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exists
    else: 
        output_d = cwd

    if args.input: input_d = os.path.join(cwd, args.input)
    else: input_d = template_dir

    if args.prefix: title_date = args.prefix
    else: title_date = datetime.datetime.now().strftime("%Y-%m-%d") 

    logger.info("Reading metadata (previously prepared by peroba)")
    metadata = df_read_genome_metadata (args.metadata, index_name = "peroba_seq_uid")
    logger.info("Reading CSV file with metadata from NORW")
    csv = df_read_genome_metadata (args.csv, index_name = "central_sample_id")
    logger.info("Reading tree file and checking if there are duplicate names")
    treestring = open(args.tree).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = ete3.Tree(treestring)
    tree_length = len([leaf.name for leaf in tree.iter_leaves()])
    tree_leaves = {str(leaf.name):leaf for leaf in tree.iter_leaves()} # dup leaves will simply overwrite node information
    if (tree_length > len(tree_leaves)):
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.debug("Found duplicate leaf names, will keep one at random")
    logger.info("%s leaves in treefile and metadata with shape %s", len(tree_leaves), str(metadata.shape))

    markdown_text = prepare_report_files (metadata, csv, tree, tree_leaves, input_d, output_d, title_date)


if __name__ == '__main__':
    main()
