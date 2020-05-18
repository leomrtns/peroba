
import logging, ete3, argparse
from utils import *

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='perobaDB %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.INFO)
logger.addHandler(stream_log)

cwd = os.getcwd()
#outdir = os.path.join(cwd, args.outdir)

def tex_formattted_string (string): # YAML header is translated into latex by pandoc; we have to tell it's already translated
    return ""\"\'" + string + "\`{=latex}\""  # returns  "`/usr/file_name.txt`{=latex}"

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

def ASR_subtrees (tree, metadata, reroot = True, asr_cols = None, method = None):
    """ ancestral state reconstruction assuming binary or not. Clusters are based on "locality": patient is from Norfolk
    *or* was sequenced here in NORW.
    One column at a time, so the n_threads doesn't help.
    """
    if method is None: method = ["DOWNPASS", "DOWNPASS"]
    if isinstance (method, str): method = [method, method]
    if asr_cols is None: asr_cols = asr_cols = ["adm2", "Postcode", "submission_org_code","date_sequenced"]

    if reroot:
        R = self.tree.get_midpoint_outgroup()
        self.tree.set_outgroup(R)

    ## subtrees will be defined by 'locality'
    yesno = (metadata["adm2"].str.contains("Norfolk", case=False)) | (metadata["submission_org_code"].str.contains("NORW", case=False))
    df = pd.DataFrame(yesno.astype(str), columns=["local"], index=self.metadata.index.copy())
    x = get_binary_trait_subtrees (self.tree, df, trait_column = "local", trait_value = "True", elements = 1, method=method[0])
    subtree, mono, result, trait_name = x  # most important is subtree
    
    ## decorate the tree with ancestral states
    df = self.metadata[asr_cols];
    df.loc[~df["submission_org_code"].str.contains("NORW", na=False), "date_sequenced"] = "nil" # only for NORW
    for col in asr_cols:
        df[col].fillna("nil", inplace=True)
    result = acr (self.tree, df, prediction_method = method[1], force_joint=False) ## annotates tree nodes with states (e.g. tre2.adm2)

    node2leaves = self.tree.get_cached_content(store_attr="name") # dict node:[leaves]
    submetadata = [self.metadata.loc[self.metadata.index.isin(node2leaves[x])] for x in subtree]
    supermetadata = [self.metadata.loc[self.metadata.index.isin(node2leaves[x.up])] for x in subtree if x.up]

    md_description="""
    Sequence clusters are based on "locality": patients from Norfolk (field `adm2`) *or* patients that were sequenced here
    (submission org = `NORW`). <br>
    """

    return subtree, submetadata, supermetadata

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
    parser.add_argument('-m', '--metadata', metavar='csv', required=True, help="metadata file formated by peroba")
    parser.add_argument('-t', '--tree', metavar='treefile', required=True,  help="single treefile in newick format")
    parser.add_argument('-i','--input', action="store", help="Input directory with pandoc templates")
    parser.add_argument('-o','--output', action="store", help="Output database directory. Default: working directory")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: output_d = os.path.join(cwd, args.output)
    else: output = cwd
    if args.input: input_d = os.path.join(cwd, args.input)
    else: input_d = cwd

    logger.info("Reading metadata (previously prepared by peroba)")
    metadata = df_read_genome_metadata (args.metadata, index_name = "peroba_seq_uid")
    logger.info("Reading tree file and checking if there are duplicate names")
    treestring = open(args.treefile).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = ete3.Tree(treestring)
    tree_length = len([leaf.name for leaf in tree.iter_leaves()])
    tree_leaves = {str(leaf.name):leaf for leaf in self.tree.iter_leaves()} # dup leaves will simply overwrite node information
    if (tree_length > len(tree_leaves)):
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.debug("Found duplicate leaf names, will keep one at random")
    logger.info("%s leaves in treefile and metadata with shape", len(tree_leaves), str(metadata.shape))

    


if __name__ == '__main__':
    main()
