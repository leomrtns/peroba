import logging, ete3, argparse, pathlib
import pandas as pd
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams

from utils import *
import phylodrawing as phylo
import statsdrawing as stdraw
import common


logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_report %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

cwd = os.getcwd()
template_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)), "data/report")

## FIXME : I must destroy UK lineage etc from CSV (since I imputed it!)

#outdir = os.path.join(cwd, args.outdir)

def tex_formattted_string (string): # YAML header is translated into latex by pandoc; we have to tell it's already translated
    return "\"`" + string + "`{=latex}\""  # returns  "`/usr/file_name.txt`{=latex}"

def plot_over_clusters (csv, tree, min_cluster_size = None, output_dir=None):
    if output_dir is None: output_dir = cwd
    if min_cluster_size is None: min_cluster_size = 5
    df = csv.copy()

    md_description = """
## Phylogenetic clusters\n
Only clusters with more than {minc_size} elements are shown.
""".format (minc_size = min_cluster_size)

    sub_csv, sub_tree, this_description, csv = phylo.ASR_subtrees (csv, tree)
    md_description += this_description

    logger.info("Decorating tree before plotting all clusters")
    colmap_dict = phylo.colormap_from_dataframe (csv)
    ts = phylo.return_treestyle_with_columns (colmap_dict)

    logger.info("Entering loop over clusters found; sit tight!")
    for i,(c,t) in enumerate(zip(sub_csv,sub_tree)):
        if len(c) > min_cluster_size:
            md_description += f"\n### Cluster {i}\n"
            md_description += phylo.plot_single_cluster (c, t, i, ts, output_dir)

            md_description += stdraw.plot_bubble_per_cluster (c, i, output_dir)
            
            md_description += stdraw.plot_time_heatmap (c, i, output_dir)
            
            md_description += stdraw.plot_postcode_map (c, i, output_dir)

    md_description += """

The complete phylogenetic tree is displayed in a separate document due to its large size.<br>

\n<br>(report generated at {today})<br>
Care must be taken when interpreting due to the phylogenetic uncertainty (many similar sequences, with ambiguous sites,
incomplete sampling, and other statistical caveats), and also due to the incomplete nature of the phylogenetic 
metadata: sequences failing sequencing quality and phylogenetic quality (excess of `N`s etc) are not
included in the metadata. We add `NORW` sequences back whenever possible. 

Furthermore this analysis uses only a small subset of all available genomic data at COG-UK: it is based on sequences
inferred to be relevant to the regional analysis, and even those may be removed (due to low quality or duplicates).

Software still under development.<br>
""".format(today=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) 
    tree.render(os.path.join(output_dir,"full_tree.pdf"), w=1200, tree_style=ts)

    return md_description, csv 

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

    remove_cols = [c for c in common.remove_from_master_cols if c in csv0.columns]  ## imputed values from previous iteration
    if len(remove_cols):
        csv0.drop (labels = remove_cols, axis=1, inplace = True) 
    csv = csv0.copy()

    matched = metadata[ metadata["central_sample_id"].isin(csv.index) ] # index of csv is "central_sample_id"
    csv["submission_org_code"] = "NORW"
    csv["submission_org"] = "Norwich"

    leaf_names = [x for x in tree_leaves.keys()] # some will be updated below to COGUK format, some are already in COGUK format
    seqs_not_in_tree = dict()  ## NORW sequences that should be in tree but are not
    # change leaf names whenever possible (i.e. they match to 'official' COGUK long names)
    for shortname, longname  in zip(matched["central_sample_id"], matched["sequence_name"]):
        if shortname in leaf_names: # leaf_names is a list; if x in y means x==z for some z in y (!= 'A' in 'ZAPPA')
            tree_leaves[str(shortname)].name = longname # leaf names E9999 --> England/E9999/2020
        else:
            if longname not in leaf_names: 
                seqs_not_in_tree[str(shortname)] = longname
    tree_leaves = {(leaf.name):leaf for leaf in tree.iter_leaves()} # dict {leaf_str_name: ete3_node}
    leaf_names = set (leaf_names + [x for x in tree_leaves.keys()]) # update with new leaf names (from COGUK)
#    norwseqnames = [x.id for x in seq_matrix]
    if len(seqs_not_in_tree):
        tbl_debug = [f"[DEBUG]\t{x}\t{y}" for x,y in seqs_not_in_tree.items()]
        logger.warning("%s\tsamples from NORW in metadata were not found on tree (excluded from analysis due to low quality?)", len(tbl_debug))
        logger.debug("And the sequences are\n%s", "\n".join(tbl_debug))

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
    csv = common.df_merge_metadata_by_index (csv, matched) 
    # replace receive leaf names in case it's NORW-E996C 
    csv["peroba_seq_uid"] = csv["peroba_seq_uid"].fillna(csv.index.to_series())
    csv["sequence_name"] = csv["sequence_name"].fillna(csv.index.to_series())
    #csv[metadata.index.names[0]] = csv[metadata.index.names[0]].fillna(csv.index.to_series()) # same as above

    ## revert index to same as global metadata ("peroba_seq_uid" usually)
    csv.reset_index (drop=False, inplace=True) ## drop=True means drop index completely, not even becomes a column
    csv.set_index (metadata.index.names, drop = True, inplace = True) # drop to avoid an extra 'peroba_seq_uid' column
    # merge with COGUK 
    csv = common.df_merge_metadata_by_index (csv, metadata) 

    csv["collection_date"] = pd.to_datetime(csv["collection_date"], infer_datetime_format=True, yearfirst=True, errors='coerce')
    #csv['days_since_Dec19'] = csv['collection_date'].map(lambda a: get_days_since_2019(a, impute = False))
    csv['days_since_Dec19'] = csv['collection_date'].map(lambda a: (a - datetime.datetime(2019,12, 1)).days) 

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
        logger.warning("Removing unmapped leaves from dictionary. List of unmapped leaves:\n%s", "\n".join(only_in_tree))
        for i in only_in_tree:
            del tree_leaves[i]
        logger.info("Will now prune these leaves from tree.")
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates

    logger.info("Comparing tree leaves to check for duplicates")
    leaf_list = [leaf.name for leaf in tree.iter_leaves()] # may have duplicates
    tree_length = len(leaf_list)
    tree_leaves = {str(leaf.name):leaf for leaf in tree.iter_leaves()} # dup leaves will simply overwrite node information
    if (tree_length > len(tree_leaves)):
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.warning("After mapping/merging, some leaves have same name -- e.g. if the same sequence was included twice in the")
        logger.warning("  phylogenetic analysis, one copy from the NORW database and one from COGUK. I will keep only one of each, at random")
        logger.warning("  some examples (whenever counter>1): %s", str(collections.Counter(leaf_list).most_common(20)))

    logger.info("Merging metadata files by sequence name, after renaming local sequences when possible")
    metadata0 = common.df_merge_metadata_by_index (csv, metadata0) 
    metadata0["collection_date"] = pd.to_datetime(metadata0["collection_date"], infer_datetime_format=True, yearfirst=True, errors='coerce')
    return metadata0, csv, tree, tree_leaves

def merge_original_with_extra_cols (csv_original, metadata):
    cols = [x for x in metadata.columns if x.startswith("peroba")] 
    cols += ["central_sample_id"]
    df =  metadata[cols]
    df = df.merge(csv_original, on="central_sample_id")
    return df  ## no csv_original were harmed in this function

def main_prepare_report_files (metadata0, csv0, tree, tree_leaves, input_dir, output_dir, title_date):
    mkd_file_name = os.path.join(output_dir,f"report_{title_date}.md")
    pdf_file_name = os.path.join(output_dir,f"report_{title_date}.pdf")
    csv_file_name = os.path.join(output_dir,f"csv_{title_date}.csv.gz")
    meta_file_name = os.path.join(output_dir,f"metadata_{title_date}.csv.gz")
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
## for landscape (use \blscape and /elscape within markdown, however it rotates layout (header, footer)
#header-includes: |
#    \usepackage{pdflscape}
#    \newcommand{\blscape}{\begin{landscape}}                                                                                                                  
#    \newcommand{\elscape}{\end{landscape}} 
#...
#\blscape 

    report_fw = open (mkd_file_name, "w")
    report_fw.write(md_description)

    ## prepare data (merging, renaming) csv0 WILL BE MODIFIED to remove columns from prev week
    metadata, csv, tree, tree_leaves = merge_metadata_with_csv (metadata0, csv0, tree, tree_leaves)
    ## start plotting 
    md_description, csv = plot_over_clusters (csv, tree, min_cluster_size = 10, output_dir=output_dir)
    metadata = phylo.save_metadata_inferred (metadata, tree)
    csv = merge_original_with_extra_cols (csv0, metadata)
    csv.to_csv (csv_file_name)
    metadata.to_csv (meta_file_name)

    report_fw.write(md_description)
    md_description = stdraw.plot_jitter_lineages (metadata, output_dir)
    report_fw.write(md_description)
    md_description = stdraw.plot_genomes_sequenced_over_time (metadata, output_dir)
    report_fw.write(md_description)

    #md_description = "\n\end{pdflscape}"
    #report_fw.write(md_description)

    report_fw.close()

    runstr = f"cd {output_dir} && pandoc {mkd_file_name} -o {pdf_file_name} --from markdown --template {pandoc_template_name} --listings"
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
    metadata = common.df_read_genome_metadata (args.metadata, index_name = "peroba_seq_uid")
    logger.info("Reading CSV file with metadata from NORW")
    csv = common.df_read_genome_metadata (args.csv, index_name = "central_sample_id")
    logger.info("Reading tree file and checking if there are duplicate names")
    treestring = open(args.tree).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = ete3.Tree(treestring)
    tree_length = len([leaf.name for leaf in tree.iter_leaves()])
    tree_leaves = {str(leaf.name):leaf for leaf in tree.iter_leaves()} # dup leaves will simply overwrite node information
    if (tree_length > len(tree_leaves)):
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.warning("Found duplicated leaf names in input treefile, will keep one at random")
    logger.info("%s leaves in treefile and metadata with shape %s", len(tree_leaves), str(metadata.shape))

    main_prepare_report_files (metadata, csv, tree, tree_leaves, input_d, output_d, title_date)

if __name__ == '__main__':
    main()
