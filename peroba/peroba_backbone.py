import logging, ete3, argparse
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams

from utils import *

logger = logging.getLogger(__name__) # https://github.com/MDU-PHL/arbow
logger.propagate = False
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba_backbone %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

current_working_dir = os.getcwd()

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

    leaf_list = [leaf.name for leaf in tree.iter_leaves()] # may have duplicates
    tree_length = len(leaf_list)
    tree_leaves = {str(leaf.name):leaf for leaf in tree.iter_leaves()} # dup leaves will simply overwrite node information
    if (tree_length > len(tree_leaves)):
        tree.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
        logger.warning("After mapping/merging, some leaves have same name -- e.g. if the same sequence was included twice in the")
        logger.warning("  phylogenetic analysis, one copy from the NORW database and one from COGUK. I will keep only one of each, at random")
        logger.warning("  some examples (whenever counter>1): %s", str(collections.Counter(leaf_list).most_common(20)))

    metadata0 = df_merge_metadata_by_index (csv, metadata0) 
    metadata0["collection_date"] = pd.to_datetime(metadata0["collection_date"], infer_datetime_format=False, errors='coerce')
    return metadata0, csv, tree, tree_leaves

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

    ## prepare data (merging, renaming)
    metadata, csv, tree, tree_leaves = merge_metadata_with_csv (metadata, csv, tree, tree_leaves)
    ## start plotting 
    md_description = plot_over_clusters (csv, tree, min_cluster_size = 2, output_dir=output_dir)
    report_fw.write(md_description)
    md_description = stdraw.plot_jitter_lineages (metadata, output_dir)
    report_fw.write(md_description)
    md_description = stdraw.plot_genomes_sequenced_over_time (metadata, output_dir)
    report_fw.write(md_description)
    
    report_fw.close()

    runstr = f"cd {output_dir} && pandoc {mkd_file_name} -o {pdf_file_name} --from markdown --template {pandoc_template_name} --listings"
    logger.info("running command:: %s", runstr)
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    logger.debug("output verbatim:\n%s",proc_run)

def read_peroba_database (f_prefix): 
    f_suffix = { ## TODO: share across files (peroba_common?)
            "metadata": ".metadata.csv.gz",
            "tree": ".tree.nhx",
            "sequences": ".sequences.fasta.bz2"
            }
    if f_prefix[-1] == ".": f_prefix = f_prefix[:-1] ## both `perobaDB.0621` and `perobaDB.0621.` are valid
    fname = f_prefix+f_suffix["metadata"]
    logger.info(f"Reading database metadata from \'{fname}\'")
    metadata = pd.read_csv (fname, compression="infer", index_col="peroba_seq_uid") 

    fname = f_prefix+f_suffix["tree"]
    logger.info(f"Reading database tree from \'{fname}\'")
    treestring = open(fname).readline().rstrip().replace("\'","").replace("\"","").replace("[&R]","")
    tree = ete3.Tree(treestring)

    fname = f_prefix+f_suffix["sequences"]
    logger.info(f"Reading database sequences from \'{fname}\'")
    sequences = read_fasta (fname, check_name = False)

    logger.info("Finished loading the database; dataframe has dimensions %s and it's assumed we have the same \
            number of sequences; the tree may be smaller", metadata.shape)
    return [metadata, sequence, tree]

def generate_backbone_dataset (database, csv, sequences, trees, prefix):


    
class ParserWithErrorHelp(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    parser = ParserWithErrorHelp(
    description="""
    peroba_backbone is the script that generates a global backbone data set (COGUK+GISAID) given a local one (NORW).
    It depends on the prefix for a perobaDB set of files (from `peroba_database`), like "perobaDB.0519".
    It's recommended that you add also a local CSV metadata, and you can furthermore add a newick file with extra 
    trees (the tree from previous run is a good choice).
    """, 
    usage='''peroba_backbone <perobaDB> [options]''')

    parser.add_argument('perobaDB')
    parser.add_argument('-d', '--debug', action="store_const", dest="loglevel", const=logging.DEBUG, default=logging.WARNING,
            help="Print debugging statements")
    parser.add_argument('-v', '--verbose', action="store_const", dest="loglevel", const=logging.INFO, help="Add verbosity")
    parser.add_argument('-i', '--input', action="store", help="Directory where perobaDB files are. Default: working directory")
    parser.add_argument('-c', '--csv', metavar='csv', help="csv table with metadata from NORW")
    parser.add_argument('-s', '--sequences', metavar='fasta.bz2', help="extra sequences from NORW")
    parser.add_argument('-t', '--trees', metavar='', help="file with trees in newick format to help produce backbone")
    parser.add_argument('-o', '--output', action="store", help="Output database directory. Default: working directory")

    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel)
    if args.output: 
        output_d = os.path.join(current_working_dir, args.output)
        pathlib.Path(output_d).mkdir(parents=True, exist_ok=True) # python 3.5+ create dir if it doesn't exist
    else: 
        output_d = current_working_dir
     prefix = os.path.join(output_d, "peroba_backbone." + datetime.datetime.now().strftime("%m%d_%H%M") + ".")

    if args.input: input_d = os.path.join(current_working_dir, args.input)
    else: input_d = current_working_dir

    logger.info("Reading metadata, sequences, and tree from peroba_database")
    database = read_peroba_database (os.path.join(input_d, args.perobaDB)) # something like "my_folder/perobaDB.0515"

    if (args.csv):
        csv = None
        fname = os.path.join(current_working_dir, args.csv)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.csv)
        if not os.path.exists(fname):
            logger.warning ("Could not find local CSV file {args.csv}; Will proceed without it")
        else:
            logger.info("Reading CSV file with metadata from NORW")
            csv = df_read_genome_metadata (fname, index_name = "central_sample_id")

    if (args.sequences):
        sequences = None
        fname = os.path.join(current_working_dir, args.sequences)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.sequences)
        if not os.path.exists(fname):
            logger.warning ("Could not find sequence file {args.sequences}; Will proceed without it")
        else:
            logger.info("Reading fasta file with sequences from NORW")
            sequences = read_fasta (fname, check_name = False)

    if (args.trees):
        trees = None
        fname = os.path.join(current_working_dir, args.trees)
        if not os.path.exists(fname):
            fname = os.path.join(input_d, args.trees)
        if not os.path.exists(fname):
            logger.warning ("Could not find tree file {args.trees}; Will proceed without it")
        else:
            logger.info("Reading file with current trees  and checking for duplicate names")
            treestring = [x.rstrip().replace("\'","").replace("\"","").replace("[&R]","") for x in open(fname)]
            trees = []
            for i,trs in enumerate (treestring): 
                tre = ete3.Tree(trs)
                tree_length = len([leaf.name for leaf in tre.iter_leaves()])
                tree_leaves = {str(leaf.name):leaf for leaf in tre.iter_leaves()} # dup leaves will simply overwrite node information
                if (tree_length > len(tree_leaves)):
                    tre.prune([node for node in tree_leaves.values()], preserve_branch_length=True) # or leafnames, but fails on duplicates
                    logger.warning(f"Found duplicated leaf names in input treefile {i}, will keep one at random")
                logger.info("%s leaves in treefile %s", len(tree_leaves), str(i))
                trees.append(tre)

    generate_backbone_dataset (database, csv, sequences, trees, prefix)


if __name__ == '__main__':
    main()
