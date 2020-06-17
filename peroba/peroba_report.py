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

def md_html_table (df):
    cols = ["central_sample_id", "biosample_source_id", "adm2", "adm2_private", "lineage", # sequence_name is same as index
            "peroba_lineage","uk_lineage", "peroba_uk_lineage", "peroba_phylotype", "acc_lineage", "del_lineage", 
            "submission_org_code", "source_age", "source_sex", "ICU_admission", 
            "collection_date", "sequencing_submission_date", "peroba_freq_acgt", "peroba_freq_n"]
    cols = [c for c in cols if c in df.columns]
    dfc = df[cols]
    return "\n\n" + dfc.to_markdown(tablefmt="pipe") + "\n\n"

def plot_over_clusters (metadata, tree, min_cluster_size = None, output_dir=None, figdir=None, pdf_report_name = None):
    if output_dir is None: output_dir = cwd
    if min_cluster_size is None: min_cluster_size = 5
    df = metadata.copy()

    md_description = """
## Phylogenetic clusters\n
Only clusters with more than {minc_size} elements are shown.\n\n
In the tables below, lineage columns starting with `peroba_` include classification inferred by us, while those without
it only have values defined upstream by COGUK.
These columns include `lineage`, `uk_lineage`, and `phylotpe`. Sometimes our method cannot infer the lineage without
error and will include several possible ones, separated by a bar ('`/`').
Some samples are left unclassified by COGUK due to uncertainty (not enough evidence for including it in a `uk_lineage`,
for instance) but our method will impute a lineage/phylotype to it nonetheless.
\n\n
""".format (minc_size = min_cluster_size)
    html_desc = pdf_desc = md_description 

    # df will be metadata with extra peroba_ columns (and sub_df is for each subtree)
    sub_df, sub_tree, md_description, df = phylo.ASR_subtrees (df, tree)
    html_desc += md_description # same text in markdown
    pdf_desc  += md_description

    logger.info("Decorating tree before plotting all clusters")
    colmap_dict = phylo.colormap_from_dataframe (df)
    ts = phylo.return_treestyle_with_columns (colmap_dict)

    logger.info("Entering loop over clusters found; sit tight!")
    for i,(c,t) in enumerate(zip(sub_df,sub_tree)):
        if len(c) > min_cluster_size:
            c = remove_imputation_from_gisaid (c) # remove imputations where it doens't make sense (uk_lineage from China)
            md_description = f"\n### Cluster {i}\n"
            html_desc += md_description
            pdf_desc  += md_description
            html_desc += md_html_table (c)

            hdesc, pdesc = phylo.plot_single_cluster (c, t, i, ts, output_dir, figdir)
            html_desc += hdesc; pdf_desc += pdesc

            hdesc, pdesc = stdraw.plot_bubble_per_cluster (c, i, output_dir, figdir)
            html_desc += hdesc; pdf_desc += pdesc
            
            hdesc, pdesc = stdraw.plot_time_heatmap (c, i, output_dir, figdir)
            html_desc += hdesc; pdf_desc += pdesc
            
            hdesc, pdesc = stdraw.plot_postcode_map (c, i, output_dir, figdir)
            html_desc += hdesc; pdf_desc += pdesc

    html_desc +=f"""
## Files for download
\n[Full tree, available for download due to its large size](./{figdir}/full_tree.pdf)\n<br>
\n[This report, in PDF format](./{pdf_report_name})\n<br>
"""
    pdf_desc +="""
The complete phylogenetic tree is displayed in a separate document due to its large size.<br>\n"""

    tree.render(os.path.join(output_dir,f"{figdir}/full_tree.pdf"), w=1200, tree_style=ts)

    return html_desc, pdf_desc, df 

def finalise_report ():
    md_description = """
\n## Timestamp: {today}<br>
Care must be taken when interpreting due to the phylogenetic uncertainty (many similar sequences, with ambiguous sites,
incomplete sampling, and other statistical caveats), and also due to the incomplete nature of the phylogenetic 
metadata: sequences failing sequencing quality and phylogenetic quality (excess of `N`s etc) are not
included in the metadata. We add `NORW` sequences back whenever possible. 

Furthermore this analysis uses only a small subset of all available genomic data at COG-UK: it is based on sequences
inferred to be relevant to the regional analysis, and even those may be removed (due to low quality or duplicates).

Software still under development.<br>
""".format(today=datetime.datetime.now().strftime("%Y-%m-%d %H:%M")) 
    return md_description

def merge_metadata_with_csv (metadata0, csv0, tree, tree_leaves):
    """
    """
    # csv will have same n_samples but columns from metadata; index is seqname if avail or central_sample_id o.w.
    metadata, csv = common.merge_global_and_local_metadata (metadata0, csv0)

    # just in case (we should not use original csv0) : csv0 has imputed values from previous iteration
    remove_cols = [c for c in common.remove_from_master_cols if c in csv0.columns]  ## imputed values from previous iteration
    if len(remove_cols):
        csv0.drop (labels = remove_cols, axis=1, inplace = True) 
    
    leaf_names = [x for x in tree_leaves.keys()] # some will be updated below to COGUK format, some are already in COGUK format
    seqs_not_in_tree = dict()  # NORW sequences that should be in tree but are not
    norw_seqs_in_tree = dict() # variable used only for debug printing
    # change leaf names whenever possible (i.e. they match to 'official' COGUK long names)
    for shortname, longname  in zip(csv["central_sample_id"], csv["sequence_name"]):
        if shortname in leaf_names: # leaf_names is a list; if x in y means x==z for some z in y (!= 'A' in 'ZAPPA')
            tree_leaves[str(shortname)].name = longname # leaf names NORW-E9999 --> England/NORW-E9999/2020
            norw_seqs_in_tree[str(shortname)] = longname 
        else:
            if longname not in leaf_names: 
                seqs_not_in_tree[str(shortname)] = longname
            else:
                norw_seqs_in_tree[str(shortname)] = longname

    tree_leaves = {(leaf.name):leaf for leaf in tree.iter_leaves()} # dict {leaf_str_name: ete3_node}
    if len(seqs_not_in_tree):
        tbl_debug = [f"[DEBUG]\t{x}\t{y}" for x,y in seqs_not_in_tree.items()]
        logger.warning("%s\tsamples in local csv metadata were not found on tree (excluded from analysis due to low quality?)", len(tbl_debug))
        logger.debug("And the sequences are\n%s", "\n".join(tbl_debug))
    logger.info ("%s\tsamples from local csv metadata were found on tree and will form the basis of the analysis", len(norw_seqs_in_tree))
    logger.debug("And the sequences are\n%s", "\n".join(norw_seqs_in_tree.keys()))
    
    # check if we have information about all leaves (o.w. we must prune the tree)
    id_meta = set(metadata.index.array) # <PandasArray> object, works like an np.array
    id_leaves = set([x for x in tree_leaves.keys()]) 
    only_in_tree = list(id_leaves - id_meta)
    if len(only_in_tree) > 0:
        logger.warning("Removing %s unmapped leaves from tree dictionary", str(len(only_in_tree)))
        logger.debug("List of unmapped leaves:\n%s", "\t".join(only_in_tree))
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
    return metadata, csv, tree, tree_leaves

def merge_original_with_extra_cols (csv_original, metadata):
    #cols =  [x for x in metadata.columns if x.startswith("peroba")] # peroba_freq_ will be twice 
    peroba_cols = ["peroba_" + c for c in common.asr_cols]
    cols =  [x for x in metadata.columns if x in peroba_cols]
    cols += ["central_sample_id"]
    df =  metadata[cols]
    df = df.merge(csv_original, on="central_sample_id")
    return df  ## no csv_original were harmed in this function

def remove_imputation_from_gisaid (metadata):
    uk_cols = [[x,"peroba_" + x] for x in metadata.columns if x in common.uk_specific_cols]
    for c1, c2 in uk_cols:
        # if not from COGUK and not estimated by COGUK pipeline (to allow for uk seqs in GISAID)
        metadata.loc[ metadata["submission_org_code"].isnull() & metadata[c1].isnull() , c2] = np.nan
    return metadata

def markdown_headers (input_dir, output_dir, title_date, figdir):
    # tex_formatted_string() protects from YAML substitution 
    titlepage_name = tex_formattted_string(os.path.join(input_dir,"titlepage-fig.pdf"))
    pagebackground_name = tex_formattted_string(os.path.join(input_dir,"letterhead-fig.pdf"))

    runstr = f"cp {template_dir}/pandoc3.css {output_dir}/pandoc.css; cp {template_dir}/QIlogo_light.png {output_dir}/{figdir}/logo.png"
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)

    html_description = """
---
title: |
  ![]({figdir}/logo.png){{width=2in}} Phylogenomic Report on SARS-CoV2 in Norfolk area

author: [Leonardo de Oliveira Martins]
date: "{title_date}"
keywords: [Markdown, SARS-CoV19]
css:
- pandoc.css
toc-depth: 2
to: html5
from: markdown
...
""".format(figdir = figdir, title_date = title_date) 

    pdf_description = """
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
""".format(title_date = title_date, titlepage=titlepage_name, pagebackground=pagebackground_name)
    ## for landscape (use \blscape and /elscape within markdown, however it rotates layout (header, footer)
    #header-includes: |
    #    \usepackage{pdflscape}
    #    \newcommand{\blscape}{\begin{landscape}}                                                                                                                  
    #    \newcommand{\elscape}{\end{landscape}} 
    #...
    #\blscape 
    return html_description, pdf_description

def main_prepare_report_files (metadata0, csv0, tree, tree_leaves, input_dir, output_dir, title_date):
    # report files
    mkd_htm_file = os.path.join(output_dir,f"_htm_{title_date}.md")
    mkd_pdf_file = os.path.join(output_dir,f"_pdf_{title_date}.md")
    pdfreportname = f"report_{title_date}.pdf"
    pdf_file_name = os.path.join(output_dir,pdfreportname)
    htm_file_name = os.path.join(output_dir,f"report_{title_date}.html")
    pandoc_template_name = os.path.join(input_dir,"eisvogel")
    # output table files
    csv_file_name = os.path.join(output_dir,f"csv_{title_date}.csv.gz")
    meta_file_name = os.path.join(output_dir,f"metadata_{title_date}.csv.gz")
   
    figdir = "figures"
    pathlib.Path(f"{output_dir}/{figdir}").mkdir(parents=True, exist_ok=True)
    html_desc, pdf_desc =  markdown_headers (input_dir, output_dir, title_date, figdir)

    fw_htm = open (mkd_htm_file, "w")
    fw_pdf = open (mkd_pdf_file, "w")
    fw_htm.write(html_desc)
    fw_pdf.write(pdf_desc)

    ## prepare data (merging, renaming) csv0 will be modified to remove columns from prev week
    metadata, csv, tree, tree_leaves = merge_metadata_with_csv (metadata0, csv0, tree, tree_leaves)
    ## start plotting (metadata will receive peroba_ imputed values) 
    html_desc, pdf_desc, metadata = plot_over_clusters (metadata, tree, min_cluster_size = 5, output_dir=output_dir, 
            figdir=figdir, pdf_report_name = pdfreportname)
    fw_htm.write(html_desc)
    fw_pdf.write(pdf_desc)

    metadata = remove_imputation_from_gisaid (metadata) # uk_lineages from Australia samples...
    metadata.to_csv (meta_file_name)
    csv = merge_original_with_extra_cols (csv, metadata)
    csv.to_csv (csv_file_name)

    html_desc, pdf_desc = stdraw.plot_jitter_lineages (metadata, output_dir, figdir)
    fw_htm.write(html_desc)
    fw_pdf.write(pdf_desc)
    
    html_desc, pdf_desc = stdraw.plot_genomes_sequenced_over_time (metadata, output_dir, figdir)
    fw_htm.write(html_desc)
    fw_pdf.write(pdf_desc)

    description = finalise_report ()
    fw_htm.write(description)
    fw_pdf.write(description)

    fw_htm.close()
    fw_pdf.close()

    # pdf pandoc (using eisvogel)
    runstr = f"cd {output_dir} && pandoc {mkd_pdf_file} -o {pdf_file_name} --from markdown --template {pandoc_template_name} --listings"
    logger.info("running command:: %s", runstr)
    proc_run = subprocess.check_output(runstr, shell=True, universal_newlines=True)
    logger.debug("output verbatim:\n%s",proc_run)
    # html pandoc (using embedded css)
    runstr = f"cd {output_dir} && pandoc -s {mkd_htm_file} -o {htm_file_name} --toc"
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
