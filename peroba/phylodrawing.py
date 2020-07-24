#!/usr/bin/env python
from peroba import common

import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import matplotlib
matplotlib.use('Agg') # first rule to prevent system of chosing X11-based
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm, colors # colormap

from pastml.acr import acr
import logging, ete3, argparse
import numpy as np, pandas as pd, seaborn as sns
from sklearn import manifold, metrics, cluster, neighbors, decomposition, preprocessing
import sys, gzip, bz2, re, glob, pickle, collections, subprocess, errno, random, itertools, pathlib


logger = logging.getLogger(__name__) 
logger.propagate = False ## otherwise it duplicates (own stream and generic stream)
stream_log = logging.StreamHandler()
log_format = logging.Formatter(fmt='peroba %(asctime)s [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M")
stream_log.setFormatter(log_format)
stream_log.setLevel(logging.DEBUG)
logger.addHandler(stream_log)

def get_ancestral_trait_subtrees (tre, csv, tiplabel_in_csv = None, elements = 1, 
        trait_column ="adm2", trait_value = "NORFOLK", n_threads = 4, method = "DOWNPASS"):
    '''
    Returns ancestral nodes predicted to have a given trait_value (e.g. "NORFOLK") for a given trait column (e.g. "adm2").
    Also returns nodes scrictly monophyletic regarding value.
    If using parsimony then it's better to store only nodes with a single state (or create a binary state?) otherwise
    we may end up chosing too close to root node...    MAP works well to find almost monophyletic nodes, but is quite slow
    max likelihood methods: pastml.ml.MPPA, pastml.ml.MAP, pastml.ml.JOINT,
    max parsimony methods: pastml.parsimony.ACCTRAN, pastml.parsimony.DELTRAN, pastml.parsimony.DOWNPASS
    '''
    if elements < 1: elements = 1 ## inferred state cardinality (how many values are allowed on matching internal node?)
    if method not in ["MPPA", "MAP", "JOINT", "ACCTRAN", "DELTRAN", "DOWNPASS"]:
        method = "DOWNPASS"
    if tiplabel_in_csv: # o.w. we assume input csv already has it
        csv_column = csv[[tiplabel_in_csv, trait_column]] # two columns now, tiplabel is removed below
        csv_column.set_index("sequence_name", drop = True, inplace = True) # acr needs index mapping ete3 leaves
    else:
        csv_column = csv[[trait_column]]  ## nodes will have e.g. n.adm2 (access through getattr(n.adm2)
    if isinstance (trait_value, list):
        trait_value = trait_value[0] ## This function can handle only one value; for grouping try the get_binary version
    ## Ancestral state reconstruction of given trait
    result = acr (tre, csv_column, prediction_method = method, force_joint=False, threads=n_threads) ## annotates tree nodes with states (e.g. tre2.adm2)
    logger.debug("Finished lowlevel ancestral_trait_subtrees()")

    ## Find all internal nodes where trait_value is possible state (b/c is seen at tips below)
    matches = filter(lambda n: not n.is_leaf() and trait_value in getattr(n,trait_column) and
            len(getattr(n,trait_column)) <= elements, tre.traverse("preorder"))
    # print ([x.__dict__ for x in matches]) dictionary of attributes; ete3 also has n.features[] with dict keys

    stored_leaves = set () # set of leaf names (created with get_cached_content)
    subtrees = [] # list of non-overlapping nodes
    node2leaves = tre.get_cached_content(store_attr="name") # set() of leaves below every node; store leaf name only
    for xnode in matches:
        if not bool (stored_leaves & node2leaves[xnode]): # both are sets; bool is just to be verbose
            stored_leaves.update (node2leaves[xnode]) # update() is append() for sets ;)
            subtrees.append(xnode)
    mono = tre.get_monophyletic (values=trait_value, target_attr = trait_column) # from ete3 
    return subtrees, mono, result

def get_binary_trait_subtrees (tre, csv,  tiplabel_in_csv = None, elements = 1, 
          trait_column ="adm2", trait_value = "NORFOLK", n_threads = 4, method = "DOWNPASS", extended_mode=0):
    '''
    Instead of reconstructing all states, we ask if ancestral state is value or not. Still we allow for `elements` > 1
    (2 in practice), so that we accept "yes and no" ancestral nodes. 
    You can group trait values into a list 
    '''
    if elements < 1: elements = 1 ## inferred state cardinality (1 makes much more sense, but you can set "2" as well)
    if method not in ["MPPA", "MAP", "JOINT", "ACCTRAN", "DELTRAN", "DOWNPASS"]:
        method = "DOWNPASS"
    if tiplabel_in_csv: # o.w. we assume input csv already has it
        csv_column = csv[[tiplabel_in_csv, trait_column]] # two columns now, tiplabel is removed below
        csv_column.set_index("sequence_name", drop = True, inplace = True) # acr needs index mapping ete3 leaves
    else:
        csv_column = csv[[trait_column]]  ## nodes will have e.g. n.adm2 (access through getattr(n.adm2)
    if isinstance (trait_value, str): # assume it's a list, below
        trait_value = [trait_value]

    # transform variable into binary
    new_trait = str(trait_column) + "_is_" + "_or_".join(trait_value)
    csv_column[new_trait] = csv_column[trait_column].map(lambda a: "yes" if a in trait_value else "no")
    print ("DEBUG::newcol:: ", collections.Counter(csv_column[new_trait]))
    tree_leaf_nodes = {leaf:leaf.name for leaf in tre.iter_leaves()} # node:leafname, not the other way around as usual...

    csv_column.drop(labels = [trait_column], axis=1, inplace = True)
    ## Ancestral state reconstruction of given trait
    result = acr (tre, csv_column, prediction_method = method, force_joint=False, threads=n_threads) ## annotates tree nodes with states (e.g. tre2.adm2)
    
    for leafnode, leafname  in tree_leaf_nodes.items(): # our tree may have dups, and pastml (correctly) replaces them by a placeholder "t123" 
        leafnode.name = leafname  # here we revert back to duplicate names (dictionary keys are nodes, and values are original names 

    logger.debug("Finished lowlevel binary_trait_subtrees()")

    stored_leaves = set () # set of leaf names (created with get_cached_content)
    subtrees = [] # list of non-overlapping nodes
    node2leaves = tre.get_cached_content(store_attr="name") # set() of leaves below every node; store leaf name only
    
    if extended_mode == 0:
        ## Find all internal nodes where trait_value is possible state (b/c is seen at tips below)
        matches = filter(lambda n: not n.is_leaf() and "yes" in getattr(n,new_trait) and # traits are sets (not lists)
            len(getattr(n,new_trait)) <= elements, tre.traverse("preorder"))
        for xnode in matches:
            if not bool (stored_leaves & node2leaves[xnode]): # both are sets; bool is just to be verbose
                stored_leaves.update (node2leaves[xnode]) # update() is append() for sets ;)
                subtrees.append(xnode)
    else:
        print ("DEBUG:: going deeper")
        matches = filter(lambda n: "yes" in getattr(n,new_trait) and len(getattr(n,new_trait)) <= (elements+1), tre.traverse("preorder"))
        for xnode in matches:
            print ("DEBUG:: match found")
            if not bool (stored_leaves & node2leaves[xnode]): # both are sets; bool is just to be verbose
                stored_leaves.update (node2leaves[xnode]) # update() is append() for sets ;)
                if extended_mode == 2 and xnode.up.up is not None:
                    subtrees.append(xnode.up.up)
                    print ("DEBUG:: up up")
                elif xnode.up is not None: # extended_mode 1 or 2
                    subtrees.append(xnode.up)
                else:
                    subtrees.append(xnode)

    mono = tre.get_monophyletic (values = "yes", target_attr = new_trait) # from ete3 
    subtrees = list(set(subtrees)) ## for cases where node->up are the same
    return subtrees, mono, result, new_trait

def colormap_from_dataframe_old (df, column_list, column_names, cmap_list = None):
    ''' returns a dictionary of lists, where key is the tip label (should be index of csv)
    column_list and column_names must have same length, mapping to names on CSV and on plot
        the default colormap are qualitative so I guess they fail if #elements>#colours...
    '''
    if cmap_list is None:
        cmap_list = ["Accent", "Dark2", "cividis", "jet", "hsv", "viridis", "plasma", "rainbow", "nipy_spectral"]
    if isinstance (cmap_list, str): # assume it's a list, below
        cmap_list = [cmap_list]
    if len(column_list) != len(column_names):
        print ("ops, list of columns (from table) must have an equivalent name (that describes it)")
        return
    ## iterate over columns (phenotypes), to generate colormaps
    d_col = {}
    for i, (c, cname) in enumerate(zip(column_list, column_names)):
        uniq = df[c].unique()
        cmap_elem = cmap_list[i%len(cmap_list)] # cycle over given colormaps
        colorlist = common.continuous_cmap (cmap_elem, len(uniq))
        d_col[cname] = {name:colors.to_hex(col) for name,col in zip(uniq, colorlist)} ## each value is another dict from csv elements to colors
    
    col_missing = colors.to_hex([1,1,1,0]) ## transparent, default to white
    ## iterate over rows (sequences)
    d_seq = {}
    d_seq_lab = {}
    for seq in df.index: # frown upon by pandas wizards
        d_seq[seq] = []
        d_seq_lab[seq] = []
        for c, cname in zip (column_list, column_names):
            if df.loc[seq,c] in d_col[cname]: # valid value for column
                d_seq[seq].append(d_col[cname][ df.loc[seq,c] ])
                d_seq_lab[seq].append(str(df.loc[seq,c]))
            else:
                d_seq[seq].append(col_missing)
                d_seq_lab[seq].append(str(" "))
    return [d_seq, d_seq_lab, d_col, column_names] 

def columnwise_color_scheme (df0):
    df = df0.copy()

    df["area_code"] = df["submission_org_code"].map(lambda a: "inter" if pd.isnull(a) else str(len(a) > 0)) # coguk 
    df.loc[df["adm2"].str.contains("folk", na=False) , "area_code"] = "folk"
    df["age_range"] = pd.qcut(pd.to_numeric(df['source_age'], errors='coerce'), 5).astype(str)
    df["collection_week"] = df["collection_date"].dt.week.astype(str)
    #df["source_age"] = df["source_age"].str.split('.').str[0] # this is before I fixed columns as_numeric
    # color bar width, label column, column title, coded labels, n1 n2 for color table
    clist = [[50, "Area", "adm2", "area_code", 0, 2]] 
    clist.append([40, "Lineages", "peroba_lineage", "peroba_lineage", 3, 3])
    clist.append([40, "UK Lineage", "peroba_uk_lineage", "peroba_uk_lineage", 15, 6])
    clist.append([40, "Hospital", "collecting_org", "collecting_org", 6, 2])
    clist.append([15, "ICU", "is_icu_patient", "is_icu_patient", 8, 2])
    #clist.append([15, "Age", "source_age",    "age_range",    11, 2])
    clist.append([20, "Age", "source_age",    "source_age",    10, 3])
    clist.append([50, "Collection Date", "collection_date", "collection_week", 13, 2])

    return df, clist

def colormap_from_dataframe (df0):
    ''' returns a dictionary of lists, where key is the tip label (should be index of csv)
    '''
    df, trio_names = columnwise_color_scheme (df0)
    width_and_names = [x[:2] for x in trio_names]
    trio_names = [x[2:] for x in trio_names]
    ## iterate over columns (phenotypes), to generate colormaps
    d_col = {} # dict of colour
    # trio: column with labels, column title, column w/ label to be mapped to colour, and colourset indices
    for i, (c, code, n1, n2) in enumerate(trio_names): 
        uniq = df[code].sort_values().unique()
        colorlist = common.list_from_custom_colorset(n1, n2, len(uniq))
        d_col[c] = {name:colors.to_hex(col) for name,col in zip(uniq, colorlist)} ## each value is another dict from csv elements to colors
    
    col_missing = "#ffffff" # white 
    ## iterate over rows (sequences)
    d_seq_color = {}
    d_seq_label = {}
    for seq in df.index: # frown upon by pandas wizards
        d_seq_color[seq] = []
        d_seq_label[seq] = []
        for c, code, n1, n2 in trio_names: # l n1 n2 not used here
            if df.loc[seq,code] in d_col[c]: # valid value for column
                d_seq_color[seq].append(d_col[c][ df.loc[seq,code] ]) # color of the code column
                d_seq_label[seq].append(str(df.loc[seq,c])) # label column
            else:
                d_seq_color[seq].append(col_missing)
                d_seq_label[seq].append(str(" "))
    return [d_seq_color, d_seq_label, width_and_names] 

def return_treestyle_with_columns (cmapvector):
    '''
    Need column names again to print header in order 
    '''
    [d_seq_color, d_seq_label, width_and_names] = cmapvector
    rect_width   = [x[0] for x in width_and_names]
    column_names = [x[1] for x in width_and_names]
    label_font_size = 7

    # default node (not used since it's lost with any customisation, so we create all node styles independently)
    ns1 = ete3.NodeStyle()
    ns1["size"] = 1;  ns1["shape"] = "square" ; ns1["fgcolor"] = "101010" 
    ns1["hz_line_type"]  = ns1["vt_line_type"]  = 0 # 0=solid, 1=dashed, 2=dotted
    ns1["hz_line_color"] = ns1["vt_line_color"] = "darkred"
    def tree_profile_layout (node):# prepare table and other node information (local function so mind the identation) 
        if "NORW" in (getattr(node, "submission_org_code")):  this_color = "darkred"
        else:     this_color = "#080816"

        node.img_style['hz_line_type']  = node.img_style['vt_line_type'] = 0 # 0=solid, 1=dashed, 2=dotted
        node.img_style['hz_line_width'] = node.img_style['vt_line_width'] = 4 
        node.img_style['hz_line_color'] = node.img_style['vt_line_color'] = this_color

        if node.is_leaf(): # the aligned leaf is "column 0", thus traits go to column+1
            node.img_style['size'] = 2; node.img_style['shape'] = "sphere";  node.img_style['fgcolor'] = this_color 
            ete3.add_face_to_node(ete3.AttrFace("name", fsize=label_font_size, text_suffix="   "), node, 0, position="aligned")
            for column, (rgb_val, lab, wdt) in enumerate(zip(d_seq_color[node.name], d_seq_label[node.name], rect_width)): 
                label = {"text": lab[:10], "color":"Black", "fontsize": label_font_size-1}
                ete3.add_face_to_node (ete3.RectFace (wdt, 12, fgcolor=rgb_val, bgcolor=rgb_val, label = label), node, 2*column+1, position="aligned")
                ete3.add_face_to_node (ete3.RectFace (2, 12, fgcolor="#ffffff", bgcolor="#ffffff", label = ""), node, 2*column+2, position="aligned")
        else:
            node.img_style['size'] = 0; 

    ts = ete3.TreeStyle()
    ts.draw_guiding_lines = True # dotted line between tip and name
    ts.guiding_lines_color = "#f4f4f4" # "#bdb76d"
    ts.guiding_lines_type = 2 #  0=solid, 1=dashed, 2=dotted
    ts.layout_fn = tree_profile_layout
    ts.branch_vertical_margin = 0
    ts.min_leaf_separation = 1 # Min separation, in pixels, between two adjacent branches
    ts.scale = 2000000 # 2e6 pixels per branch length unit (i.e. brlen=1 should be how many pixels?)
    ts.show_scale = False
    show_branch_length = True
    ts.show_leaf_name = False # we handle this in the layout function

    ## STILL dont know how to do it
    #ts.legend.add_face(CircleFace(10, "red"), column=0)
    #ts.legend.add_face(TextFace("0.5 support"), column=1)
    #ts.legend_position = 3 #  TopLeft corner if 1, TopRight if 2, BottomLeft if 3, BottomRight if 4
    for col, label in enumerate(column_names): # the first are tip labels
        labelFace = ete3.TextFace(label, fsize=9, fgcolor="DimGray") # fsize controls interval betweel columns
        labelFace.rotation = 270
        labelFace.vt_align = 1  # 0 top, 1 center, 2 bottom
        labelFace.hz_align = 1  # 0 left, 1 center, 2 right 
        ts.aligned_header.add_face(labelFace, 2 * col + 1)
    return ts

def ASR_subtrees (metadata0, tree, extended_mode = 0, reroot = True, method = None):
    """ ancestral state reconstruction assuming binary or not. Clusters are based on "locality": patient is from Norfolk
    *or* was sequenced here in NORW. One column at a time, so the n_threads doesn't help.
    """
    if method is None: method = ["DOWNPASS", "ACCTRAN"]
    if isinstance (method, str): method = [method, method]

    metadata = metadata0.copy()  # work with copies
    csv_cols = [x for x in common.asr_cols if x in metadata.columns]
    if reroot: ## this should be A,B 
        R = tree.get_midpoint_outgroup()
        tree.set_outgroup(R)

    md_description="""
Sequence clusters are based on **locality**: patients from Norfolk (field `adm2`) *or* patients that were sequenced here
(submission org = `NORW`).
This definition is not equivalent to the `UK lineages`, which relies on sequence properties and is estimated at the
national level. 
Our method does not take the `UK lineages` explicitly into account and as result we may find clusters spanning several
lineages, and we may also see samples from the same lineage scattered across clusters (reflecting their locality).
The **locality** allows us to focus on the local scale, by "zooming in" into geographycally connected lineages.
<br>
    """
    if extended_mode == 0:
        yesno = ( (metadata["adm2"].str.contains("Norfolk", case=False, na=False))
                | (metadata["submission_org_code"].str.contains("NORW", case=False, na=False)))
    else:   
        md_description = "Extended mode of **peroba**, for non-COGUK analyses<br><br>"
        yesno = (metadata["submission_org_code"].str.contains("NORW", case=False, na=False))

    logger.info("Start estimating ancestral states by %s for locality and %s for others", method[0], method[1])
    df = pd.DataFrame(list(yesno.astype(str)), columns=["local"], index=metadata.index.copy())
    print ("DEBUG::yesno and local:: ", collections.Counter(yesno.astype(str)), collections.Counter(df["local"]))
    x = get_binary_trait_subtrees (tree, df, trait_column = "local", trait_value = "True", elements = 1,
            method=method[0], extended_mode = extended_mode)
    subtree, mono, result, trait_name = x  # most important is subtree
    logger.info("Finished estimating ancestral state for 'locality', which defines clusters")
    
    ## decorate the tree with ancestral states
    csv, csv_cols = prepare_csv_columns_for_asr (metadata, csv_cols) # csv is used only in pastml (imputation go to tree)
    if (csv_cols):
        logger.info("Will now estimate ancestral states for %s", " ".join(csv_cols))
        tree_leaf_nodes = {leaf:leaf.name for leaf in tree.iter_leaves()} # in case we have duplicated names 
        result = acr (tree, csv, prediction_method = method[1], force_joint=False) ## annotates tree nodes with states (e.g. tre2.adm2)
        for leafnode, leafname  in tree_leaf_nodes.items(): # pastml (correctly) replaces duplicated names by a placeholder like "t123" 
            leafnode.name = leafname  # reverts back to duplicate names 
    # adds new peroba_ columns with imputed and original values:
    metadata = save_metadata_inferred (metadata, tree, csv_cols)

    node2leaves = tree.get_cached_content(store_attr="name") # dict node:[leaves]
    submetadata = [metadata.loc[metadata.index.isin(node2leaves[x])] for x in subtree]
    # not currently used
    supermetadata = [metadata.loc[metadata.index.isin(node2leaves[x.up])] for x in subtree if x.up]

    return submetadata, subtree, md_description, metadata  ## metadata now has tip reconstruction

def prepare_csv_columns_for_asr (csv, csv_cols=None):
    if csv_cols is None:  csv_cols = common.asr_cols # csv_cols = estimate tips, add "peroba_" to name and export
    csv_cols = [x for x in csv_cols if x in csv.columns]
    csv = csv[csv_cols + ["submission_org_code"]]; ## dataframe csv will also have submission_org_code column (for tree colouring)
    #csv["ICU_admission"] = csv["ICU_admission"].replace("Unknown", "", regex=True)

    csv.loc[~csv["submission_org_code"].str.contains("NORW", na=False), "date_sequenced"] = "nil" # only for NORW

    csv_cols = [x for x in csv_cols if x in csv.columns]
    for col in csv_cols:
        csv[col].fillna("", inplace=True) ## to estimate tip values
        csv[col] = csv[col].astype(str)

    logger.info ("Follow-up columns found in CSV: %s", " ".join(csv_cols))
    return csv, csv_cols

def plot_single_cluster (csv, tree, counter, ts = None, output_dir=None, figdir=None):
    if output_dir is None: output_dir = cwd
    exclude=["Postcode","Repeat Sample ID", "acc_lineage", "submission_org_code", "subsample_omit","swab_site","title","url","virus"]
    exclude = [x for x in csv.columns if x in exclude]

    df = csv.drop (labels = exclude, axis=1)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    caption = f"Tree of cluster {counter}, with red branches connecting samples sequenced in NORW."

    fname = f"{figdir}/tree{counter}.pdf"
    tree.render(os.path.join(output_dir,fname), w=800, tree_style=ts)
    pdf_desc = f"\n![{caption}]({fname})\n\n"
    fname = f"{figdir}/tree{counter}.png"
    tree.render(os.path.join(output_dir,fname), w=1400, tree_style=ts)
    html_desc = f"\n![{caption}]({fname}){{ width=100% }}\n\n<br>"

    # report_md += df.to_html(max_rows=4, max_cols=20)
    # report_md += "<hr>\n"

    return html_desc, pdf_desc 

def save_metadata_inferred (df, tree, csv_cols = None):
    ## write only a copy of the original column, since here we also mix imputed values
    if csv_cols is None: csv_cols = common.asr_cols
    prefix = "peroba_"
    for col in csv_cols:
        df[prefix+col] = df[col] # same information as original 

    for leaf in tree.iter_leaves():
        for col in csv_cols:
            x = getattr(leaf,col)
            df.loc[str(leaf.name), prefix+col] = "/".join([str(i) for j,i in enumerate(x) if j < 3])
#            if pd.isnull(df.loc[str(leaf.name), col]): ## just impute if it was nan
#                df.loc[str(leaf.name), col] = "/".join([str(i) for i in x]) # this is original column
    return df

#cmap_list = ["Accent", "Dark2", "cividis", "jet", "hsv", "viridis", "plasma", "rainbow", "nipy_spectral"]

