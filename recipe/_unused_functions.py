
## the commented_out loop tries to fix raw COGUK sequence names, that do not follow same pattern : 
# seq names are a mess, it's better to map using metadata as key (however we have keys as "D02" etc
# COG-UK/SHEF-D34A0/SHEF:20200423_1347_X1_FAN43269_d941f767|SHEF-D34A0|... 
# COG-UK/OXON-AF08B/OXON
# EPI_ISL_422016|hCoV-19/Wales/PHWC-26796/2020|Wales|WALES|2020-03-28|PHWC|...
# hCoV-19/England/201140442/2020|PENDING|England|GREATER_LONDON|2020-03-10|PHEC|2
# hCoV-19/England/CAMB-74D3D/2020||England|CAMBRIDGESHIRE|2020-03-19|  
# also hCoV Hcov hcoV etc.

def merge_data_sequence (self): # keep only intersection
    if self.sequences is None or self.metadata is None: return
    id_meta = self.metadata.index.array # <PandasArray> object, works like an np.array
    id_seqs = [x for x in self.sequences.keys()] # copy needed since I'll delete elements (modifying dict)

    ## remove sequences absent from the metadata
    idx_sequence_only = set (id_seqs) - set (id_meta) # names unique to sequence, not found in metadata
    idx_matched = set(id_seqs) & set(id_meta)
    df_matched = self.metadata.loc[idx_matched,] ## these rows have well-behaved sequences, with same name 

    unique_id = "central_sample_id" # this is an alternative, unique ID, that is always (?) found within seq names
    _this_is_a_test_ = False  ## problem happens with "raw" elan.consensus.fasta but not with (1 week older) phylogenetics/alignment
    if _this_is_a_test_ and unique_id in self.metadata.columns: # somes thing is wrong if this column is missing from COGUK metadata
        idx_metadata_only = set (id_meta) - set (id_seqs) # names unique to metadata, not found in sequence
        ## dataframe with unique_id of rows without sequences: (notice the list [] of one element, o.w. returns Series
        df_u = self.metadata.loc[idx_metadata_only,[unique_id]] 
        found_indices = []
        # map current seqname to a string likely to have central_sample_id
        description_leaves = {sn:self.sequences[sn].description.replace("EPI_ISL","").replace("COG-UK/","")[:50] for sn in idx_sequence_only} 

        print ("DEBUG:: start test, length:", df_u.shape) # this block can be safely removed (bugged), loop below is slow
        desc_dict = {self.sequences[sn].description.replace("EPI_ISL","").replace("COG-UK/","")[:50]:sn for sn in idx_sequence_only} 
        df1 = pd.DataFrame.from_dict(desc_dict, orient='index', columns=["description"])
        df_u = df_u[df_u[unique_id].isin(df1.index)] ## FIXME: length zero (none found)
        print ("DEBUG:: end test, length:", df_u.shape)

        print ("DEBUG:: entering loop")
        for s_id, s_desc in description_leaves.items(): ## try to find sequence through "central_sample_id", a short code in metadata
            for index, row in df_u.iterrows(): # loop over rows of unmatched metadata
                if index not in found_indices and str(row[unique_id]) in s_desc: # unique_id was found within the sequence long name
                    self.sequences[index] = self.sequences[s_id]
                    self.sequences[index].id = index
                    #SeqRecord(Seq.Seq(str(xq.seq),Alphabet.IUPAC.ambiguous_dna),id=str(index),description=xq.description)
                    # update dataframes to speed up next iteration
                    found_indices.append(index)
                    if len(found_indices)%10 == 0: print (".", end="")
                    break # next seqname
        print ("DEBUG:: left loop")
        df_matched = df_matched.append(self.metadata.loc[found_indices,]) ## these rows have well-behaved sequences, with same name 

    # sequences found will be new SeqRecords, old ones can be removed
    for seqname in idx_sequence_only:
        del self.sequences[seqname] # delete this dict element
    self.metadata = df_matched
    print_redblack ("LOG:: new dataframe size, after removing missing sequences:", self.metadata.shape)
    # remove table rows without sequence information
    samples_found_in_sequences = self.metadata.index.isin([x for x in self.sequences.keys()])
    self.metadata = self.metadata[samples_found_in_sequences]
    print_redblack ("LOG:: new table size is (after excluding missing sequences):", sum(samples_found_in_sequences))

    self.add_sequence_counts_to_metadata (from_scratch = False) ## add column if missing
