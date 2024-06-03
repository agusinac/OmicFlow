# Specify filenames ------------------------------------------------------------
ps <- phyloseq_from_biom(biom_data = data_00$biomData, 
                         meta_data = data_00$metaData, 
                         tree_data = data_00$treeData, 
                         refseq_data = data_00$refseqData)
