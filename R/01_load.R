# Construct phyloseq -----------------------------------------------------------
ps <- phyloseq_from_biom(biom_data = data_00$biomData, 
                         meta_data = data_00$metaData, 
                         tree_data = data_00$treeData, 
                         refseq_data = data_00$refseqData)


# Fetch headers ----------------------------------------------------------------
meta_tab <- get_meta(ps)
meta_tab.colnames <- colnames(meta_tab)

# Saves columns as dataframes
RANKSTAT_data <- select(meta_tab, meta_tab.colnames[grepl("RANKSTAT_", meta_tab.colnames)])
CORRELATION_data <- select(meta_tab, meta_tab.colnames[grepl("CORRELATION_", meta_tab.colnames)])
PAIREDGROUPBY_data <- select(meta_tab, meta_tab.colnames[grepl("PAIREDGROUPBY_", meta_tab.colnames)])

# Fetch additional files from preprocessing pipeline ---------------------------
preprocess_stats <- readr::read_delim("automated-omics-analysis/data/read_stats.txt")
shannon_file <- readr::read_csv(file = "automated-omics-analysis/data/shannon.csv")
