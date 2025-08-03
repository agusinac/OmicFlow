sparse_to_dtable <- function(sparsemat) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(sparsemat, "sparseMatrix"))
    cli::cli_abort("sparsemat must be a sparseMatrix.")

  ## MAIN
  #--------------------------------------------------------------------#

  return(data.table::data.table(as.matrix(sparsemat)))
}

read_tsv_matrix <- function(filepath) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filepath))
    cli::cli_abort("{filepath} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  distmat <- check_input(filepath)

  rownames <- distmat[[1]]
  colnames <- names(distmat)[-1]
  values <- as.matrix(distmat[, -1, with = FALSE])

  sparse_matrix <- sparseMatrix(
    i = rep(1:nrow(values), ncol(values)),
    j = rep(1:ncol(values), each = nrow(values)),
    x = as.vector(values)
  )
  rownames(sparse_matrix) <- rownames
  colnames(sparse_matrix) <- colnames

  return(sparse_matrix)
}

read_rarefraction_qiime <- function(filepath) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filepath))
    cli::cli_abort("{filepath} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  df_shannon <- data.table::fread(filepath)

  # Pivot into long table
  shannon_long <- data.table::melt(data = df_shannon,
                                   measure.vars = colnames(df_shannon)[grepl("depth-", colnames(df_shannon))],
                                   variable.name = "iters",
                                   variable.factor = FALSE,
                                   value.name = "alpha_div")
  # Corrects colnames
  colnames(shannon_long) <- c("SAMPLE_ID", "iters", "alpha_div")

  return(shannon_long)
}

column_exists <- function(column, table) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.character(column) && length(column) != 1)
    cli::cli_abort("{column} needs to contain characters with length of 1.")

  if (!inherits(table, "data.frame") && !inherits(table, "data.table"))
    cli::cli_abort("table must be a data.frame or data.table.")

  ## MAIN
  #--------------------------------------------------------------------#

  valid_columns <- column[column %in% colnames(table)]

  if (length(valid_columns) == 0 ) {
    return(FALSE)
  }

  # For each existing column, check if it's *not entirely NA*
  columns_empty <- all(sapply(valid_columns, function(col) {
    any(!is.na(table[[col]]))
  }))

  return ( length(valid_columns) == length(column) && columns_empty )
}

parse_commandline <- function() {
  option_list <- list (
    optparse::make_option("--omics",
                          action = "store",
                          type = "character",
                          default = "metagenomics",
                          help="tab seperated file"),
    optparse::make_option(c("-m", "--metadata"),
                          action = "store",
                          type = "character",
                          help="tab seperated file"),
    optparse::make_option(c("-b", "--biom"),
                          action = "store",
                          type = "character",
                          help="biom format file"),
    optparse::make_option(c("-t", "--tree"),
                          action = "store",
                          type = "character",
                          help="Phylogenetic tree in newick format"),
    optparse::make_option(c("-c", "--cpus"),
                          action = "store",
                          type = "numeric",
                          help="Number of cores",
                          default = 4),
    optparse::make_option(c("-o", "--outdir"),
                          action = "store",
                          type = "character",
                          help="Output directory",
                          default = normalizePath(getwd())),
    optparse::make_option(c("-f", "--filename"),
                          action = "store",
                          type = "character",
                          help="Name of the HTML report",
                          default = "report.html"),
    optparse::make_option(c("--i-beta-div"),
                          action = "store",
                          type = "character",
                          help="custom beta diversity from qiime2"),
    optparse::make_option(c("--i-alpha-div"),
                          action = "store",
                          type = "character",
                          help="custom alpha diversity with rarefraction from qiime2")
  )

  parser <- optparse::OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  return(arguments$options)
}

update_citations_md <- function() {
  # Get imported packages from DESCRIPTION
  imports <- desc::desc_get_deps("DESCRIPTION")
  imported_pkgs <- imports$package[imports$type == "Imports"]

  # Open output file
  outfile <- file("CITATION.md", "w")
  writeLines("# Citations for Imported Packages\n", outfile)

  for (pkg in imported_pkgs) {
    writeLines(paste0("## ", pkg, "\n"), outfile)

    # Try to get citations
    cites <- tryCatch(utils::citation(pkg), error = function(e) NA)

    if (is.na(cites)) {
      writeLines("No citation found.\n", outfile)
      next
    }

    # Process each citation entry using BibTeX
    bibtex_entry <- tryCatch(utils::toBibtex(cites), error = function(e) NA)

    if (all(is.na(bibtex_entry))) {
      writeLines("Could not generate BibTeX citation.\n", outfile)
      next
    }

    # Write in a code block
    writeLines("```", outfile)
    writeLines(as.character(bibtex_entry), outfile)
    writeLines("```", outfile)
    writeLines("\n", outfile)
  }

  # Close file
  close(outfile)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.character(x)) {
    return(FALSE)
  } else {
    abs(x - round(x)) < tol
  }
}

combine_conditions <- function(condition1, condition2) {
  if (!is.null(condition1) && !is.null(condition2)) {
    if (!inherits(condition1, "data.frame") && !inherits(condition1, "data.table"))
      cli::cli_abort("condition1 must be a data.frame or data.table.")

    if (!inherits(condition2, "data.frame") && !inherits(condition2, "data.table"))
      cli::cli_abort("condition2 must be a data.frame or data.table.")
  }

  # Combine to strings for easy comparison
  cond1_str <- paste(
    pmin(condition1$group1, condition1$group2),
    pmax(condition1$group1, condition1$group2), 
    sep = "_")

  cond2_str <- paste(
    pmin(condition2$group1, condition2$group2),
    pmax(condition2$group1, condition2$group2), 
    sep = "_")

  # Find which condition2 are NOT already in condition1
  new_pairs_idx <- !cond2_str %in% cond1_str

  if (any(new_pairs_idx)) {
    # There are new pairs in condition2 not in condition1;
    # append only the new ones
    new_rows <- condition2[new_pairs_idx, ]
    updated_conditions <- rbind(condition1, new_rows)
  } else {
    # All pairs in condition2 are already included
    updated_conditions <- condition1
  }

  return(updated_conditions)
}