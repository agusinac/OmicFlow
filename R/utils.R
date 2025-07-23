sparse_to_dtable <- function(sparsemat) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(sparsemat, "sparseMatrix"))
    cli::cli_abort("sparsemat must be a sparseMatrix.")

  ## MAIN
  #--------------------------------------------------------------------#

  return(data.table::data.table(as.matrix(sparsemat)))
}

read_tsv_matrix <- function(filename) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filename))
    cli::cli_abort("{filename} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  distmat <- data.table::fread(filename,
                               header = TRUE)

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

read_rarefraction_qiime <- function(filename) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filename))
    cli::cli_abort("{filename} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  df_shannon <- data.table::fread(filename)

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

read_sparseTable <- function(filename) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filename))
    cli::cli_abort("{filename} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  # Read text file, supports csv, excel and tsv formats
  dt <- data.table::fread(filename)
  dt[, (names(dt)) := lapply(.SD, function(x) {
    x <- gsub("\\s+", "", x)                      # Removes spaces between strings
    x <- gsub("^[A-Za-z]*", "", x)                # Removes letters
    })]

  # Convert to matrix format
  mat_1 <- as.matrix(dt,
                     rownames = rownames(dt),
                     colnames = colnames(dt))

  # Change character values to numeric
  mat_2 <- matrix(data = as.numeric(mat_1),
                  ncol = ncol(dt))
  colnames(mat_2) <- colnames(dt)

  mat_2[is.na(mat_2) | mat_2 == ""] <- 0          # Empty strings from cleaning step

  # Return sparseMatrix
  return(as(mat_2, "sparseMatrix"))
}

column_exists <- function(column, table) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.character(column) && length(column) != 1)
    cli::cli_abort("{column} needs to contain characters with length of 1.")

  if (!inherits(table, "data.frame") || !inherits(table, "data.table"))
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
  option_list <- list (optparse::make_option(c("-m", "--metadata"),
                                             action = "store",
                                             help="tab seperated file"),
                       optparse::make_option(c("-b", "--biom"),
                                             action = "store",
                                             help="biom format file"),
                       optparse::make_option(c("-t", "--tree"),
                                             action = "store",
                                             help="Phylogenetic tree in newick format"),
                       optparse::make_option(c("-c", "--cpus"),
                                             action = "store",
                                             help="Number of cores",
                                             default = 4),
                       optparse::make_option(c("-o", "--outdir"),
                                             action = "store",
                                             help="Output directory",
                                             default = normalizePath(getwd())),
                       optparse::make_option(c("--i-beta-div"),
                                             action = "store",
                                             help="custom beta diversity from qiime2"),
                       optparse::make_option(c("--i-alpha-div"),
                                             action = "store",
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
