### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------

## Global character vector to hold default names for circular sequences.
DEFAULT_CIRC_SEQS = c("chrM","MT","mit","2micron","2-micron");


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DB related.
###
### Most of this stuff was copy/pasted from AnnotationDbi (trying to avoid
### depending on AnnotationDbi for now).
### It is NOT exported (to avoid collisions with AnnotationDbi).
###

### Environment for storing run-time objects
RTobjs <- new.env(hash=TRUE, parent=emptyenv())

assign("debugSQL", FALSE, envir=RTobjs)

debugSQL <- function()
{
    debugSQL <- !get("debugSQL", envir=RTobjs)
    assign("debugSQL", debugSQL, envir=RTobjs)
    debugSQL
}


### Use dbQuery(conn, SQL, 1) instead of dbQuery(conn, SQL)[[1]],
### it's much safer!
dbEasyQuery <- function(conn, SQL, j0=NA)
{
    if (get("debugSQL", envir=RTobjs)) {
        if (!is.character(SQL) || length(SQL) != 1L || is.na(SQL))
            stop("[debugSQL] 'SQL' must be a single string")
        cat("[debugSQL] SQL query: ", SQL, "\n", sep="")
        st <- system.time(data0 <- dbGetQuery(conn, SQL))
        cat("[debugSQL]      time: ", st["user.self"], " seconds\n", sep="")
    } else {
        data0 <- dbGetQuery(conn, SQL)
    }
    if (is.na(j0))
        return(data0)
    ## Needed to deal properly with data frame with 0 column ("NULL data
    ## frames with 0 rows") returned by RSQLite when the result of a SELECT
    ## query has 0 row
    if (nrow(data0) == 0L)
        character(0)
    else
        data0[[j0]]
}

dbEasyPreparedQuery <- function(conn, SQL, bind.data)
{
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the early bail out.
    if (nrow(bind.data) == 0L)
        return()
    if (get("debugSQL", envir=RTobjs)) {
        if (!is.character(SQL) || length(SQL) != 1L || is.na(SQL))
            stop("[debugSQL] 'SQL' must be a single string")
        cat("[debugSQL] SQL prepared query: ", SQL, "\n", sep="")
        cat("[debugSQL]     dim(bind.data): ",
            paste(dim(bind.data), collapse=" x "), "\n", sep="")
        st <- system.time({
                  dbBeginTransaction(conn)
                  dbGetPreparedQuery(conn, SQL, bind.data)
                  dbCommit(conn)})
        cat("[debugSQL]               time: ", st["user.self"],
            " seconds\n", sep="")
    } else {
        dbBeginTransaction(conn)
        dbGetPreparedQuery(conn, SQL, bind.data)
        dbCommit(conn)
    }
}

.dbFileConnect <- function(dbfile)
{
    if (!file.exists(dbfile))
        stop("DB file '", dbfile, "' not found")
    dbConnect(SQLite(), dbname=dbfile, cache_size=64000, synchronous=0,
              flags=SQLITE_RO)
}

get_dbfile <- function(libname, pkgname)
{
    filename <- paste(pkgname, ".sqlite", sep="")
    system.file("extdata", filename, package=pkgname, lib.loc=libname)
}

get_dbconn <- function(libname, pkgname)
{
    dbfile <- get_dbfile(libname, pkgname)
    .dbFileConnect(dbfile)
}

get_cached_dbfile <- function(datacache)
{
    if (!exists("dbfile", envir=datacache))
        stop("symbol \"dbfile\" not found in 'datacache'")
    get("dbfile", envir=datacache)
}

get_cached_dbconn <- function(datacache)
{
    if (!exists("dbconn", envir=datacache))
        stop("symbol \"dbconn\" not found in 'datacache'")
    get("dbconn", envir=datacache)
}

get_dbtable <- function(tablename, datacache)
{
    if (!exists(tablename, envir=datacache)) {
        dbtable <- dbReadTable(get_cached_dbconn(datacache), tablename, row.names=NULL)
        assign(tablename, dbtable, envir=datacache)
    }
    get(tablename, envir=datacache)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data frame related (NOT exported).
###
### TODO: Find a better home for these low-level data.frame utils.
###

### Not data.frame specific. Would work on any matrix-like object.
hasCol <- function(x, colnames) {colnames %in% colnames(x)}

### Sets the class of (all or some of) the columns of a data.frame.
### Typical use:
###   x <- setDataFrameColClass(x, c(colA="integer", colB="factor"))
### Note that if 'x' has more than one "colA" col, then *all* of them are
### coerced to integer.
setDataFrameColClass <- function(x, col2class, drop.extra.cols=FALSE)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    if (!is.character(col2class) || is.null(names(col2class)))
        stop("'col2class' must be a named character vector")
    if (!all(names(col2class) %in% colnames(x)))
        stop("'col2class' has invalid names")
    if (!isTRUEorFALSE(drop.extra.cols))
        stop("'drop.extra.cols' must be TRUE or FALSE")
    if (drop.extra.cols) {
        col_idx <- which(colnames(x) %in% names(col2class))
    } else {
        col_idx <- seq_len(ncol(x))
    }
    tmp <- lapply(col_idx,
                  function(j)
                  {
                      col <- x[[j]]
                      colname <- colnames(x)[j]
                      if (!(colname %in% names(col2class)))
                          return(col)
                      class <- col2class[[colname]]
                      FUNname <- paste("as", class, sep=".")
                      if (exists(FUNname) && is.function(FUN <- get(FUNname)))
                          return(FUN(col))
                      as(col, class)
                  })
    names(tmp) <- colnames(x)[col_idx]
    return(data.frame(tmp, check.names=FALSE, stringsAsFactors=FALSE))
}

### Acts like an SQL *inner* join.
### 'x' must be a data frame. 'name2val' must be a named atomic vector or
### a named factor. 'join_colname' must be the name of the col in 'x' whose
### values are matched against 'names(name2val)'. 'vals_colname' must be the
### name of the col that will be populated with the appropriate 'name2val'
### vals and bound to 'x'.
### Note that this acts like an SQL *inner* join, not a *left* join, i.e.
### rows in 'x' that can't be mapped to a value in 'name2val' are dropped.
joinDataFrameWithName2Val <- function(x, join_colname, name2val, vals_colname)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    if (!isSingleString(join_colname) || !hasCol(x, join_colname))
        stop("'join_colname' must be a valid colname for 'x'")
    if (!is.vector(name2val) && !is.factor(name2val))
        stop("'name2val' must be a vector (or factor)")
    if (!is.atomic(name2val) || is.null(names(name2val)))
        stop("'name2val' must be atomic and have names")
    if (!isSingleString(vals_colname))
        stop("invalid 'vals_colname'")
    join_col <- as.character(x[[join_colname]])
    common_names <- intersect(join_col, names(name2val))
    name2val <- name2val[names(name2val) %in% common_names]
    x <- x[join_col %in% common_names, ]
    tmp <- split(as.vector(name2val), names(name2val))
    ## as.character() is required below just because 'x[[join_colname]]'
    ## could be a factor (subsetting by a factor is equivalent to subsetting
    ## by an integer vector but this is not what we want here).
    tmp <- tmp[as.character(x[[join_colname]])]
    x <- x[rep.int(seq_len(nrow(x)), elementNROWS(tmp)), ]
    row.names(x) <- NULL
    if (nrow(x) == 0L)
        vals <- name2val[FALSE]
    else if (is.factor(name2val))
        vals <- factor(unname(unlist(tmp)), levels=levels(name2val))
    else
        vals <- unname(unlist(tmp))
    x[[vals_colname]] <- vals
    x
}

### Returns the vector of ids such that 'unique(x)[ids, ]' is identical
### to 'x' (in the same way that 'levels(f)[f]' is identical to
### 'as.vector(f)' when 'f' is a character factor).
### This unambiguously defines 'ids'. In particular, it's not Locale
### specific, despite the fact that the current implementation uses a
### sorting approach.
makeIdsForUniqueDataFrameRows <- function(x)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    x_order <- do.call(order, x)
    x_dups <- duplicated(x)
    ## First we make "provisory" ids. Those ids *are* Locale specific.
    prov_ids <- integer(nrow(x))
    prov_ids[x_order] <- cumsum(!x_dups[x_order])
    ## Convert the "provisory" ids into the final ids. The final ids are
    ## *not* Locale specific anymore.
    as.integer(factor(prov_ids, levels=unique(prov_ids)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous (NOT exported).
###

### AFAIK UCSC doesn't flag circular sequences.
### As of Sep 21, 2010 (Ensembl release 59), Ensembl was still not flagging
### circular sequences in their db (see this thread for the details
### http://lists.ensembl.org/pipermail/dev/2010-September/000139.html),
### This just takes the list of things that users are calling circular in the
### circ_seqs argument and then marks those things as being circular and
### returns the vector all marked up

### TODO: still need to get the new parameter passed along to where this is called...  :P
matchCircularity <- function(seqnames, circ_seqs)
{
    ## shorten and put to lowercase (for simplicity in subsequent comparisons)
    seqs <- tolower(seqnames)
    circs <- tolower(circ_seqs)
    ## checks
    if(length(intersect(seqs,circs))<1 && length(circs)>0){
      warning("None of the strings in your circ_seqs argument match your seqnames.")  
    }
    int <- intersect(seqs,circs)
    is_circular <- rep.int(FALSE, length(seqs))
    if(length(int)>0){
      for(i in seq_len(length(int))){
        idx <- grep(int[i], seqs)
        is_circular[idx] <- TRUE
      }
    }
    is_circular
}

### 'exon_count' must be a vector of positive integers and 'tx_strand' a
### character vector with "+" or "-" values. Both vectors must have the
### same length.
makeExonRankCol <- function(exon_count, tx_strand)
{
    ans <- lapply(seq_len(length(exon_count)),
        function(i)
        {
            if (tx_strand[i] == "+")
                seq_len(exon_count[i])
            else
                (exon_count[i]):1L
        }
    )
    unlist(ans)
}



#' Read the cached result of an expression from a locally cached file if it exists,
#' else evaluate the expression, cache it, and return result.
#' 
#' @param expression The expression to be evaluated and returned if the cached result does not exist.
#' @param local_cache_path The path that stores the cached objects. If NULL, no caching is performed.
#' @param object_name The name of the object used to read/create the cache file.
#' @return The result of evaluating the expression, either read from cache or from actually evaluating it.
#' 
#' @note If evaluating the expression has side-effects,
#' they will not happen if the result is returned from cache.
#' 
#' @examples
#' read_or_update_local_cache(2+2, temp_dir(), "TwoPlusTwo")
read_or_update_local_cache <- function(expression, local_cache_path, object_name)
{
  if (is.null(local_cache_path))
  {
    return (expression)
  }
  else
  {
    cached_filepath = paste0(local_cache_path, "/", object_name, ".rds")
    if (file.exists(cached_filepath))
    {
      return (readRDS(cached_filepath))
    }
    else
    {
      if (!dir.exists(local_cache_path)) { dir.create(local_cache_path, recursive=TRUE) }
      result = expression
      saveRDS(result, file=cached_filepath)
      return (result)
    }
  }
}

read_or_update_local_cacheDb <- function(expression, local_cache_path, object_name)
{
  if (is.null(local_cache_path))
  {
    return (expression)
  }
  else
  {
    cached_filepath = paste0(local_cache_path, "/", object_name, ".sqlite")
    if (file.exists(cached_filepath))
    {
      return (loadDb(cached_filepath))
    }
    else
    {
      if (!dir.exists(local_cache_path)) { dir.create(local_cache_path, recursive=TRUE) }
      result = expression
      saveDb(result, file=cached_filepath)
      return (result)
    }
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Expectation: is the object equal to a reference value stored in a file?
#'
#' This expectation is equivalent to [expect_equal()], except that the
#' expected value is stored in an RDS file instead of being specified literally.
#' This can be helpful when the value is necessarily complex. If the file does
#' not exist then it will be created using the value of the specified object,
#' and subsequent tests will check for consistency against that generated value.
#' The test can be reset by deleting the RDS file.
#'
#' It is important to initialize the reference RDS file within the source
#' package, most likely in the `tests/testthat` directory. Testing spawned
#' by [devtools::test()/test_file()], for example, will accomplish this. But note
#' that testing spawned by `R CMD check` and [devtools::check()] will NOT.
#' In the latter cases, the package source is copied to an external location
#' before tests are run. The resulting RDS file will not make its way back into
#' the package source and will not be available for subsequent comparisons.
#'
#' @inheritParams expect_that
#' @param file The file name used to store the object. Should have an "rds"
#'   extension.
#' @param label For the full form, a label for the expected object, which is
#'   used in error messages. Useful to override the default (which is based on
#'   the file name), when doing tests in a loop. For the short-cut form, the
#'   object label, which is computed from the deparsed object by default.
#' @param expected.label Equivalent of `label` for shortcut form.
#' @param on.update If non-NULL and if the reference file does not exist,
#'  [on.update(a, b)] will be called with 2 parameters:
#'  1) the object itself, and
#'  2) a name composed of the reference file name and the object name.
#'  If on.update is NULL and the reference file does not exist, a short message
#'  will be printed to alert the tester that the reference file was created instead of tested.
#' @param ... other values passed to [expect_equal()]
#' @family expectations
#' @import testthat
#' @examples
#' \dontrun{
#' expect_equal_to_reference(1, "one.rds")
#' }
expect_equal_to_reference = function(object, file, ..., info=NULL, label=NULL, expected.label=NULL,
                                     on.update=NULL) {
  lab_act <- testthat::make_label(object, label)
  lab_exp <- expected.label %||% paste0("reference from `", file, "`")
  
  if (!file.exists(file)) {
    # first time always succeeds
    saveRDS(object, file)

    if (is.null(on.update))
      print(paste0("Updated reference file \"", file, "\" for object ", deparse(substitute(object))))
    else
      on.update(object, paste0(file, ":", deparse(substitute(object))))

    succeed()
  } else {
    reference <- readRDS(file)
    
    comp <- testthat::compare(object, reference, ...)
    expect(
      comp$equal,
      sprintf("%s not equal to %s.\n%s", lab_act, lab_exp, comp$message),
      info = info
    )
  }
  
  invisible(object)
}


# @return full path to this script
#' current script file (in full path)
#' @return The path of the currently executing script or code snippet.
#' @examples
#' current_script_file()
#' @note Works with Rscript, source() or in RStudio Run selection
current_script_file <- function() {
  # http://stackoverflow.com/a/32016824/2292993
  cmdArgs = commandArgs(trailingOnly = FALSE)
  needle = "--file="
  match = grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript via command line
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    ls_vars = ls(sys.frames()[[1]])
    if ("fileName" %in% ls_vars) {
      # Source'd via RStudio
      return(normalizePath(sys.frames()[[1]]$fileName)) 
    } else {
      if (!is.null(sys.frames()[[1]]$ofile)) {
        # Source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
      } else if (rstudioapi::isAvailable() && nchar(rstudioapi::getActiveDocumentContext()$path) > 0) {
        # RStudio Run Selection
        # http://stackoverflow.com/a/35842176/2292993  
        return(normalizePath(rstudioapi::getActiveDocumentContext()$path))
      } else {
        return("")
      }
    }
  }
}
