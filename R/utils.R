### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------

## Global character vector to hold default names for circular sequences.
DEFAULT_CIRC_SEQS = c("chrM","MT","mit","2micron","2-micron");


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


.complements <- c("A"="T","T"="A",
                  "G"="C","C"="G",
                  "M"="K","K"="M",
                  "R"="Y","Y"="R",
                  "W"="W","S"="S",
                  "D"="H","H"="D",
                  "B"="V","V"="B")

.fastComplement = function(base) .complements[base]

#' Pure R function to turn a character vector of nucleotides into a reverse complement vector.
#' When running on many strings, this is much faster than Biostrings::reverseComplement(DNAString).
#'
#' @param naString a character vector of nucleotides; IUPAC ambiguous nucleotides are supported
#'
#' @return the reverse complement of naString
#' @export
#'
#' @examples
#' fastComplement(c("GAT", "TACA")) # "CTA"  "ATGT"
#' fastComplement("MAW") # "KTW"
fastComplement = function(naString)
{
  sapply(lapply(strsplit(naString, ""), .fastComplement), paste0, collapse="")
}


#' Get URL to download protein sequence FASTA from UCSC genome browser for a given dbkey.
#'
#' @param dbkey The UCSC dbkey to get protein sequences for, e.g. hg19, hg38, mm10.
#'
#' @return A URL which can be downloaded with \code{\link{download.file}}
#' @export
#'
#' @examples
#' getProteinFastaUrlFromUCSC("hg38")
getProteinFastaUrlFromUCSC = function(dbkey)
{
    refseqTrack = ifelse(dbkey=="hg38", "refSeqComposite", "refGene")
    paste0("http://genome.ucsc.edu/cgi-bin/hgTables?db=", dbkey,
           "&hgta_geneSeqType=protein",
           "&hgta_doGenePredSequence=submit",
           "&hgta_table=refGene",
           "&hgta_track=", refseqTrack)
}

#' Get URL to download coding sequence FASTA from UCSC genome browser for a given dbkey.
#'
#' @param dbkey The UCSC dbkey to get coding sequences for, e.g. hg19, hg38, mm10.
#'
#' @return A URL which can be downloaded with \code{\link{download.file}}
#' @export
#'
#' @examples
#' getCodingFastaUrlFromUCSC("hg38")
getCodingFastaUrlFromUCSC = function(dbkey)
{
    paste0("http://genome.ucsc.edu/cgi-bin/hgTables?db=", dbkey,
           "&hgSeq.cdsExon=on&hgSeq.granularity=gene",
           "&hgSeq.casing=exon&hgSeq.repMasking=lower",
           "&hgta_doGenomicDna=get+sequence",
           "&hgta_group=genes",
           "&hgta_track=refGene",
           "&hgta_table=refGene",
           "&hgta_regionType=genome")
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
#' \dontrun{
#' read_or_update_local_cache(2+2, tempdir(), "TwoPlusTwo")
#' }
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
#' @param object object to test
#' @param info extra information to be included in the message (useful when
#'   writing tests in loops).
#' @param file The file name used to store the object. Should have an "rds"
#'   extension.
#' @param label For the full form, a label for the expected object, which is
#'   used in error messages. Useful to override the default (which is based on
#'   the file name), when doing tests in a loop. For the short-cut form, the
#'   object label, which is computed from the deparsed object by default.
#' @param expected.label Equivalent of `label` for shortcut form.
#' @param on.update If non-NULL and if the reference file does not exist,
#'  [on.update(o, t)] will be called with 2 parameters:
#'  \enumerate{
#'  \item the new object
#'  \item a title string composed of the reference file name and the object name
#'  }
#'  If on.update is NULL and the reference file does not exist, a short message
#'  will be printed to alert the tester that the reference file was created instead of tested.
#' @param on.fail If non-NULL and the tested object does not match the reference object,
#'  [on.fail(ref, new, title)] will be called with 3 parameters:
#'  \enumerate{
#'  \item the reference object
#'  \item the new object
#'  \item a title string composed of the reference file name and the object name
#'  }
#'  If on.fail is NULL, normal testthat reporting will still occur.
#' @param ... other values passed to [expect_equal()]
#' @family expectations
#' @examples
#' \dontrun{
#' expect_equal_to_reference(1, "one.rds")
#' }
expect_equal_to_reference = function(object, file, ..., info=NULL, label=NULL, expected.label=NULL,
                                     on.update=getOption("testthat.on.update"),
                                     on.fail=getOption("testthat.on.fail")) {
  lab_act <- testthat:::quasi_label(rlang::enquo(object), label)
  lab_exp <- expected.label %||% paste0("reference from `", file, "`")
  
  # check for reference files either in current directory or in ./tests/testthat directory
  if (!file.exists(file) && file.exists(file.path("tests", "testthat", file)))
    file = file.path("tests", "testthat", file)
  
  if (!file.exists(file)) {
    # first time always succeeds
    saveRDS(object, file)

    if (is.null(on.update))
      print(paste0("Updated reference file \"", file, "\" for object ", deparse(substitute(object))))
    else
      on.update(object, paste0(file, ":", deparse(substitute(object))))

    testthat::succeed()
  } else {
    reference <- readRDS(file)
    
    if (class(reference)[1] == "GRanges")
    {
        reference = GenomicRanges::updateObject(reference)
        object = GenomicRanges::updateObject(object)
    }
    
    comp <- testthat::compare(object, reference, ...)
    
    if (!comp$equal && !is.null(on.fail))
        on.fail(reference, object, paste0(file, ":", deparse(substitute(object))))
    
    testthat::expect(
      comp$equal,
      sprintf("%s not equal to %s.\n%s", lab_act, lab_exp, comp$message),
      info = info
    )
  }
  
  invisible(object)
}

# stock on.update function that will open the new data with View (suitable for rectangular data)
on.update.view = function(new_ref, t) { utils::View(new_ref, t) }

# stock on.update function that will open the new data in a browser window
on.update.edit = function(new_ref, t)
{
  old = options(max.print = 1000000)
  on.exit(options(old), add = TRUE)

  diff = as.character(diffobj::diffPrint("Reference file missing", new_ref, interactive=FALSE, mode="unified", format="html", style=list(html.output="diff.w.style"), disp.width=2000))
  if (!"package:knitr" %in% search())
  {
    diffFile = tempfile(gsub("[^-\\w^&'@{},$=!#().%+~ ]", "_", t, perl=TRUE), fileext=".html")
    cat(diff, file=diffFile)
    utils::browseURL(diffFile)
  }
  else
    cat(diff)
}

# stock on.fail function that will open a diff between the reference and new data in a browser window
on.fail.diff = function(reference, new, t)
{
  old = options(max.print = 1000000)
  on.exit(options(old), add = TRUE)

  diff = as.character(diffobj::diffPrint(reference, new, interactive=FALSE, mode="unified", format="html", style=list(html.output="diff.w.style"), disp.width=2000))
  if (!"package:knitr" %in% search())
  {
    diffFile = tempfile(gsub("[^-\\w^&'@{},$=!#().%+~ ]", "_", t, perl=TRUE), fileext=".html")
    cat(diff, file=diffFile)
    utils::browseURL(diffFile)
  }
  else
    cat(diff)
}

# @return full path to this script
#' current script file (in full path)
#' @return The path of the currently executing script or code snippet.
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
      #} else if (rstudioapi::isAvailable() && nchar(rstudioapi::getActiveDocumentContext()$path) > 0) {
      #  # RStudio Run Selection
      #  # http://stackoverflow.com/a/35842176/2292993  
      #  return(normalizePath(rstudioapi::getActiveDocumentContext()$path))
      } else {
        return("")
      }
    }
  }
}

# Downloaded from: https://gist.github.com/xhdong-umd/6429e7f96735142fa467f3b1daa91a2c
# To decompress zip, gz, bzip2, xz into temp file, run function then remove temp file.
.temp_unzip <- function(filename, fun, ...){
  BFR.SIZE <- 1e7
  if (!file.exists(filename)) {
    stop("No such file: ", filename);
  }
  if (!is.function(fun)) {
    stop(sprintf("Argument 'fun' is not a function: %s", mode(fun)));
  }
  temp_dir <- tempdir()
  # test if it's zip
  files_in_zip <- try(utils::unzip(filename, list = TRUE)$Name, silent = TRUE)
  if (class(files_in_zip) == "character") {
    # hidden files can be ignored: starting with ., ending with $, __MACOSX folder 
    visible_files <- files_in_zip[!grepl("((^__MACOSX\\/.*)|(^\\..*)|(^.*\\$$))", 
                                         files_in_zip)]
    # will not continue for multiple non-hidden files since behavior is not well defined.
    if(length(visible_files)>1) { 
      stop(paste0("Zip file contains multiple visible files:\n", 
                  paste0("    ", visible_files, collapse = "\n")))
    }
    if(length(visible_files) == 0) { stop("\n  No visible file found in Zip file")}
    # proceed with single non-hidden file
    utils::unzip(filename, files = visible_files[1], exdir = temp_dir, overwrite = TRUE)
    dest_file <- file.path(temp_dir, visible_files[1])
  } else {
    dest_file <- tempfile()
    # Setup input and output connections
    inn <- gzfile(filename, open = "rb")
    out <- file(description = dest_file, open = "wb")
    # Process
    nbytes <- 0
    repeat {
      bfr <- readBin(inn, what=raw(0L), size=1L, n=BFR.SIZE)
      n <- length(bfr)
      if (n == 0L) break;
      nbytes <- nbytes + n
      writeBin(bfr, con=out, size=1L)
      bfr <- NULL  # Not needed anymore
    }
    close(inn)
    close(out)
  }
  on.exit(file.remove(dest_file))
  # call fun with temp file
  res <- fun(dest_file, ...)
  return(res)
}
