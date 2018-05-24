#' @title CustomProDB
#' @description Database search is the most widely used approach for peptide and protein identification in
#' mass-spectrometry-based proteomics studies. Sample-specific protein databases derived from RNA-Seq data can 
#' better approximate the real protein pools in the samples and thus improve protein identification. 
#' More importantly, single nucleotide variations, short insertion and deletions and novel junctions 
#' identified from RNA-Seq data make protein databases more complete and sample-specific. CustomProDB 
#' enables the easy generation of customized databases from RNA-Seq data for proteomics search. It bridges 
#' genomics and proteomics studies and facilitates cross-omics data integration.
#' @docType package
#' @name customProDB
#' @importFrom methods is
#' @importFrom stats quantile
#' @importFrom utils download.file read.table setTxtProgressBar txtProgressBar write.table
#' @importFrom AnnotationDbi saveDb loadDb
#' @importFrom data.table data.table rbindlist setkey setDT 
#' @importFrom GetoptLong qq
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom RCurl getURL
#' @importFrom rtracklayer import.gff
#' @import S4Vectors GenomicRanges GenomicFeatures Biostrings
NULL
