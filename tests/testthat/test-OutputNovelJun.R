library(testthat)
library(AnnotationDbi)
library(customProDB)
suppressPackageStartupMessages(library('BSgenome.Hsapiens.UCSC.hg19'))

context("OutputNovelJun")

load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite", package="customProDB"))
bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
jun <-  Bed2Range(bedfile,skip=1,covfilter=5)

test_that("JunctionType sets the junctions to the correct types", {
  junction_type <- JunctionType(jun, splicemax, txdb, ids)
  junctionCountByType = table(junction_type[, 'jun_type'])
  expect_equivalent(junctionCountByType["connect a known exon and a region overlap with known exon"], 1)
  expect_equivalent(junctionCountByType["connect two non-exon region"], 9)
  expect_equivalent(junctionCountByType["known junction"], 46)
})

test_that("OutputNovelJun creates 3 RData files", {
  junction_type <- JunctionType(jun, splicemax, txdb, ids)
  chrom <- paste('chr',c(1:22,'X','Y','M'), sep='')
  junction_type <- subset(junction_type, seqnames %in% chrom)
  outf_junc <- paste(tempdir(), '/test_junc.fasta', sep='')
  #outf_junc_coding <- paste(tempdir(), '/test_junc_coding.fasta', sep='')
  expect_warning(OutputNovelJun(junction_type, Hsapiens, outf_junc, proteinseq))
  
  load(paste(outf_junc, '_jun_anno.RData', sep=''))
  if(!file.exists("hg19_jun_anno.rds")) { View(jun_anno) }
  expect_equal_to_reference(jun_anno, "hg19_jun_anno.rds")
  
  load(paste(outf_junc, '_coding.RData', sep=''))
  if(!file.exists("hg19_junpepcoding.rds")) { View(junpepcoding) }
  expect_equal_to_reference(junpepcoding, "hg19_junpepcoding.rds")
  
  load(paste(outf_junc, '_junpep.RData', sep=''))
  if(!file.exists("hg19_junpep.rds")) { View(junpep) }
  expect_equal_to_reference(junpep, "hg19_junpep.rds")
})
