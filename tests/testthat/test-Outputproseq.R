library(testthat)
library(customProDB)

context("OutputProseq")

load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))

test_that("Outputproseq with cutoff=0 creates a FASTA file with all proteins", {
  bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
  RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
  
  outfile <- paste(tempdir(), '/test_rpkm0.fasta', sep='')
  Outputproseq(RPKM, cutoff=0, proteinseq, outfile, ids)
  
  fasta = readLines(outfile)
  if (!file.exists('test_rpkm0.fasta.rds')) { utils::edit(fasta) }
  expect_equal_to_reference(fasta, 'test_rpkm0.fasta.rds')
})

test_that("Outputproseq with cutoff=1000 creates a FASTA file with 2 proteins", {
  bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
  RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
  
  outfile <- paste(tempdir(), '/test_rpkm1000.fasta', sep='')
  Outputproseq(RPKM, cutoff=1000, proteinseq, outfile, ids)
  
  fasta = readLines(outfile)
  if (!file.exists('test_rpkm1000.fasta.rds')) { utils::edit(fasta) }
  expect_equal_to_reference(fasta, 'test_rpkm1000.fasta.rds')
})