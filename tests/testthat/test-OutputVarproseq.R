library(testthat)
library(customProDB)
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library('BSgenome.Hsapiens.UCSC.hg19'))

context("OutputVarproseq")

load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")

test_that("INDELs are separated from SNVs", {
  vcf <- InputVcf(vcffile)
  indelCounts = table(values(vcf[[1]])[['INDEL']])
  expect_equivalent(indelCounts["FALSE"], 54)
  expect_equivalent(indelCounts["TRUE"], 7)
})

test_that("OutputVarproseq creates a FASTA file and returns a data.frame", {
  vcf <- InputVcf(vcffile)
  index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
  SNVvcf <- vcf[[1]][index]
  postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding)
  txlist <- unlist(unique(postable_snv[, 'txid']))
  codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist, ]
  stopifnot(nrow(postable_snv) > 0)
  stopifnot(nrow(codingseq) > 0)
  
  mtab <- aaVariation(postable_snv, codingseq, show_progress=F)
  stopifnot(nrow(mtab) > 0)
  
  outfile <- paste(tempdir(), '/test_snv.fasta', sep='')
  snvproseq <- OutputVarproseq(mtab, proteinseq, outfile, ids, lablersid=TRUE, RPKM=NULL)
  
  fasta = readLines(outfile)
  expect_equal_to_reference(fasta, 'test_snv.fasta.rds')
  expect_equal_to_reference(snvproseq, 'test_snv.fasta_snvproseq.rds')
})
