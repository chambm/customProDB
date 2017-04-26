library(testthat)
library(customProDB)
suppressPackageStartupMessages(library(VariantAnnotation))

context("OutputVarprocodingseq")

load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")

test_that("OutputVarprocodingseq returns correct coding sequences", {
  vcf <- InputVcf(vcffile)
  index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
  SNVvcf <- vcf[[1]][index]
  postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding)
  firstTxId = postable_snv[, .(txid=min(txid)), txname]
  codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% firstTxId$txid, ]
  stopifnot(nrow(postable_snv) > 0)
  stopifnot(nrow(codingseq) > 0)
  
  mtab <- aaVariation(postable_snv, codingseq, show_progress=F)
  stopifnot(nrow(mtab) > 0)
  
  snvprocodingseq = OutputVarprocodingseq(mtab, codingseq, ids, lablersid=TRUE)
  expect_equal_to_reference(snvprocodingseq, "test_snv.fasta_snvprocodingseq.rds", on.update=on.update.view, on.fail=on.fail.diff)
})

