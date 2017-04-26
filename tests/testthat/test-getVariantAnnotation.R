library(testthat)
library(customProDB)

context("getVariantAnnotation")

test_that("aaVariation works", {
  
  load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
  load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
  load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
  load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
  
  variantAnnotation = getVariantAnnotation(system.file("extdata/vcfs", "test1.vcf", package="customProDB"),
                                           ids, exon,
                                           proteinseq, procodingseq,
                                           dbsnpinCoding, cosmic)
  
  expect_equal_to_reference(variantAnnotation$variantTable,
                            "test1-hg19-aaVariation.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
})


test_that("proteinseq and procodingseq return correctly", {
  
  load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
  load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
  load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
  load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
  
  variantAnnotation = getVariantAnnotation(system.file("extdata/vcfs", "test1.vcf", package="customProDB"),
                                           ids, exon,
                                           proteinseq, procodingseq,
                                           dbsnpinCoding, cosmic)
  
  expect_equal_to_reference(variantAnnotation$snvprocoding,
                            "test1-hg19-snvprocoding.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$snvproseq,
                            "test1-hg19-snvproseq.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$indelprocoding,
                            "test1-hg19-indelprocoding.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$indelproseq,
                            "test1-hg19-indelproseq.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
})


test_that("variantType returns correct variant types", {
  expect_equal(variantType("A", "C"), "snp")
  expect_equal(variantType("A", "C,G"), "snp")
  expect_equal(variantType("AA", "CT"), "mnp")
  expect_equal(variantType("AT", "A"), "del")
  expect_equal(variantType("AT", "T"), "del")
  expect_equal(variantType("A", "AT"), "ins")
  expect_equal(variantType("A", "TA"), "ins")
  expect_equal(variantType("AA", "C"), "complex")
  expect_equal(variantType("AA", "CAT"), "complex")
  expect_equal(variantType("AA", "TAAC"), "complex")
  expect_equal(variantType("ATTA", "ATA"), "del")
  expect_equal(variantType("AAA", "AATA"), "ins")
})
