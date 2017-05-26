library(testthat)
library(AnnotationDbi)
library(customProDB)

context("easyRun")

bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")


test_that("easyRun works for just SNVs", {
  outfile_path <- tempdir()

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=FALSE, nov_junction=FALSE,
          lablersid=FALSE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
})


test_that("easyRun works for SNVs with dbSNP and COSMIC annotation", {
  outfile_path <- tempdir()
  
  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=FALSE, nov_junction=FALSE,
          lablersid=TRUE, COSMIC=TRUE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv_dbsnp.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
})


test_that("easyRun works for SNVs and INDELs", {
  outfile_path <- tempdir()

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=TRUE, nov_junction=FALSE,
          lablersid=TRUE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv_dbsnp.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
  
  fasta = read.delim(file.path(outfile_path, "test_indel.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_indel.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
})


test_that("easyRun works for SNVs, INDELs, and novel junctions", {
  outfile_path <- tempdir()
  bedFile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=TRUE, nov_junction=TRUE,
          bedFile=bedFile, genome=BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
          lablersid=FALSE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
  
  fasta = read.delim(file.path(outfile_path, "test_indel.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_indel.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
  
  fasta = read.delim(file.path(outfile_path, "test_junc.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_junc.fasta.rds", on.update=on.update.edit, on.fail=on.fail.diff)
})

