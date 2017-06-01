library(testthat)
library(AnnotationDbi)
library(customProDB)

context("easyRun")
options(testthat.on.update = on.update.edit)
options(testthat.on.fail = on.fail.diff)


## test hg19
bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")


test_that("easyRun works for just SNVs", {
  outfile_path <- tempdir()

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=FALSE, nov_junction=FALSE,
          lablersid=FALSE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv.fasta.rds")
})


test_that("easyRun works for SNVs with dbSNP and COSMIC annotation", {
  outfile_path <- tempdir()
  
  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=FALSE, nov_junction=FALSE,
          lablersid=TRUE, COSMIC=TRUE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv_dbsnp.fasta.rds")
})


test_that("easyRun works for SNVs and INDELs", {
  outfile_path <- tempdir()

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=TRUE, nov_junction=FALSE,
          lablersid=TRUE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv_dbsnp.fasta.rds")
  
  fasta = read.delim(file.path(outfile_path, "test_indel.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_indel.fasta.rds")
})


test_that("easyRun works for SNVs, INDELs, and novel junctions", {
  outfile_path <- tempdir()
  bedFile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")

  easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "test",
          rpkm_cutoff=1, INDEL=TRUE, nov_junction=TRUE,
          bedFile=bedFile, genome=BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
          lablersid=FALSE, COSMIC=FALSE)
  
  fasta = read.delim(file.path(outfile_path, "test_snv.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_snv.fasta.rds")
  
  fasta = read.delim(file.path(outfile_path, "test_indel.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_indel.fasta.rds")
  
  fasta = read.delim(file.path(outfile_path, "test_junc.fasta"), header=FALSE)
  customProDB:::expect_equal_to_reference(fasta, "test-easyRun-test_junc.fasta.rds")
})


## test mouse Ensembl
bamFile <- system.file("extdata/mmusculus_gene_ensembl_87", "test.bam", package="customProDB")
annotation_path <- system.file("extdata/mmusculus_gene_ensembl_87", package="customProDB")

test_that("easyRun works with Ensembl annotation and NULL vcfFile", {
    outfile_path <- tempdir()
    
    easyRun(bamFile, RPKM=NULL, NULL, annotation_path, tempdir(), "mmusculus_gene_ensembl_87",
            rpkm_cutoff=0, INDEL=TRUE, nov_junction=FALSE,
            lablersid=TRUE, COSMIC=FALSE)
    
    fasta = read.delim(file.path(outfile_path, "mmusculus_gene_ensembl_87_rpkm.fasta"), header=FALSE)
    customProDB:::expect_equal_to_reference(fasta, "test-easyRun-mmusculus_gene_ensembl_87_rpkm.fasta.rds")
})

vcffile <- system.file("extdata/mmusculus_gene_ensembl_87", "test.vcf", package="customProDB")
test_that("easyRun works for SNVs with Ensembl annotation", {
    outfile_path <- tempdir()
    
    easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, tempdir(), "mmusculus_gene_ensembl_87",
            rpkm_cutoff=0, INDEL=FALSE, nov_junction=FALSE,
            lablersid=TRUE, COSMIC=FALSE)
    
    fasta = read.delim(file.path(outfile_path, "mmusculus_gene_ensembl_87_rpkm.fasta"), header=FALSE)
    customProDB:::expect_equal_to_reference(fasta, "test-easyRun-mmusculus_gene_ensembl_87_rpkm.fasta.rds")
    
    fasta = read.delim(file.path(outfile_path, "mmusculus_gene_ensembl_87_snv.fasta"), header=FALSE)
    customProDB:::expect_equal_to_reference(fasta, "test-easyRun-mmusculus_gene_ensembl_87_snv_dbsnp.fasta.rds")
})
