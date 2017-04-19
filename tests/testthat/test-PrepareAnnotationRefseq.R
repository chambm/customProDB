#Sys.setenv("R_TESTS" = "")
library(testthat)
library(customProDB)
library(GetoptLong)

context("PrepareAnnotationRefseq")

testthat_path = current_script_file()
if (!nzchar(testthat_path)) {
  testthat_path = getwd()
  if (!dir.exists(paste0(testthat_path, "/../../inst/extdata"))) {
    extdata_path = system.file("extdata", package="customProDB", mustWork=TRUE)
  } else {
    extdata_path = paste0(testthat_path, "/../../inst/extdata")
  }
} else {
  testthat_path = dirname(testthat_path)
  extdata_path = paste0(testthat_path, "/../../inst/extdata")
}
stopifnot(dir.exists(extdata_path))

pepfasta = paste0(extdata_path, "/refseq_pro_seq.fasta")
CDSfasta = paste0(extdata_path, "/refseq_coding_seq.fasta")
local_cache_path = paste0(extdata_path, "/cache")
annotation_path = tempdir()
stopifnot(file.exists(pepfasta))
stopifnot(file.exists(CDSfasta))
stopifnot(dir.exists(local_cache_path))
stopifnot(dir.exists(annotation_path))


load_annotations = function(annotation_path, envir, dbsnp=FALSE, cosmic=FALSE) {
  
  stopifnot(file.exists(paste0(annotation_path, "/exon_anno.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/ids.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/procodingseq.RData")))
  stopifnot(file.exists(paste0(annotation_path, "/proseq.RData")))
  
  load(paste0(annotation_path, "/exon_anno.RData"), envir=envir)
  load(paste0(annotation_path, "/ids.RData"), envir=envir)
  load(paste0(annotation_path, "/procodingseq.RData"), envir=envir)
  load(paste0(annotation_path, "/proseq.RData"), envir=envir)
  
  expect_is(envir$exon, "data.frame")
  expect_is(envir$ids, "data.frame")
  expect_is(envir$proteinseq, "data.frame")
  expect_is(envir$procodingseq, "data.frame")
  
  if (dbsnp) {
    load(paste0(annotation_path, "/dbsnpinCoding.RData"), envir=envir)
    expect_is(envir$dbsnpinCoding, "GRanges")
  }
  
  if (cosmic) {
    load(paste0(annotation_path, "/cosmic.RData"), envir=envir)
    expect_is(envir$cosmic, "GRanges")
  }
}


## hg19/snp146
genome = "hg19"
dbsnp = "snp146"
transcript_ids = c("NM_001126112", "NM_033360", "NR_073499", "NM_004448",
                   "NM_000179", "NR_029605", "NM_004333", "NM_001127511")

test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 146 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)

  env = new.env()
  load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{dbsnp}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 146 and COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=TRUE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{dbsnp}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=TRUE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}.rds"), on.update=on.update.view)
})

## hg19/snp141
dbsnp = "snp141"
test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 141 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{dbsnp}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 141 and COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=TRUE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{dbsnp}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=TRUE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}.rds"), on.update=on.update.view)
})


## mm10/snp
genome = "mm10"
dbsnp = "snp142"
pepfasta = paste0(extdata_path, "/mm10/refseq_pro_seq.fasta")
CDSfasta = paste0(extdata_path, "/mm10/refseq_coding_seq.fasta")
transcript_ids = c("NM_007462", "NR_029825", "NM_021284", "NM_001003817",
                   "NM_139294", "NM_011317", "NM_010830", "NM_001127233")

test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 142 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{dbsnp}.rds"), on.update=on.update.view)
})
