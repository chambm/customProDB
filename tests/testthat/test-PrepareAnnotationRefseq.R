library(testthat)
library(customProDB)
library(GetoptLong)

context("PrepareAnnotationRefseq")

transcript_ids = c("NM_001126112", "NM_033360", "NR_073499", "NM_004448",
                   "NM_000179", "NR_029605", "NM_004333", "NM_001127511")
testthat_path = current_script_file()
if (nchar(testthat_path) == 0) {
  testthat_path = getwd()
} else {
  testthat_path = dirname(testthat_path)
}

extdata_path = paste0(testthat_path, "/../../inst/extdata")
pepfasta = paste0(extdata_path, "/refseq_pro_seq.fasta")
CDSfasta = paste0(extdata_path, "/refseq_coding_seq.fasta")
local_cache_path = paste0(extdata_path, "/cache")
annotation_path = tempdir()

test_hg19_annotations = function(annotation_path) {
  load(paste0(annotation_path, "/exon_anno.RData"))
  load(paste0(annotation_path, "/ids.RData"))
  load(paste0(annotation_path, "/procodingseq.RData"))
  load(paste0(annotation_path, "/proseq.RData"))
  expect_is(exon, "data.frame")
  expect_is(ids, "data.frame")
  expect_is(proteinseq, "data.frame")
  expect_is(procodingseq, "data.frame")
  expect_equal(nrow(exon), 98)
  expect_equal(nrow(ids), 8)
  expect_equal(nrow(proteinseq), 6)
  expect_equal(nrow(procodingseq), 6)
}

test_cosmic = function(annotation_path) {
  load(paste0(annotation_path, "/cosmic.RData"))
  expect_is(cosmic, "GRanges")
  expect_equal(length(cosmic), 8667)
}


## hg19/snp146
genome = "hg19"
dbsnp = "snp146"
test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 146 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  
  load(paste0(annotation_path, "/dbsnpinCoding.RData"))
  expect_is(dbsnpinCoding, "GRanges")
  expect_equal(length(dbsnpinCoding), 4670)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 146 and COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  test_cosmic(annotation_path)
  
  load(paste0(annotation_path, "/dbsnpinCoding.RData"))
  expect_is(dbsnpinCoding, "GRanges")
  expect_equal(length(dbsnpinCoding), 4670)
})

test_that(qq("Downloading @{genome} RefSeq annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  test_cosmic(annotation_path)
})

## hg19/snp141
dbsnp = "snp141"
test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 141 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  
  load(paste0(annotation_path, "/dbsnpinCoding.RData"))
  expect_is(dbsnpinCoding, "GRanges")
  expect_equal(length(dbsnpinCoding), 1446)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 141 and COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  test_cosmic(annotation_path)
  
  load(paste0(annotation_path, "/dbsnpinCoding.RData"))
  expect_is(dbsnpinCoding, "GRanges")
  expect_equal(length(dbsnpinCoding), 1446)
})

test_that(qq("Downloading @{genome} RefSeq annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_cosmic")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=TRUE,
                          local_cache_path=local_cache_path)
  test_hg19_annotations(annotation_path)
  test_cosmic(annotation_path)
})


## mm10/snp
genome = "mm10"
dbsnp = "snp142"
pepfasta = paste0(extdata_path, "/mm10/refseq_pro_seq.fasta")
CDSfasta = paste0(extdata_path, "/mm10/refseq_coding_seq.fasta")
transcript_ids = c("NM_007462", "NR_029825", "NM_021284", "NM_001003817",
                   "NM_139294", "NM_011317", "NM_010830", "NM_001127233")

test_mm10_annotations = function(annotation_path) {
  load(paste0(annotation_path, "/exon_anno.RData"))
  load(paste0(annotation_path, "/ids.RData"))
  load(paste0(annotation_path, "/procodingseq.RData"))
  load(paste0(annotation_path, "/proseq.RData"))
  expect_is(exon, "data.frame")
  expect_is(ids, "data.frame")
  expect_is(proteinseq, "data.frame")
  expect_is(procodingseq, "data.frame")
  expect_equal(nrow(exon), 101)
  expect_equal(nrow(ids), 8)
  expect_equal(nrow(proteinseq), 7)
  expect_equal(nrow(procodingseq), 7)
}

test_that(qq("Downloading basic @{genome} RefSeq annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=NULL, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_mm10_annotations(annotation_path)
})

test_that(qq("Downloading @{genome} RefSeq annotations with dbSNP 142 for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{dbsnp}")
  PrepareAnnotationRefseq(genome=genome, CDSfasta, pepfasta, annotation_path,
                          dbsnp=dbsnp, transcript_ids=transcript_ids,
                          splice_matrix=FALSE, COSMIC=FALSE,
                          local_cache_path=local_cache_path)
  test_mm10_annotations(annotation_path)
  
  load(paste0(annotation_path, "/dbsnpinCoding.RData"))
  expect_is(dbsnpinCoding, "GRanges")
  expect_equal(length(dbsnpinCoding), 304)
})
