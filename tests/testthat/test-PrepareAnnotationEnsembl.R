library(testthat)
library(customProDB)
library(GetoptLong)
library(biomaRt)

context("PrepareAnnotationEnsembl")

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
    expect_is(envir$dbsnpInCoding, "GRanges")
  }

  if (cosmic) {
    load(paste0(annotation_path, "/cosmic.RData"), envir=envir)
    expect_is(envir$cosmic, "GRanges")
  }
}


on.update.view = function(o, t) { View(o, t) }


## Ensembl human 82
genome = "hsapiens_gene_ensembl"
dbsnp = "snp146"
transcript_ids = c("ENST00000234420", "ENST00000269305", "ENST00000445888", 
                   "ENST00000257430", "ENST00000508376", "ENST00000288602",
                   "ENST00000269571", "ENST00000256078", "ENST00000384871")
ensembl_version = list(number=82, host="sep2015.archive.ensembl.org")
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=genome, host=ensembl_version$host)

test_that(qq("Downloading basic @{genome} Ensembl @{ensembl_version$number} annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}")
  
  PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                           splice_matrix=FALSE, dbsnp=NULL, COSMIC=FALSE,
                           local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
})

## Ensembl human 87 + COSMIC
ensembl_version = list(number=87, host="dec2016.archive.ensembl.org")
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=genome, host=ensembl_version$host)
test_that(qq("Downloading @{genome} Ensembl @{ensembl_version$number} annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}_cosmic")
  
  PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                           splice_matrix=FALSE, dbsnp=NULL, COSMIC=TRUE,
                           local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=TRUE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
})


## Ensembl 87 mouse
genome = "mmusculus_gene_ensembl"
transcript_ids = c("ENSMUST00000083500", "ENSMUST00000156486", "ENSMUST00000005503",
                   "ENSMUST00000066257", "ENSMUST00000108657", "ENSMUST00000066133",
                   "ENSMUST00000058295", "ENSMUST00000002487")
ensembl_version = list(number=87, host="dec2016.archive.ensembl.org")
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=genome, host=ensembl_version$host)

test_that(qq("Downloading basic @{genome} Ensembl @{ensembl_version$number} annotations for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}")
  
  PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                           splice_matrix=FALSE, dbsnp=NULL, COSMIC=FALSE,
                           local_cache_path=local_cache_path)
  
  env = new.env()
  load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"), on.update=on.update.view)
})