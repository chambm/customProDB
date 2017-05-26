#Sys.setenv("R_TESTS" = "")
library(testthat)
library(customProDB)
library(GetoptLong)
library(biomaRt)

context("PrepareAnnotationEnsembl")
options(testthat.on.update = on.update.view)
options(testthat.on.fail = on.fail.diff)

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

local_cache_path = paste0(extdata_path, "/cache")
annotation_path = tempdir()
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


## Ensembl human 82
genome = "hsapiens_gene_ensembl"
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
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
})

## Ensembl human 82 + dbSNP
dbsnp = "snp146"
test_that(qq("Downloading @{genome} Ensembl @{ensembl_version$number} annotations with dbSNP for a few transcripts"), {
    annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}_@{dbsnp}")
    
    PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                             splice_matrix=FALSE, dbsnp=dbsnp, COSMIC=FALSE,
                             local_cache_path=local_cache_path)
    
    env = new.env()
    load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=FALSE)
    expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{ensembl_version$number}_@{dbsnp}.rds"))
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
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}_@{ensembl_version$number}.rds"))
})


## Ensembl human 87 + dbSNP + COSMIC
ensembl_version = list(number=87, host="dec2016.archive.ensembl.org")
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=genome, host=ensembl_version$host)
test_that(qq("Downloading @{genome} Ensembl @{ensembl_version$number} annotations with dbSNP and COSMIC for a few transcripts"), {
    annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}_@{dbsnp}_cosmic")
    
    PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                             splice_matrix=FALSE, dbsnp=dbsnp, COSMIC=TRUE,
                             local_cache_path=local_cache_path)
    
    env = new.env()
    load_annotations(annotation_path, env, dbsnp=TRUE, cosmic=TRUE)
    expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$dbsnpinCoding, qq("dbsnpinCoding_@{genome}_@{ensembl_version$number}_@{dbsnp}.rds"))
    expect_equal_to_reference(env$cosmic, qq("cosmic_@{genome}_@{ensembl_version$number}.rds"))
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
  expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
  expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
})


## Ensembl 87 dog
genome = "cfamiliaris_gene_ensembl"
transcript_ids = c("ENSCAFT00000028439", "ENSCAFT00000011866", "ENSCAFT00000047757",
                   "ENSCAFT00000036929", "ENSCAFT00000040033")
ensembl_version = list(number=87, host="dec2016.archive.ensembl.org")
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=genome, host=ensembl_version$host)

test_that(qq("Downloading basic @{genome} Ensembl @{ensembl_version$number} annotations for a few transcripts"), {
    annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}")
    
    PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                             splice_matrix=FALSE, dbsnp=NULL, COSMIC=FALSE,
                             local_cache_path=local_cache_path)
    
    env = new.env()
    load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
    expect_equal_to_reference(env$exon, qq("exon_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$ids, qq("ids_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}_@{ensembl_version$number}.rds"))
    expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}_@{ensembl_version$number}.rds"))
})

