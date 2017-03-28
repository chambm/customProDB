library(testthat)
library(customProDB)
library(GetoptLong)

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


## Ensembl 82 / snp146
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
                           local_cache_path=)
})

test_that(qq("Downloading @{genome} Ensembl @{ensembl_version$number} annotations with COSMIC for a few transcripts"), {
  annotation_path = qq("@{annotation_path}/@{genome}_@{ensembl_version$number}_cosmic")
  
  PrepareAnnotationEnsembl(mart=ensembl_mart, annotation_path=annotation_path, transcript_ids=transcript_ids,
                           splice_matrix=FALSE, dbsnp=NULL, COSMIC=TRUE)
})
