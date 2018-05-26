library(testthat)
library(customProDB)
library(GetoptLong)

context("PrepareAnnotationGCF")
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

genome = "Tur_tru_v1"
test_that(qq("Annotating subset of Tur_tru_v1 GCF"), {
    annotation_path = qq("@{annotation_path}/@{genome}")
    PrepareAnnotationGCF(qq("@{extdata_path}/gff"), annotation_path)
    env = new.env()
    load_annotations(annotation_path, env, dbsnp=FALSE, cosmic=FALSE)
    expect_equal_to_reference(env$exon, qq("exon_@{genome}.rds"))
    expect_equal_to_reference(env$ids, qq("ids_@{genome}.rds"))
    expect_equal_to_reference(env$proteinseq, qq("proseq_@{genome}.rds"))
    expect_equal_to_reference(env$procodingseq, qq("procodingseq_@{genome}.rds"))
})
