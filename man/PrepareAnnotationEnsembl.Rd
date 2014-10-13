\name{PrepareAnnotationEnsembl}
\alias{PrepareAnnotationEnsembl}
\title{prepare annotation from ENSEMBL}
\usage{
  PrepareAnnotationEnsembl(mart, annotation_path,
    splice_matrix = FALSE, dbsnp = NULL,
    transcript_ids = NULL, COSMIC = FALSE, ...)
}
\arguments{
  \item{mart}{which version of ENSEMBL dataset to use. see
  useMart from package biomaRt for more detail.}

  \item{annotation_path}{specify a folder to store all the
  annotations}

  \item{dbsnp}{specify a snp dataset you want to use for
  the SNP annotation, default is NULL.}

  \item{transcript_ids}{optionally, only retrieve
  transcript annotation data for the specified set of
  transcript ids}

  \item{splice_matrix}{whether generate a known exon splice
  matrix from the annotation. this is not necessary if you
  don't want to analyse junction results, default is
  FALSE.}

  \item{COSMIC}{whether to download COSMIC data, default is
  FALSE.}

  \item{...}{additional arguments}
}
\value{
  several .RData file containing annotations needed for
  following analysis.
}
\description{
  prepare the annotation from ENSEMBL through biomaRt.
}
\details{
  this function automaticlly prepares all annotation
  infromation needed in the following analysis.
}
\examples{
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
host="feb2012.archive.ensembl.org", path="/biomart/martservice",
archive=FALSE)

annotation_path <- tempdir()
transcript_ids <- c("ENST00000234420", "ENST00000269305", "ENST00000445888",
    "ENST00000257430", "ENST00000457016", "ENST00000288602",
    "ENST00000269571", "ENST00000256078", "ENST00000384871")

PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path,
    splice_matrix=FALSE, dbsnp=NULL, transcript_ids=transcript_ids,
    COSMIC=FALSE)
}
\author{
  Xiaojing Wang
}

