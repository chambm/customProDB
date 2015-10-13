\name{calculateRPKM}
\alias{calculateRPKM}
\title{Caculate RPKM for each transcripts based on exon read counts.}
\usage{
  calculateRPKM(bamFile, exon, proteincodingonly = TRUE,
    ids = NULL, ...)
}
\arguments{
  \item{bamFile}{a the input BAM file name.}

  \item{exon}{a dataframe of exon annotations.}

  \item{proteincodingonly}{if TRUE only output RPKMs for
  protein coding transcripts, the name of output vector
  will be protein id. if FALSE, output the RPKM for all
  transcripts.}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information.}

  \item{...}{additional arguments}
}
\value{
  RPKM value for all transcripts or protein coding
  transcripts.
}
\description{
  Normalized expression level based on exon read counts.
  The default output is a vector containing RPKMs for each
  transcript. vector name is the transcript name. calculate
  the RPKMs by chromosome. If proteincodingonly=TRUE, vetor
  name will be set to protein name, and only output RPKMs
  for the protein coding transcripts.
}
\details{
  caculate RPKM from a BAM file based on exon read counts
}
\examples{
##test1.bam file is part of the whole bam file.
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
RPKM <- calculateRPKM(bamFile,exon,proteincodingonly=TRUE,ids)
}
\author{
  Xiaojing Wang
}

