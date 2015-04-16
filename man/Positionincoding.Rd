\name{Positionincoding}
\alias{Positionincoding}
\title{Find the position in coding sequence for each variation.}
\usage{
  Positionincoding(Vars, exon, dbsnp = NULL, COSMIC = NULL,
    ...)
}
\arguments{
  \item{Vars}{a GRanges object of variations}

  \item{exon}{a dataframe of exon annotations for protein
  coding transcripts.}

  \item{dbsnp}{provide a GRanges object of known dbsnp
  information to include dbsnp evidence into the output
  table, default is NULL.}

  \item{COSMIC}{provide a GRanges object of known COSMIC
  information to include COSMIC evidence into the output
  table, default is NULL.}

  \item{...}{additional arguments}
}
\value{
  a data frame containing the position in coding sequence
  for each variation
}
\description{
  For those variations labeled with "Coding",
  positionincoding() function computes the position of
  variation in the coding sequence of each transcript.
}
\details{
  this function prepares input data frame for
  aaVariation().
}
\examples{
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
table(values(vcf[[1]])[['INDEL']])
index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
indelvcf <- vcf[[1]][index]

index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
SNVvcf <- vcf[[1]][index]
load(system.file("extdata/refseq", "exon_anno.RData",
    package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData",
    package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData",
    package="customProDB"))
load(system.file("extdata/refseq", "cosmic.RData",
    package="customProDB"))
postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
}
\author{
  Xiaojing Wang
}

