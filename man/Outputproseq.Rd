\name{Outputproseq}
\alias{Outputproseq}
\title{output FASTA format file contains proteins that have expression level above the cutoff}
\usage{
  Outputproseq(rpkm, cutoff = "30\%", proteinseq, outfile,
    ids, ...)
}
\arguments{
  \item{rpkm}{a numeric vector containing RPKM for each
  protein}

  \item{cutoff}{cutoff of RPKM value. Two options are
  available, percentage format or RPKM. By default we use
  "30%" or the RPKM value of 1. "30%" means we keep top 70%
  proteins according to their RPKMs.}

  \item{proteinseq}{a dataframe containing protein ids and
  protein sequences.}

  \item{outfile}{output file name.}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information.}

  \item{...}{additional arguments}
}
\value{
  FASTA file contains proteins with RPKM above the cutoff.
}
\description{
  Get the FASTA file of proteins that pass RPKM cutoff. the
  FASTA ID line contains protein ID, gene ID, HGNC symbol
  and description
}
\details{
  by taking the RPKM value as input, the function outputs
  sequences of the proteins that pass the cutoff.
}
\examples{
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
bamFile <- system.file("extdata/bams", "test1_sort.bam",
    package="customProDB")
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
outf1 <- paste(tempdir(), '/test_rpkm.fasta', sep='')
Outputproseq(RPKM, 1, proteinseq, outf1, ids)
}
\author{
  Xiaojing Wang
}

