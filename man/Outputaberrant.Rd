\name{Outputaberrant}
\alias{Outputaberrant}
\title{generate FASTA file containing short INDEL}
\usage{
  Outputaberrant(positiontab, outfile, coding, proteinseq,
    ids, RPKM = NULL, ...)
}
\arguments{
  \item{positiontab}{a data frame which is the output of
  function Positionincoding() for INDELs.}

  \item{outfile}{output file name}

  \item{coding}{a data frame cotaining coding sequence for
  each protein.}

  \item{proteinseq}{a data frame cotaining amino acid
  sequence for each protein.}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information.}

  \item{RPKM}{if includes the RPKM value in the header of
  each sequence, default is NULL.}

  \item{...}{Additional arguments.}
}
\value{
  FASTA file containing aberrant proteins.
}
\description{
  Short insertion/deletion may lead to aberrant proteins in
  cells. We provide a function to generate FASTA file
  containing this kind of proteins.
}
\details{
  the function applys the INDEL into the coding sequence,
  then translates them into protein sequence, terminated by
  stop codon. Remove the sequences the same as normal ones
  or as part of normal ones.
}
\examples{
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
table(values(vcf[[1]])[['INDEL']])
index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
indelvcf <- vcf[[1]][index]

load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData",
        package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData",
        package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
postable_indel <- Positionincoding(indelvcf, exon)
txlist_indel <- unique(postable_indel[, 'txid'])
codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] \%in\% txlist_indel, ]
outfile <-  paste(tempdir(), '/test_indel.fasta', sep='')
Outputaberrant(postable_indel, coding=codingseq_indel,
proteinseq=proteinseq, outfile=outfile, ids=ids)

}
\author{
  Xiaojing Wang
}

