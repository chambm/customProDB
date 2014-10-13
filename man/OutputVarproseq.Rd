\name{OutputVarproseq}
\alias{OutputVarproseq}
\title{Output the variant(SNVs) protein sequences into FASTA format}
\usage{
  OutputVarproseq(vartable, proteinseq, outfile, ids,
    lablersid = FALSE, RPKM = NULL, ...)
}
\arguments{
  \item{vartable}{A data frame which is the output of
  aaVariation().}

  \item{proteinseq}{A dataframe containing protein ids and
  the protein sequence.}

  \item{outfile}{Output file name.}

  \item{ids}{A dataframe containing gene/transcript/protein
  id mapping information.}

  \item{lablersid}{If includes the dbSNP rsid in the header
  of each sequence, default is FALSE.  Must provide dbSNP
  information in function Positionincoding() if put TRUE
  here.}

  \item{RPKM}{If includes the RPKM value in the header of
  each sequence, default is NULL.}

  \item{...}{Additional arguments}
}
\value{
  FASTA file containing proteins with single nucleotide
  variation.
}
\description{
  Output the non-synonymous SNVs into FASTA file.
}
\details{
  This function uses the output of aaVariation() as input,
  introduces the nonsynonymous variation into the protein
  database.
}
\examples{
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
table(values(vcf[[1]])[['INDEL']])
index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
SNVvcf <- vcf[[1]][index]
load(system.file("extdata/refseq", "exon_anno.RData",
package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData",
    package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData",
    package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding)
txlist <- unique(postable_snv[, 'txid'])
codingseq <- procodingseq[procodingseq[, 'tx_id'] \%in\% txlist, ]
mtab <- aaVariation (postable_snv, codingseq)
outfile <- paste(tempdir(), '/test_snv.fasta',sep='')
OutputVarproseq(mtab, proteinseq, outfile, ids, lablersid=TRUE)
}
\author{
  Xiaojing Wang
}

