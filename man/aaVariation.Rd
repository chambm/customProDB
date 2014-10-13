\name{aaVariation}
\alias{aaVariation}
\title{get the functional consequencece of SNVs located in coding region}
\usage{
  aaVariation(position_tab, coding, ...)
}
\arguments{
  \item{position_tab}{a data frame from Positionincoding()}

  \item{coding}{a data frame cotaining coding sequence for
  each protein.}

  \item{...}{Additional arguments}
}
\value{
  a data frame containing consequence for each variations.
}
\description{
  Variations can be divided into SNVs and INDELs. By taking
  the output of positionincoding() as input, aaVariation()
  function predicts the consequences of SNVs in the
  harbored transcript, such as synonymous or
  non-synonymous.
}
\details{
  this function predicts the consequence for SNVs. for
  INDELs, use Outputabberrant().
}
\examples{
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
table(values(vcf[[1]])[['INDEL']])

index <- which(values(vcf[[1]])[['INDEL']]==FALSE)
SNVvcf <- vcf[[1]][index]
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
postable_snv <- Positionincoding(SNVvcf,exon,dbsnpinCoding)
txlist <- unique(postable_snv[,'txid'])
codingseq <- procodingseq[procodingseq[,'tx_id'] \%in\% txlist,]
mtab <- aaVariation (postable_snv,codingseq)
mtab[1:3,]
}
\author{
  Xiaojing Wang
}

