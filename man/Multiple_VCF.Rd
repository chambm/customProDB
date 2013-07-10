\name{Multiple_VCF}
\alias{Multiple_VCF}
\title{Generate shared variation dataset from multiple VCF files}
\usage{
  Multiple_VCF(vcfs, share_num, ...)
}
\arguments{
  \item{vcfs}{a list of GRanges object which input from
  multiple VCF files using function InputVcf.}

  \item{share_num}{Two options, percentage format or sample
  number.}

  \item{...}{additional arguments}
}
\value{
  a GRange object that contains the shared variations
}
\description{
  Load multiple vcf files and output a GRange object with
  SNVs present in multiple samples.
}
\details{
  This function allows to limit SNVs that are present in at
  least m out of n VCF files.
}
\examples{
path <- system.file("extdata/vcfs", package="customProDB")
vcfFiles<- paste(path, '/', list.files(path, pattern="*vcf$"), sep='')
vcfs <- lapply(vcfFiles, function(x) InputVcf(x))
shared <- Multiple_VCF(vcfs, share_num=2)
}
\author{
  Xiaojing Wang
}

