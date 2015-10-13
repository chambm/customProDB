\name{InputVcf}
\alias{InputVcf}
\title{Generate a list of GRanges objects from a VCF file.}
\usage{
  InputVcf(vcfFile, ...)
}
\arguments{
  \item{vcfFile}{a character contains the path and name of
  a VCF file}

  \item{...}{additional arguments}
}
\value{
  a list of GRanges object containing a representation of
  data from the VCF file
}
\description{
  The InputVcf() function generates a list of GRanges
  object from a single VCF file.
}
\details{
  Read all fields in a VCF file into GRanges object.
}
\examples{
## multiple samples in one VCF file

vcffile <- system.file("extdata", "test_mul.vcf", package="customProDB")
vcfs <- InputVcf(vcffile)
length(vcfs)

## single sample

vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)
length(vcf)
}
\author{
  Xiaojing Wang
}

