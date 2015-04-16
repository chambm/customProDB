\name{SharedJunc}
\alias{SharedJunc}
\title{Generate shared junctions dataset from multiple BED files}
\usage{
  SharedJunc(juns, share_num = 2, ext_up = 100,
    ext_down = 100, ...)
}
\arguments{
  \item{juns}{a list of GRanges object which input from
  multiple VCF files using function InputVcf.}

  \item{share_num}{Junctions must occurs in this number of
  samples to be consider. Two options, percentage format or
  sample number.}

  \item{ext_up}{upstream extension from the junction point,
  Default is 100nt.}

  \item{ext_down}{downstream extension from the junction
  point, Default is 100nt.}

  \item{...}{additional arguments}
}
\value{
  a GRange object that contains the shared junctions
}
\description{
  Load multiple BED files and output a GRange object with
  junctions present in multiple samples.
}
\details{
  This function allows to limit junctions that are present
  in at least m out of n BED files.
}
\examples{
path <- system.file("extdata/beds", package="customProDB")
bedFiles<- paste(path, '/', list.files(path, pattern="*bed$"), sep='')
juncs <- lapply(bedFiles, function(x) Bed2Range(x, skip=1, covfilter=5))
shared <- SharedJunc(juncs, share_num=2, ext_up=100, ext_down=100)
shared
}
\author{
  Xiaojing Wang
}

