\name{Bed2Range}
\alias{Bed2Range}
\title{Generate a GRanges objects from BED file.}
\usage{
  Bed2Range(bedfile, skip = 1, covfilter = 5, ...)
}
\arguments{
  \item{bedfile}{a character contains the path and name of
  a BED file.}

  \item{skip}{the number of lines of the BED file to skip
  before beginning to read data, default 1.}

  \item{covfilter}{the number of minimum coverage for the
  candidate junction, default 5.}

  \item{...}{additional arguments}
}
\value{
  a GRanges object containing all candidate junctions from
  the BED file.
}
\description{
  Read BED file into a GRanges object. This function
  requires complete BED file. Go to
  https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for
  more information about BED format.
}
\details{
  Read BED file contain junctions into a GRanges object.
}
\examples{
bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
jun <-  Bed2Range(bedfile, skip=1,covfilter=5)
length(jun)
}
\author{
  Xiaojing Wang
}

