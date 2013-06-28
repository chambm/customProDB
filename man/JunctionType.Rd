\name{JunctionType}
\alias{JunctionType}
\title{Annotates the junctions in a bed file.}
\usage{
  JunctionType(bedfile, skip = 1, covfilter = 5, splicemax,
    txdb, ids, ...)
}
\arguments{
  \item{bedfile}{bed file path}

  \item{skip}{default 1. skip the annotation rows in the
  bed file.}

  \item{covfilter}{filter the junctions by reads coverage,
  default 5.}

  \item{splicemax}{a known exon splice matrix from the
  annotation.}

  \item{txdb}{a TranscriptDb object.}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information.}

  \item{...}{additional arguments}
}
\value{
  a data frame of type and source for each junction.
}
\description{
  For identified splice junctions from RNA-Seq, the
  JunctionType() function finds the genomic locations for
  each entry according to the given annotation. four types
  of junctions are classified. find more details in the
  tutorial. This function requires complete BED file.
}
\details{
  Go to https://genome.ucsc.edu/FAQ/FAQformat.html#format1
  for more information about BED format.
}
\examples{
bedfile <- system.file("extdata", "junctions.bed", package="customProDB")
load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite",
            package="customProDB"))
junction_type <- JunctionType(bedfile, skip=1, covfilter=5, splicemax, txdb,
                    ids)
table(junction_type[, 'jun_type'])
}
\author{
  Xiaojing Wang
}

