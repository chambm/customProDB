\name{JunctionType}
\alias{JunctionType}
\title{Annotates the junctions in a bed file.}
\usage{
  JunctionType(jun, splicemax, txdb, ids, ...)
}
\arguments{
  \item{jun}{a GRange object for junctions, the output of
  function Bed2Range.}

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
  For identified splice junctions from RNA-Seq, this
  function finds the junction types for each entry
  according to the given annotation. Six types of junctions
  are classified. find more details in the tutorial.
}
\details{
  Go to https://genome.ucsc.edu/FAQ/FAQformat.html#format1
  for more information about BED format.
}
\examples{
bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
jun <-  Bed2Range(bedfile,skip=1,covfilter=5)
load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite",
            package="customProDB"))
junction_type <- JunctionType(jun, splicemax, txdb, ids)
table(junction_type[, 'jun_type'])
}
\author{
  Xiaojing Wang
}

