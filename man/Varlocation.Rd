\name{Varlocation}
\alias{Varlocation}
\title{Annotates the variations with genomic location.}
\usage{
  Varlocation(Vars, txdb, ids, ...)
}
\arguments{
  \item{Vars}{a GRange object of variations}

  \item{txdb}{a TranscriptDb object.}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information}

  \item{...}{additional arguments}
}
\value{
  a data frame of locations for each variation
}
\description{
  For a given GRange object of variations, the
  Varlocation() function finds the genomic locations for
  each entry according to the given annotation. Seven
  labels are used to describe the location (intergenic,
  intro_nonProcoding, exon_nonProcoding, intron, 5utr, 3utr
  and coding). details of the definition can be found in
  the tutorial.
}
\details{
  see 'introduction' for more details
}
\examples{
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
vcf <- InputVcf(vcffile)

table(values(vcf[[1]])[['INDEL']])
index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
indelvcf <- vcf[[1]][index]

index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
SNVvcf <- vcf[[1]][index]

txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
SNVloc <- Varlocation(SNVvcf,txdb,ids)
indelloc <- Varlocation(indelvcf,txdb,ids)
table(SNVloc[,'location'])
}
\author{
  Xiaojing Wang
}

