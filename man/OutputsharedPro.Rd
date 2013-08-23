\name{OutputsharedPro}
\alias{OutputsharedPro}
\title{Output the sequences of proteins with high expressions in multiple samples.}
\usage{
  OutputsharedPro(RPKMs, cutoff = "30\%",
    share_sample = "50\%", proteinseq, outfile, ids, ...)
}
\arguments{
  \item{RPKMs}{RPKM matrix; row name (protein name) is
  required.}

  \item{cutoff}{a percentage format cutoff (e.g. '30%'), or
  a vector with each element as a vlaue cutoff referring to
  one sample}

  \item{share_sample}{the minimum share sample numbers for
  proteins which pass the cutoff.}

  \item{proteinseq}{a dataframe containing protein ids and
  protein sequences}

  \item{outfile}{output file name}

  \item{ids}{a dataframe containing gene/transcript/protein
  id mapping information.}

  \item{...}{additional arguments}
}
\value{
  a FASTA file containing proteins with RPKM above the
  cutoff in at least certain number of samples
}
\description{
  Output a FASTA file containing shared proteins with
  expression above cutoff in multiple samples
}
\details{
  this function takes RPKM matrix as input, users can set
  two paramteters,cutoff and shared, to generated a
  consensus expressed database
}
\examples{
path <- system.file("extdata/bams", package="customProDB")
load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
bamFile<- paste(path, '/', list.files(path, pattern="*bam$"), sep='')
rpkms <- sapply(bamFile,function(x)
            calculateRPKM(x, exon, proteincodingonly=TRUE, ids))
outfile <- paste(tempdir(), '/test_rpkm_share.fasta', sep='')
OutputsharedPro(rpkms, cutoff=1, share_sample=2, proteinseq,
            outfile, ids)
}
\author{
  Xiaojing Wang
}

