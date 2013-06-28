\name{OutputNovelJun}
\alias{OutputNovelJun}
\title{generate peptide FASTA file that contains novel junctions.}
\usage{
  OutputNovelJun(junction_type, genome, outfile,
    proteinseq, ...)
}
\arguments{
  \item{junction_type}{a data frame which is the output of
  function JunctionType()}

  \item{genome}{a BSgenome object. (e.g. Hsapiens)}

  \item{outfile}{output file name}

  \item{proteinseq}{a data frame cotaining amino acid
  sequence for each protein.}

  \item{...}{Additional arguments.}
}
\value{
  FASTA file that contains novel junction peptides.
}
\description{
  Three-frame translation of novel junctions. And remove
  those could be found in normal protein sequences. This
  function requires a genome built by BSgenome package.
}
\examples{
bedfile <- system.file("extdata", "junctions.bed", package="customProDB")
load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite",
    package="customProDB"))
junction_type <- JunctionType(bedfile, skip=1, covfilter=5, splicemax, txdb, ids)
table(junction_type[, 'jun_type'])
outf_junc <- paste(tempdir(), '/test_junc.fasta', sep='')
load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
library('BSgenome.Hsapiens.UCSC.hg19')
OutputNovelJun <- OutputNovelJun(junction_type, Hsapiens, outf_junc,
            proteinseq)
}
\author{
  Xiaojing Wang
}

