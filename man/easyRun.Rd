\name{easyRun}
\alias{easyRun}
\title{An integrated function to generate customized protein database for a single sample}
\usage{
  easyRun(bamFile, RPKM = NULL, vcfFile, annotation_path,
    outfile_path, outfile_name, rpkm_cutoff = 1,
    INDEL = FALSE, lablersid = FALSE, COSMIC = FALSE,
    nov_junction = FALSE, bedFile = NULL, genome = NULL,
    ...)
}
\arguments{
  \item{bamFile}{Input BAM file name}

  \item{RPKM}{Alternative to bamFile,default NULL, a vector
  containing expression level for proteins. (e.g. FPKMs
  from cufflinks)}

  \item{vcfFile}{Input VCF file name.}

  \item{outfile_path}{Folder path for the output FASTA
  files.}

  \item{outfile_name}{Output FASTA file name.}

  \item{annotation_path}{The path of saved annotation.}

  \item{rpkm_cutoff}{The cutoff of RPKM value. see 'cutoff'
  in function Outputproseq for more detail.}

  \item{INDEL}{If the vcfFile contains the short
  insertion/deletion. Default is FALSE.}

  \item{lablersid}{If includes the dbSNP rsid in the header
  of each sequence, default is FALSE.  Users should provide
  dbSNP information when running function
  Positionincoding() if put TRUE here.}

  \item{COSMIC}{If output the cosmic ids in the variation
  table.Default is FALSE. If choose TRUE, there must have
  cosmic.RData in the annotation folder.}

  \item{nov_junction}{If output the peptides that cover
  novel junction into the database. if TRUE, there should
  be splicemax.RData in the annotation folder.}

  \item{bedFile}{The path of bed file which contains the
  splice junctions identified in RNA-Seq.}

  \item{genome}{A BSgenome object(e.g. Hsapiens). Default
  is NULL.}

  \item{...}{Additional arguments}
}
\value{
  A table file contains detailed variation information and
  several FASTA files.
}
\description{
  Generate a customized protein database for a single
  sample.
}
\details{
  The function gives a more convenient way for proteomics
  researchers to generate customized database for a single
  sample.
}
\examples{
bamFile <- system.file("extdata/bams", "test1_sort.bam",
            package="customProDB")
vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")
outfile_path <- tempdir()
outfile_name <- 'test'

easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, outfile_path,
        outfile_name, rpkm_cutoff=1, INDEL=TRUE, lablersid=TRUE,
        COSMIC=TRUE, nov_junction=FALSE)
}
\author{
  Xiaojing Wang
}

