\name{easyRun_mul}
\alias{easyRun_mul}
\title{An integrated function to generate consensus protein database from multiple samples}
\usage{
  easyRun_mul(bamFile_path, RPKM_mtx = NULL, vcfFile_path,
    annotation_path, rpkm_cutoff, share_num = 2,
    var_shar_num = 2, outfile_path, outfile_name,
    INDEL = FALSE, lablersid = FALSE, COSMIC = FALSE,
    nov_junction = FALSE, bedFile_path = NULL,
    genome = NULL, junc_shar_num = 2, ...)
}
\arguments{
  \item{bamFile_path}{The path of BAM files}

  \item{RPKM_mtx}{Alternative to bamFile_path,default NULL,
  a matrix containing expression level for proteins in each
  sample. (e.g. FPKMs from cufflinks)}

  \item{vcfFile_path}{The path of VCF files}

  \item{annotation_path}{The path of already saved
  annotation, which will be used in the function}

  \item{rpkm_cutoff}{Cutoffs of RPKM values. see 'cutoff'
  in function OutputsharedPro for more information}

  \item{share_num}{The minimum share sample numbers for
  proteins which pass the cutoff.}

  \item{var_shar_num}{Minimum sample number of recurrent
  variations.}

  \item{outfile_path}{The path of output FASTA file}

  \item{outfile_name}{The name prefix of output FASTA file}

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

  \item{bedFile_path}{The path of BED files which contains
  the splice junctions identified in RNA-Seq.}

  \item{genome}{A BSgenome object(e.g. Hsapiens). Default
  is NULL. Required if nov_junction==TRUE.}

  \item{junc_shar_num}{Minimum sample number of recurrent
  splicing junctions.}

  \item{...}{Additional arguments}
}
\value{
  A table file contains detailed variation information and
  several FASTA files.
}
\description{
  Generate consensus protein database for multiple samples
  in a single function.
}
\details{
  The function give a more convenient way for proteinomics
  researchers to generate customized database of multiple
  samples.
}
\examples{
bampath <- system.file("extdata/bams", package="customProDB")
vcfFile_path <- system.file("extdata/vcfs", package="customProDB")
annotation_path <- system.file("extdata/refseq", package="customProDB")
outfile_path <- tempdir()
outfile_name <- 'mult'

easyRun_mul(bampath, RPKM_mtx=NULL, vcfFile_path, annotation_path, rpkm_cutoff=1,
            share_num=2, var_shar_num=2, outfile_path, outfile_name, INDEL=TRUE,
            lablersid=TRUE, COSMIC=TRUE, nov_junction=FALSE)
}
\author{
  Xiaojing Wang
}

