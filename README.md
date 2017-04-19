# customProDB

[![TravisCI Build Status](https://travis-ci.org/chambm/customProDB.svg?branch=master)](https://travis-ci.org/chambm/customProDB)
[![Coverage Status](https://img.shields.io/codecov/c/github/chambm/customProDB/master.svg)](https://codecov.io/github/chambm/customProDB?branch=master)
[![Platforms version](http://bioconductor.org/shields/availability/3.4/customProDB.svg)](http://bioconductor.org/packages/devel/bioc/html/customProDB.html)
<!---
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/hadley/testthat?branch=master&svg=true)](https://ci.appveyor.com/project/hadley/testthat)
-->

Package Homepage: http://bioconductor.org/packages/devel/bioc/html/customProDB.html
Bug Reports: https://support.bioconductor.org/p/new/post/?tag_val=customProDB

CustomProDB enables the easy generation of customized databases from RNA-Seq data for proteomics search. It bridges 
genomics and proteomics studies and facilitates cross-omics data integration.

Database search is the most widely used approach for peptide and protein identification in mass-
spectrometry-based proteomics studies. Sample-specific protein databases derived from RNA-Seq data can 
better approximate the real protein pools in the samples and thus improve protein identification. 
More importantly, single nucleotide variations, short insertion and deletions and novel junctions 
identified from RNA-Seq data make protein databases more complete and sample-specific. 

## Installation

To install this package, start R and enter:

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("customProDB")
```

Alternatively, you can install the latest version from github using devtools:

```R
# install.packages("devtools")
devtools::install_github("chambm/customProDB")
```

## Documentation

To view documentation for the version of this package installed in your system, start R and enter:
```R
browseVignettes("customProDB")
```