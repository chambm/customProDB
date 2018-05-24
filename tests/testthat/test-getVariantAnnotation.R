library(testthat)
library(customProDB)
library(data.table)

context("getVariantAnnotation")


test_that("aaVariation works", {

    codingseq = data.table(read.table(header=TRUE, stringsAsFactors=FALSE,text=
"tx_id pro_name tx_name                      coding
     1     Pro1     Tx1 ATGCAGCGCCGGGACGATCCTGCCTAG
     2     Pro2     Tx2 ATGCAGCGCCGGGACGATCCTGCCTAG
     3     Pro3     Tx3 ATGCAGCGCCGGGACGATCCTGCCTAG
     4     Pro4     Tx4 ATGCAGCGCCGGGACGATCCTGCCTAG
     5     Pro5     Tx5 ATGCAGCGCCGGGACGATCCTGCCTAG
     6     Pro6     Tx6 ATGCAGCGCCGGGACGATCCTGCCTAG
"), key="tx_id")

    # 123 456 789 10  13  16  19  22  25 
    # ATG CAG CGC CGG GAC GAT CCT GCC TAG
    #  M   Q   R   R   D   D   P   A   *

    postable = data.table(read.table(header=TRUE, stringsAsFactors=FALSE, flush=TRUE, text=
"genename txname txid proname chr strand pos refbase varbase pincoding rsid 
    Gene1    Tx1    1    Pro1   1      + 123       C       G         4   r1  # snp
    Gene2    Tx2    2    Pro2   2      + 234       A     G,T         5   r2  # ambiguous snp
    Gene2    Tx2    2    Pro2   2      + 234       C     A,T        15   NA  # ambiguous snp, half synonymous
    Gene3    Tx3    3    Pro3   3      - 345       G       T         6   r3  # reverse complement snp

    Gene4    Tx4    4    Pro4   4      + 456       C       A         7   NA  # 2 variants on 1 codon
    Gene4    Tx4    4    Pro4   4      + 456       G       A         8   r4

    Gene5    Tx5    5    Pro5   5      - 567       C       A         7   r5  # 2 codons with 1 variant
    Gene5    Tx5    5    Pro5   5      - 567       C       T        10   NA

    Gene6    Tx6    6    Pro6   6      + 789       C       T         4   r6  # stop gained
"), key="txid")

    variantTable = aaVariation(postable, codingseq)
    
    new_columns = data.table(read.table(header=TRUE, stringsAsFactors=FALSE, text=
"rsid CodonStart CodonVariants RefCodon varcode        vartype aaref aapos aavar
   r1          4           4:G      CAG     GAG non-synonymous     Q     2     E
   r2          4           5:K      CAG     CKG non-synonymous     Q     2     R,L
   NA         13          15:W      GAC     GAW non-synonymous     D     5     E
   r3          4           6:T      CAG     CAA     synonymous     Q     2     Q

   r4          7       7:A:8:A      CGC     AAC non-synonymous     R     3     N

   r5          7           7:A      CGC     TGC non-synonymous     R     3     C
   NA         10          10:T      CGG     AGG     synonymous     R     4     R

   r6          4           4:T      CAG     TAG non-synonymous     Q     2     *
"))
    expected_table = data.table(variantTable[, names(new_columns)])
    
    expect_equal(new_columns, expected_table)
})


test_that("getVariantAnnotation works", {
  
  load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
  load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
  load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
  load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
  
  variantAnnotation = getVariantAnnotation(system.file("extdata/vcfs", "test1.vcf", package="customProDB"),
                                           ids, exon,
                                           proteinseq, procodingseq,
                                           dbsnpinCoding, cosmic)
  
  expect_equal_to_reference(variantAnnotation$variantTable,
                            "test1-hg19-aaVariation.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
})


test_that("snvproseq, snvprocodoing, indelproseq, indelprocoding are created correctly", {
  
  load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
  load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
  load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
  load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
  load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
  
  variantAnnotation = getVariantAnnotation(system.file("extdata/vcfs", "test1.vcf", package="customProDB"),
                                           ids, exon,
                                           proteinseq, procodingseq,
                                           dbsnpinCoding, cosmic)
  
  expect_equal_to_reference(variantAnnotation$snvprocoding,
                            "test1-hg19-snvprocoding.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$snvproseq,
                            "test1-hg19-snvproseq.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$indelprocoding,
                            "test1-hg19-indelprocoding.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
  
  expect_equal_to_reference(variantAnnotation$indelproseq,
                            "test1-hg19-indelproseq.rds",
                            on.update=on.update.view, on.fail=on.fail.diff)
})


test_that("variantType returns correct variant types", {
  expect_equal(variantType("A", "C"), "snp")
  expect_equal(variantType("A", "C,G"), "snp")
  expect_equal(variantType("AA", "CT"), "mnp")
  expect_equal(variantType("AT", "A"), "del")
  expect_equal(variantType("AT", "T"), "del")
  expect_equal(variantType("A", "AT"), "ins")
  expect_equal(variantType("A", "TA"), "ins")
  expect_equal(variantType("AA", "C"), "mix")
  expect_equal(variantType("AA", "CAT"), "mix")
  expect_equal(variantType("AA", "TAAC"), "mix")
  expect_equal(variantType("ATTA", "ATA"), "del")
  expect_equal(variantType("AAA", "AATA"), "ins")
  expect_equal(variantType("A", "T,AT"), "mix")
  expect_equal(variantType("TAA", "TAA,TAAAA"), "ins")
  expect_equal(variantType("c", "cT"), "ins")
})
