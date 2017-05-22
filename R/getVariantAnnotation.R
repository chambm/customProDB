#' Get genome annotations for variants called in one or more VCF files. Each sample is called separately
#' so that variants on the same codon are handled properly. This is a convenience function which calls
#' OutputVarprocodingseq, OutputVarproseq, Outputaberrant, and OutputNovelJun.
#'
#' @title Get genome annotations for variants called in one or more VCF files.
#' @param vcfFilepaths one or more VCF files to create variant proteins
#' @param ids a dataframe mapping between gene/transcript/protein ids
#' @param exon_anno a dataframe of exon annotations for each protein coding transcript
#' @param proteinseq a dataframe of reference amino acid sequences for each protein
#' @param procodingseq a dataframe of reference coding sequences for each protein coding transcript
#' @param dbsnpinCoding a GRanges object of dbSNP variants to include in the output tables; default is NULL
#' @param cosmic a GRanges object of COSMIC variants to include in the output tables; default is NULL
#'
#' @return A list of annotation tables with named members:
#' \describe{
#' \item{variantTable}{describes the SNP variants and their coding consequences}
#' \item{snvprocoding}{the CDS sequences for the variant proteins from SNPs}
#' \item{snvproseq}{the translated sequences for the variant proteins from SNPs}
#' \item{indelprocoding}{the CDS sequences for the variant proteins from INDELs}
#' \item{indelproseq}{the translated sequences for the variant proteins from INDELs}
#' }
#' @export
#'
#' @examples
#' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
#' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
#' load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
#' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
#' load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
#' load(system.file("extdata/refseq", "cosmic.RData", package="customProDB"))
#' vcfFilepath = system.file("extdata/vcfs", "test1.vcf", package="customProDB")
#' 
#' variantAnnotation = getVariantAnnotation(vcfFilepath, ids, exon,
#'                                          proteinseq, procodingseq,
#'                                          dbsnpinCoding, cosmic)
#' head(variantAnnotation$variantTable)
#' head(variantAnnotation$snvprocoding)
#' head(variantAnnotation$snvproseq)
#' head(variantAnnotation$indelprocoding)
#' head(variantAnnotation$indelproseq)
getVariantAnnotation <- function(vcfFilepaths,
                                 ids, exon_anno,
                                 proteinseq, procodingseq,
                                 dbsnpinCoding = NULL,
                                 cosmic = NULL)
{
  stopifnot(all(file.exists(vcfFilepaths)))
  
  variantTable = vector('list', length(vcfFilepaths))
  snvprocoding = vector('list', length(vcfFilepaths))
  snvproseq = vector('list', length(vcfFilepaths))
  indelprocoding = vector('list', length(vcfFilepaths))
  indelproseq = vector('list', length(vcfFilepaths))
  
  label_rsid = !is.null(dbsnpinCoding)
  
  for(i in 1:length(vcfFilepaths))
  {
    vcfFilepath = vcfFilepaths[i]
    readvcf <- VariantAnnotation::readVcf(vcfFilepath, row.names=FALSE,
                                          param=VariantAnnotation::ScanVcfParam(geno=NA, info=NA))
    
    # read REF and ALT columns with the super-fast data.table::fread
    headerLinesToSkip = grep("#CHROM", readLines(vcfFilepath, n=500))
    vcftable = .temp_unzip(vcfFilepath, data.table::fread,
                           skip=headerLinesToSkip-1, sep="\t", select=c("#CHROM", "POS", "REF", "ALT"))
    ref = toupper(vcftable$REF)
    alt = toupper(vcftable$ALT)
    variantTypes = variantType(ref, alt)
    
    vcfRanges = SummarizedExperiment::rowRanges(readvcf)
    GenomicRanges::mcols(vcfRanges)$REF = ref
    GenomicRanges::mcols(vcfRanges)$ALT = alt
    
    snpVariants = vcfRanges[which(variantTypes == "snp")]
    if (length(snpVariants) > 0)
    {
      postable_snv = Positionincoding(snpVariants, exon_anno, dbsnp=dbsnpinCoding, COSMIC=cosmic)
      if (nrow(postable_snv) > 0)
      {
        firstTxId_snp = postable_snv[, .(txid=min(txid)), txname]
        codingseq_snp = procodingseq[procodingseq$tx_id %in% firstTxId_snp$txid, ]
        
        variantTable[[i]] = aaVariation(postable_snv, codingseq_snp)
        snvprocoding[[i]] = OutputVarprocodingseq(variantTable[[i]], codingseq_snp, ids, lablersid=label_rsid)
        snvproseq[[i]] = OutputVarproseq(variantTable[[i]], proteinseq,
                                         file.path(tempdir(), "snv.fasta"),
                                         ids, lablersid=label_rsid)
      }
    }
    
    indelVariants = vcfRanges[which(variantTypes == "ins" | variantTypes == "del")]
    if (length(indelVariants) > 0)
    {
      postable_indel = Positionincoding(indelVariants, exon_anno, dbsnp=dbsnpinCoding, COSMIC=cosmic)
      if (nrow(postable_indel) > 0)
      {
        firstTxId_indel = postable_indel[, .(txid=min(txid)), txname]
        codingseq_indel = procodingseq[procodingseq$tx_id %in% firstTxId_indel$txid, ]
        indelvariants = Outputaberrant(postable_indel, file.path(tempdir(), "indel.fasta"),
                                       codingseq_indel, proteinseq, ids)
        indelprocoding[[i]] = indelvariants$indelprocoding
        indelproseq[[i]] = indelvariants$indelproseq
      }
    }
  }
  
  list(variantTable = unique(do.call("rbind", variantTable)),
       snvprocoding = unique(do.call("rbind", snvprocoding)),
       snvproseq = unique(do.call("rbind", snvproseq)),
       indelprocoding = unique(do.call("rbind", indelprocoding)),
       indelproseq = unique(do.call("rbind", indelproseq)))
}

.variantType = function(ref, alt)
{
  refLength = stri_length(ref)
  altLength = stri_length(alt)
  
  alternates = c("A"="^[TGC](,[TGC])*$",
                 "T"="^[AGC](,[AGC])*$",
                 "G"="^[ATC](,[ATC])*$",
                 "C"="^[ATG](,[ATG])*$")
  
  if(refLength==1 && grepl(alternates[ref], alt))
    "snp"
  else if(refLength>1 && refLength==altLength)
    "mnp"
  else if(refLength > altLength && (stri_startswith_fixed(ref, alt) || stri_endswith_fixed(ref, alt)))
    "del"
  else if(refLength < altLength && (stri_startswith_fixed(alt, ref) || stri_endswith_fixed(alt, ref)))
    "ins"
  else
  {
    lcs = Biobase::lcSuffix(c(ref, alt))
    if(stri_length(lcs) > 0)
    {
        ref = stri_sub(ref, 1, -(stri_length(lcs)+1))
        alt = stri_sub(alt, 1, -(stri_length(lcs)+1))
        refLength = stri_length(ref)
        altLength = stri_length(alt)
        if(refLength > altLength && (stri_startswith_fixed(ref, alt) || stri_endswith_fixed(ref, alt)))
          "del"
        else if(refLength < altLength && (stri_startswith_fixed(alt, ref) || stri_endswith_fixed(alt, ref)))
          "ins"
        else
          "mix"
    }
    else
      "mix"
  }
}

#' Get the type of variant given one or more reference alleles and the corresponding variant alleles.
#'
#' @param ref the reference alleles
#' @param alt the alternate (variant) alleles
#'
#' @return variant type for each ref/alt pair, one of:
#' \enumerate{
#' \item snp (single nucleotide polymorphism)
#' \item mnp (multiple nucleotide polymorphism)
#' \item ins (insertion)
#' \item del (deletion)
#' \item mix (some combination of the above)
#' }
#' @export
#' @importFrom stringi stri_length stri_sub stri_startswith_fixed stri_endswith_fixed
#'
#' @examples
#' ref = c("A", "A",   "AA", "AT", "A",  "AA",  "AA",  "AA",  "ATTA", "AAA",  "A",    "TAAA")
#' alt = c("C", "C,G", "CT", "A",  "AT", "C",  "CAT", "TAAC", "ATA",  "AATA", "T,AT", "TAA,TAAAA")
#' print(cbind(ref, alt, variantType(ref, alt)))
variantType = Vectorize(.variantType, USE.NAMES=FALSE)
