#' Generate a customized protein database for a single sample.
#'
#' This function provides a convenient way to generate customized databases for a single sample based on transcript
#' expression, SNV and INDEL variants, and novel junctions detected during RNA-Seq mapping.
#' @title An integrated function to generate customized protein databases for a single sample
#' @param bamFile Input BAM file name for calculating transcript FPKMs.
#' @param RPKM Alternative to bamFile. If non-NULL, must be a vector containing expression level for proteins (e.g. FPKMs from cufflinks).
#' @param vcfFile Input VCF file name.
#' @param outfile_path Folder path for the output FASTA files.
#' @param outfile_name Output FASTA file prefix.
#' @param annotation_path The path to read annotation files from (e.g. ids.RData, exon_anno.RData, etc.).
#' @param rpkm_cutoff The cutoff for RPKM values. See 'cutoff' in function Outputproseq for more detail.
#' @param INDEL Set to TRUE to create a FASTA of the variations due to short insertions and deletions in the VCF (default is FALSE).
#' @param lablersid Set to TRUE to include dbSNP rsid as a prefix to each SNP in the FASTA headers (default is FALSE).
#'   If TRUE, the annotation folder must have the dbsnpinCoding.RData file.
#' @param COSMIC Set to TRUE to include COSMIC ids in the variation table (default is FALSE).
#'   If TRUE, the annotation folder must have the cosmic.RData file.
#' @param nov_junction Set to TRUE to create a FASTA of the peptides that cover novel junctions.
#'   If TRUE, the annotation folder must have the splicemax.RData file, and the bedFile and genome 
#'   parameters must be provided.
#' @param bedFile Path to a .bed file which contains the splice junctions identified in RNA-Seq (e.g. junctions.bed from TopHat).
#' @param genome A BSgenome object (e.g. BSgenome.Hsapiens.UCSC.hg19::Hsapiens). Default is NULL.
#' @param ... Additional arguments
#' @return NULL. It creates several FASTA files in <outfile_path>; always creates <outfile_name>_rpkm.fasta
#'   and <outfile_name>_snv.fasta. Optionally creates <outfile_name>_indel.fasta and <outfile_name>_junc.fasta.
#' @author Xiaojing Wang
#' @export
#' @examples
#'
#' bamFile <- system.file("extdata/bams", "test1_sort.bam",
#'             package="customProDB")
#' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
#' annotation_path <- system.file("extdata/refseq", package="customProDB")
#' outfile_path <- tempdir()
#' outfile_name <- 'test'
#'
#' easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, outfile_path,
#'         outfile_name, rpkm_cutoff=1, INDEL=TRUE, lablersid=TRUE,
#'         COSMIC=TRUE, nov_junction=FALSE)
#'
#'
easyRun <- function(bamFile, RPKM=NULL, vcfFile, annotation_path, outfile_path, 
            outfile_name, rpkm_cutoff=1, INDEL=FALSE, lablersid=FALSE, 
            COSMIC=FALSE, nov_junction=FALSE, bedFile=NULL, genome=NULL, 
 ...) {
    if(missing(annotation_path)) {
        stop("must specify the path of annotation files")
    }
    if(nov_junction == TRUE&(is.null(bedFile)|is.null(genome))){
        stop("must supply BED formatted junction file and genome 
            (BSgenome format) if you want to include novel junctions")
    }
    exon <- ''
    ids <- ''
    proteinseq <- ''
    procodingseq <- ''
    
    message("Calculate RPKMs and Output proteins pass the cutoff into FASTA file ... ", 
            appendLF=FALSE)
    load(paste(annotation_path, '/exon_anno.RData', sep=''))
    load(paste(annotation_path, '/ids.RData', sep=''))
    if(is.null(RPKM)){
        RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
    }else RPKM <- RPKM
    load(paste(annotation_path, '/proseq.RData', sep=''))
    outf_rpkm <- paste(outfile_path, '/', outfile_name, '_rpkm.fasta', sep='')
    Outputproseq(RPKM, cutoff=rpkm_cutoff, proteinseq, outf_rpkm, ids)
    packageStartupMessage(" done")
    
    
    load(paste(annotation_path, '/procodingseq.RData', sep=''))

    readvcf <- VariantAnnotation::readVcf(vcfFile, row.names=FALSE,
                                          param=VariantAnnotation::ScanVcfParam(geno=NA, info=NA))
    
    # read REF and ALT columns with the super-fast data.table::fread
    headerLinesToSkip = grep("#CHROM", readLines(vcfFile, n=1000))
    vcftable = .temp_unzip(vcfFile, data.table::fread,
                           skip=headerLinesToSkip-1, sep="\t",
                           select=c("#CHROM", "POS", "REF", "ALT"),
                           showProgress=FALSE)
    ref = toupper(vcftable$REF)
    alt = toupper(vcftable$ALT)
    variantTypes = variantType(ref, alt)
    
    vcfRanges = SummarizedExperiment::rowRanges(readvcf)
    GenomicRanges::mcols(vcfRanges)$REF = ref
    GenomicRanges::mcols(vcfRanges)$ALT = alt
    
    if(lablersid) {
      dbsnpinCoding <- ''
      load(paste0(annotation_path, '/dbsnpinCoding.RData'))
    }else{
      dbsnpinCoding <- NULL
    }
    
    if(COSMIC){
      cosmic <- ''
      load(paste0(annotation_path, '/cosmic.RData'))
    }
    else {
      cosmic <- NULL
    }

    snpVariants = vcfRanges[which(variantTypes == "snp")]
    if (length(snpVariants) > 0)
    {
      postable_snv = Positionincoding(snpVariants, exon, dbsnp=dbsnpinCoding, COSMIC=cosmic)
      if (nrow(postable_snv) > 0)
      {
        message("Output variation table and variant protein sequence caused by SNVs... ", appendLF=FALSE)
        firstTxId_snp = postable_snv[, .(txid=min(txid)), txname]
        codingseq_snp = procodingseq[procodingseq$tx_id %in% firstTxId_snp$txid, ]
        
        variantTable = aaVariation(postable_snv, codingseq_snp)
        OutputVarproseq(variantTable, proteinseq,
                        paste0(outfile_path, '/', outfile_name, "_snv.fasta"),
                        ids, lablersid=lablersid, RPKM=RPKM)
        packageStartupMessage(" done")
      }
    }
    
    indelVariants = vcfRanges[which(variantTypes == "ins" | variantTypes == "del")]
    if (length(indelVariants) > 0)
    {
      postable_indel = Positionincoding(indelVariants, exon, dbsnp=dbsnpinCoding, COSMIC=cosmic)
      if (nrow(postable_indel) > 0)
      {
        message("Output abberant protein FASTA file caused by short INDEL... ", appendLF=FALSE)
        firstTxId_indel = postable_indel[, .(txid=min(txid)), txname]
        codingseq_indel = procodingseq[procodingseq$tx_id %in% firstTxId_indel$txid, ]
        indelvariants = Outputaberrant(postable_indel,
                                       paste0(outfile_path, '/', outfile_name, '_indel.fasta'),
                                       codingseq_indel,
                                       proteinseq,
                                       ids, RPKM=RPKM)
        packageStartupMessage(" done")
      }
    }

    if(nov_junction == TRUE&!is.null(bedFile)&!is.null(genome)){
        message("Output novel junction peptides... ", appendLF=FALSE)
        splicemax <- ''
        load(paste(annotation_path, '/splicemax.RData', sep=''))
        txdb <- loadDb(paste(annotation_path, '/txdb.sqlite', sep=''))
        jun <-  Bed2Range(bedFile, skip=1, covfilter=5)
        junction_type <- JunctionType(jun, splicemax, txdb, ids)
        outf_junc <- paste(outfile_path, '/', outfile_name, 
                        '_junc.fasta', sep='')
        OutputNovelJun(junction_type, genome, outf_junc, 
                        proteinseq)
        packageStartupMessage(" done")
    }
    
    #message("Combine database... ", appendLF=FALSE)
    #files<- paste(outfile_path, 
    #            list.files(outfile_path, '/', pattern="*fasta$"), sep='')
    #pos_files <- c(outf_rpkm, outf_indel, outf_snv, outf_junc)
    #files_idx <- which(files %in% pos_files)
    #files <- files[files_idx]
    #cmd1 <- paste("cat ", paste(files,collapse=' '), ' > ', outfile_path, 
    #           outfile_name, '_combined.fasta', sep='')
    #p <- pipe(cmd1, "r")
    #close(p)
    #packageStartupMessage(" done")
    
 }

