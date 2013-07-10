##' Generate consensus protein database for multiple samples in a single function.
##'
##' The function give a more convenient way for proteinomics researchers to generate customized database of multiple samples.
##' @title An integrated function to generate consensus protein database from multiple samples
##' @param bamFile_path The path of BAM files
##' @param RPKM_mtx Alternative to bamFile_path,default NULL, a matrix containing expression level for proteins in each sample. (e.g. FPKMs from cufflinks) 
##' @param vcfFile_path The path of VCF files
##' @param annotation_path The path of already saved annotation, which will be used in the function
##' @param rpkm_cutoff Cutoffs of RPKM values. see 'cutoff' in function OutputsharedPro for more information
##' @param share_num Which snp dataset you want to use for the SNP annotation
##' @param var_shar_num Which snp dataset you want to use for the SNP annotation
##' @param outfile_path The path of output FASTA file
##' @param outfile_name The name of output FASTA file
##' @param INDEL If the vcfFile contains the short insertion/deletion. Default is FALSE.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Users should provide dbSNP information when running function Positionincoding() if put TRUE here.
##' @param COSMIC If output the cosmic ids in the variation table.Default is FALSE. If choose TRUE, there must have cosmic.RData in the annotation folder. 
##' @param ... Additional arguments
##' @return A table file contains detailed variation information and several FASTA files.
##' @author Xiaojing Wang
##' @examples
##' 
##' bampath <- system.file("extdata/bams", package="customProDB")
##' vcfFile_path <- system.file("extdata/vcfs", package="customProDB")
##' annotation_path <- system.file("extdata/refseq", package="customProDB")
##' outfile_path <- tempdir()    
##' outfile_name <- 'mult'
##' 
##' easyRun_mul(bampath,RPKM_mtx=NULL,vcfFile_path,annotation_path,rpkm_cutoff=1,
##'             share_num=2,var_shar_num=2,outfile_path, outfile_name,INDEL=TRUE,
##'             lablersid=TRUE,COSMIC=TRUE)
##' 
##' 
##'
 
easyRun_mul <- function(bamFile_path, RPKM_mtx=NULL, vcfFile_path, 
    annotation_path, rpkm_cutoff, share_num, var_shar_num, outfile_path, 
    outfile_name, INDEL=FALSE, lablersid=FALSE, COSMIC=FALSE, 
 ...) {
    if(missing(annotation_path)) {
        stop("must specify the path of annotation files")
    }
    exon <- ''
    ids <- ''
    proteinseq <- ''
    procodingseq <- ''
    dbsnpinCoding <- ''
    message("Calculate RPKMs and Output proteins pass the cutoff in multiple samples into FASTA file ... ", 
            appendLF=FALSE)
    load(paste(annotation_path, '/exon_anno.RData', sep=''))
    load(paste(annotation_path, '/ids.RData', sep=''))
    if(is.null(RPKM_mtx)){
        bamFile<- paste(bamFile_path, '/', 
            list.files(bamFile_path,pattern="*bam$"), sep='')
        rpkms <- sapply(bamFile, 
            function(x) calculateRPKM(x,exon,proteincodingonly=TRUE, ids))
    }else rpkms <- RPKM_mtx
    load(paste(annotation_path, '/proseq.RData', sep=''))
    outf_rpkm <- paste(outfile_path,'/', outfile_name, 
            '_rpkm_shared.fasta', sep='')
    OutputsharedPro(rpkms, cutoff=rpkm_cutoff, share_sample=share_num, 
        proteinseq, outf_rpkm, ids)
    packageStartupMessage(" done")
    
    meanRPKM <- apply(rpkms, 1, mean)
    load(paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
    load(paste(annotation_path, '/procodingseq.RData', sep=''))

    vcfFiles<- paste(vcfFile_path, '/', list.files(vcfFile_path, pattern="*vcf$"), 
                sep='')
    vcfs <- lapply(vcfFiles, function(x) InputVcf(x))
    vcf <- Multiple_VCF(vcfs, var_shar_num)
    
    if(INDEL){
        
        idx_snv <- which(values(vcf)[['INDEL']] == FALSE)
        SNVvcf <- vcf[idx_snv]
        idx_indel <- which(values(vcf)[['INDEL']] == TRUE)
        indelvcf <- vcf[idx_indel]
        
        message("Output abberant protein FASTA file caused by short INDEL... ",
            appendLF=FALSE)
        postable_indel <- Positionincoding(indelvcf, exon)
        outf_indel <- paste(outfile_path, '/', outfile_name, 
                '_indel.fasta', sep='')
        if(!is.null(postable_indel)){
            txlist_indel <- unique(postable_indel[, 'txid'])
            codingseq_indel <- procodingseq[procodingseq[,'tx_id'] %in% 
                    txlist_indel, ]
            Outputaberrant(postable_indel, coding=codingseq_indel, 
                        proteinseq=proteinseq, outfile=outf_indel, ids=ids, 
                        RPKM=meanRPKM)
        }
        
        packageStartupMessage(" done")
    }else{
        SNVvcf <- vcf
    }
    
    message("Output Output variation table and variant protein sequence caused by SNVs... ", 
            appendLF=FALSE)
    if(lablersid){
        dbsnpinCoding <- ''
        load(paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
        if(COSMIC){
            cosmic <- ''
            load(paste(annotation_path, '/cosmic.RData', sep=''))
            postable_snv <- Positionincoding(SNVvcf, exon, dbsnp=dbsnpinCoding, 
                    COSMIC=cosmic)
        }else{
            postable_snv <- Positionincoding(SNVvcf, exon, dbsnp=dbsnpinCoding)
        }
    }else{
        if(COSMIC){
            cosmic <- ''
            load(paste(annotation_path, '/cosmic.RData', sep=''))
            postable_snv <- Positionincoding(SNVvcf, exon, dbsnp=NULL, 
                        COSMIC=cosmic)
        }else{
            postable_snv <- Positionincoding(SNVvcf, exon)
        }
    }
    
    txlist <- unique(postable_snv[, 'txid'])
    codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist, ]
    mtab <- aaVariation (postable_snv, codingseq)
    
    outf_mtab <- paste(outfile_path,'/', outfile_name, '_snv.tab', sep='')
    write.table(mtab, file=outf_mtab, sep='\t', quote=F, row.names=F)
    
    outf_snv <- paste(outfile_path,'/', outfile_name, '_snv.fasta', sep='')
    OutputVarproseq(mtab, proteinseq, outf_snv, ids, RPKM=meanRPKM)
    
    packageStartupMessage(" done")



    #message("Combine database... ", appendLF=FALSE)
    #files<- paste(outfile_path,'/', 
    #        list.files(outfile_path, pattern="*fasta$"),sep='')
    #pos_files <- c(outf_rpkm, outf_indel, outf_snv)
    #files_idx <- which(files %in% pos_files)
    #files <- files[files_idx]
    #cmd1 <- paste("cat ", paste(files,collapse=' '), ' > ', outfile_path, '/', 
    #               outfile_name, '_combined.fasta', sep='')
    #p <- pipe(cmd1, "r")
    #close(p)
    #packageStartupMessage(" done")

 }




