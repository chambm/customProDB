##' Generate consensus protein database for multiple samples in a single function.
##'
##' The function give a more convenient way for proteinomics researchers to generate customized database of multiple samples.
##' @title An integrated function to generate consensus protein database from multiple samples
##' @param bamFile_path The path of BAM files
##' @param RPKM_mtx Alternative to bamFile_path,default NULL, a matrix containing expression level for proteins in each sample. (e.g. FPKMs from cufflinks) 
##' @param vcfFile_path The path of VCF files
##' @param annotation_path The path of already saved annotation, which will be used in the function
##' @param rpkm_cutoff Cutoffs of RPKM values. see 'cutoff' in function OutputsharedPro for more information
##' @param share_num The minimum share sample numbers for proteins which pass the cutoff.
##' @param var_shar_num Minimum sample number of recurrent variations.
##' @param outfile_path The path of output FASTA file
##' @param outfile_name The name prefix of output FASTA file
##' @param INDEL If the vcfFile contains the short insertion/deletion. Default is FALSE.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Users should provide dbSNP information when running function Positionincoding() if put TRUE here.
##' @param COSMIC If output the cosmic ids in the variation table.Default is FALSE. If choose TRUE, there must have cosmic.RData in the annotation folder. 
##' @param nov_junction If output the peptides that cover novel junction into the database. if TRUE, there should be splicemax.RData in the annotation folder.
##' @param bedFile_path The path of BED files which contains the splice junctions identified in RNA-Seq. 
##' @param genome A BSgenome object(e.g. Hsapiens). Default is NULL. Required if nov_junction==TRUE.
##' @param junc_shar_num Minimum sample number of recurrent splicing junctions. 
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
##' easyRun_mul(bampath, RPKM_mtx=NULL, vcfFile_path, annotation_path, rpkm_cutoff=1,
##'             share_num=2, var_shar_num=2, outfile_path, outfile_name, INDEL=TRUE,
##'             lablersid=TRUE, COSMIC=TRUE, nov_junction=FALSE)
##' 
##' 
##'
 
easyRun_mul <- function(bamFile_path, RPKM_mtx=NULL, vcfFile_path, 
    annotation_path, rpkm_cutoff, share_num=2, var_shar_num=2, outfile_path, 
    outfile_name, INDEL=FALSE, lablersid=FALSE, COSMIC=FALSE, nov_junction=FALSE, 
    bedFile_path=NULL, genome=NULL, junc_shar_num=2,
 ...) {
    if(missing(annotation_path)) {
        stop("must specify the path of annotation files")
    }
    if(nov_junction == TRUE&(is.null(bedFile_path)|is.null(genome))){
        stop("must supply BED formatted junction file and genome 
            (BSgenome format) if you want to include novel junctions")
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
            chrlist <- paste('chr',c(seq(1:22),'X','Y'),sep='')
            indexchr <-which(postable_indel[,'chr'] %in% chrlist)
            postable_indel <- postable_indel[indexchr,]

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
    OutputVarproseq_single(mtab, proteinseq, outf_snv, ids, RPKM=meanRPKM)
    
    packageStartupMessage(" done")

    if(nov_junction == TRUE&!is.null(bedFile_path)&!is.null(genome)){
        message("Output novel junction peptides... ", appendLF=FALSE)
        splicemax <- ''
        load(paste(annotation_path, '/splicemax.RData', sep=''))
        txdb <- loadDb(paste(annotation_path, '/txdb.sqlite', sep=''))
        bedFiles<- paste(bedFile_path, '/', list.files(bedFile_path, pattern="*bed$"), 
                    sep='')
        juncs <- lapply(bedFiles, function(x) Bed2Range(x, skip=1, covfilter=5))
        shared <- SharedJunc(juncs, junc_shar_num=2, ext_up=100, ext_down=100)
        junction_type <- JunctionType(shared, splicemax, txdb, ids)
        outf_junc <- paste(outfile_path, '/', outfile_name, 
                        '_junc.fasta', sep='')
        OutputNovelJun(junction_type, genome, outf_junc, 
                        proteinseq)
        packageStartupMessage(" done")
    }    

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




