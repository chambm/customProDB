##' Generate a customized protein database for a single sample.
##'
##' The function gives a more convenient way for proteomics researchers to generate customized database for a single sample.
##' @title An integrated function to generate customized protein database for a single sample
##' @param bamFile Input BAM file name
##' @param RPKM Alternative to bamFile,default NULL, a vector containing expression level for proteins. (e.g. FPKMs from cufflinks) 
##' @param vcfFile Input VCF file name.
##' @param outfile_path Folder path for the output FASTA files.
##' @param outfile_name Output FASTA file name.
##' @param annotation_path The path of saved annotation.
##' @param rpkm_cutoff The cutoff of RPKM value. see 'cutoff' in function Outputproseq for more detail.
##' @param INDEL If the vcfFile contains the short insertion/deletion. Default is FALSE.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Users should provide dbSNP information when running function Positionincoding() if put TRUE here.
##' @param COSMIC If output the cosmic ids in the variation table.Default is FALSE. If choose TRUE, there must have cosmic.RData in the annotation folder. 
##' @param nov_junction If output the peptides that cover novel junction into the database. if TRUE, there should be splicemax.RData in the annotation folder.
##' @param bedFile The path of bed file which contains the splice junctions identified in RNA-Seq.
##' @param genome A BSgenome object(e.g. Hsapiens). Default is NULL.
##' @param ... Additional arguments
##' @return A table file contains detailed variation information and several FASTA files. 
##' @author Xiaojing Wang
##' @examples
##' \dontrun{
##' bamFile <- system.file("extdata/bams", "test1_sort.bam", 
##'             package="customProDB")
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' annotation_path <- system.file("extdata/refseq", package="customProDB")
##' outfile_path <- tempdir()
##' outfile_name <- 'test'
##' 
##' easyRun(bamFile, RPKM=NULL, vcffile, annotation_path, outfile_path, 
##'         outfile_name, rpkm_cutoff=1, INDEL=TRUE, lablersid=TRUE, 
##'         COSMIC=TRUE, nov_junction=FALSE) 
##'    }
##'
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

    vcf <- InputVcf(vcfFile)
    #table(values(vcf[[1]])[['INDEL']])
    if(INDEL){
        
        idx_snv <- which(values(vcf[[1]])[['INDEL']] == FALSE)
        SNVvcf <- vcf[[1]][idx_snv]
        idx_indel <- which(values(vcf[[1]])[['INDEL']] == TRUE)
        indelvcf <- vcf[[1]][idx_indel]
        
        message("Output abberant protein FASTA file caused by short INDEL... ",
            appendLF=FALSE)
        postable_indel <- Positionincoding(indelvcf, exon)
        outf_indel <- paste(outfile_path, '/', outfile_name, '_indel.fasta', 
                        sep='')
        if(!is.null(postable_indel)){
            txlist_indel <- unique(postable_indel[, 'txid'])
            codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% 
                        txlist_indel, ]
            Outputaberrant(postable_indel, coding=codingseq_indel, 
                proteinseq=proteinseq, outfile=outf_indel, ids=ids, RPKM=RPKM)
        }
        
        packageStartupMessage(" done")
    }else{
        SNVvcf <- vcf[[1]]
    }
    
    message("Output variation table and variant protein sequence caused by SNVs... ", 
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
    
    outf_mtab <- paste(outfile_path, '/', outfile_name, '_snv.tab', sep='')
    write.table(mtab, file=outf_mtab, sep='\t', quote=F, row.names=F)
    
    outf_snv <- paste(outfile_path, '/', outfile_name, '_snv.fasta', sep='')
    OutputVarproseq(mtab, proteinseq, outf_snv, ids, lablersid=lablersid, 
                    RPKM=RPKM)
    
    packageStartupMessage(" done")
    
    if(nov_junction == TRUE&!is.null(bedFile)&!is.null(genome)){
        message("Output novel junction peptides... ", appendLF=FALSE)
        splicemax <- ''
        load(paste(annotation_path, '/splicemax.RData', sep=''))
        txdb <- loadDb(paste(annotation_path, '/txdb.sqlite', sep=''))
        junction_type <- JunctionType(bedFile, skip=1, covfilter=5, 
                        splicemax, txdb, ids)
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

