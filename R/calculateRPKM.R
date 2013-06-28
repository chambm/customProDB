##' Normalized expression level based on exon read counts. The default output is a vector containing RPKMs for each transcript. vector name is the transcript name.
##' calculate the RPKMs by chromosome. If proteincodingonly=TRUE, vetor name will be set to protein name, and only output RPKMs for the protein coding transcripts. 
##'
##' caculate RPKM from a BAM file based on exon read counts
##' @title Caculate RPKM for each transcripts based on exon read counts.
##' @param bamFile a the input BAM file name.
##' @param exon a dataframe of exon annotations.
##' @param proteincodingonly if TRUE only output RPKMs for protein coding transcripts, the name of output vector will be protein id. if FALSE, output the RPKM for all transcripts.
##' @param ids a dataframe containing gene/transcript/protein id mapping information.
##' @param ... additional arguments
##' @return RPKM value for all transcripts or protein coding transcripts.
##' @author Xiaojing Wang
##' @examples
##' \dontrun{
##' ##test1.bam file is part of the whole bam file. 
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' bamFile <- system.file("extdata/bams", "test1_sort.bam", package="customProDB")
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' RPKM <- calculateRPKM(bamFile,exon,proteincodingonly=TRUE,ids)
##' }
##' 

calculateRPKM <- function(bamFile,exon, proteincodingonly = TRUE,ids=NULL,...)
    {
        #aln <- readGappedAlignments(bamFile)
        #galn <- granges(aln)
        if(proteincodingonly&&is.null(ids))
            stop("must supply ids mapping information if you choose proteincodingonly= TRUE")
        anno <- GRanges(seqnames = exon$chromosome_name,ranges = IRanges(start=exon$exon_chrom_start,end=exon$exon_chrom_end), strand = exon$strand,
             tr_name = exon$tx_name)
        
        
        targets <- scanBamHeader(bamFile)[[1]][["targets"]]
        which <- GRanges(names(targets), IRanges(1, unname(targets)))
        all_tr <- c()
        readbychr <- c()
        
        for (i in seq_along(which)){
            param <- ScanBamParam(which=which[i], what=character())
            aln <- readGAlignments(bamFile, param=param)
            galn <- granges(aln)
            #seqlevels(galn)
        #    if(TRUE%in% grep('chr',seqlevels(galn)) > 0 ) {
        #        rchar <- sub('chr','',seqlevels(galn))
        #        names(rchar) <- seqlevels(galn)
        #        galn <- renameSeqlevels(galn, rchar) }
        #    if('M'%in%seqlevels(galn)) galn <- renameSeqlevels(galn, c( M='MT'))
            
            anno <- keepSeqlevels(anno,seqlevels(galn))
            galn <- keepSeqlevels(galn,seqlevels(anno))
            if(length(galn)>0){
                anno_1 <- anno[seqnames(anno)==seqnames(galn)[1]]
                exon_len <- as.data.frame(cbind(values(anno_1)[,'tr_name'],width(anno_1)))
                exonlenByTrans <- tapply(as.numeric(as.character(exon_len$V2)),exon_len$V1,sum)

                
                readcount <- countOverlaps(anno_1,galn)
                names(readcount) <-  values(anno_1)[,'tr_name']
                countbypro <- tapply(readcount, names(readcount), sum)
                RPK <- countbypro/(exonlenByTrans/1000)
                
                all_tr <- c(all_tr,RPK)    
                readbychr <- c(readbychr,sum(readcount))    
            }
            #print(as.character(seqnames(which[i])))
        }
        
        
        totalReads <- sum(readbychr)
        RPKM <- all_tr/(totalReads/1e+06)
        #full.rpkm <- cbind(RPKM,countbypro,exonlenByTrans)
        #ttt <- as.numeric(RPKM)
        #names(ttt) <- names(RPKM)
        #ttt
        
        options(stringsAsFactors=FALSE)
        
        if(proteincodingonly==TRUE){
            proex <- RPKM
            names(proex) <- ids[match(names(RPKM),ids[,'tx_name']),'pro_name']    
            proex <- proex[-which((names(proex)==''))]
            
            proex
        }else{
            RPKM
        }
    }


