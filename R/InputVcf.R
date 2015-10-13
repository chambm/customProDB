##' The InputVcf() function generates a list of GRanges object from a single VCF file. 
##'
##' Read all fields in a VCF file into GRanges object.
##' @title Generate a list of GRanges objects from a VCF file.
##' @param vcfFile a character contains the path and name of a VCF file 
##' @param ... additional arguments
##' @return a list of GRanges object containing a representation of data from the VCF file
##' @author Xiaojing Wang
##' @examples
##' ## multiple samples in one VCF file
##' 
##' vcffile <- system.file("extdata", "test_mul.vcf", package="customProDB")
##' vcfs <- InputVcf(vcffile)
##' length(vcfs)
##' 
##' ## single sample
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' length(vcf)
##' 


InputVcf <- function(vcfFile, ...)
    {
        options(stringsAsFactors=FALSE)
        vcf <- scanVcf(vcfFile)
        #vcf_header <-scanBcfHeader(vcfFile)
        vcf_header <-scanVcfHeader(vcfFile)
        #vcf_header
        #samname <- vcf_header[[1]]$Sample
        samname <- samples(vcf_header)
        samnum <- length(samname)
        
        ### remove position have multiple ALT
        #index_mul <- grep(',',vcf[[1]]$ALT, fixed=T)
        #if(length(index_mul)!=0L){
        #    vcf_new <- lapply(vcf[[1]], function(x) x[-index_mul])
        #}else vcf_new <- vcf[[1]]

        #vcf_unpack <- VariantAnnotation:::.unpackVcf(vcf[[1]], vcf_header)
        #vcf_unpack <- unpackVcf(vcf,vcfFile,info=TRUE, geno=TRUE)
        #unpackVCF <- vcf_unpack[[1]]

        #info <- unpackVCF[['INFO']]
        #info <- vcf_unpack[['INFO']]
        
        info <- vcf[[1]]$INFO
        index_ar <- which(lapply(info, class) == "array")
        if (length(index_ar) != 0)
            for (i in index_ar) {
                coln <- rep(names(info[i]), dim(info[[i]])[3])
                info[[i]] <- matrix(info[[i]],nrow=dim(info[[i]])[1], byrow=F, 
                                dimnames=list(NULL,coln))
            }

        index_li <- which(lapply(info, class) == "list")
        if (length(index_li) != 0)
            for (i in index_li) {

                info[[i]] <- unlist(info[[i]])
            }

        info_fr <- data.frame(info)
        info_df <- DataFrame(info_fr)
        
        #if(length(index_mul)!=0L){
        #    info_df <- info_df[-index_mul,]
        #}else info_df <- info_df


        
        #vcf_unpack <- unpackVcf(vcf,geno=TRUE)
        #geno  <- vcf_unpack[[1]][['GENO']]
        #geno  <- vcf_unpack[['GENO']]
        geno  <- vcf[[1]]$GENO
        
        index_ar <- which(lapply(geno, class) == "array")
        if(length(index_ar) > 0)geno <- geno[-index_ar]
        
        #geno_fr <- lapply(1:samnum, function(x)  data.frame(lapply(geno, function(y) y[,x]),stringsAsFactors =F))
        #names(geno_fr) <- samname
        #geno_df <- lapply(geno_fr, DataFrame)
        
        #if(length(index_mul)!=0L){
        #    geno <- lapply(geno, function(x) as.matrix(x[-index_mul,]))
        #}else geno <- geno

        

        geno_fr <- lapply(1:samnum, function(z){
            tmp <- lapply(geno, function(y)  y[, z])
            index_li <- which(lapply(tmp, mode) == "list")
            
            if (length(index_li) != 0){
                for (i in 1:length(index_li)) {
                    tmp[[index_li[i]]][which(lapply(tmp[[index_li[i]]],length)==0L)] <- NA
                }
            }
                  
            tmp1 <- lapply(tmp, function(x) 
                            if(mode(x)=='list') do.call(rbind, x) else x)
            tmp2 <- data.frame(do.call(cbind, tmp1))
            
            if (length(index_li) != 0){
                coln <- c()
                for (i in 1:length(tmp)) {
                    if(i %in% index_li){
                        coln <- c(coln, rep(names(tmp)[i], dim(tmp1[[i]])[2]))
                    }else  coln <- c(coln, names(tmp)[i])
                }
                
                names(tmp2) <- coln
            }
            tmp2
        } )
        names(geno_fr) <- samname
        geno_df <- lapply(geno_fr, DataFrame)
        
        
        
        ALT_new <- lapply(vcf[[1]]$ALT, function(x) paste(x, collapse=","))
        partA <- DataFrame(REF=as.character(vcf[[1]]$REF), 
                        ALT=as.character(ALT_new), QUAL=vcf[[1]]$QUAL, 
                        FILTER=vcf[[1]]$FILTER)
        #partA <- DataFrame(REF=as.character(vcf_new$REF), 
        #                ALT=as.character(vcf_new$ALT), QUAL=vcf_new$QUAL, 
        #                FILTER=vcf_new$FILTER)
        vcf_granges <- vcf[[1]]$rowRanges
        partAll <- lapply(geno_df, function(x) 
                    cbind(values(vcf_granges), DataFrame(partA, info_df, x) ))
        vcfs <- lapply(partAll, function(x) GRanges(seqnames=seqnames(vcf_granges),
                            ranges=ranges(vcf_granges), strand='*', x) )
        #vcfsGR <- GRangesList(vcfs)

        #names(vcfs) <- samname
        vcfs

    }
