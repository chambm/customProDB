##' Load multiple vcf files and output a GRange object with SNVs present in multiple samples. 
##'
##' This function allows to limit SNVs that are present in at least m out of n VCF files.
##' @title Generate shared variation dataset from multiple VCF files
##' @param vcfs a list of GRanges object which input from multiple VCF files using function InputVcf.
##' @param share_num Two options, percentage format or sample number.
##' @param ... additional arguments
##' @return a GRange object that contains the shared variations
##' @author Xiaojing Wang
##' @examples
##' 
##' path <- system.file("extdata/vcfs", package="customProDB")
##' vcfFiles<- paste(path, '/', list.files(path, pattern="*vcf$"), sep='')
##' vcfs <- lapply(vcfFiles, function(x) InputVcf(x))
##' shared <- Multiple_VCF(vcfs, share_num=2)
##' 
##'

Multiple_VCF <-  function(vcfs, share_num, ...)
    {

        allsample <-  as.list(unlist(vcfs))
        allsample_basic <- lapply(allsample, function(x){
                if('INDEL' %in% names(values(x))){
                    x[, c('REF', 'ALT', 'INDEL')]
                }else x[, c('REF', 'ALT')]
                })
                
        if(grepl('%', share_num)){
            share_num <- round(as.numeric(gsub('%', '', share_num)) * 
                                length(allsample) / 100) 
        }else share_num <- as.numeric(share_num)
        vcfgls <- GRangesList(allsample_basic)
        total <- unlist(vcfgls)
        if('INDEL' %in% names(values(total))){
            uniquetotal <- unique(total[, c('REF', 'ALT', 'INDEL')])
        }else uniquetotal <- unique(total[, c('REF', 'ALT')])

        ctab <- do.call('cbind', lapply(allsample_basic, function(x) 
                            countOverlaps(uniquetotal, x)))
        sumcount <- apply(ctab, 1, sum)
        index <- which(sumcount >= share_num)
        vcf <- uniquetotal[index]

    }




