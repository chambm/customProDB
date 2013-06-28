##' Load multiple vcf files and output a GRange object with SNVs present in multiple samples. 
##'
##' This function allows to limit SNVs that are present in at least m out of n VCF files.
##' @title Generate shared variation dataset from multiple VCF files
##' @param path a directory of VCF files. read all VCF files in this directory.
##' @param share_num Two options, percentage format or sample number.
##' @param ... additional arguments
##' @return a GRange object that contains the shared variations
##' @author Xiaojing Wang
##' @examples
##' 
##' path <- system.file("extdata/vcfs", package="customProDB")
##' 
##' shared <- Multiple_VCF(path, share_num=2)
##' 
##'

Multiple_VCF <-  function(path, share_num, ...)
    {
        vcfFile<- paste(path, '/', list.files(path, pattern="*vcf$"), sep='')

        vcfs <- lapply(vcfFile, function(x){ 
            InputVcf(x)}
        )
        allsample <-  as.list(unlist(vcfs))

        if(grepl('%', share_num)){
            share_num <- round(as.numeric(gsub('%', '', share_num)) * 
                                length(allsample) / 100) 
        }else share_num <- as.numeric(share_num)
        vcfgls <- GRangesList(allsample)
        total <- unlist(vcfgls)
        if('INDEL' %in% names(values(total))){
            uniquetotal <- unique(total[, c('REF', 'ALT', 'INDEL')])
        }else uniquetotal <- unique(total[, c('REF', 'ALT')])

        ctab <- do.call('cbind', lapply(allsample, function(x) 
                            countOverlaps(uniquetotal, x)))
        sumcount <- apply(ctab, 1, sum)
        index <- which(sumcount >= share_num)
        vcf <- uniquetotal[index]

    }




