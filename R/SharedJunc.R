##' Load multiple BED files and output a GRange object with junctions present in multiple samples. 
##'
##' This function allows to limit junctions that are present in at least m out of n BED files.
##' @title Generate shared junctions dataset from multiple BED files
##' @param juns a list of GRanges object which input from multiple VCF files using function InputVcf.
##' @param share_num Junctions must occurs in this number of samples to be consider. Two options, percentage format or sample number.
##' @param ext_up upstream extension from the junction point, Default is 100nt.
##' @param ext_down downstream extension from the junction point, Default is 100nt.
##' @param ... additional arguments
##' @return a GRange object that contains the shared junctions
##' @author Xiaojing Wang
##' @examples
##' 
##' path <- system.file("extdata/beds", package="customProDB")
##' bedFiles<- paste(path, '/', list.files(path, pattern="*bed$"), sep='')
##' juncs <- lapply(bedFiles, function(x) Bed2Range(x, skip=1, covfilter=5))
##' shared <- SharedJunc(juncs, share_num=2, ext_up=100, ext_down=100)
##' shared
##'

SharedJunc <-  function(juns, share_num=2, ext_up=100, ext_down=100, ...)
    {
   
        jungls <- GRangesList(juns)
        
        jungls_basic <- lapply(jungls, function(x){
                c(x, ignore.mcols=TRUE)
                })
        jungls_basic <- GRangesList(jungls_basic)
        total <- unlist(jungls_basic)       
        if(grepl('%', share_num)){
            share_num <- round(as.numeric(gsub('%', '', share_num)) * 
                                length(allsample) / 100) 
        }else share_num <- as.numeric(share_num)
        uniquetotal <- unique(total)

        ctab <- do.call('cbind', lapply(jungls, function(x) 
                            countOverlaps(uniquetotal, x, type='equal')))
        sumcount <- apply(ctab, 1, sum)
        index <- which(sumcount >= share_num)
        sharejunc <- uniquetotal[index]
        
        junRange <- GRanges(seqnames=seqnames(sharejunc), ranges=ranges(sharejunc), 
                strand=strand(sharejunc), id=paste('JUNC', 1:length(index), sep=''), 
                cov=sumcount[index], part1_len=ext_up + 1, part2_len=ext_down + 1, 
                part1_sta=start(sharejunc)-ext_up, part1_end=start(sharejunc), 
                part2_sta=end(sharejunc), part2_end=end(sharejunc)+ext_down)
 
}