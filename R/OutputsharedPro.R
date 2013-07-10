##' Output a FASTA file containing shared proteins with expression above cutoff in multiple samples
##'
##' this function takes RPKM matrix as input, users can set two paramteters,cutoff and shared, to generated a consensus expressed database 
##' @title Output the sequences of proteins with high expressions in multiple samples.
##' @param RPKMs RPKM matrix; row name (protein name) is required.
##' @param cutoff a percentage format cutoff (e.g. '30%'), or a vector with each element as a vlaue cutoff referring to one sample
##' @param share_sample the minimum share sample numbers for proteins which pass the cutoff.
##' @param proteinseq a dataframe containing protein ids and protein sequences
##' @param outfile output file name
##' @param ids a dataframe containing gene/transcript/protein id mapping information.
##' @param ... additional arguments
##' @return a FASTA file containing proteins with RPKM above the cutoff in at least certain number of samples
##' @author Xiaojing Wang
##' @examples
##' 
##' path <- system.file("extdata/bams", package="customProDB")
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' bamFile<- paste(path, '/', list.files(path, pattern="*bam$"), sep='')
##' rpkms <- sapply(bamFile,function(x) 
##'             calculateRPKM(x, exon, proteincodingonly=TRUE, ids))
##' outfile <- paste(tempdir(), '/test_rpkm_share.fasta', sep='')
##' pro <- OutputsharedPro(rpkms, cutoff=1, share_sample=2, proteinseq, 
##'             outfile, ids)
##' 

OutputsharedPro <- function(RPKMs, cutoff="30%",share_sample='50%', proteinseq, 
            outfile, ids,...)
    {

        if(length(cutoff) > 1){
            cutoffs <- cutoff 
        }else{
            cutoffs <- apply(RPKMs, 2, function(x) ifelse(grepl('%', cutoff), 
                    quantile(x, as.numeric(gsub('%', '', cutoff))/100), cutoff))
        }
        #s<-RPKMs[RPKMs>cutoff]
        rpkms_new <- RPKMs
        for(i in 1:dim(RPKMs)[2]) {
            rpkms_new[which(rpkms_new[, i] <= cutoffs[i]), i] <- 0
            }

        test <- apply(rpkms_new, 1, function(x) sum(x > 0))
        if(grepl('%', share_sample)){
            share_sample <- round(as.numeric(gsub('%', '', share_sample))*dim(RPKMs)[2]/100) 
        }else share_sample <- as.numeric(share_sample)
        rpkms_res <- rpkms_new[which(test >= share_sample), ]
        proid <- rownames(rpkms_res)

        seqs <- proteinseq[proteinseq[, 'pro_name'] %in% proid,]
        mean_rpkm <- unlist(lapply(seqs[, 'pro_name'], function(x) 
                        round(mean(RPKMs[x,]),4)))
        seqs <- cbind(seqs, mean_rpkm)

        

        ftab <- merge(ids, seqs, by.x='pro_name', by.y='pro_name', all=F, 
                    stringsAsFactors = FALSE)

    
        tmp <- paste('>', ftab[, 'pro_name'], " |", ftab[, 'mean_rpkm'], "|", 
                ftab[, 'tx_name.x'],"|",ftab[, 'gene_name'], "|", 
                ftab[, 'description'], '\n', 
                unlist(strsplit(ftab[, 'peptide'], '\\*'))[1], sep='')
        write(tmp,file=outfile)

    }




