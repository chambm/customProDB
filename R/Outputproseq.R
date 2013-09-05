##' Get the FASTA file of proteins that pass RPKM cutoff. the FASTA ID line contains protein ID, gene ID, HGNC symbol and description
##'
##' by taking the RPKM value as input, the function outputs sequences of the proteins that pass the cutoff.
##' @title output FASTA format file contains proteins that have expression level above the cutoff
##' @param rpkm a numeric vector containing RPKM for each protein
##' @param cutoff cutoff of RPKM value. Two options are available, percentage format or RPKM. By default we use "30%" or the RPKM value of 1. "30%" means we keep top 70% proteins according to their RPKMs.
##' @param proteinseq a dataframe containing  protein ids and protein sequences.
##' @param outfile output file name.
##' @param ids a dataframe containing gene/transcript/protein id mapping information.
##' @param ... additional arguments
##' @return FASTA file contains proteins with RPKM above the cutoff.
##' @author Xiaojing Wang
##' @examples
##' 
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
##' bamFile <- system.file("extdata/bams", "test1_sort.bam", 
##'     package="customProDB")
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' RPKM <- calculateRPKM(bamFile, exon, proteincodingonly=TRUE, ids)
##' outf1 <- paste(tempdir(), '/test_rpkm.fasta', sep='')
##' Outputproseq(RPKM, 1, proteinseq, outf1, ids)
##' 
##' 

Outputproseq <- function(rpkm, cutoff="30%", proteinseq, outfile, ids, ...)
    {
        if(grepl('%', cutoff)){
            cutoff <- quantile(rpkm, as.numeric(gsub('%', '', cutoff))/100)
            }else cutoff <- as.numeric(cutoff)
        s<-rpkm[rpkm >= cutoff]
        seqs <- proteinseq[proteinseq[, 'pro_name'] %in% names(s), ]
        v <- s[seqs[, 'pro_name']]
        seqs <- cbind(seqs, v)

        ftab <- merge(ids,seqs, by.x='pro_name', by.y='pro_name', all=FALSE, 
                    stringsAsFactors=FALSE)
        ftab <- ftab[order(ftab[, 'v'], decreasing=TRUE), ]

        tmp <- apply(ftab, 1, function(x) 
                paste('>', x['pro_name'], " |", round(x['v'], 4), "|", 
                x['tx_name.x'], "|", x['gene_name'], "|", x['description'], '\n', 
                unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
        
        write(tmp, file=outfile)

    }
