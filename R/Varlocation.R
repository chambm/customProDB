##' For a given GRange object of variations, the Varlocation() function finds the genomic locations for each entry according to the given annotation.
##' Seven labels are used to describe the location (intergenic, intro_nonProcoding,  exon_nonProcoding, intron, 5utr, 3utr and coding).
##' details of the definition can be found in the tutorial.
##'
##' see 'introduction' for more details
##' @title Annotates the variations with genomic location.
##' @param Vars a GRange object of variations
##' @param txdb a TranscriptDb object.
##' @param ids a dataframe containing gene/transcript/protein id mapping information
##' @param ... additional arguments
##' @return a data frame of locations for each variation
##' @author Xiaojing Wang
##' @examples
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##'
##' table(values(vcf[[1]])[['INDEL']])
##' index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
##' indelvcf <- vcf[[1]][index]
##'
##' index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
##' SNVvcf <- vcf[[1]][index]
##'    
##' txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite", package="customProDB"))
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' SNVloc <- Varlocation(SNVvcf,txdb,ids)
##' indelloc <- Varlocation(indelvcf,txdb,ids)
##' table(SNVloc[,'location'])
##' 
##'

Varlocation <- function(Vars, txdb, ids, ...)
    {

        #if(TRUE%in% grep('chr',seqlevels(Vars)) > 0 ) {
        #    rchar <- sub('chr','',seqlevels(Vars))
        #    names(rchar) <- seqlevels(Vars)
        #    Vars <- renameSeqlevels(Vars, rchar) }
        #if('M'%in%seqlevels(Vars)) Vars <- renameSeqlevels(Vars, c( M='MT'))
        tr_coding <- subset(ids, pro_name!="")
        tr_noncoding <- subset(ids, pro_name == "")
     
        trans <- transcripts(txdb)
        tr_coding_id <- values(trans)['tx_id'][, 1][which(values(trans)['tx_name'][, 1] %in% tr_coding[, 'tx_name'])]
        tr_noncoding_id <- values(trans)['tx_id'][, 1][which(values(trans)['tx_name'][, 1] %in% tr_noncoding[, 'tx_name'])]
        
        
        exonsByTx <- exonsBy(txdb, "tx", use.names=F)
        intronsByTx <- intronsByTranscript(txdb, use.names=F)
        cdsByTx <- cdsBy(txdb, "tx", use.names=F)
        fiveutrByTx <- fiveUTRsByTranscript(txdb, use.names=F)
        threeutrByTx <- threeUTRsByTranscript(txdb, use.names=F)
        
        # regions for nocoding trans
        exonsByTx_noncoding <- exonsByTx[which(names(exonsByTx) %in% 
                        tr_noncoding_id)]
        intronsByTx_noncoding <- intronsByTx[which(names(intronsByTx) %in% 
                        tr_noncoding_id)]

        # regions for coding trans
        cdsByTx <- cdsByTx[which(names(cdsByTx) %in% tr_coding_id)]
        intronsByTx <- intronsByTx[which(names(intronsByTx) %in% tr_coding_id)]
        fiveutrByTx <- fiveutrByTx[which(names(fiveutrByTx) %in% tr_coding_id)]
        threeutrByTx <- threeutrByTx[which(names(threeutrByTx) %in% 
                    tr_coding_id)]


        
        

        pos <- start(ranges(Vars))
        chr <- as.character(seqnames(Vars))
        #consensusQ <- values(Vars)[["consensusQuality"]]
        #snpQ <- values(Vars)[["QUAL"]]
        #maxmappQ <- values(Vars)[["maxMappingQuality"]]
        #depth <- values(Vars)[["DP"]]
        refbase <- as.character(values(Vars)[["REF"]])
        varbase <- as.character(values(Vars)[["ALT"]])

        location <- rep("Intergenic", length=length(Vars))
        
        #exonsByTx_noncoding <- keepSeqlevels(exonsByTx_noncoding,Vars)
        unknown <- which(!chr %in% c(seqlevels(exonsByTx), 
                seqlevels(intronsByTx)))
        location[unknown] <- "Unknown"
        
        matchexon_noncoding <- queryHits(findOverlaps(Vars, 
                    exonsByTx_noncoding))
        matchintron_noncoding <- queryHits(findOverlaps(Vars, 
                    intronsByTx_noncoding))

        location[matchexon_noncoding] <- "Exon_nonprocoding"
        location[matchintron_noncoding] <- "Intron_nonprocoding"


        matchcoding <- queryHits(findOverlaps(Vars, cdsByTx))
        match3utr <- queryHits(findOverlaps(Vars, threeutrByTx))
        match5utr <- queryHits(findOverlaps(Vars, fiveutrByTx))
        matchintron <- queryHits(findOverlaps(Vars, intronsByTx))

        location[matchintron] <- "Intron"
        location[matchcoding] <- "Coding"
        location[match3utr] <- "3'UTR"
        location[match5utr] <- "5'UTR"

        #res <- data.frame(chr=character(), pos=integer(), refbase=character(), 
        #           varbase=character(), location=character(), snpQ=integer(), 
        #            depth=integer())
        #res <- data.frame(chr,pos,refbase,varbase,location,snpQ,depth)
        res <- data.frame(chr=character(), pos=integer(), refbase=character(), 
                        varbase=character(), location=character())
        res <- data.frame(chr, pos, refbase, varbase, location)

    }


