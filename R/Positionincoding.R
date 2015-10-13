##' For those variations labeled with "Coding", positionincoding() function computes the position of variation in the coding sequence of each transcript.
##'
##' this function prepares input data frame for aaVariation().
##' @title Find the position in coding sequence for each variation.
##' @param Vars a GRanges object of variations
##' @param exon a dataframe of exon annotations for protein coding transcripts.
##' @param dbsnp provide a GRanges object of known dbsnp information to include dbsnp evidence into the output table, default is NULL. 
##' @param COSMIC provide a GRanges object of known COSMIC information to include COSMIC evidence into the output table, default is NULL. 
##' @param ... additional arguments
##' @return a data frame containing the position in coding sequence for each variation
##' @author Xiaojing Wang
##' @examples
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' table(values(vcf[[1]])[['INDEL']])
##' index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
##' indelvcf <- vcf[[1]][index]
##' 
##' index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
##' SNVvcf <- vcf[[1]][index]
##' load(system.file("extdata/refseq", "exon_anno.RData", 
##'     package="customProDB"))
##' load(system.file("extdata/refseq", "dbsnpinCoding.RData", 
##'     package="customProDB"))
##' load(system.file("extdata/refseq", "procodingseq.RData", 
##'     package="customProDB"))
##' load(system.file("extdata/refseq", "cosmic.RData", 
##'     package="customProDB"))
##' postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
##' 


Positionincoding<- function(Vars, exon, dbsnp=NULL, COSMIC=NULL,...)
    {
        options(stringsAsFactors=FALSE)

        exon_cds <- subset(exon, !is.na(cds_start))
        exon_cds <- subset(exon_cds, !is.na(cds_end))

    
        cdscombine <- exon_cds
        
        #cdsanno <- GRanges(seqnames=cdscombine$chromosome_name, ranges=IRanges(start=cdscombine$cds_chr_start,end=cdscombine$cds_chr_end), strand=cdscombine$strand,
        #        trans_name=cdscombine$ensembl_transcript_id, gene_name=cdscombine$ensembl_gene_id, pro_name=cdscombine$ensembl_peptide_id, cds_s2coding=cdscombine$cds_start,cds_e2coding=cdscombine$cds_end)

        cdsanno <- GRanges(seqnames=cdscombine$chromosome_name, 
            ranges=IRanges(start=cdscombine$cds_chr_start, 
            end=cdscombine$cds_chr_end), strand=cdscombine$strand, 
            trans_id=cdscombine$tx_id, trans_name=cdscombine$tx_name, 
            pro_name=cdscombine$pro_name, gene_name=cdscombine$gene_name, 
            cds_s2coding=cdscombine$cds_start, cds_e2coding=cdscombine$cds_end)

        #if(TRUE%in% grep('chr',seqlevels(Vars)) > 0 ) {
        #    rchar <- sub('chr','',seqlevels(Vars))
        #    names(rchar) <- seqlevels(Vars)
        #    Vars <- renameSeqlevels(Vars, rchar) }
        #if('M'%in%seqlevels(Vars)) Vars <- renameSeqlevels(Vars, c( M='MT'))

        
        #cdsanno <- keepSeqlevels(cdsanno, Vars)
        #Vars <- keepSeqlevels(Vars, cdsanno)
        candiincds <- subsetByOverlaps(Vars, cdsanno)
        candimtch <- findOverlaps(candiincds, cdsanno)
        queryindex <- queryHits(candimtch)
        hitindex <- subjectHits(candimtch)

        pos <- start(ranges(candiincds))[queryindex]
        chr <- as.character(seqnames(candiincds))[queryindex]
        #consensusQ <- values(candiincds)[["consensusQuality"]][queryindex]
        #snpQ <- values(candiincds)[["snpQuality"]][queryindex]
        #maxmappQ <- values(candiincds)[["maxMappingQuality"]][queryindex]
        #depth <- values(candiincds)[["coverage"]][queryindex]
        refbase <- as.character(values(candiincds)[["REF"]][queryindex])
        varbase <- as.character(values(candiincds)[["ALT"]][queryindex])

        txname <- values(cdsanno)["trans_name"][hitindex, ]
        txid <- values(cdsanno)["trans_id"][hitindex, ]
        proname <- values(cdsanno)["pro_name"][hitindex, ]
        genename <- values(cdsanno)["gene_name"][hitindex, ]
        strand <- as.character(strand(cdsanno))[hitindex]
        cds_s2coding <- values(cdsanno)["cds_s2coding"][hitindex, ]
        cds_s <- start(ranges(cdsanno))[hitindex]
        cds_e <- end(ranges(cdsanno))[hitindex]
        pincoding <- ifelse(strand == '+', cds_s2coding+pos-cds_s, 
                        cds_s2coding+cds_e-pos-nchar(refbase)+1)
        
        
        #res <- data.frame(genename=character(),txname=character(),proname=character(),chr=character(),strand=character(),
        #            pos=integer(),refbase=character(),varbase=character(),pincoding=integer(), consensusQ=integer(),snpQ=integer(),maxmappQ=integer(),depth=integer())
        #res <- data.frame(genename,txname,proname,chr,strand,pos,refbase,varbase,pincoding,consensusQ,snpQ,maxmappQ,depth)
        res <- data.frame(genename=character(), txname=character(), 
                txid=character(), proname=character(), chr=character(), 
                strand=character(), pos=integer(), refbase=character(), 
                varbase=character(), pincoding=integer())
        res <- data.frame(genename, txname,txid, proname, chr, strand, pos, 
                    refbase, varbase, pincoding)
        
        #dbsnp <- subset(dbsnp,transcript %in% unique(txname))
        if(!is.null(dbsnp)){
            dbsnp <- as.data.frame(dbsnp)
            dbsnp_tmp <- subset(dbsnp,select=c(seqnames:end, rsid))
            dbsnp_unique <- unique(dbsnp_tmp)        
            
            chrstart <- paste(dbsnp_unique$seqnames, (dbsnp_unique$start)+1, sep=':')
            dbsnp_id <- data.frame(unique(cbind(name=as.character(dbsnp_unique$rsid), chrstart)))

            chrpos <- paste(res$chr, res$pos,sep=':')
            res_tmp <- cbind(res, chrpos)

            #tmptab <- merge(nonsy,dbsnp_id,by.x='chrpos',by.y='chrstart',all.x=T)
            rsid <- dbsnp_id[match(res_tmp[, 'chrpos'], dbsnp_id[, 'chrstart']), 'name']
            res <- cbind(res, rsid)
        }
        
        if(!is.null(COSMIC)){
            COSMIC <- as.data.frame(COSMIC)
            cos_chrstart <- paste(COSMIC$seqnames, (COSMIC$start)+1, sep=':')
            cos_id <- data.frame(unique(cbind(name=as.character(COSMIC$cosid), cos_chrstart)))

            chrpos <- paste(res$chr, res$pos, sep=':')
            res_tmp <- cbind(res, chrpos)

            #tmptab <- merge(nonsy,dbsnp_id,by.x='chrpos',by.y='chrstart',all.x=T)
            COSMIC_id <- cos_id[match(res_tmp[, 'chrpos'],cos_id[, 'cos_chrstart']), 'name']
            res <- cbind(res, COSMIC_id)
        }
        
        res
        
        index <- which(res[, 'pincoding'] > 0)
        if(length(index)>0)  res <- res[index, ]
        
    }
