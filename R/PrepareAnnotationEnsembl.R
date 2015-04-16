##' prepare the annotation from ENSEMBL through biomaRt.
##'
##' this function automaticlly prepares all  annotation infromation needed in the following analysis.
##' @title prepare annotation from ENSEMBL
##' @param mart which version of ENSEMBL dataset to use. see useMart from package biomaRt for more detail.
##' @param annotation_path specify a folder to store all the annotations
##' @param dbsnp specify a snp dataset you want to use for the SNP annotation, default is NULL.
##' @param transcript_ids optionally, only retrieve transcript annotation data for the specified set of transcript ids
##' @param splice_matrix whether generate a known exon splice matrix from the annotation. this is not necessary if you don't want to analyse junction results, default is FALSE. 
##' @param COSMIC whether to download COSMIC data, default is FALSE.
##' @param ... additional arguments
##' @return several .RData file containing annotations needed for following analysis.
##' @author Xiaojing Wang
##' @examples
##' 
##' ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
##' host="feb2012.archive.ensembl.org", path="/biomart/martservice", 
##' archive=FALSE)
##' 
##' annotation_path <- tempdir()
##' transcript_ids <- c("ENST00000234420", "ENST00000269305", "ENST00000445888", 
##'     "ENST00000257430", "ENST00000457016", "ENST00000288602", 
##'     "ENST00000269571", "ENST00000256078", "ENST00000384871")
##' 
##' PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path, 
##'     splice_matrix=FALSE, dbsnp=NULL, transcript_ids=transcript_ids, 
##'     COSMIC=FALSE)
##' 
##' 



PrepareAnnotationEnsembl <- function(mart, annotation_path, splice_matrix=FALSE, 
                dbsnp=NULL, transcript_ids=NULL, COSMIC=FALSE,
 ...) {
    options(stringsAsFactors=FALSE)
    dataset <- biomaRt:::martDataset(mart)
    biomart <- biomaRt:::martBM(mart)
    host <- strsplit(strsplit(biomaRt:::martHost(mart), ':')[[1]][2], '//')[[1]][2]
    if (!is.null(dbsnp)) {
        session  <- browserSession()
        if(dataset == 'hsapiens_gene_ensembl') {
            if(host == 'may2009.archive.ensembl.org'){
                genome(session) <- 'hg18'
                dbsnps <- 'snp130'
            }else{
                genome(session) <- 'hg19'
                dbsnps <- trackNames(session)[grep('snp', trackNames(session), fixed=T)]
            }
        }
        
        if(dataset == 'mmusculus_gene_ensembl') {
            if(host == 'jan2013.archive.ensembl.org'||host == 'oct2012.archive.ensembl.org'){
                genome(session) <- 'mm10'
                dbsnps <- 'snp137'
            }else{
                genome(session) <- 'mm9'
                dbsnps <- 'snp128'
            }
        }
        
        dbsnp <- pmatch(dbsnp, dbsnps)
        if (is.na(dbsnp)) 
                stop("invalid dbsnp name for specified genome")
        if (dbsnp == -1) 
                stop("ambiguous dbsnp name")
    }
    
    message("Prepare gene/transcript/protein id mapping information (ids.RData) ... ", 
            appendLF=FALSE)
    
    if(is.null(transcript_ids)){ 
        transcript_ids <- getBM(attributes=c("ensembl_transcript_id"), mart=mart)[,1]
        
    }
    attributes.id <- c("ensembl_gene_id", "hgnc_symbol", "description") 
    idstab <- getBM(attributes=attributes.id, mart=mart, 
                filters='ensembl_transcript_id', values=transcript_ids)
    colnames(idstab) <- c("ensembl_gene_id", "hgnc_symbol", "description") 
    
            idssum <- ddply(idstab, .(ensembl_gene_id), function(x) {
             new.x <- x[1, ]
             new.x$hgnc_symbol <- paste(x$hgnc_symbol, collapse=",")
             new.x
            })
    
    
    attributes.tr <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id")
    tr <- getBM(attributes=attributes.tr, mart=mart, 
                filters='ensembl_transcript_id', values=transcript_ids)
    colnames(tr) <- c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id")
    ids <- merge(tr, idssum, by='ensembl_gene_id')
    description <- paste(ids[, 'hgnc_symbol'], ids[, 'description'], sep='|')
    ids <- cbind(ids[, 1:3], description)
    colnames(ids) <- c('gene_name', 'tx_name', 'pro_name', 'description')
    save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    
    packageStartupMessage(" done")
    
    message("Build TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)
    tr_coding <- subset(ids, pro_name != "")
    tr_noncoding <- subset(ids, pro_name == "")
    
    txdb<- makeTranscriptDbFromBiomart_archive(biomart=biomart, dataset=dataset, 
            host=host, path="/biomart/martservice", archive=FALSE, 
            transcript_ids=transcript_ids)
    #saveFeatures(txdb, file=paste(annotation_path,'/txdb.sqlite',sep=''))
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    #txdb_coding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #   dataset=dataset, host=host, path="/biomart/martservice", archive=FALSE, 
    #   transcript_ids=tr_coding[, "tx_name"])
    #saveFeatures(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    #saveDb(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))

    #txdb_noncoding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #    dataset =dataset, host=host, path="/biomart/martservice", 
    #    archive=FALSE, transcript_ids=tr_noncoding[, "tx_name"])
    #saveFeatures(txdb_noncoding, file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))
    #saveDb(txdb_noncoding, file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))

    transGrange <- transcripts(txdb)
    transintxdb <- IRanges::as.data.frame(transGrange)[, c('tx_id', 'tx_name')]
    
    
    message("Prepare exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)
    
    attributes.exon <- c("ensembl_exon_id", "ensembl_peptide_id", 'ensembl_gene_id', 
        "chromosome_name", "start_position", "end_position", "exon_chrom_start", 
        "exon_chrom_end", "strand", "5_utr_start", "5_utr_end", "3_utr_start", 
        "3_utr_end", "cds_start", "cds_end", "rank", "ensembl_transcript_id")
    exon <- getBM(attributes=attributes.exon, mart=mart, 
            filters='ensembl_transcript_id', values=transcript_ids)
    colnames(exon) <- attributes.exon
    exon <- merge(exon, transintxdb, by.x="ensembl_transcript_id", by.y="tx_name")
    
    colnames(exon) <- c("tx_name", "exon_name", "pro_name", "gene_name", 
        "chromosome_name", "start_position", "end_position", "exon_chrom_start", 
        "exon_chrom_end", "strand", "5_utr_start", "5_utr_end", "3_utr_start", 
        "3_utr_end", "cds_start", "cds_end", "rank", "tx_id")
    
    
    cdsByTx <- cdsBy(txdb, "tx", use.names=T)
    cdss <-  IRanges::as.data.frame(cdsByTx)
    cds_chr_p <- data.frame(tx_name=cdss[, "group"], cds_chr_start=cdss[, "start"], 
                    cds_chr_end=cdss[, "end"], rank=cdss[, "exon_rank"])
    
    cds_chr_p_coding <- subset(cds_chr_p, tx_name %in% tr_coding[, 'tx_name'])
    
    exon <- merge(exon, cds_chr_p_coding, by.y=c("tx_name", "rank"), 
            by.x=c("tx_name", "rank"), all.x=T)
        
    save(exon,file=paste(annotation_path, '/exon_anno.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein coding sequence (procodingseq.RData)... ", 
            appendLF=FALSE)
    attributes.codingseq <- c("coding", "ensembl_peptide_id", 
            "ensembl_transcript_id") 
    if(length(tr_coding[, 'pro_name']<10000)){
        coding <- getBM(attributes=attributes.codingseq, 
                    filters="ensembl_peptide_id", values=tr_coding[, 'pro_name'], 
                    mart=mart)
    }else{ 
        index <- floor(length(tr_coding[, 'pro_name'])/10000)
        coding <- c()
        for(i in 1:index) {
            st <- (i-1)*10000+1
            ed <- i*10000
            tmp <- getBM(attributes=attributes.codingseq, filters="ensembl_peptide_id", 
                        values=tr_coding[st:ed, 'pro_name'], mart=mart)
            coding <- rbind(coding, tmp)
            #print(i)
        }
        tmp <- getBM(attributes=attributes.codingseq, filters="ensembl_peptide_id", 
                values=tr_coding[ed+1:length(tr_coding[, 'pro_name']), 'pro_name'], 
                mart=mart)
        coding <- rbind(coding, tmp)
    }
    colnames(coding) <- attributes.codingseq 
    tx_id <- transintxdb[match(coding[, 'ensembl_transcript_id'], transintxdb[, 'tx_name']), 
                'tx_id']
    procodingseq <- cbind(coding, tx_id)
    colnames(procodingseq) <- c("coding", "pro_name", "tx_name", "tx_id")
    save(procodingseq,file=paste(annotation_path, '/procodingseq.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare protein sequence (proseq.RData) ... ", appendLF=FALSE)
    attributes.proseq <- c("peptide", "ensembl_peptide_id", "ensembl_transcript_id") 
    if(length(tr_coding[, 'pro_name']<10000)){
        proteinseq <- getBM(attributes=attributes.proseq, 
            filters="ensembl_peptide_id", values=tr_coding[,'pro_name'], mart=mart)
    }else{ 
        index <- floor(length(tr_coding[, 'pro_name'])/10000)
        proteinseq <- c()
        for(i in 1:index) {
            st <- (i-1)*10000+1
            ed <- i*10000
            tmp <- getBM(attributes=attributes.proseq, filters="ensembl_peptide_id", 
                    values= tr_coding[st:ed, 'pro_name'], mart=mart)
            proteinseq <- rbind(proteinseq, tmp)
            #print(i)
        }
        tmp <- getBM(attributes=attributes.proseq, filters="ensembl_peptide_id", 
            values=tr_coding[ed+1:length(tr_coding[, 'pro_name']), 'pro_name'], 
            mart=mart)
        proteinseq <- rbind(proteinseq, tmp)
    }
    colnames(proteinseq) <- c("peptide", "pro_name", "tx_name")
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    packageStartupMessage(" done")
    
    
    if (!is.null(dbsnp)) {
        
        message("Prepare dbSNP information (dbsnpinCoding.RData) ... ", 
            appendLF=FALSE)
        
        if(length(dbsnps) == 1&&dbsnps == 'snp128'){
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], table='snp128')
        }else{
            dbsnp_query <- ucscTableQuery(session, dbsnps[dbsnp], 
                    table=paste(dbsnps[dbsnp], 'CodingDbSnp', sep=''))
        }
        snpCodingTab <- getTable(dbsnp_query)
        snpCodingTab[, 'chrom'] <- gsub('chr', '', snpCodingTab[, 'chrom'])
        chrlist <- paste(c(seq(1:22),'X','Y'))
        snpCoding <- subset(snpCodingTab, chrom %in% chrlist ,select=c(chrom:name, alleleCount, alleles))
        snpCoding <- unique(snpCoding)
        #snpCoding[, 'chrom'] <- gsub('chrM', 'MT', snpCoding[, 'chrom'])
        #
        
        #save(snpCoding,file=paste(annotation_path,'/snpcoding.RData',sep=''))
        snpCoding <- GRanges(seqnames=snpCoding[, 'chrom'], 
                    ranges=IRanges(start=snpCoding[, 'chromStart'], 
                    end=snpCoding[, 'chromEnd']), strand='*', 
                    rsid=snpCoding[, 'name'], alleleCount=snpCoding[, 'alleleCount'], 
                    alleles=snpCoding[, 'alleles'])
        
        #seqlevels(snpCoding)
        
        #if(TRUE%in% grep('chr',seqlevels(snpCoding)) > 0 ) {
        #    rchar <- sub('chr','',seqlevels(snpCoding))
        #    names(rchar) <- seqlevels(snpCoding)
        #    snpCoding <- renameSeqlevels(snpCoding, rchar) }
        #if('M'%in%seqlevels(snpCoding)) snpCoding <- renameSeqlevels(snpCoding, c( M='MT'))
        #chrlist <- paste(c(seq(1:22),'X','Y'))
        transGrange_snp <- transGrange
        #transGrange_snp <- keepSeqlevels(transGrange_snp, snpCoding)
        #snpCoding <- keepSeqlevels(snpCoding, transGrange_snp)
        
        #snpCoding <- keepSeqlevels(snpCoding, transGrange)
        
        dbsnpinCoding <- subsetByOverlaps(snpCoding,transGrange_snp)
        
        save(dbsnpinCoding,file=paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
        packageStartupMessage(" done")
    
    }
    
    if (COSMIC) {
        message("Prepare COSMIC information (cosmic.RData) ... ", 
                appendLF=FALSE)
    
        cosmic_query <- ucscTableQuery(session, 'cosmic', table='cosmic')
        cosmicTab <- getTable(cosmic_query)
        cosmicTab[,'chrom'] <- gsub('chrM', 'MT', cosmicTab[, 'chrom'])
        cosmicTab[,'chrom'] <- gsub('chr', '', cosmicTab[, 'chrom']) 
        chrlist <- paste(c(seq(1:22),'X','Y','MT')) 
        cosmicTab <- subset(cosmicTab, chrom %in% chrlist)
        cosmic <- GRanges(seqnames=cosmicTab[, 'chrom'], 
                ranges=IRanges(start=cosmicTab[, 'chromStart'], 
                end=cosmicTab[, 'chromEnd']), strand = '*', 
                cosid=cosmicTab[, 'name'])   
        
        transGrange_cosmic <- transGrange 
        #transGrange_cosmic <- keepSeqlevels(transGrange_cosmic, cosmic)        
        #cosmic <- keepSeqlevels(cosmic, transGrange_cosmic)
        cosmic <- subsetByOverlaps(cosmic, transGrange_cosmic)
        
        save(cosmic,file=paste(annotation_path, '/cosmic.RData', sep=''))    
        packageStartupMessage(" done")        
    }
    
    
    if(splice_matrix){
        message("Prepare exon splice information (splicemax.RData) ... ", 
                appendLF=FALSE)
        exonByTx <- exonsBy(txdb, "tx", use.names=F)
        index <- which(elementLengths(exonByTx)==1)
        exonByTx_mul <- exonByTx[-index]
        exons_mul <- IRanges::as.data.frame(exonByTx_mul)
        exonslist <- split(exons_mul, exons_mul$group)
       
        splicemax_list <- lapply(exonslist, .gen_splicmatrix)
        splicemax <- do.call(rbind, splicemax_list)
        save(splicemax, file=paste(annotation_path, '/splicemax.RData', sep=''))
        packageStartupMessage(" done")    
    }
    #splice <- paste(splicemax[, 1], splicemax[, 2], sep='-')
    
}

###########################################################################################
#generate splicing matrix
###########################################################
.gen_splicmatrix <- function(x, 
     ...) {           
         mystrand=x[1,'strand']
         a=x[,'exon_rank']
         b=x[,'exon_id']
         n=length(a)
         if (mystrand=='+'){
            tmp=order(a)
            
        }else{
            tmp=order(a,decreasing=T)
            
        }
        mm=cbind(b[tmp[1:(n-1)]], b[tmp[2:n]])
        mm
    }

