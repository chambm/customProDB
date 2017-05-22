#' prepare the annotation for Refseq through UCSC table browser.
#'
#' @title prepare annotation for Refseq
#' @param genome specify the UCSC DB identifier (e.g. "hg19")
#' @param CDSfasta path to the fasta file of coding sequence.
#' @param pepfasta path to the fasta file of protein sequence, check 'introduction' for more detail.
#' @param annotation_path specify a folder to store all the annotations.
#' @param dbsnp specify a snp dataset to be used for the SNP annotation, default is NULL. (e.g. "snp135")
#' @param transcript_ids optionally, only retrieve transcript annotation data for the specified set of transcript ids. Default is NULL.
#' @param splice_matrix whether generate a known exon splice matrix from the annotation. this is not necessary if you don't want to analyse junction results, default is FALSE.
#' @param COSMIC whether to download COSMIC data, default is FALSE.
#' @param local_cache_path if non-NULL, refers to a directory where previously downloaded resources
#' (like protein coding sequences and COSMIC data) are cached so that the function can be re-run without
#' needing to download identical data again
#' @param ... additional arguments
#' @return several .RData file containing annotations needed for further analysis.
#' @author Xiaojing Wang
#' @importFrom rtracklayer browserSession ucscTableQuery tableNames getTable trackNames ucscSchema genome<-
#' @export
#' @examples
#'
#' transcript_ids <- c("NM_001126112", "NM_033360", "NR_073499", "NM_004448",
#'                     "NM_000179", "NR_029605", "NM_004333", "NM_001127511")
#' pepfasta <- system.file("extdata", "hg19/hg19_protein.fasta", package="customProDB")
#' CDSfasta <- system.file("extdata", "hg19/hg19_coding.fasta", package="customProDB")
#' cache_path <- system.file("extdata", "cache", package="customProDB")
#' annotation_path <- tempdir()
#' PrepareAnnotationRefseq(genome='hg19', CDSfasta, pepfasta, annotation_path,
#'             dbsnp=NULL, transcript_ids=transcript_ids,
#'             splice_matrix=FALSE, COSMIC=FALSE, local_cache_path=cache_path)
#'
#'\dontrun{
#' dbkey = "hg38"
#' tempLocalCache = tempdir()
#' refseqTrack = ifelse(dbkey=="hg38", "refSeqComposite", "refGene")
#' codingFastaFilepath = paste0(tempLocalCache, "/", dbkey, ".cds.fa")
#' proteinFastaFilepath = paste0(tempLocalCache, "/", dbkey, ".protein.fa")
#' 
#' options(timeout=3600)
#' if (!file.exists(codingFastaFilepath)) {
#'   cat(paste("Downloading coding FASTA from:", ucscTableCodingFastaURL, "\n"))
#'   download.file(getCodingFastaUrlFromUCSC(dbkey), codingFastaFilepath, quiet=T, mode='wb')
#' }
#' 
#' if (!file.exists(proteinFastaFilepath)) {
#'   cat(paste("Downloading protein FASTA from:", ucscTableProteinFastaURL, "\n"))
#'   download.file(getProteinFastaUrlFromUCSC(dbkey), proteinFastaFilepath, quiet=T, mode='wb')
#' }
#'
#' cat(paste("Preparing Refseq annotation files\n"))
#' customProDB::PrepareAnnotationRefseq(dbkey, codingFastaFilepath, proteinFastaFilepath,
#'                                      annotation_path=".",
#'                                      dbsnp="snp146", COSMIC=FALSE,
#'                                      local_cache_path=tempLocalCache)
#'                                      
#'}


PrepareAnnotationRefseq <- function(genome='hg19', CDSfasta, pepfasta, 
                annotation_path, dbsnp=NULL, transcript_ids=NULL, 
                splice_matrix=FALSE, COSMIC=FALSE, local_cache_path=NULL, ...) {
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
  
    message("Creating UCSC browser session ... ", appendLF=TRUE)
    session  <- browserSession()
    genome(session) <- genome
    tablename <- 'refGene'
    
    if (!dir.exists(annotation_path) && !dir.create(annotation_path, recursive=TRUE)) {
      stop("error creating annotation_path: ", annotation_path)
    }
    
    if (!is.null(local_cache_path)) {
      local_cache_path = paste0(local_cache_path, "/", genome)
    }

    message("Building TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)
    txdb <- read_or_update_local_cacheDb(makeTxDbFromUCSC(genome=genome, tablename=tablename,
                                                        transcript_ids=transcript_ids),
                                         local_cache_path, "txdb")
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    
    message("Prepare gene/transcript/protein id mapping information (ids.RData) ... ", 
            appendLF=FALSE)

    # UCSC has a different RefSeq structure for hg38
    refgeneTrack = ifelse(genome=="hg38", "NCBI RefSeq", "refGene")
    refGene <- read_or_update_local_cache(getTable(ucscTableQuery(session, refgeneTrack, table="refGene",
                                                                  names=transcript_ids)),
                                          local_cache_path, "refGene")
    stopifnot(nrow(refGene) > 0)
    
    # refLink is not supported with hg38, but it's not genome specific, so always request it from hg19
    if (genome == "hg38") genome(session) <- "hg19"
    reflink <- read_or_update_local_cache(getTable(ucscTableQuery(session, "refGene", table="hgFixed.refLink",
                                                                  names=refGene$name2)),
                                          local_cache_path, "reflink")
    stopifnot(nrow(reflink) > 0)
    if (genome == "hg38") genome(session) <- genome
    
    ids <- subset(reflink, mrnaAcc %in% refGene$name, select = name:protAcc)
    stopifnot(nrow(ids) > 0)
    colnames(ids) <- c('gene_name', 'description', 'tx_name', 'pro_name')
    save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    packageStartupMessage(" done")
    
    #tr_coding <- subset(ids,pro_name!="")
    #tr_noncoding <- subset(ids,pro_name == "")
    #txdb_coding <- makeTranscriptDbFromUCSC(genome=genome, tablename=tablename, 
    #                transcript_ids=tr_coding$tx_name )
    #saveDb(txdb_coding, 
    #       file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    
    #txdb_noncoding <- makeTranscriptDbFromUCSC(genome=genome, 
    #        tablename=tablename, transcript_ids=tr_noncoding$tx_name )
    #saveDb(txdb_noncoding, 
    #        file=paste(annotation_path, '/txdb_noncoding.sqlite', sep=''))
    
    message("Preparing exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)
    
    transGrange <- transcripts(txdb)
    tr <- IRanges::as.data.frame(transGrange)
    stopifnot(nrow(tr) > 0)
    
    cdsByTx <- cdsBy(txdb, "tx", use.names=F)
    exonByTx <- exonsBy(txdb, "tx", use.names=F)
    fiveutrByTx <- fiveUTRsByTranscript(txdb, use.names=F)
    threeutrByTx <- threeUTRsByTranscript(txdb, use.names=F)
    
    #####################################################
    # get error in most recent version of IRanges package
    # previous: rownames after unlist 1.1 1.2 1.3
    # now: .1 .2 .3 ... were removed, so there are some rows with same rownames
    ######################################################
    #cdss <-  IRanges::as.data.frame(IRanges::unlist(cdsByTx))     
    #exons <- IRanges::as.data.frame(IRanges::unlist(exonByTx))
    #fiveutrs <- IRanges::as.data.frame(IRanges::unlist(fiveutrByTx))
    #threeutrs <- IRanges::as.data.frame(IRanges::unlist(threeutrByTx))

    cdss <-  IRanges::as.data.frame(cdsByTx)
    exons <- IRanges::as.data.frame(exonByTx)
    fiveutrs <- IRanges::as.data.frame(fiveutrByTx)
    threeutrs <- IRanges::as.data.frame(threeutrByTx)
    
    #txid <- matrix(unlist(strsplit(rownames(exons), '\\.')), ncol = 2, byrow =T)[, 1]
    #txid <- gsub('=','\\.', txid)
    exon_p <- data.frame(txid=exons$group_name, chr=exons$seqnames, 
                exon_s=exons$start, exon_e=exons$end, 
                exon_rank=exons$exon_rank)
    exon2tr <-  merge(exon_p, tr,by.y="tx_id", by.x="txid")
    exon2tr <- exon2tr[, -which(names(exon2tr) %in% c("seqnames", "width"))]

    #txid <- matrix(unlist(strsplit(rownames(cdss), '\\.')), ncol = 2, 
    #       byrow =T)[, 1]
    #txid <- gsub('=','\\.',txid)
    cds_p <- data.frame(txid=cdss$group_name, cds_s=cdss$start, 
                cds_e=cdss$end, exon_rank=cdss$exon_rank, 
                width=cdss$width)
    ttt <- split(cds_p, cds_p$txid)
    
    cds_p_new <- as.data.frame(data.table::rbindlist(lapply(ttt, function(x){
        #len <- x$cds_e-x$cds_s+1
        #cum <- cumsum(len)
        cum <- cumsum(x$width)
        rdis <- cbind(c(1, cum[1:length(cum)-1]+1), cum)
        colnames(rdis) <- c('cds_start', 'cds_end')
        tmp <- cbind(x, rdis)
        tmp
        })))
    

    #cds_p_new <- do.call(rbind, cds_p_new_list)
    cds_p_new <- cds_p_new[, -which(names(cds_p_new) %in% c("width"))]
    
    #for(i in 1:length(ttt)) {
    #    print(i)
    #    ttt1 <- ttt[[i]]
    #    len <- ttt1$cds_e-ttt1$cds_s+1
    #    cum <- cumsum(len)
    #    rdis <- cbind(c(1,cum[1:length(cum)-1]+1),cum)
    #    colnames(rdis) <- c('cds_start','cds_end')
    #    tmp <- cbind(ttt1,rdis)
    #    cds_p_new <- rbind(cds_p_new,tmp)
    #}

    cds2exon <- merge(exon2tr, cds_p_new, by.x=c("txid", "exon_rank"), 
                        by.y=c("txid", "exon_rank"), all.x = TRUE)
    #txid <- matrix(unlist(strsplit(rownames(fiveutrs), '\\.')), ncol = 2, 
    #               byrow =T)[, 1]
    #txid <- gsub('=','\\.', txid)
    if(dim(fiveutrs)[1] > 0){
        fiveutr_p <- data.frame(txid=fiveutrs$group_name, 
                                fiveutr_s=fiveutrs$start, 
                                fiveutr_e=fiveutrs$end, 
                                exon_rank=fiveutrs$exon_rank)
        fiveutr2exon <- merge(cds2exon, fiveutr_p, by.x=c("txid", "exon_rank"), 
                    by.y =c("txid", "exon_rank"), all.x = TRUE)
                 
        #txid <- matrix(unlist(strsplit(rownames(threeutrs),'\\.')), 
        #s    ncol = 2,byrow =T)[, 1]
        #txid <- gsub('=','\\.', txid)
        threeutr_p <- data.frame(txid=threeutrs$group_name, 
                                threeutr_s=threeutrs$start, 
                                threeutr_e=threeutrs$end, 
                                exon_rank=threeutrs$exon_rank)
        threeutr2exon <- merge(fiveutr2exon, threeutr_p, 
                    by.x=c("txid", "exon_rank"),
                    by.y=c("txid", "exon_rank"), all.x = TRUE)
        
        #exon <- merge(threeutr2exon,ids,by.x=c("tx_name"),by.y="tx_name",all.x = TRUE)
        exon <- threeutr2exon[order(threeutr2exon$txid, threeutr2exon$exon_rank), ]
        
        colnames(exon) <- c("tx_id","rank", "chromosome_name", 
            "exon_chrom_start", "exon_chrom_end", "start_position", 
            "end_position", "strand", "tx_name", "cds_chr_start", 
            "cds_chr_end", "cds_start", "cds_end", "5_utr_start", 
            "5_utr_end", "3_utr_start", "3_utr_end")
    }else{
        exon <- cds2exon
        colnames(exon) <- c("tx_id","rank", "chromosome_name", 
            "exon_chrom_start", "exon_chrom_end", "start_position", 
            "end_position", "strand", "tx_name", "cds_chr_start", 
            "cds_chr_end", "cds_start", "cds_end")
    }
    pro_name <- ids[match(exon$tx_name, ids$tx_name), 'pro_name']
    gene_name <- ids[match(exon$tx_name, ids$tx_name), 'gene_name']
    exon <- cbind(exon, pro_name, gene_name)
    stopifnot(nrow(exon) > 0)
    
    save(exon, file=paste(annotation_path, '/exon_anno.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Preparing protein sequence (proseq.RData) ... ", appendLF=FALSE)
    pro_seqs <- readAAStringSet(pepfasta, format= 'fasta')
    pro_name_v <- names(pro_seqs)
    pro_name <- unlist(lapply(pro_name_v, function(x) strsplit(x, '\\.')[[1]][1]))
    tx_name <- ids[match(pro_name, ids$pro_name), 'tx_name']
    proteinseq <- as.data.frame(pro_seqs)
    proteinseq <- cbind(proteinseq, pro_name_v, pro_name, tx_name)
    colnames(proteinseq) <- c("peptide", "pro_name_v", "pro_name", "tx_name")
    proteinseq <- subset(proteinseq, tx_name %in% refGene$name)
    stopifnot(nrow(proteinseq) > 0)
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    packageStartupMessage(" done")
    
    message("Preparing protein coding sequence (procodingseq.RData)... ", 
                appendLF=FALSE)
    cds_seqs <- readDNAStringSet(CDSfasta, format= 'fasta')
    tx_name_tmp <- unlist(lapply(names(cds_seqs), function(x) strsplit(x, ' ')[[1]][1]))
    tx_name <- unlist(lapply(tx_name_tmp, function(x) paste(strsplit(x, '_')[[1]][3:4], collapse='_')))
    tx_range_tmp <- unlist(lapply(names(cds_seqs), function(x) strsplit(x, ' ')[[1]][2]))
    tx_chr <- unlist(lapply(tx_range_tmp, function(x) substring(strsplit(x, ':')[[1]][1],7)))
    tx_cds_sta <- unlist(lapply(tx_range_tmp, function(x) strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][[1]]))
    tx_cds_end <- unlist(lapply(tx_range_tmp, function(x) strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][[2]]))
    
    tx_name_cds <- refGene[match(paste(tx_name, tx_chr, tx_cds_sta, tx_cds_end, sep=' '), 
                    paste(refGene$name, refGene$chrom, 
                    refGene$cdsStart+1, refGene$cdsEnd, sep=' ')),
                    c('name','chrom','txStart','txEnd')]
    
    tx_id <- tr[match(paste(tx_name_cds$name, tx_name_cds$chrom, 
            tx_name_cds$txStart+1, tx_name_cds$txEnd, sep=' '), 
            paste(tr$tx_name, tr$seqnames, tr$start, 
            tr$end, sep=' ')), 'tx_id']
    
    pro_name <- ids[match(tx_name,ids$tx_name), 'pro_name']
    procodingseq <- as.data.frame(cds_seqs)
    procodingseq <- cbind(procodingseq, names(cds_seqs), pro_name, tx_name, tx_id)
    colnames(procodingseq) <- c("coding", "tx_name_full", "pro_name", "tx_name", "tx_id")
    procodingseq <- subset(procodingseq, tx_name %in% refGene$name)
    stopifnot(nrow(procodingseq) > 0)
    save(procodingseq, file=paste(annotation_path, '/procodingseq.RData', sep=''))
    
    packageStartupMessage(" done")
    
    if (!is.null(dbsnp) && nzchar(dbsnp)) {
      
        if (!is.null(local_cache_path))
            dbsnp_cache_path = paste0(local_cache_path, "/", dbsnp)
        else
            dbsnp_cache_path = NULL
        
        message("Preparing dbSNP information (dbsnpinCoding.RData) ... ", appendLF=FALSE)
        dbsnps <- trackNames(session)[grep('snp', trackNames(session), fixed=T)]
        dbsnpName <- dbsnp
        dbsnp <- pmatch(dbsnp, dbsnps)
        if (is.na(dbsnp)) 
            stop("invalid dbsnp name for specified genome")
        if (dbsnp == -1) 
            stop("ambiguous dbsnp name")
        
        # for a small set of transcripts, query UCSC's MySQL table directly
        if (!is.null(transcript_ids) && length(transcript_ids) < 1000) {
            getSnpTable = function(genome, transcripts) {
                snpTable = paste0(dbsnps[dbsnp], 'CodingDbSnp')
                transcriptSet = paste0('("', paste(transcript_ids, collapse='", "'), '")', collapse="")
                sql = qq('SELECT snp.chrom, chromStart, chromEnd, snp.name, snp.transcript, alleleCount, alleles
    FROM @{snpTable} snp
    JOIN (SELECT chrom, txStart, txEnd FROM refGene WHERE name IN @{transcriptSet}) txInfo ON snp.chrom=txInfo.chrom
    AND snp.chromStart BETWEEN txInfo.txStart AND txInfo.txEnd
    WHERE snp.transcript IN @{transcriptSet}')
                ucscDb = RMySQL::dbConnect(RMySQL::MySQL(), host="genome-mysql.soe.ucsc.edu", user="genome", dbname=genome)
                suppressWarnings(result = DBI::dbGetQuery(ucscDb, sql)) # Unsigned INTEGER in col 1 imported as numeric
                RMySQL::dbDisconnect(ucscDb)
                return (result)
            }
          
            snpCodingTab = read_or_update_local_cache(getSnpTable(genome, transcript_ids), dbsnp_cache_path, "snpCodingTab")
        } else {
            # for the full dataset, download from their FTP
            
            dbSnpFile = qq("@{dbsnpName}CodingDbSnp.txt.gz")
            dbSnpURL = qq("http://hgdownload.soe.ucsc.edu/goldenPath/@{genome}/database/@{dbSnpFile}")
            if (!is.null(local_cache_path)) {
                dbSnpFile = file.path(local_cache_path, dbSnpFile)
                if (!file.exists(dbSnpFile))
                    download.file(dbSnpURL, dbSnpFile, quiet=T, mode='wb')
            } else
                download.file(dbSnpURL, dbSnpFile, quiet=T, mode='wb')
            snpCodingTab = .temp_unzip(dbSnpFile, data.table::fread, showProgress=FALSE,
                                       select=c(2:6, 8, 10),
                                       col.names=c("chrom", "chromStart", "chromEnd", "name",
                                                   "transcript", "alleleCount","alleles"))
        }
        snpCoding <- subset(snpCodingTab, transcript %in% refGene$name, 
                        select=c(chrom:name, alleleCount, alleles))
        snpCoding <- unique(snpCoding)
        #save(snpCoding,file=paste(annotation_path, '/snpcoding.RData', sep=''))
        dbsnpinCoding <- GRanges(seqnames=snpCoding$chrom, 
            ranges=IRanges(start=snpCoding$chromStart, 
            end=snpCoding$chromEnd), strand='*', 
            rsid=snpCoding$name, alleleCount=snpCoding$alleleCount, 
            alleles=snpCoding$alleles)    
        save(dbsnpinCoding, file=paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
        packageStartupMessage(" done")
    }
    
    if (COSMIC) {
        #cosmic <- trackNames(session)[grep('cosmic',trackNames(session), fixed=T)]
        message("Preparing COSMIC information (cosmic.RData) ... ", appendLF=FALSE)
        
        getCOSMIC = function() {
          cosmicTab = getTable(ucscTableQuery(session, 'cosmic', table='cosmic'))
          cosmic <- GRanges(seqnames=cosmicTab$chrom,
                            ranges=IRanges(start=cosmicTab$chromStart,
                                           end=cosmicTab$chromEnd), 
                            strand = '*',
                            cosid=cosmicTab$name)    
          
          #cosmic <- keepSeqlevels(cosmic,transGrange)
          cosmic <- subsetByOverlaps(cosmic, transGrange)
          cosmic
        }
        cosmic <- read_or_update_local_cache(getCOSMIC(), local_cache_path, "cosmic")
        save(cosmic, file=paste(annotation_path, '/cosmic.RData', sep=''))
        packageStartupMessage(" done")        
    }
    
    if(splice_matrix){
        message("Preparing exon splice information (splicemax.RData) ... ", appendLF=FALSE)
        index <- which(elementNROWS(exonByTx)==1)
        exonByTx_mul <- exonByTx[-index]
        exons_mul <- IRanges::as.data.frame(exonByTx_mul)
        exonslist <- split(exons_mul, exons_mul$group_name)
        #system.time( exonByTx <- exonsBy(txdb,"tx", use.names=F))
        splicemax_list <- lapply(exonslist, .gen_splicmatrix)
        splicemax <- do.call(rbind, splicemax_list)
        save(splicemax, file=paste(annotation_path, '/splicemax.RData', sep=''))
        packageStartupMessage(" done")        
    }
                
}
    
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

