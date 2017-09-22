#' prepare the annotation from ENSEMBL through biomaRt.
#'
#' this function automaticlly prepares all  annotation infromation needed in the following analysis.
#' @title prepare annotation from ENSEMBL
#' @param mart which version of ENSEMBL dataset to use. see useMart from package biomaRt for more detail.
#' @param annotation_path specify a folder to store all the annotations
#' @param dbsnp specify a snp dataset you want to use for the SNP annotation, default is NULL.
#' @param transcript_ids optionally, only retrieve transcript annotation data for the specified set of transcript ids
#' @param splice_matrix whether generate a known exon splice matrix from the annotation; not necessary if you don't want to analyse junction results, default is FALSE.
#' @param COSMIC whether to download COSMIC data, default is FALSE.
#' @param local_cache_path if non-NULL, refers to a directory where previously downloaded resources
#' (like protein coding sequences and COSMIC data) are cached so that the function can be re-run without
#' needing to download identical data again
#' @param ensembl_to_UCSC_genome_map a named list of named lists used to look up the UCSC dbkey for a given biomart; only used for downloading dbSNPs;
#' if DEFAULT_ENSEMBL_UCSC_GENOME_MAP does not contain an up-to-date mapping, pass a new mapping like
#' \code{list("<species>_gene_ensembl" = list("<month>.archive.ensembl.org" = "<ucsc_dbkey>"))}
#' @param ... additional arguments, currently unused
#' @return several .RData file containing annotations needed for following analysis.
#' @author Xiaojing Wang
#' @importFrom rtracklayer browserSession ucscTableQuery tableNames getTable trackNames ucscSchema genome<-
#' @importFrom plyr ddply .
#' @import biomaRt
#' @export
#' @examples
#'
#' ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
#'                             dataset="hsapiens_gene_ensembl",
#'                             host="sep2015.archive.ensembl.org")
#'
#' cache_path <- system.file("extdata", "cache", package="customProDB")
#' annotation_path <- tempdir()
#' transcript_ids <- c("ENST00000234420", "ENST00000269305", "ENST00000445888",
#'                     "ENST00000257430", "ENST00000508376", "ENST00000288602",
#'                     "ENST00000269571", "ENST00000256078", "ENST00000384871")
#'
#' PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path,
#'     splice_matrix=FALSE, dbsnp=NULL, transcript_ids=transcript_ids,
#'     COSMIC=FALSE, local_cache_path=cache_path)
#'
#'\dontrun{
#' # full annotation tests
#' 
#' test_datasets = c("hsapiens", "mmusculus", "cfamiliaris", "scerevisiae", "ggorilla")
#' test_releases = c("mar2017", "may2009", "may2012", "may2017")
#' for (d in 1:length(test_datasets))
#'   for (r in 1:length(test_releases)) {
#'     dataset = paste0(test_datasets[d], "_gene_ensembl")
#'     host = paste0(test_releases[r], ".archive.ensembl.org")
#'     mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset, host)
#'     PrepareAnnotationEnsembl(mart, annotation_path, splice_matrix=FALSE, dbsnp=NULL, COSMIC=FALSE,
#'                              local_cache_path=file.path(annotation_path, "cache"))
#'   }
#' 
#'}
#'



PrepareAnnotationEnsembl <- function(mart, annotation_path, splice_matrix=FALSE, 
                                     dbsnp=NULL, transcript_ids=NULL, COSMIC=FALSE, local_cache_path=NULL,
                                     ensembl_to_UCSC_genome_map = DEFAULT_ENSEMBL_UCSC_GENOME_MAP, 
                                     dbsnp_and_cosmic_only=FALSE, ...) {
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
  
    dataset <- mart@dataset
    biomart <- mart@biomart
    version <- sub("Ensembl (?:Genes )?(\\d+)", "\\1", listEnsembl(mart)[listEnsembl(mart)["biomart"]=="ensembl", 2])
    host <- strsplit(strsplit(mart@host, ':')[[1]][2], '//')[[1]][2]
  
    if (!dir.exists(annotation_path) && !dir.create(annotation_path, recursive=TRUE)) {
      stop("error creating annotation_path: ", annotation_path)
    }
    
    if (!is.null(local_cache_path)) {
      local_cache_path = qq("@{local_cache_path}/@{dataset}_@{version}")
    }
  
    if (!is.null(dbsnp) && nzchar(dbsnp)) {
        genome = ensembl_to_UCSC_genome_map[[dataset]][[host]]
        if (is.null(genome))
            stop(qq("no UCSC dbkey in ensembl_to_UCSC_genome_map for genome @{dataset} and host @{host}"))
      
        if (!is.null(local_cache_path))
            dbsnp_cache_path = paste0(local_cache_path, "/", dbsnp)
        else
            dbsnp_cache_path = NULL
    }
    
    message("Prepare gene/transcript/protein id mapping information (ids.RData) ... ", appendLF=FALSE)
    
    original_transcript_ids = transcript_ids
    if(is.null(transcript_ids)){ 
        transcript_ids <- read_or_update_local_cache(getBM(attributes=c("ensembl_transcript_id"), mart=mart)[,1],
                                                     local_cache_path, "transcript_ids")
    }
    
    external_gene_name = ifelse(version < 77, "external_gene_id", "external_gene_name")
    attributes.id <- c("ensembl_gene_id", external_gene_name, "description") 
    idstab <- read_or_update_local_cache(getBM(attributes=attributes.id, mart=mart, 
                                               filters='ensembl_transcript_id', values=transcript_ids),
                                         local_cache_path, "idstab")
    stopifnot(nrow(idstab) > 0)
    colnames(idstab) <- c("ensembl_gene_id", external_gene_name, "description") 
    
            idssum <- ddply(idstab, .(ensembl_gene_id), function(x) {
             new.x <- x[1, ]
             new.x$external_gene_name <- paste(x$external_gene_name, collapse=",")
             new.x
            })
    
    
    attributes.tr <- c("ensembl_gene_id", "ensembl_transcript_id",
        "ensembl_peptide_id")
    tr <- read_or_update_local_cache(getBM(attributes=attributes.tr, mart=mart,
                                           filters='ensembl_transcript_id', values=transcript_ids),
                                     local_cache_path, "tr")
    stopifnot(nrow(tr) > 0)
    colnames(tr) <- c("ensembl_gene_id", "ensembl_transcript_id", 
        "ensembl_peptide_id")
    ids <- merge(tr, idssum, by='ensembl_gene_id')
    description <- paste(ids$external_gene_name, ids$description, sep='|')
    ids <- cbind(ids[, 1:3], description)
    colnames(ids) <- c('gene_name', 'tx_name', 'pro_name', 'description')
    if (!dbsnp_and_cosmic_only) {
        save(ids, file=paste(annotation_path, '/ids.RData', sep=''))
    }
    packageStartupMessage(" done")
    
    message("Build TranscriptDB object (txdb.sqlite) ... ", appendLF=TRUE)
    tr_coding <- subset(ids, pro_name != "")
    tr_noncoding <- subset(ids, pro_name == "")
    
    txdb<- read_or_update_local_cacheDb(makeTxDbFromBiomart(biomart=biomart, dataset=dataset, host=host,  
                                                            transcript_ids=original_transcript_ids),
            local_cache_path, "txdb")
    #saveFeatures(txdb, file=paste(annotation_path,'/txdb.sqlite',sep=''))
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    #txdb_coding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #   dataset=dataset, host=host, path="/biomart/martservice", archive=FALSE, 
    #   transcript_ids=tr_coding$tx_name)
    #saveFeatures(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))
    #saveDb(txdb_coding, file=paste(annotation_path, '/txdb_coding.sqlite', sep=''))

    #txdb_noncoding <- makeTranscriptDbFromBiomart_archive(biomart=biomart, 
    #    dataset =dataset, host=host, path="/biomart/martservice", 
    #    archive=FALSE, transcript_ids=tr_noncoding$tx_name)
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
    exon <- read_or_update_local_cache(getBM(attributes=attributes.exon, mart=mart, 
                                             filters='ensembl_transcript_id', values=transcript_ids),
                                       local_cache_path, "exon")
    stopifnot(nrow(exon) > 0)
    colnames(exon) <- attributes.exon
    exon <- merge(exon, transintxdb, by.x="ensembl_transcript_id", by.y="tx_name")
    
    colnames(exon) <- c("tx_name", "exon_name", "pro_name", "gene_name", 
        "chromosome_name", "start_position", "end_position", "exon_chrom_start", 
        "exon_chrom_end", "strand", "5_utr_start", "5_utr_end", "3_utr_start", 
        "3_utr_end", "cds_start", "cds_end", "rank", "tx_id")
    
    
    cdsByTx <- cdsBy(txdb, "tx", use.names=FALSE)
    cdss <-  IRanges::as.data.frame(cdsByTx)
    cds_chr_p <- data.frame(tx_id=cdss$group_name,
                    cds_chr_start=cdss$start, 
                    cds_chr_end=cdss$end, rank=cdss$exon_rank)
    
    
    cds_chr_p_coding <- subset(cds_chr_p, tx_id %in% exon[which(exon$pro_name != ''), 'tx_id'])
    
    exon <- merge(exon, cds_chr_p_coding, by.y=c("tx_id", "rank"), 
            by.x=c("tx_id", "rank"), all.x=T)
    ###Ensembl use 1 & -1, change it to +/-
	exon$strand <- unlist(lapply(exon$strand, function(x) 
			ifelse(x=='1', '+', '-')))
    if (!dbsnp_and_cosmic_only) {
        save(exon,file=paste(annotation_path, '/exon_anno.RData', sep=''))
    }
    packageStartupMessage(" done")
    
    message("Prepare protein coding sequence (procodingseq.RData)... ", appendLF=FALSE)
    attributes.codingseq <- c("coding", "ensembl_peptide_id", 
            "ensembl_transcript_id") 

    getCoding = function() {
      if(length(tr_coding$pro_name)<1000){
        getBM(attributes=attributes.codingseq, filters="ensembl_peptide_id", 
              values=tr_coding$pro_name, mart=mart)
      }else{ 
        index <- floor(length(tr_coding$pro_name)/1000)
        coding <- c()
        ptm <- proc.time()
        for(i in 1:index) {
          st <- (i-1)*1000+1
          message(st)
          ed <- i*1000
          tmp <- getBM(attributes=attributes.codingseq, filters="ensembl_peptide_id", 
                       values=tr_coding[st:ed, 'pro_name'], mart=mart)
          coding <- rbind(coding, tmp)
          message(round(object.size(coding) / (proc.time()-ptm)[[3]], 3), " B/s")
          #print(i)
        }
        tmp <- getBM(attributes=attributes.codingseq, filters="ensembl_peptide_id", 
                     values=tr_coding[ed+1:length(tr_coding$pro_name), 'pro_name'], 
                     mart=mart)
        coding <- rbind(coding, tmp)
      }
    }
    coding <- read_or_update_local_cache(getCoding(), local_cache_path, "coding")
    stopifnot(nrow(coding) > 0)
    colnames(coding) <- attributes.codingseq 
    tx_id <- transintxdb[match(coding$ensembl_transcript_id, 
                transintxdb$tx_name), 'tx_id']
    procodingseq <- cbind(coding, tx_id)
    colnames(procodingseq) <- c("coding", "pro_name", "tx_name", "tx_id")
    if (!dbsnp_and_cosmic_only) {
        save(procodingseq,file=paste(annotation_path, '/procodingseq.RData', sep=''))
    }
    packageStartupMessage(" done")
    
    message("Prepare protein sequence (proseq.RData) ... ", appendLF=FALSE)
    attributes.proseq <- c("peptide", "ensembl_peptide_id", "ensembl_transcript_id") 
    getProteinseq = function() {
      if(length(tr_coding$pro_name)<1000){
        getBM(attributes=attributes.proseq, filters="ensembl_peptide_id",
              values=tr_coding$pro_name, mart=mart)
      }else{ 
        index <- floor(length(tr_coding$pro_name)/1000)
        proteinseq <- c()
        for(i in 1:index) {
          st <- (i-1)*1000+1
          message(st)
          ed <- i*1000
          tmp <- getBM(attributes=attributes.proseq, filters="ensembl_peptide_id", 
                       values= tr_coding[st:ed, 'pro_name'], mart=mart)
          proteinseq <- rbind(proteinseq, tmp)
          #print(i)
        }
        tmp <- getBM(attributes=attributes.proseq, filters="ensembl_peptide_id", 
                     values=tr_coding[ed+1:length(tr_coding$pro_name), 'pro_name'], 
                     mart=mart)
        proteinseq <- rbind(proteinseq, tmp)
      }
    }
    proteinseq <- read_or_update_local_cache(getProteinseq(), local_cache_path, "proteinseq")
    stopifnot(nrow(proteinseq) > 0)
    colnames(proteinseq) <- c("peptide", "pro_name", "tx_name")
    if (!dbsnp_and_cosmic_only) {
        save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    }
    packageStartupMessage(" done")
    
    
    if (!is.null(dbsnp) && nzchar(dbsnp)) {
        
        message("Prepare dbSNP information (dbsnpinCoding.RData) ... ", appendLF=FALSE)
        
        getSnpTable = function(genome, transGrange) {
            # can't easily query UCSC MySQL with Ensembl transcripts, so always download full coding dbSNP file
            dbSnpFile = qq("@{dbsnp}CodingDbSnp.txt.gz")
            dbSnpURL = qq("ftp://hgdownload.soe.ucsc.edu/goldenPath/@{genome}/database/@{dbSnpFile}")
            download_error = function(e) {stop(qq("error downloading @{dbSnpURL}; either @{dbsnp} is not available for @{dataset} or UCSC's website is down"))}
            if (!is.null(local_cache_path) && length(transcript_ids) > 1000)
                dbSnpFile = file.path(local_cache_path, dbSnpFile)
            if (!file.exists(dbSnpFile))
                tryCatch({download.file(dbSnpURL, dbSnpFile, quiet=T, mode='wb')}, error=download_error, warning=download_error)
            snpCodingTab = .temp_unzip(dbSnpFile, data.table::fread, showProgress=FALSE,
                                       select=c(2:5, 8, 10), sep="\t",
                                       col.names=c("chrom", "chromStart", "chromEnd", "name",
                                                   "alleleCount", "alleles"))
            snpCodingTab$chrom <- gsub('chr', '', snpCodingTab$chrom)
            snpCodingTab <- unique(snpCodingTab)
            snpCoding <- GRanges(seqnames=snpCodingTab$chrom, 
                                 ranges=IRanges(start=snpCodingTab$chromStart, 
                                                end=snpCodingTab$chromEnd), strand='*', 
                                 rsid=snpCodingTab$name, alleleCount=snpCodingTab$alleleCount, 
                                 alleles=snpCodingTab$alleles)
            return(subsetByOverlaps(snpCoding, transGrange))
        }
        dbsnpinCoding = read_or_update_local_cache(getSnpTable(genome, transGrange),
                                                  dbsnp_cache_path, "dbsnpinCoding")

        # snpCodingTab$chrom <- gsub('chr', '', snpCodingTab$chrom)
        # chrlist <- paste(c(seq(1:22),'X','Y'))
        # snpCoding <- subset(snpCodingTab, chrom %in% chrlist ,select=c(chrom:name, alleleCount, alleles))
        # snpCoding <- unique(snpCoding)
        # #snpCoding$chrom <- gsub('chrM', 'MT', snpCoding$chrom)
        # #
        # 
        # #save(snpCoding,file=paste(annotation_path,'/snpcoding.RData',sep=''))
        # snpCoding <- GRanges(seqnames=snpCoding$chrom, 
        #             ranges=IRanges(start=snpCoding$chromStart, 
        #             end=snpCoding$chromEnd), strand='*', 
        #             rsid=snpCoding$name, alleleCount=snpCoding$alleleCount, 
        #             alleles=snpCoding$alleles)
        # 
        # #seqlevels(snpCoding)
        # 
        # #if(TRUE%in% grep('chr',seqlevels(snpCoding)) > 0 ) {
        # #    rchar <- sub('chr','',seqlevels(snpCoding))
        # #    names(rchar) <- seqlevels(snpCoding)
        # #    snpCoding <- renameSeqlevels(snpCoding, rchar) }
        # #if('M'%in%seqlevels(snpCoding)) snpCoding <- renameSeqlevels(snpCoding, c( M='MT'))
        # #chrlist <- paste(c(seq(1:22),'X','Y'))
        # transGrange_snp <- transGrange
        # #transGrange_snp <- keepSeqlevels(transGrange_snp, snpCoding)
        # #snpCoding <- keepSeqlevels(snpCoding, transGrange_snp)
        # 
        # #snpCoding <- keepSeqlevels(snpCoding, transGrange)
        # 
        # dbsnpinCoding <- subsetByOverlaps(snpCoding,transGrange_snp)
        
        save(dbsnpinCoding,file=paste(annotation_path, '/dbsnpinCoding.RData', sep=''))
        packageStartupMessage(" done")
    
    }
    
    if (COSMIC) {
        message("Prepare COSMIC information (cosmic.RData) ... ", appendLF=FALSE)
      
        getCOSMIC = function() {
            # hsapiens_snp_som.default.snp.refsnp_id
            # hsapiens_snp_som.default.snp.chr_name
            # hsapiens_snp_som.default.snp.chrom_start
            # hsapiens_snp_som.default.snp.chrom_end
            #FILTERS=hsapiens_snp_som.default.filters.variation_source."COSMIC"
            varmart = useMart("ENSEMBL_MART_SNP", host=host, dataset="hsapiens_snp_som")
            cosmicAttributes = c("refsnp_id", "chr_name", "chrom_start")
            if (length(unique(exon$gene_name)) < 500) {
                cosmicTab = getBM(attributes=cosmicAttributes,
                                  filters=c("variation_source", "ensembl_gene"),
                                  values=list("COSMIC", unique(exon$gene_name)),
                                  mart=varmart)
            } else {
                cosmicTab = getBM(attributes=cosmicAttributes,
                                  filters="variation_source",
                                  values="COSMIC",
                                  mart=varmart)
            }
            
            chrlist <- paste(c(seq(1:22),'X','Y','MT')) 
            cosmicTab <- subset(cosmicTab, chr_name %in% chrlist)
            cosmic <- GRanges(seqnames=cosmicTab$chr_name, 
                              ranges=IRanges(start=cosmicTab$chrom_start, 
                                             end=cosmicTab$chrom_start), strand = '*', 
                              cosid=cosmicTab$refsnp_id)   
            
            transGrange_cosmic <- transGrange 
            #transGrange_cosmic <- keepSeqlevels(transGrange_cosmic, cosmic)        
            #cosmic <- keepSeqlevels(cosmic, transGrange_cosmic)
            cosmic <- subsetByOverlaps(cosmic, transGrange_cosmic)
            cosmic
        }
        cosmic <- read_or_update_local_cache(getCOSMIC(), local_cache_path, "cosmic")

        save(cosmic, file=paste(annotation_path, '/cosmic.RData', sep=''))    
        packageStartupMessage(" done")        
    }
    
    
    if(!dbsnp_and_cosmic_only && splice_matrix){
        message("Prepare exon splice information (splicemax.RData) ... ", 
                appendLF=FALSE)
        exonByTx <- exonsBy(txdb, "tx", use.names=F)
        index <- which(elementNROWS(exonByTx)==1)
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


if (BiocInstaller::biocVersion() != "3.5")
{
    # GenomicFeature's .Ensembl_getMySQLCoreDir is currently broken for mouse because Ensembl has multiple
    # strains; this version handles it properly: mmusculus_gene_ensembl maps to mus_musculus_core_xx_x;
    # TODO: remove this hack when GenomicFeatures is fixed
    .Ensembl_getMySQLCoreDir <- function(dataset, release=NA, url=NA,
                                         use.grch37=FALSE)
    {
        if (is.na(url))
            url <- GenomicFeatures:::ftp_url_to_Ensembl_mysql(release, use.grch37=use.grch37)
        core_dirs <- GenomicFeatures:::.Ensembl_listMySQLCoreDirs(release=release, url=url,
                                                                  use.grch37=use.grch37)
        shortnames <- shortnames <- sub("(\\w)\\w*?_(\\w+?)_core_\\S+", "\\1\\2", core_dirs, perl=TRUE)
        if (dataset == "mfuro_gene_ensembl") {
            shortname0 <- "mputorius_furo"
        } else {
            shortname0 <- strsplit(dataset, "_", fixed=TRUE)[[1L]][1L]
        }
        core_dir <- core_dirs[shortnames == shortname0]
        if (length(core_dir) != 1L)
            stop("found 0 or more than 1 subdir for \"", dataset,
                 "\" dataset at ", url)
        core_dir
    }
    
    ls_ftp_url <- function (url) 
    {
        doc <- RCurl::getURL(url)
        listing <- strsplit(doc, "\n", fixed = TRUE)[[1L]]
        listing <- listing[stringi::stri_startswith_fixed(listing, "d")]
        pattern <- paste(c("^", rep.int("[^[:space:]]+[[:space:]]+", 
                                        8L)), collapse = "")
        listing <- sub(pattern, "", listing)
        sub("[[:space:]].*$", "", listing)
    }
    
    .extractEnsemblReleaseFromDbVersion <- function (db_version) 
    {
        db_version <- tolower(db_version)
        sub("^ensembl(?: genes)? ([0-9]+).*$", "\\1", db_version, perl=TRUE)
    }
    
    assignInNamespace(".Ensembl_getMySQLCoreDir", .Ensembl_getMySQLCoreDir, "GenomicFeatures")
    assignInNamespace("ls_ftp_url", ls_ftp_url, "GenomicFeatures")
    assignInNamespace(".extractEnsemblReleaseFromDbVersion", .extractEnsemblReleaseFromDbVersion, "GenomicFeatures")
}

# convenient data structure for mapping Ensembl genome and archive hostname to a UCSC dbkey;
# only needed for dbSNP genomes (human and mouse currently); users can override this mapping in case a new
# Ensembl archive hostname or new genome assembly becomes available
DEFAULT_ENSEMBL_UCSC_GENOME_MAP = list("hsapiens_gene_ensembl" = list("may2009.archive.ensembl.org" = "hg18",
                                                                      "may2012.archive.ensembl.org" = "hg19",
                                                                      "dec2013.archive.ensembl.org" = "hg19",
                                                                      "feb2014.archive.ensembl.org" = "hg19",
                                                                      "aug2014.archive.ensembl.org" = "hg38",
                                                                      "oct2014.archive.ensembl.org" = "hg38",
                                                                      "dec2014.archive.ensembl.org" = "hg38",
                                                                      "mar2015.archive.ensembl.org" = "hg38",
                                                                      "may2015.archive.ensembl.org" = "hg38",
                                                                      "jul2015.archive.ensembl.org" = "hg38",
                                                                      "sep2015.archive.ensembl.org" = "hg38",
                                                                      "dec2015.archive.ensembl.org" = "hg38",
                                                                      "mar2016.archive.ensembl.org" = "hg38",
                                                                      "jul2016.archive.ensembl.org" = "hg38",
                                                                      "oct2016.archive.ensembl.org" = "hg38",
                                                                      "dec2016.archive.ensembl.org" = "hg38",
                                                                      "mar2017.archive.ensembl.org" = "hg38",
                                                                      "may2017.archive.ensembl.org" = "hg38",
                                                                      "www.ensembl.org" = "hg38"),
                                       "mmusculus_gene_ensembl" = list("may2009.archive.ensembl.org" = "mm9",
                                                                       "may2012.archive.ensembl.org" = "mm9",
                                                                       "dec2013.archive.ensembl.org" = "mm10",
                                                                       "feb2014.archive.ensembl.org" = "mm10",
                                                                       "aug2014.archive.ensembl.org" = "mm10",
                                                                       "oct2014.archive.ensembl.org" = "mm10",
                                                                       "dec2014.archive.ensembl.org" = "mm10",
                                                                       "mar2015.archive.ensembl.org" = "mm10",
                                                                       "may2015.archive.ensembl.org" = "mm10",
                                                                       "jul2015.archive.ensembl.org" = "mm10",
                                                                       "sep2015.archive.ensembl.org" = "mm10",
                                                                       "dec2015.archive.ensembl.org" = "mm10",
                                                                       "mar2016.archive.ensembl.org" = "mm10",
                                                                       "jul2016.archive.ensembl.org" = "mm10",
                                                                       "oct2016.archive.ensembl.org" = "mm10",
                                                                       "dec2016.archive.ensembl.org" = "mm10",
                                                                       "mar2017.archive.ensembl.org" = "mm10",
                                                                       "may2017.archive.ensembl.org" = "mm10",
                                                                       "www.ensembl.org" = "mm10"))
