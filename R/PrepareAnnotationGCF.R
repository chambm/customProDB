#' Prepare an annotation set from an NCBI GCF assembly.
#'
#' @title Prepare an annotation set from an NCBI GCF assembly.
#' @param ncbiFtpUrl the location of _cds_from_genomic.fna, _translated_from_cds.faa, and _genomic.gff files
#' @param annotation_path specify a folder to store all the annotations.
#' @param ... additional arguments
#' @return No return value. It creates several .RData files in annotation_path containing annotations needed for further analysis.
#' @author Matt Chambers
#'
#' @export
#' 
#' @examples
#' annotation_path <- tempdir()
#' PrepareAnnotationGCF("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/922/835/GCF_001922835.1_NIST_Tur_tru_v1/", annotation_path)
#'
PrepareAnnotationGCF <- function(ncbiFtpUrl, annotation_path, ...) {
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
    
    if (!stringi::stri_endswith_fixed(ncbiFtpUrl, "/"))
        ncbiFtpUrl = paste0(ncbiFtpUrl, "/")
    
    if (!stringi::stri_endswith_fixed(annotation_path, "/"))
        annotation_path = paste0(annotation_path, "/")
    
    files = strsplit(getURL(ncbiFtpUrl, dirlistonly = TRUE), "\r*\n")[[1]]
    assembly_report = files[grep(".*_assembly_report.txt", files)]
    cds_fna = files[grep(".*_cds_from_genomic.fna.gz", files)]
    translated_faa = files[grep(".*_translated_cds.faa.gz", files)]
    genomic_gff = files[grep(".*_genomic.gff.gz", files)]
    stopifnot(stringr::str_length(assembly_report) > 0,
              stringr::str_length(cds_fna) > 0,
              stringr::str_length(translated_faa) > 0,
              stringr::str_length(genomic_gff) > 0)
    
    download.filename(ncbiFtpUrl, assembly_report, annotation_path, quiet=TRUE)
    organism_name = sub("# Organism name:\\s*(.*?)\\s*\\(.*\\)", "\\1", readLines(paste0(annotation_path, assembly_report), n=2)[2])
    
    download.filename(ncbiFtpUrl, cds_fna, annotation_path)
    download.filename(ncbiFtpUrl, translated_faa, annotation_path)
    download.filename(ncbiFtpUrl, genomic_gff, annotation_path)
    
# # Organism name:  Tursiops truncatus (bottlenose dolphin)
    message(paste0("Build TranscriptDB object for '", organism_name, "' (txdb.sqlite) ..."), appendLF=TRUE)    
    txdb <- makeTxDbFromGFF(file=paste0(annotation_path, genomic_gff),
                            format="gff3",
                            dataSource=ncbiFtpUrl,
                            organism=organism_name)
    
    
    saveDb(txdb, file=paste(annotation_path, '/txdb.sqlite', sep=''))
    packageStartupMessage(" done")
    
    tx_by_protein = sqldf("SELECT t._tx_id AS tx_id, tx_name, cds_name AS pro_name 
                           FROM splicing s, transcript t, cds c 
                           WHERE t._tx_id=s._tx_id 
                             AND c._cds_id=s._cds_id 
                             AND cds_name IS NOT NULL 
                           GROUP BY cds_name",
                          dbname = paste0(annotation_path, "/txdb.sqlite"))
    setDT(tx_by_protein)
    setkey(tx_by_protein, pro_name)
    
    protein_by_tx = data.table::copy(tx_by_protein)
    setkey(protein_by_tx, tx_name)
    
    
    message(paste("Prepare gene/transcript/protein id mapping information (ids.RData) ... "), appendLF=FALSE)
    
    gff <- import.gff(paste0(annotation_path, genomic_gff), colnames = c("transcript_id", "gene", "model_evidence", "product"), feature.type="mRNA")
    ids <- data.frame(gff$gene, transcript_id=gff$transcript_id, gff$transcript_id,
                      gff$gene, gff$product, gff$transcript_id,
                      gff$model_evidence, gff$transcript_id,
                      gff$transcript_id, protein_by_tx[gff$transcript_id]$pro_name)

    colnames(ids) <- c('gene_name', 'tx_name', 'gene_type', 
                       'gene_status', 'description', 'transcript_type', 
                       'transcript_status', 'transcript_name', 
                       'level', 'pro_name')
    
    save(ids, file=paste0(annotation_path, '/ids.RData'))
    packageStartupMessage(" done")
    
    message("Prepare exon annotation information (exon_anno.RData) ... ", 
            appendLF=FALSE)
    
    transGrange <- transcripts(txdb)
    tr <- IRanges::as.data.frame(transGrange)
    cdsByTx <- cdsBy(txdb, "tx", use.names=FALSE)
    exonByTx <- exonsBy(txdb, "tx", use.names=FALSE)
    fiveutrByTx <- fiveUTRsByTranscript(txdb, use.names=FALSE)
    threeutrByTx <- threeUTRsByTranscript(txdb, use.names=FALSE)
    
    cdss <-  IRanges::as.data.frame(cdsByTx)
    exons <- IRanges::as.data.frame(exonByTx)
    fiveutrs <- IRanges::as.data.frame(fiveutrByTx)
    threeutrs <- IRanges::as.data.frame(threeutrByTx)
    
    #txid <- matrix(unlist(strsplit(rownames(exons), '\\.')), ncol = 2, 
    #    byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    exon_p <- data.frame(txid=exons[, "group_name"], chr=exons[, "seqnames"], 
                         exon_s=exons[, "start"], exon_e=exons[, "end"], 
                         exon_rank=exons[, "exon_rank"])
    exon2tr <-  merge(exon_p, tr, by.y="tx_id", by.x="txid")
    exon2tr <- exon2tr[, -which(names(exon2tr) %in% c("seqnames", "width"))]
    
    #txid <- matrix(unlist(strsplit(rownames(cdss), '\\.')), ncol = 2, 
    #       byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.',txid)
    cds_p <- data.frame(txid=cdss[, "group_name"], cds_s=cdss[, "start"], 
                        cds_e=cdss[, "end"], exon_rank=cdss[, "exon_rank"], 
                        width=cdss[, "width"])
    ttt <- split(cds_p, cds_p$txid)
    
    cds_p_new_list <-lapply(ttt, function(x){
        #len <- x[,'cds_e']-x[,'cds_s']+1
        #cum <- cumsum(len)
        cum <- cumsum(x[, 'width'])
        rdis <- cbind(c(1, cum[1:length(cum)-1]+1), cum)
        colnames(rdis) <- c('cds_start', 'cds_end')
        tmp <- cbind(x, rdis)
        tmp
    })
    
    
    cds_p_new <- do.call(rbind, cds_p_new_list)
    cds_p_new <- cds_p_new[, -which(names(cds_p_new) %in% c("width"))]
    
    #for(i in 1:length(ttt)) {
    #    print(i)
    #    ttt1 <- ttt[[i]]
    #    len <- ttt1[,'cds_e']-ttt1[,'cds_s']+1
    #    cum <- cumsum(len)
    #    rdis <- cbind(c(1,cum[1:length(cum)-1]+1),cum)
    #    colnames(rdis) <- c('cds_start','cds_end')
    #    tmp <- cbind(ttt1,rdis)
    #    cds_p_new <- rbind(cds_p_new,tmp)
    #}
    
    cds2exon <- merge(exon2tr, cds_p_new, by.x=c("txid", "exon_rank"), 
                      by.y=c("txid", "exon_rank"), all.x = TRUE)
    #txid <- matrix(unlist(strsplit(rownames(fiveutrs), '\\.')), ncol = 2, 
    #               byrow=TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    fiveutr_p <- data.frame(txid=fiveutrs[, "group_name"], 
                            fiveutr_s=fiveutrs[, "start"], 
                            fiveutr_e=fiveutrs[, "end"], 
                            exon_rank=fiveutrs[, "exon_rank"])
    fiveutr2exon <- merge(cds2exon, fiveutr_p, by.x=c("txid", "exon_rank"), 
                          by.y =c("txid", "exon_rank"), all.x = TRUE)
    
    #txid <- matrix(unlist(strsplit(rownames(threeutrs),'\\.')), ncol = 2, 
    #         byrow =TRUE)[, 1]
    #txid <- gsub('=','\\.', txid)
    threeutr_p <- data.frame(txid=threeutrs[, "group_name"], 
                             threeutr_s=threeutrs[, "start"], 
                             threeutr_e=threeutrs[, "end"], 
                             exon_rank=threeutrs[, "exon_rank"])
    threeutr2exon <- merge(fiveutr2exon, threeutr_p, 
                           by.x=c("txid", "exon_rank"),
                           by.y=c("txid", "exon_rank"), all.x = TRUE)
    
    #exon <- merge(threeutr2exon,ids,by.x=c("tx_name"), by.y="tx_name", 
    #         all.x = TRUE)
    exon <- threeutr2exon[order(threeutr2exon$txid, threeutr2exon$exon_rank), ]
    
    colnames(exon) <- c("tx_id","rank", "chromosome_name", "exon_chrom_start", 
                        "exon_chrom_end", "start_position", "end_position", "strand", "tx_name", 
                        "cds_chr_start", "cds_chr_end", "cds_start", "cds_end", "5_utr_start", 
                        "5_utr_end", "3_utr_start", "3_utr_end")
    
    pro_name <- ids[match(exon$tx_name, ids$tx_name), 'pro_name']
    gene_name <- ids[match(exon$tx_name, ids$tx_name), 'gene_name']
    exon <- cbind(exon, pro_name, gene_name)
    
    save(exon, file=paste0(annotation_path, '/exon_anno.RData'))
    packageStartupMessage(" done")

    message("Prepare protein and coding sequences (proseq.RData and procodingseq.RData) ... ", appendLF=FALSE)
    pro_seqs <- readAAStringSet(paste0(annotation_path, translated_faa), format='fasta')
    pro_name_v <- names(pro_seqs)
    
    cds_seqs <- readDNAStringSet(paste0(annotation_path, cds_fna), format='fasta')
    tx_name_v <- names(cds_seqs)
    
    stopifnot(length(pro_seqs) == length(cds_seqs))
    
    valid_pro_ids = grepl("protein_id=", pro_name_v, fixed=TRUE)
    valid_cds_ids = grepl("protein_id=", tx_name_v, fixed=TRUE)
    
    stopifnot(sum(valid_pro_ids) == sum(valid_cds_ids))

    pro_name <- sub("lcl\\|\\S+\\.\\d_prot_([NXY]P_\\d+\\.\\d+).*", "\\1", pro_name_v, perl=TRUE)
    tx_name <- tx_by_protein[pro_name]$tx_name
    
    proteinseq <- as.data.frame(pro_seqs[valid_pro_ids, ])
    proteinseq <- cbind(proteinseq, pro_name_v[valid_pro_ids], pro_name[valid_pro_ids], tx_name[valid_pro_ids])
    colnames(proteinseq) <- c("peptide", "pro_name_v", "pro_name", "tx_name")
    save(proteinseq, file=paste(annotation_path, '/proseq.RData', sep=''))
    
    procodingseq <- as.data.frame(cds_seqs[valid_pro_ids, ])
    procodingseq <- cbind(procodingseq, tx_name_v[valid_pro_ids], pro_name[valid_pro_ids], tx_name[valid_pro_ids], tx_by_protein[pro_name[valid_pro_ids]]$tx_id)
    colnames(procodingseq) <- c("coding", "tx_name_full", "pro_name", "tx_name", "tx_id")
    save(procodingseq, file=paste0(annotation_path, '/procodingseq.RData'))
    
    packageStartupMessage(" done")
}


download.filename = function(url, filename, directory, overwrite=FALSE, ...)
{
    if (file.exists(paste0(directory, filename)))
    {
        if (overwrite) file.remove(paste0(directory, filename))
        else return()
    }
    download.file(paste0(url, filename), paste0(directory, filename), ...)
}