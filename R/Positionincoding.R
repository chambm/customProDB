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
##' @export
##' @examples
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' table(GenomicRanges::values(vcf[[1]])[['INDEL']])
##' index <- which(GenomicRanges::values(vcf[[1]])[['INDEL']] == TRUE)
##' indelvcf <- vcf[[1]][index]
##' 
##' index <- which(GenomicRanges::values(vcf[[1]])[['INDEL']] == FALSE)
##' SNVvcf <- vcf[[1]][index]
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "cosmic.RData",package="customProDB"))
##' postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding, COSMIC=cosmic)
##' 


Positionincoding<- function(Vars, exon, dbsnp=NULL, COSMIC=NULL,...)
{
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  
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
  
  seqlevelsInCommon <- intersect(seqlevels(Vars), seqlevels(cdsanno))
  
  # TODO: remove this when customProDB is officially updated to BioC 3.5
  if (packageVersion("GenomeInfoDb") >= "1.11.6") {
      cdsanno <- keepSeqlevels(cdsanno, seqlevelsInCommon, pruning.mode="coarse")
      Vars <- keepSeqlevels(Vars, seqlevelsInCommon, pruning.mode="coarse")
  } else {
      cdsanno <- keepSeqlevels(cdsanno, seqlevelsInCommon)
      Vars <- keepSeqlevels(Vars, seqlevelsInCommon)
  }
  
  #if(TRUE%in% grep('chr',seqlevels(Vars)) > 0 ) {
  #    rchar <- sub('chr','',seqlevels(Vars))
  #    names(rchar) <- seqlevels(Vars)
  #    Vars <- renameSeqlevels(Vars, rchar) }
  #if('M'%in%seqlevels(Vars)) Vars <- renameSeqlevels(Vars, c( M='MT'))

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
  refbase <- values(candiincds)[["REF"]][queryindex]
  varbase <- values(candiincds)[["ALT"]][queryindex]
  
  txname <- values(cdsanno)["trans_name"][hitindex, ]
  txid <- as.integer(values(cdsanno)["trans_id"][hitindex, ])
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
  #res <- data.frame(genename=character(), txname=character(), 
  #                  txid=character(), proname=character(), chr=character(), 
  #                  strand=character(), pos=integer(), refbase=character(), 
  #                  varbase=character(), pincoding=integer())
  res <- data.table(genename, txname, txid, proname, chr, strand, pos, 
                    refbase, varbase, pincoding)
  setkey(res, chr, pos)
  
  #dbsnp <- subset(dbsnp,transcript %in% unique(txname))
  if(!is.null(dbsnp)){
    dbsnp = data.table(seqnames=as.character(seqnames(dbsnp)),
                       startplus1=as.integer(start(dbsnp)+1),
                       rsid=dbsnp$rsid)
    setkey(dbsnp, seqnames, startplus1) # res$pos matches to dbsnp$start+1
    
    #dbsnp_tmp <- subset(dbsnp,select=c(seqnames:end, rsid))
    #dbsnp_unique <- unique(dbsnp_tmp)        
    
    #chrstart <- paste(dbsnp_unique$seqnames, (dbsnp_unique$start)+1, sep=':')
    #dbsnp_id <- data.frame(unique(cbind(name=as.character(dbsnp_unique$rsid), chrstart)))
    
    #chrpos <- paste(res$chr, res$pos,sep=':')
    #res_tmp <- cbind(res, chrpos)
    
    #tmptab <- merge(nonsy,dbsnp_id,by.x='chrpos',by.y='chrstart',all.x=T)
    #rsid <- dbsnp_id[match(res_tmp$chrpos, dbsnp_id$chrstart), 'name']
    #res <- cbind(res, rsid)
    
    res = dbsnp[res, .(genename, txname, txid, proname, chr, strand, pos, refbase, varbase, pincoding, rsid), mult="first"]
    setkey(res, chr, pos)
  }
  
  if(!is.null(COSMIC)){
    COSMIC <- data.table(seqnames=as.character(seqnames(COSMIC)),
                         startplus1=as.integer(start(COSMIC)+1),
                         cosid=as.character(COSMIC$cosid))
    setkey(COSMIC, seqnames, startplus1) # res$pos matches to dbsnp$start+1
    #cos_chrstart <- paste(COSMIC$seqnames, (COSMIC$start)+1, sep=':')
    #cos_id <- data.frame(unique(cbind(name=as.character(COSMIC$cosid), cos_chrstart)))
    
    #chrpos <- paste(res$chr, res$pos, sep=':')
    #res_tmp <- cbind(res, chrpos)
    
    #tmptab <- merge(nonsy,dbsnp_id,by.x='chrpos',by.y='chrstart',all.x=T)
    #COSMIC_id <- cos_id[match(res_tmp$chrpos,cos_id$cos_chrstart), 'name']
    #res <- cbind(res, COSMIC_id)
    
    if(is.null(dbsnp))
      res = COSMIC[res, .(genename, txname, txid, proname, chr, strand, pos, refbase, varbase, pincoding, cosid), mult="first"]
    else
      res = COSMIC[res, .(genename, txname, txid, proname, chr, strand, pos, refbase, varbase, pincoding, rsid, cosid), mult="first"]
  }
  
  res[pincoding > 0]
}
