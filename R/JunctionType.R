##' For identified splice junctions from RNA-Seq, this function finds the junction types for each entry according to the given annotation.
##' Six types of junctions are classified. find more details in the tutorial.
##' 
##' Go to https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for more information about BED format.
##' @title Annotates the junctions in a bed file.
##' @param jun a GRange object for junctions, the output of function Bed2Range.
##' @param splicemax a known exon splice matrix from the annotation.
##' @param txdb a TxDb object.
##' @param ids a dataframe containing gene/transcript/protein id mapping information.
##' @param ... additional arguments
##' @return a data frame of type and source for each junction.
##' @author Xiaojing Wang
##' @examples
##' 
##' bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
##' jun <-  Bed2Range(bedfile,skip=1,covfilter=5)
##' load(system.file("extdata/refseq", "splicemax.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' txdb <- loadDb(system.file("extdata/refseq", "txdb.sqlite", 
##'             package="customProDB"))
##' junction_type <- JunctionType(jun, splicemax, txdb, ids)
##' table(junction_type[, 'jun_type'])
##' 

             
             

JunctionType <- function(jun, splicemax, txdb, ids, ...)
    {
        options(stringsAsFactors=FALSE)
        #jun <- read.table(bedfile, sep='\t', header=F, quote = "\"", 
        #            stringsAsFactors = F, skip=skip)
        #jun5 <- subset(jun, V5 > covfilter)
        
        #part1_len <- as.numeric(as.data.frame(strsplit(jun5[, 'V11'], ','))[1, ])
        #part2_len <- as.numeric(as.data.frame(strsplit(jun5[, 'V11'], ','))[2, ])
        #gap_len <- as.numeric(as.data.frame(strsplit(jun5[, 'V12'], ','))[2, ])
        #part1_sta <- as.numeric(jun5[, 'V2']) + 1
        #part1_end <- part1_sta + part1_len - 1
        #part2_sta <- part1_sta + gap_len
        #part2_end <- as.numeric(jun5[, 'V3'])
        #junction <- as.data.frame(jun)
        #colnames(junction) <- c('chr', 'start', 'end', 'width', 'strand', 'id', 'cov')
        #junction <- data.frame(chr=jun5[, 'V1'], id=jun5[, 'V4'], 
        #            start=jun5[, 'V2'], end=jun5[, 'V3'], cov=jun5[, 'V5'], 
        #            strand=jun5[, 'V6'], part1_len, part2_len, part1_sta, 
        #            part1_end, part2_sta, part2_end)
        #if('chrM' %in% junction$chr){
        #    junction <- junction[-which(junction$chr == 'chrM'), ]
        #}
        #if('MT' %in% junction$chr){
        #    junction <- junction[-which(junction$chr == 'MT'), ]
        #}
        junRange1 <- GRanges(seqnames=seqnames(jun), 
                ranges=IRanges(start=values(jun)[['part1_sta']], end=values(jun)[['part1_end']]), 
                strand=strand(jun), junction_id=values(jun)[['id']])
        junRange2 <- GRanges(seqnames=seqnames(jun), 
                ranges=IRanges(start=values(jun)[['part2_sta']], end=values(jun)[['part2_end']]), 
                strand=strand(jun), junction_id=values(jun)[['id']]) 
        
        ## map to exon
        splice <- paste(splicemax[, 1], splicemax[, 2], sep='-')
        exons <- exons(txdb)
        
        jun_type <- rep("connect two non-exon region", length=length(junRange1))
        part1_type <- rep("non-exon region", length=length(junRange1))
        part2_type <- rep("non-exon region", length=length(junRange2))
        part1_exon <- rep(NA, length=length(junRange1))
        part2_exon <- rep(NA, length=length(junRange2))
        
        match1_any <- findOverlaps(junRange1, exons)
        match2_any <- findOverlaps(junRange2, exons)

        part1_type[queryHits(match1_any)] <- "overlap with known exon"
        part2_type[queryHits(match2_any)] <- "overlap with known exon"
        part1_exon[queryHits(match1_any)] <- 
            values(exons)["exon_id"][subjectHits(match1_any), ]
        part2_exon[queryHits(match2_any)] <- 
            values(exons)["exon_id"][subjectHits(match2_any), ]


        match1 <- findOverlaps(junRange1, exons, type='end')
        match2 <- findOverlaps(junRange2, exons, type='start')

        part1_type[queryHits(match1)] <- "known exon (same end)"
        part2_type[queryHits(match2)] <- "known exon (same start)"
        part1_exon[queryHits(match1)] <- values(exons)["exon_id"][subjectHits(match1), ]
        part2_exon[queryHits(match2)] <- values(exons)["exon_id"][subjectHits(match2), ]
        
        ########################junction type
        # the order below is matters
        ##################################
        jun_type[intersect(unique(queryHits(match1_any)), 
            unique(queryHits(match2_any)))] <- 
                'connect two regions overlaped with known exons'
        jun_type[intersect(unique(queryHits(match1)),
            unique(queryHits(match2)))] <- 'connect two known exon'
        jun_type[intersect(setdiff(unique(queryHits(match1)), 
            unique(queryHits(match2))),unique(queryHits(match2_any)))] <- 
            'connect a known exon and a region overlap with known exon'
        jun_type[intersect(setdiff(unique(queryHits(match2)), 
            unique(queryHits(match1))),unique(queryHits(match1_any)))] <- 
            'connect a known exon and a region overlap with known exon'
      
        jun_type[setdiff(unique(queryHits(match1_any)), 
            unique(queryHits(match2_any)))] <- 
            'connect a region overlap with known exon and a non-exon region'
        jun_type[setdiff(unique(queryHits(match2_any)), 
            unique(queryHits(match1_any)))] <- 
            'connect a region overlap with known exon and a non-exon region'

        jun_type[intersect(setdiff(unique(queryHits(match1)), unique(queryHits(match2))), 
            setdiff(1:length(junRange2),unique(queryHits(match2_any))))] <- 
            'connect a known exon and a non-exon region'
        jun_type[intersect(setdiff(unique(queryHits(match2)), unique(queryHits(match1))), 
            setdiff(1:length(junRange1),unique(queryHits(match1_any))))] <- 
            'connect a known exon and a non-exon region'
             
        #######################find transcript
        trans <- transcripts(txdb)
        
        tx_part1 <- .map2trans(junRange1, trans,ids)
        colnames(tx_part1) <- c("tx_id_part1", "tx_name_part1", "ge_name_part1")
        tx_part2 <- .map2trans(junRange2 ,trans,ids)  
        colnames(tx_part2) <- c("tx_id_part2", "tx_name_part2", "ge_name_part2")
        
        ########################junction type for 'connect two known exon'
        matchid <- merge(as.matrix(match1), as.matrix(match2), by='queryHits', 
                          all=T)
        index_NA <- which(apply(matchid, 1, function(x) any(is.na(x))) == TRUE)
        if(length(index_NA) >= 0) matchid_new <- matchid[-index_NA, ]
        match1_exon <- values(exons[matchid_new[, 'subjectHits.x'], ])['exon_id'][, 1]
        match2_exon <- values(exons[matchid_new[, 'subjectHits.y'], ])['exon_id'][, 1]
        matchexon <- cbind(matchid_new, match1_exon, match2_exon)

        exonmatrix <- paste(matchexon[, 4], matchexon[, 5], sep='-')

        index_know <- which(exonmatrix %in% splice)
        index_unknown <- which(!exonmatrix %in% splice)

        jun_type[matchid_new[index_unknown, 'queryHits']] <- 
            'novel alternative splicing junction'
        jun_type[matchid_new[index_know, 'queryHits']] <- 'known junction'
        
        index_diff_ge <- which(mapply(function(x, y) 
            length(intersect(unlist(strsplit(x, '\\,')), 
            unlist(strsplit(y, '\\,')))), 
           tx_part1[, 'ge_name_part1'], tx_part2[, 'ge_name_part2']) == 0)
        
        index_fu <- intersect(matchid_new[index_unknown, 'queryHits'],index_diff_ge)
        if(length(index_fu) > 0) jun_type[index_fu] <- 'gene fusion'
        
        junction_type <- cbind(as.data.frame(jun), part1_type, 
                    part2_type, part1_exon, part2_exon, jun_type, tx_part1, tx_part2)

        
    }

    

.map2trans <- function(junRan, txs,ids)
    {
        tx_id_part <- rep(NA, length=length(junRan))
        tx_name_part <- rep(NA, length=length(junRan))
        ge_name_part <- rep(NA, length=length(junRan))
        match_tx <- findOverlaps(junRan, txs)
        
        tx_id_match <- values(txs)["tx_id"][subjectHits(match_tx), ]
        tx_name_match <- values(txs)["tx_name"][subjectHits(match_tx), ]
        gene_match <- ids[match(tx_name_match,ids[,'tx_name']),'gene_name']
        tx_part <- cbind(as.data.frame(match_tx), tx_id_match, tx_name_match, 
                    gene_match)
        ttt <- split(tx_part, tx_part$queryHits)
        tx_id_combine <- unlist(lapply(ttt, function(x) 
                                paste(x$tx_id_match, collapse=",")))
        tx_name_combine <- unlist(lapply(ttt, function(x) 
                                paste(x$tx_name_match, collapse=",")))
        ge_name_combine <- unlist(lapply(ttt, function(x) 
                                paste(x$gene_match, collapse=",")))                        
        
        tx_id_part[as.numeric(names(ttt))] <- tx_id_combine
        tx_name_part[as.numeric(names(ttt))] <- tx_name_combine
        ge_name_part[as.numeric(names(ttt))] <- ge_name_combine
        cbind(tx_id_part, tx_name_part,ge_name_part)
    }
