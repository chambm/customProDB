##' Three-frame translation of novel junctions. And remove those could be found in normal protein sequences.
##' This function requires a genome built by BSgenome package.
##' 
##' @title generate peptide FASTA file that contains novel junctions.
##' @param junction_type a data frame which is the output of function JunctionType()
##' @param genome a BSgenome object. (e.g. Hsapiens)
##' @param outfile output file name
##' @param proteinseq a data frame cotaining amino acid sequence for each protein.
##' @param ... Additional arguments.
##' @return FASTA file that contains novel junction peptides.
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
##' chrom <- paste('chr',c(1:22,'X','Y','M'),sep='')
##' junction_type <- subset(junction_type, seqnames %in% chrom)
##' outf_junc <- paste(tempdir(), '/test_junc.fasta', sep='')
##' outf_junc_coding <- paste(tempdir(), '/test_junc_coding.fasta', sep='')
##' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
##' library('BSgenome.Hsapiens.UCSC.hg19')
##' OutputNovelJun <- OutputNovelJun(junction_type, Hsapiens, outf_junc, 
##'             outf_junc_coding, proteinseq)
##' 


OutputNovelJun <- function(junction_type, genome, outfile, outfile_c, 
    proteinseq, ...)
    {
        options(stringsAsFactors=FALSE)
        #ids <- subset(ids,pro_name!='')
        
        #trans <- transcripts(txdb)
        #index <- which(values(trans)[['tx_name']] %in% ids[,'tx_name'])
        #pro_trans <- trans[index]
        novel_junc <- subset(junction_type, jun_type != 'known junction')
        if(!length(grep('chr', novel_junc[, 'seqnames'], fixed=T))>0) {
            novel_junc[, 'seqnames'] <- paste('chr', novel_junc[, 'seqnames'], sep='')
            idx <- which(novel_junc[, 'seqnames'] %in% seqnames(genome))
            novel_junc <- novel_junc[idx, ]
        }
        
        ###remove abnormal junctions
        idx_abn <- union(which(novel_junc[, 'start'] < 0),  
                    which(novel_junc[, 'end'] < 0))
        if(length(idx_abn > 0)) novel_junc <- novel_junc[-idx_abn, ]
        
        junRange1 <- GRanges(seqnames=novel_junc$seqnames, 
                ranges=IRanges(start=novel_junc$part1_sta, 
                end=novel_junc$part1_end), strand=novel_junc$strand, 
                junction_id=novel_junc$id)
        
        junRange2 <- GRanges(seqnames=novel_junc$seqnames, 
                ranges=IRanges(start=novel_junc$part2_sta, 
                end=novel_junc$part2_end), strand=novel_junc$strand, 
                junction_id=novel_junc$id)
        
        #match1_protx <- findOverlaps(junRange1,pro_trans)
        #match2_protx <- findOverlaps(junRange2,pro_trans)
        
        #juntransRange1 <- junRange1[unique(queryHits(match1_protx))]
        #juntransRange2 <- junRange2[unique(queryHits(match2_protx))]
        
        #junseq1 <- getSeq(genome,'chr1',start=1000,end=2000,as.character=TRUE)
        junseq1 <- getSeq(genome, junRange1)
        junseq2 <- getSeq(genome, junRange2)
        
        junseq_cat <- DNAStringSet(mapply(function(x, y) 
                paste(x, y, sep=''), as.data.frame(junseq1), 
                as.data.frame(junseq2))[, 1])
        
        index_plus <- which(strand(junRange1) == '+')
        index_minus <- which(strand(junRange1) == '-')
        seqs_plus <- junseq_cat[index_plus]
        seqs_minus <- reverseComplement(junseq_cat[index_minus])
        seqs <- c(seqs_plus, seqs_minus)
        
        novel_junc_new <- rbind(novel_junc[index_plus, ], 
                                        novel_junc[index_minus, ])
        
        ##Remove sequences contains NNN
        Nindx <- grep('N', seqs)
        if(length(Nindx) > 0){
            seqs <- seqs[-Nindx]
            novel_junc_new <- novel_junc_new[-Nindx, ]
        }
        
        seqs_name <- paste(paste(novel_junc_new[, 'id'], '_', 
                    novel_junc_new[, 'seqnames'], ':', 
                    novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'], 
                    sep=''), novel_junc_new[, 'cov'],sep='|')
                    
        junpepcoding <- data.frame('pro_name'=seqs_name, 
                                    'coding'=as.data.frame(seqs)[, 1])
        save(junpepcoding, file=outfile_c)
        
        peptides_r1 <- translate(seqs)
        peptides_r2 <- translate(subseq(seqs, start=2))
        peptides_r3 <- translate(subseq(seqs, start=3))
        
        junpos_rna_p1 <- ifelse(novel_junc_new[, 'strand'] == '+', 
                        as.numeric(novel_junc_new[, 'part1_len']), 
                        as.numeric(novel_junc_new[, 'part2_len']))
        junpos_rna_p2 <- ifelse(novel_junc_new[, 'strand'] == '+', 
                        as.numeric(novel_junc_new[, 'part1_len'])+1, 
                        as.numeric(novel_junc_new[, 'part2_len'])+1)
        
        junpos_r1_p1 <- ceiling(junpos_rna_p1/3)
        junpos_r1_p2 <- ceiling(junpos_rna_p2/3)
        
        junpos_r2_p1 <- ceiling((junpos_rna_p1-1)/3)
        junpos_r2_p2 <- ceiling((junpos_rna_p2-1)/3)
        
        junpos_r3_p1 <- ceiling((junpos_rna_p1-2)/3)
        junpos_r3_p2 <- ceiling((junpos_rna_p2-2)/3)
        
        name_r1 <- paste(paste(novel_junc_new[, 'id'], '_', novel_junc_new[, 'seqnames'], 
                    ':', novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'], 
                    sep=''), novel_junc_new[, 'cov'], 'ORF1 ', 
                    paste('Junpos:', junpos_r1_p1, '-', junpos_r1_p2, sep=''), 
                    novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'], 
                    novel_junc_new[, 'tx_name_part2'], 
                    novel_junc_new[, 'jun_type'], sep='|')
        name_r2 <- paste(paste(novel_junc_new[, 'id'],'_', novel_junc_new[, 'seqnames'], 
                    ':',novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'], 
                    sep=''), novel_junc_new[, 'cov'], 'ORF2 ', 
                    paste('Junpos:', junpos_r2_p1, '-', junpos_r2_p2, sep=''), 
                    novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'], 
                    novel_junc_new[, 'tx_name_part2'], 
                    novel_junc_new[, 'jun_type'], sep='|')
        name_r3 <- paste(paste(novel_junc_new[, 'id'],'_', novel_junc_new[, 'seqnames'], 
                    ':',novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'], 
                    sep=''), novel_junc_new[, 'cov'], 'ORF3 ', 
                    paste('Junpos:', junpos_r3_p1, '-', junpos_r3_p2, sep=''), 
                    novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'], 
                    novel_junc_new[, 'tx_name_part2'], 
                    novel_junc_new[, 'jun_type'], sep='|')
        
        all_pep<- rbind(cbind(name_r1, as.data.frame(peptides_r1)[, 1]), 
                    cbind(name_r2, as.data.frame(peptides_r2)[, 1]), 
                    cbind(name_r3, as.data.frame(peptides_r3)[, 1]))
        
        ### remove peptide contain stop codon
        index_stop <- grep('*', all_pep[, 2], fixed=T)
        if(length(index_stop) > 0){
            all_pep_rmstop <- all_pep[-index_stop, ]
        }else all_pep_rmstop <- all_pep
        ### check if any peptides can be found in the normal database, remove those

        index_nor <- lapply(all_pep_rmstop[, 2], function(x) 
                            grep(x, proteinseq[, 'peptide'], fixed=T))
        index_nor <- which(unlist(lapply(index_nor, length)) > 0)
        
        if(length(index_nor) > 0){
            all_pep_new <- all_pep_rmstop[-index_nor, ] 
        }else all_pep_new <- all_pep_rmstop
        tmp <- paste('>', all_pep_new[, 1], '\n', all_pep_new[, 2], sep='')
        write(tmp, file=outfile)
        
    }
