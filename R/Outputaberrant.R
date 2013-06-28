##' Short insertion/deletion may lead to aberrant proteins in cells. We provide a function to generate FASTA file containing this kind of proteins.
##'
##'
##' the function applys the INDEL into the coding sequence, then translates them into protein sequence, terminated by stop codon. Remove the sequences the same as normal ones or as part of normal ones.
##' @title generate FASTA file containing short INDEL
##' @param positiontab a data frame which is the output of function Positionincoding() for INDELs.
##' @param outfile output file name
##' @param coding a data frame cotaining coding sequence for each protein.
##' @param proteinseq a data frame cotaining amino acid sequence for each protein.
##' @param ids a dataframe containing gene/transcript/protein id mapping information.
##' @param RPKM if includes the RPKM value in the header of each sequence, default is NULL.
##' @param ... Additional arguments.
##' @return FASTA file containing aberrant proteins.
##' @author Xiaojing Wang
##' @examples
##' \dontrun{
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' table(values(vcf[[1]])[['INDEL']])
##' index <- which(values(vcf[[1]])[['INDEL']] == TRUE)
##' indelvcf <- vcf[[1]][index]
##' 
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "dbsnpinCoding.RData", 
##'         package="customProDB"))
##' load(system.file("extdata/refseq", "procodingseq.RData", 
##'         package="customProDB"))
##' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' postable_indel <- Positionincoding(indelvcf, exon)
##' txlist_indel <- unique(postable_indel[, 'txid'])
##' codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% txlist_indel, ]
##' outfile <-  paste(tempdir(), '/test_indel.fasta', sep='')
##' Outputaberrant(postable_indel, coding=codingseq_indel, 
##' proteinseq=proteinseq, outfile=outfile, ids=ids)
##'                 
##' }
##' 
 

Outputaberrant <- function(positiontab, outfile, coding, proteinseq, ids, 
                    RPKM=NULL, ...)
    {
        options(stringsAsFactors=FALSE)
        idx <- grep(',', positiontab[, 'varbase'], fixed=T)
        if(length(idx) > 0) {
            muti <- positiontab[idx, ]
            muti_new <- c()
            for(i in 1:dim(muti)[1]){
                tmp <- cbind(muti[i, 1:8], 
                        unlist(strsplit(muti[i, 'varbase'], ',', fixed=T)), 
                        muti[i, 10])
                muti_new <- rbind(muti_new, tmp)
            }
            colnames(muti_new) <- colnames(positiontab)
            positiontab <- positiontab[-idx, ]
            positiontab <- rbind(positiontab, muti_new)
        }


        mtable <- merge(positiontab, coding, by.x='txid', by.y='tx_id', all.x=T, 
                        stringsAsFactors = FALSE)
        #mtable <- cbind(positiontab,codingseq)

        mtab_plus <- subset(mtable, strand == '+')
        mtab_minus <- subset(mtable, strand == '-')

        str_sub(mtab_plus[, 'coding'], mtab_plus[, 'pincoding'], 
        nchar(mtab_plus[, 'refbase']) + mtab_plus[, 'pincoding'] - 1) <- mtab_plus[, 'varbase']

        str_sub(mtab_minus[, 'coding'], mtab_minus[, 'pincoding'], 
        mtab_minus[, 'pincoding'] + nchar(mtab_minus[, 'refbase']) - 1) <- 
        as.character(reverseComplement(DNAStringSet(mtab_minus[, 'varbase'])))

        total <- rbind(mtab_plus,mtab_minus)

        ## the translate() function cann't deal with base which is not 'ATGC', so replace 'N' with 'A' and lable it in the header
        lable <- unlist(lapply(total[,'coding'], function(x){
                         if(grepl ('N',x,fixed=T)){
                            paste('with N in pos ', paste(as.character(
                            gregexpr('N', x, fixed=T)[[1]]), collapse=' '), 
                            sep='')
                         }else  ""}))
        total[,'coding'] <- gsub('N', 'A', total[, 'coding'], fixed=T)
        aa <- as.character(translate(DNAStringSet(total[, 'coding'])))

        total <- cbind(total, aa, lable)
        pep <- apply(total, 1, function(x) unlist(strsplit(x['aa'], '\\*'))[1])
        total <- cbind(total, pep)

        #########remove the entry identical to normal sequence
        plist <- unique(total[, 'proname'])
        proteinseq <- subset(proteinseq, pro_name %in% plist)
        mseq <- merge(total,proteinseq, by.x='proname', by.y='pro_name', 
                    all=FALSE, stringsAsFactors = FALSE)
        #index <- apply(mseq, 1, function(x)  grep(x['pep'],x['peptide'],fixed=T))
        #grep(mseq[4,'pep'],mseq[4,'peptide'])

        index <-c()
        for(i in 1:dim(mseq)[1]){
            index <- c(index, ifelse(grep(mseq[i, 'pep'], 
                            mseq[i, 'peptide'], fixed=T) == 1, i, '' ))
            #print(i)
        }
        if(length(index) > 0) mseq_r <- mseq[-index, ] else mseq_r <- mseq

        mseq_f <- merge(mseq_r, ids, by.x='pro_name', by.y='pro_name', 
                    all=FALSE, stringsAsFactors = FALSE)
        outformat <- paste('>', mseq_f[, 'pro_name'], "_", mseq_f[, 'pincoding'], 
                    ":", mseq_f[, 'refbase'], '>', mseq_f[, 'varbase'], " |", 
                    mseq_f[, 'tx_name'], "|", mseq_f[, 'gene_name'], "|", 
                    mseq_f[, 'description'], "|", mseq_f[, 'lable'], '\n', 
                    mseq_f[, 'pep'], sep='')
        
        if(!is.null(RPKM)){
            v <- unlist(lapply(mseq_f[, 'pro_name'], function(x) 
                ifelse(x %in% names(RPKM), paste(round(RPKM[x], 4)), paste('NA'))))
            mseq_f <- cbind(mseq_f, v)
            mseq_f <- mseq_f[order(as.numeric(mseq_f[, 'v']), decreasing=T),]
            outformat <- paste('>', mseq_f[, 'pro_name'], "_", 
                        mseq_f[, 'pincoding'], ":", mseq_f[, 'refbase'], '>', 
                        mseq_f[, 'varbase'], " |", mseq_f[, 'v'], "|", 
                        mseq_f[, 'tx_name'], "|", mseq_f[, 'gene_name'], "|", 
                        mseq_f[, 'description'], "|", mseq_f[, 'lable'], '\n', 
                        mseq_f[, 'pep'], sep='')
    
        }
                
        write(outformat, file=outfile)


    }
