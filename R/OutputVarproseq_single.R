##' Output the non-synonymous SNVs into FASTA file, one SNV per sequence.
##'
##' This function uses the output of aaVariation() as input, introduces the nonsynonymous variation into the protein database.
##' If a protein have more than one SNVs, introduce one SNV each time, end up with equal number of sequences.
##' @title Output the variant(SNVs) protein sequences into FASTA format
##' @param vartable A data frame which is the output of aaVariation().
##' @param proteinseq A dataframe containing protein ids and the protein sequence.
##' @param outfile Output file name.
##' @param ids A dataframe containing gene/transcript/protein id mapping information.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Must provide dbSNP information in function Positionincoding() if put TRUE here.
##' @param RPKM If includes the RPKM value in the header of each sequence. default is NULL.
##' @param ... Additional arguments
##' @return FASTA file containing proteins with single nucleotide variation.
##' @author Xiaojing Wang
##' @examples
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' table(values(vcf[[1]])[['INDEL']])
##' index <- which(values(vcf[[1]])[['INDEL']] == FALSE)
##' SNVvcf <- vcf[[1]][index]
##' load(system.file("extdata/refseq", "exon_anno.RData", 
##' package="customProDB"))
##' load(system.file("extdata/refseq", "dbsnpinCoding.RData", 
##'     package="customProDB"))
##' load(system.file("extdata/refseq", "procodingseq.RData", 
##'     package="customProDB"))
##' load(system.file("extdata/refseq", "ids.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "proseq.RData", package="customProDB"))
##' postable_snv <- Positionincoding(SNVvcf, exon, dbsnpinCoding)
##' txlist <- unique(postable_snv[, 'txid'])
##' codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist, ]
##' mtab <- aaVariation (postable_snv, codingseq)
##' outfile <- paste(tempdir(), '/test_snv_single.fasta',sep='')
##' OutputVarproseq_single(mtab, proteinseq, outfile, ids, lablersid=TRUE)
##' 


OutputVarproseq_single <- function(vartable, proteinseq, outfile, ids, 
            lablersid=FALSE, RPKM=NULL, ...)
    {
        options(stringsAsFactors=FALSE)
        nonsy <- vartable[vartable[, 'vartype'] == "non-synonymous", ]
        
        if(lablersid){
            aavar2pro <- subset(nonsy, select=c(genename, txname, proname, 
                        aaref, aapos, aavar, rsid))
        }else{
            aavar2pro <- subset(nonsy,select=c(genename, txname, proname, 
                        aaref, aapos, aavar))
        }
        aavar2pro <- aavar2pro[aavar2pro[, 'aaref']!="*", ]
        #aavar2pro <- aavar2pro[aavar2pro[, 'aavar']!="*", ]
        aavar2pro <- unique(aavar2pro)
        
        pep_all<- c()
        for(i in 1:dim(aavar2pro)[1]){
            pvar <- aavar2pro[i, ]
            pep_nor <- proteinseq[proteinseq[, 'pro_name'] == as.character(pvar['proname']), ]
            pepseq_nor <- as.character(pep_nor['peptide'])
            pepseq_var <- pepseq_nor
            
            substr(pepseq_var, as.integer(pvar['aapos']), 
                as.integer(pvar['aapos'])) <- substr(pvar['aavar'], 1, 1)
            if(pepseq_var != pepseq_nor){
                if(lablersid){
                    var_name <- ifelse(is.na(pvar['rsid']), paste(pvar['aaref'], 
                        pvar['aapos'], pvar['aavar'], sep=""), paste(pvar['rsid'], 
                        ":", pvar['aaref'], pvar['aapos'], pvar['aavar'], sep=""))
                }else{
                    var_name <- paste(pvar['aaref'], pvar['aapos'], pvar['aavar'], 
                                sep="")
                }
                pepseq_name <- cbind(peptide=pepseq_var, 
                        pro_name=pep_nor['pro_name'], var_name=toString(var_name))
                pep_all <- rbind(pep_all, pepseq_name)

            }
        }
       
        ftab <- merge(pep_all, ids, by.x='pro_name', by.y='pro_name', all=F, 
                    stringsAsFactors = FALSE)
        outformat <- apply(ftab, 1, function(x) 
                    paste('>', x['pro_name'], "_", x['var_name'], " |", 
                    x['tx_name'], "|", x['gene_name'], "|", x['description'], 
                    '\n', unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
        

        if(!is.null(RPKM)){
            v <- unlist(lapply(ftab[, 'pro_name'], function(x) 
                        ifelse(x %in% names(RPKM), 
                        paste(round(RPKM[x], 4), sep=';'), paste('NA'))))
            ftab <- cbind(ftab, v)            
            ftab <- ftab[order(as.numeric(ftab[, 'v']), decreasing=T), ]
            outformat <- apply(ftab, 1, function(x) paste('>', x['pro_name'], 
                        "_", x['var_name'], " |", x['v'], "|", x['tx_name'], 
                        "|", x['gene_name'], "|", x['description'],'\n', 
                        unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
        
        }    
        
            
        write(outformat, file=outfile)
        
    }
