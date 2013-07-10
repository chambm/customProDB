##' Output the non-synonymous SNVs into FASTA file.
##'
##' This function uses the output of aaVariation() as input, introduces the nonsynonymous variation into the protein database.
##' @title Output the variant(SNVs) protein sequences into FASTA format
##' @param vartable A data frame which is the output of aaVariation().
##' @param proteinseq A dataframe containing protein ids and the protein sequence.
##' @param outfile Output file name.
##' @param ids A dataframe containing gene/transcript/protein id mapping information.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Must provide dbSNP information in function Positionincoding() if put TRUE here.
##' @param RPKM If includes the RPKM value in the header of each sequence, default is NULL.
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
##' outfile <- paste(tempdir(), '/test_snv.fasta',sep='')
##' OutputVarproseq(mtab, proteinseq, outfile, ids, lablersid=TRUE)
##' 


OutputVarproseq <- function(vartable, proteinseq, outfile, ids, lablersid=FALSE, 
            RPKM=NULL, ...)
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
        
        plist <- unique(aavar2pro[, 'proname'])
        pep <- proteinseq[proteinseq[, 'pro_name'] %in% plist, ]

        pep_var <- pep
        pep_all<- c()
        test <- c()
        for(i in 1:dim(pep_var)[1]){
            #print(i)
            pvar <-subset(aavar2pro,aavar2pro[, 'proname'] == pep_var[i, 'pro_name'])
            pvar <- pvar[order(as.numeric(pvar[, 'aapos'])), ]
            for(j in 1:dim(pvar)[1]){
                substr(pep_var[i, 'peptide'], as.integer(pvar[j, 'aapos']), 
                   as.integer(pvar[j, 'aapos'])) <- substr(pvar[j, 'aavar'], 1, 1)
            }
            if(pep_var[i, 'peptide']!=pep[i, 'peptide']){
                if(lablersid){
                    var_name <- apply(pvar, 1, function(x) ifelse(is.na(x['rsid']),
                            paste(x['aaref'], x['aapos'], x['aavar'], sep=""), 
                            paste(x['rsid'], ":", x['aaref'], x['aapos'], x['aavar'], sep="")))
                }else{
                    var_name <- apply(pvar, 1, function(x) 
                            paste(x['aaref'], x['aapos'], x['aavar'], sep=""))
                }
                pep_name <- cbind(pep_var[i,], 
                            var_name=gsub(" ", "", toString(var_name)))
                pep_all <- rbind(pep_all, pep_name)

            }else{
                test <- c(test,pep_var[i, 'pro_name'])
            }
        }
                
        ftab <- merge(pep_all, ids, by.x='pro_name', by.y='pro_name', all=F, 
                    stringsAsFactors = FALSE)
        outformat <- apply(ftab, 1, function(x) 
                    paste('>', x['pro_name'], "_", x['var_name'], " |", 
                    x['tx_name.x'], "|", x['gene_name'], "|", x['description'], 
                    '\n', unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
        

        if(!is.null(RPKM)){
            v <- unlist(lapply(ftab[, 'pro_name'], function(x) 
                        ifelse(x %in% names(RPKM), 
                        paste(round(RPKM[x], 4), sep=';'), paste('NA'))))
            ftab <- cbind(ftab, v)            
            ftab <- ftab[order(as.numeric(ftab[, 'v']), decreasing=T), ]
            outformat <- apply(ftab, 1, function(x) paste('>', x['pro_name'], 
                        "_", x['var_name'], " |", x['v'], "|", x['tx_name.x'], 
                        "|", x['gene_name'], "|", x['description'],'\n', 
                        unlist(strsplit(x['peptide'], '\\*'))[1], sep=''))
        
        }    
        
            
        write(outformat, file=outfile)
        
    }
