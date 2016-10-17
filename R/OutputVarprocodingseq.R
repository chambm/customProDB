##' Output 'snvprocoding'
##'
##' This function uses the output of aaVariation() as input, introduces the nonsynonymous variation into the protein database.
##' @title Output the variant(SNVs) protein coding sequences
##' @param vartable A data frame which is the output of aaVariation().
##' @param procodingseq A dataframe containing protein ids and coding sequence for the protein.
##' @param ids A dataframe containing gene/transcript/protein id mapping information.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Must provide dbSNP information in function Positionincoding() if put TRUE here.
##' @param ... Additional arguments
##' @return a data frame containing protein coding sequence proteins with single nucleotide variation.
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
##' OutputVarprocodingseq(mtab, codingseq, ids, lablersid=TRUE)
##' 

OutputVarprocodingseq <- function(vartable, procodingseq, ids, 
            lablersid=FALSE, ...)
    {
        options(stringsAsFactors=FALSE)
        nonsy <- vartable[vartable[, 'vartype'] == "non-synonymous", ]
        
        aavar2pro <- nonsy
        
        refbase <- mapply(function(x, y) 
                ifelse(y=='+', x, toString(reverseComplement(DNAStringSet(x)))), 
                    aavar2pro[, 'refbase'], aavar2pro[, 'strand'])
        varbase <- mapply(function(x,y) 
                ifelse(y=='+', x, toString(reverseComplement(DNAStringSet(x)))), 
                    aavar2pro[, 'varbase'], aavar2pro[, 'strand'])
        aavar2pro[, 'refbase'] <- refbase
        aavar2pro[, 'varbase'] <- varbase
        aavar2pro <- aavar2pro[aavar2pro[, 'aaref']!="*", ]
        #aavar2pro <- aavar2pro[aavar2pro[, 'aavar']!="*", ]
        aavar2pro <- unique(aavar2pro)

        
        plist <- unique(aavar2pro[, 'proname'])
        pepcoding <- procodingseq[procodingseq[, 'pro_name'] %in% plist, ]

        pep_var <- pepcoding
        pep_all<- c()
        test <- c()
        for(i in 1:dim(pep_var)[1]){
            #print(i)
            pvar <-subset(aavar2pro,aavar2pro[, 'proname'] == pep_var[i, 'pro_name'])
            pvar <- pvar[order(as.numeric(pvar[, 'aapos'])), ]
            for(j in 1:dim(pvar)[1]){
                substr(pep_var[i, 'coding'],as.integer(pvar[j, 'pincoding']), 
                    as.integer(pvar[j, 'pincoding'])) <- 
                                substr(pvar[j, 'varbase'], 1, 1)
            }
            if(pep_var[i, 'coding']!=pepcoding[i, 'coding']){
                if(lablersid){
                    var_name <- apply(pvar, 1, function(x) ifelse(is.na(x['rsid']),
                            paste(x['aaref'], x['aapos'], x['aavar'], sep=""), 
                            paste(x['rsid'], ":", x['aaref'], x['aapos'], x['aavar'], 
                            sep="")))
                }else{
                    var_name <- apply(pvar, 1, function(x) 
                            paste(x['aaref'], x['aapos'], x['aavar'], sep=""))
                }
                pep_name <- cbind(pep_var[i,], 
                            var_name=gsub(" ", "", toString(var_name)))
                pep_all <- rbind(pep_all, pep_name)

            }else{
                test <- c(test, pep_var[i, 'pro_name'])
            }
        }
                
        ftab <- merge(pep_all, ids, by.x='pro_name', by.y='pro_name', all=F, 
                    stringsAsFactors = FALSE)
        snvprocoding <-  ftab
        snvprocoding[, 'pro_name'] <- paste(snvprocoding[, 'pro_name'], 
                                        "_", snvprocoding[, 'var_name'], sep='')
        
        
        snvprocoding <- snvprocoding[, c(1:4, 6, 8)]
        colnames(snvprocoding) <- c('pro_name', 'coding', 'tx_name', 'tx_id',
                                    'gene_name', 'description')
        snvprocoding    
        #write(outformat, file=outfile)
        
    }
