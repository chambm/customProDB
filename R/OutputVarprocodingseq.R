##' Output 'snvprocoding'
##'
##' This function uses the output of aaVariation() as input, introduces the nonsynonymous variation into the protein database.
##' @title Output the variant(SNVs) protein coding sequences
##' @param vartable A data frame which is the output of aaVariation().
##' @param procodingseq A dataframe containing protein ids and coding sequence for the protein.
##' @param ids A dataframe containing gene/transcript/protein id mapping information.
##' @param lablersid If includes the dbSNP rsid in the header of each sequence, default is FALSE. 
##'             Must provide dbSNP information in function Positionincoding() if put TRUE here.
##' @param show_progress If true, a progress bar will be shown.
##' @param ... Additional arguments
##' @return a data frame containing protein coding sequence proteins with single nucleotide variation.
##' @author Xiaojing Wang
##' @import Biostrings
##' @importFrom data.table data.table rbindlist setkey setDT set
##' @export
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
##' txlist <- unique(postable_snv$txid)
##' codingseq <- procodingseq[procodingseq$tx_id %in% txlist, ]
##' mtab <- aaVariation (postable_snv, codingseq)
##' OutputVarprocodingseq(mtab, codingseq, ids, lablersid=TRUE)
##' 

OutputVarprocodingseq <- function(vartable, procodingseq, ids, lablersid=FALSE, show_progress=FALSE, ...)
    {
        options(stringsAsFactors=FALSE)
        nonsy <- vartable[vartable$vartype == "non-synonymous", ]
        
        aavar2pro <- nonsy
        
        complements <- c("A"="T","T"="A",
                         "G"="C","C"="G",
                         "M"="K","K"="M",
                         "R"="Y","Y"="R",
                         "W"="W","S"="S",
                         "D"="H","H"="D",
                         "B"="V","V"="B")
        
        fastComplement = function(base) complements[[base]]
        refbase <- mapply(function(base, strand) ifelse(strand=='+', base, fastComplement(base)), aavar2pro$refbase, aavar2pro$strand)
        varbase <- mapply(function(base, strand) ifelse(strand=='+', base, fastComplement(base)), aavar2pro$varbase, aavar2pro$strand)

        aavar2pro$refbase <- refbase
        aavar2pro$varbase <- varbase
        aavar2pro <- aavar2pro[aavar2pro$aaref != "*", ]
        #aavar2pro <- aavar2pro[aavar2pro$aavar!="*", ]
        aavar2pro <- unique(aavar2pro)

        setDT(aavar2pro, key=c("proname", "aapos"))
        
        plist <- unique(aavar2pro$proname)
        pepcoding <- procodingseq[procodingseq$pro_name %in% plist, ]
        pep_vars_by_name = lapply(plist, function(x) aavar2pro[x, .(varbase, aaref, aavar=unlist(aavar), aapos, pincoding, rsid)])
        names(pep_vars_by_name) = plist
        
        coding_index = which(colnames(pepcoding)=="coding")
        complement_index = which(colnames(pepcoding)=="dna_complement")
        
        #pep_var <- pepcoding
        var_names <- vector('character', length(plist))
        if (show_progress) { pb = txtProgressBar(0, nrow(pepcoding), style=3) }
        for(i in 1:nrow(pepcoding)){
            if (show_progress) { setTxtProgressBar(pb, i) }
            #print(i)
            #pvar <- aavar2pro[pep_var$pro_name[i, 'pro_name'])
            #pvar <- pvar[order(as.numeric(pvar$aapos)), ]
            pro_name = pepcoding$pro_name[[i]]
            varcoding = pepcoding$coding[[i]]
            varcomplement = pepcoding$dna_complement[[i]]
            
            pvar = pep_vars_by_name[[pro_name]]
            pincoding = pvar$pincoding
            aaref = pvar$aaref
            aapos = pvar$aapos
            aavar = pvar$aavar
            rsid = pvar$rsid
            varbase = pvar$varbase
            
            var_names_each = vector('character', nrow(pvar))
            
            # apply the SNPs to the reference coding sequence
            for(j in 1:nrow(pvar)){
                p = pincoding[[j]]
                var = substr(varbase[[j]], 1, 1)
                substr(varcoding, p, p) <- var
                substr(varcomplement, p, p) <- fastComplement(var)
                
                if(lablersid && !is.na(rsid[[j]])){
                  var_names_each[[j]] <- paste0(rsid[[j]], ":", aaref[[j]], aapos[[j]], aavar[[j]])
                }else{
                  var_names_each[[j]] <- paste0(aaref[[j]], aapos[[j]], aavar[[j]])
                }
            }
            # assign by reference
            set(pepcoding, i, coding_index, varcoding)
            set(pepcoding, i, complement_index, varcomplement)

            var_names[[i]] = paste0(var_names_each, collapse=",")

        }
        if (show_progress) { close(pb) }
        
        pep_all = data.table(pro_name=pepcoding$pro_name,
                             tx_name=pepcoding$tx_name,
                             tx_id=pepcoding$tx_id,
                             coding=pepcoding$coding,
                             dna_complement=pepcoding$dna_complement,
                             var_name=var_names,
                             key="pro_name")

        setDT(ids, key="pro_name")
        snvprocoding = pep_all[ids, nomatch=0]
        snvprocoding$pro_name <- paste0(snvprocoding$pro_name, "_", snvprocoding$var_name)
        
        snvprocoding[, .(pro_name, tx_name, tx_id,
                         coding, dna_complement,
                         gene_name, description)]
        #write(outformat, file=outfile)
    }
