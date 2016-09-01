##' Variations can be divided into SNVs and INDELs.
##' By taking the output of positionincoding() as input, aaVariation() function predicts the consequences of SNVs in the harbored transcript, such as synonymous or non-synonymous.
##'
##' this function predicts the consequence for SNVs. for INDELs, use Outputabberrant().
##' @title get the functional consequencece of SNVs located in coding region
##' @param position_tab a data frame from Positionincoding()
##' @param coding a data frame cotaining coding sequence for each protein.
##' @param ... Additional arguments
##' @return a data frame containing consequence for each variations.
##' @author Xiaojing Wang
##' @examples
##' 
##' vcffile <- system.file("extdata/vcfs", "test1.vcf", package="customProDB")
##' vcf <- InputVcf(vcffile)
##' table(values(vcf[[1]])[['INDEL']])
##' 
##' index <- which(values(vcf[[1]])[['INDEL']]==FALSE)
##' SNVvcf <- vcf[[1]][index]
##' load(system.file("extdata/refseq", "exon_anno.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "dbsnpinCoding.RData", package="customProDB"))
##' load(system.file("extdata/refseq", "procodingseq.RData", package="customProDB"))
##' postable_snv <- Positionincoding(SNVvcf,exon,dbsnpinCoding)
##' txlist <- unique(postable_snv[,'txid'])
##' codingseq <- procodingseq[procodingseq[,'tx_id'] %in% txlist,]
##' mtab <- aaVariation (postable_snv,codingseq)
##' mtab[1:3,]
##' 
##' 
##' 

aaVariation <-  function(position_tab, coding, ...)
{
  options(stringsAsFactors=FALSE)
  
  setkey(position_tab, txid)
  setDT(coding, key="tx_id")
  mtable <- position_tab[coding]
  
  iub_mul <- list("A,C"="M", "C,A"="M",
                  "A,G"="R", "G,A"="R",
                  "A,T"="W", "T,A"="W",
                  "C,G"="S", "G,C"="S",
                  "C,T"="Y", "T,C"="Y",
                  "G,T"="K", "T,G"="K",
                  "C,G,T"="B", "C,T,G"="B", "G,C,T"="B", "G,T,C"="B", "T,C,G"="B", "T,G,C"="B",
                  "A,G,T"="D", "A,T,G"="D", "G,A,T"="D", "G,T,A"="D", "T,A,G"="D", "T,G,A"="D",
                  "A,C,T"="H", "A,T,C"="H", "C,A,T"="H", "C,T,A"="H", "T,A,C"="H", "T,C,A"="H",
                  "A,C,G"="V", "A,G,C"="V", "C,A,G"="V", "C,G,A"="V", "G,A,C"="V", "G,C,A"="V",
                  "A"="A", "T"="T", "G"="G", "C"="C")
  iub <- list("M"=c("A","C"),"R"=c("A","G"),"W"=c("A","T"),"S"=c("C","G"),"Y"=c("C","T"),"K"=c("G","T"),
              "B"=c("G","T","C"),"D"=c("G","T","A"),"H"=c("A","T","C"),"V"=c("G","A","C"),
              "A"="A","T"="T","G"="G","C"="C")
  forward <- c("A","T","G","C","M","R","W","S","Y","K","B","D","H","V")
  reverse <- c("T","A","C","G","K","Y","W","S","R","M","V","H","D","B")
  nucle <- cbind(forward, reverse)
  
  index <- which(nchar(mtable$varbase) > 1)
  var_new <- unlist(lapply(mtable$varbase, function(x) iub_mul[[x]]))
  
  #vars <- mtable$varbase
  #vars[index] <- var_new
  #class(mtable[, varbase]) <- 'character'
  
  mtable$varbase <- var_new
  
  strand <- mtable$strand
  pincoding <- mtable$pincoding
  #refbase <- mtable[,refbase]
  refbase <- mapply(function(x,y) ifelse(y=='+', x, nucle[nucle[, 1]== x, 2]), mtable$refbase, strand)
  varbase <- mapply(function(x,y) ifelse(y=='+', x, nucle[nucle[, 1]== x, 2]), mtable$varbase, strand)
  codeindex <- ceiling(pincoding/3)
  code_s <- (codeindex-1)*3+1
  code_e <-  codeindex*3
  refcode <- substr(mtable$coding, code_s, code_e)
  
  pcode <- ifelse(pincoding%%3==0, 3, pincoding%%3)
  #substr(refcode,pcode,pcode) <- refbase
  varcode <- refcode
  substr(varcode, pcode, pcode) <- varbase
  
  fastTranslate = function(codon) tryCatch(GENETIC_CODE[[codon]], error=function(e) "X")
  
  pb = txtProgressBar(style=3, min=1, max=length(varcode))
  
  varaa = vector('character', length(varcode))
  vartype <- vector('character', length(varcode))
  aaref <- vector('character', length(varcode))
  aapos <- vector('integer', length(varcode))
  aavar <- vector('character', length(varcode))
  for (i in 1:length(varcode))
  {
    setTxtProgressBar(pb, i)
    
    # if there are no ambiguous bases, simply do a quick translation
    if (!stringi::stri_detect_regex(varcode[[i]], "[^ACGT]")) {
      varaa[[i]] = fastTranslate(varcode[[i]])
    } else {
      # expand the ambiguous bases into combinations of unambiguous bases
      tt = iub[unlist(strsplit(varcode[[i]], split = ""))]
      combine = expand.grid(tt, KEEP.OUT.ATTRS=F, stringsAsFactors=F)
      vcodes = apply(combine, 1, paste0, collapse='')
      varaa[[i]] = unique(lapply(vcodes, fastTranslate))
    }

    aaref[[i]] = fastTranslate(refcode[[i]])
    aapos[[i]] = codeindex[[i]]

    cur_aaref = aaref[[i]]
    cur_varaa = varaa[[i]]
    
    if(is.na(match(cur_aaref, cur_varaa))) {
      vartype[[i]] = 'non-synonymous'
      if (length(cur_varaa) > 1) aavar[[i]] = paste(cur_varaa, collapse='')
      else aavar[[i]] = cur_varaa
    } else {
      varaaunique <- cur_varaa[-match(cur_aaref, cur_varaa)]
      if(length(varaaunique)==0) {
        vartype[[i]] = 'synonymous'
        aavar[[i]] = unique(unlist(cur_varaa))
      } else {
        vartype[[i]] = 'non-synonymous'
        aavar[[i]] = paste(varaaunique, collapse='')
      }
    }
  }
  close(pb)

  # return input table with new columns added
  mtable$refcode = unlist(refcode)
  mtable$varcode = unlist(varcode)
  mtable$vartype = vartype
  mtable$aaref = aaref
  mtable$aapos = aapos
  mtable$aavar = aavar
  as.data.frame(mtable)
}

