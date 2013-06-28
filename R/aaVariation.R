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

        mtable <- merge(position_tab,coding,by.x='txid',by.y='tx_id',all=F,stringsAsFactors = FALSE)

        iub <- list("M"=c("A","C"),"R"=c("A","G"),"W"=c("A","T"),"S"=c("C","G"),"Y"=c("C","T"),"K"=c("G","T"),
            "B"=c("G","T","C"),"D"=c("G","T","A"),"H"=c("A","T","C"),"V"=c("G","A","C"),
            "A"="A","T"="T","G"="G","C"="C")
        forward <- c("A","T","G","C","M","R","W","S","Y","K","B","D","H","V")
        reverse <- c("T","A","C","G","K","Y","W","S","R","M","V","H","D","B")
        nucle <- cbind(forward,reverse)

        index <- which(nchar(as.character(mtable[,'varbase'])) > 1)
        var_mul <- strsplit(as.character(mtable[index,'varbase']),',')
        var_new <- unlist(lapply(var_mul, function(x)
              names(which(lapply(iub, function(y) setequal(x,y))==TRUE))
        ))

        vars <- as.character(mtable[,'varbase'])
        vars[index] <- var_new
        class(mtable[,'varbase']) <- 'character'
        mtable[,'varbase'] <- vars

        strand <- as.character(mtable[,'strand'])
        pincoding <- mtable[,'pincoding']
        #refbase <- mtable[,'refbase']
        refbase <- mapply(function(x,y) ifelse(y=='+',x,nucle[nucle[,1]== x,2]),as.character(mtable[,'refbase']),strand)
        varbase <- mapply(function(x,y) ifelse(y=='+',x,nucle[nucle[,1]== x,2]),as.character(mtable[,'varbase']),strand)
        codeindex <- ceiling(pincoding/3)
        code_s <- (codeindex-1)*3+1
        code_e <-  codeindex*3
        refcode <- substr(mtable[,'coding'],code_s,code_e)
        
        pcode <- ifelse(pincoding%%3==0,3,pincoding%%3)
        #substr(refcode,pcode,pcode) <- refbase
        varcode <- refcode
        substr(varcode,pcode,pcode) <- varbase

        tt <- lapply(varcode, function(x) c(iub[substr(x,1,1)],iub[substr(x,2,2)],iub[substr(x,3,3)]))
        combine <- lapply(tt, function(x) expand.grid(x))
        vcodes <- lapply(combine, function(y) apply(y,1, function(x) paste(x[1],x[2],x[3],sep="")))
        vcodesDNAS <- lapply(vcodes,function(x) DNAStringSet(x))
        refcodesDNAS <- lapply(refcode, function(x) DNAString(x))
        
        #refaa <- lapply(refcodesDNAS,function(x) translate(x))
        refaa <- lapply(refcodesDNAS,function(x){
                         if(grepl ('N',x,fixed=T)) AAString('X')
                         else  translate(x)})
        
        #varaa <- lapply(vcodesDNAS,function(x) unique(translate(x),package='Biostrings'))
        varaa <- lapply(vcodesDNAS,function(x) unique(translate(x)))
        
        
        mtype <- mapply(function(x,y,z)
            if(is.na(match(as.character(x),as.character(y)))){
                    mtype <- paste('non-synonymous',x,z, as.character(unlist(y)),sep="\t")
                }else {
                    varaaunique <- as.character(y)[-match(as.character(x),as.character(y))]
                    if(length(varaaunique)==0){ 
                        mtype <- paste("synonymous",x,z,unique(as.character(y)),sep="\t")
                    }else{
                        mtype <- paste('non-synonymous',x,z,varaaunique,sep="\t")
                    }
                },
        refaa,varaa,codeindex
        )

        #mtab <- data.frame(genename=character(),txname=character(),proname=character(),chr=character(),strand=character(),
        #                    pos=integer(),refbase=character(),varbase=character(),pincoding=integer(), consensusQ=integer(),snpQ=integer(),maxmappQ=integer(),depth=integer(),
        #                    refcode=character(),varcode=character(),mtype=character())
        #mtab <- cbind(mtable[,-dim(mtable)[2]],refcode=unlist(refcode),varcode=unlist(varcode),mtype)
        #mtab <- data.frame(genename=character(),txname=character(),proname=character(),chr=character(),strand=character(),
        #                    pos=integer(),refbase=character(),varbase=character(),pincoding=integer(),
        #                    refcode=character(),varcode=character(),mtype=character(),stringsAsFactors = FALSE)
        aapos <- do.call('rbind',strsplit(mtype,'\t'))[,3]
        vartype <- do.call('rbind',strsplit(mtype,'\t'))[,1]
        aaref <- do.call('rbind',strsplit(mtype,'\t'))[,2]
        aavar <- do.call('rbind',strsplit(mtype,'\t'))[,4]
        mtab <- cbind(mtable[,1:dim(position_tab)[2]],refcode=unlist(refcode),varcode=unlist(varcode),vartype,aaref,aapos,aavar)
                
    }

