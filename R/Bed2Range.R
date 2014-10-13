##' Read BED file into a GRanges object. 
##' This function requires complete BED file. 
##' Go to https://genome.ucsc.edu/FAQ/FAQformat.html#format1 for more information about BED format.
##'
##' Read BED file contain junctions into a GRanges object.
##' @title Generate a GRanges objects from BED file.
##' @param bedfile a character contains the path and name of a BED file. 
##' @param skip the number of lines of the BED file to skip before beginning to read data, default 1.
##' @param covfilter the number of minimum coverage for the candidate junction, default 5. 
##' @param ... additional arguments
##' @return a GRanges object containing all candidate junctions from the BED file.
##' @author Xiaojing Wang
##' @examples
##' 
##' bedfile <- system.file("extdata/beds", "junctions1.bed", package="customProDB")
##' jun <-  Bed2Range(bedfile, skip=1,covfilter=5)
##' length(jun)
##' 

Bed2Range <- function(bedfile, skip=1, covfilter=5, ...)
    {
		options(stringsAsFactors=FALSE)
		jun <- read.table(bedfile, sep='\t', header=F, quote = "\"", stringsAsFactors=F, 
                skip=skip)
		jun5 <- subset(jun,V5 > covfilter)
		
		part1_len <- as.numeric(as.data.frame(strsplit(jun5[,'V11'],','))[1,])
		part2_len <- as.numeric(as.data.frame(strsplit(jun5[,'V11'],','))[2,])
		gap_len <- as.numeric(as.data.frame(strsplit(jun5[,'V12'],','))[2,])
		part1_sta <- as.numeric(jun5[,'V2'])+1
		part1_end <- part1_sta+part1_len-1
		part2_sta <- part1_sta+gap_len
		part2_end <- as.numeric(jun5[,'V3'])
		
		junction <- data.frame(chr=jun5[, 'V1'], id=jun5[, 'V4'], start=jun5[, 'V2'], 
            end=jun5[,'V3'], cov=jun5[, 'V5'], strand=jun5[, 'V6'], part1_len, 
            part2_len, part1_sta, part1_end, part2_sta, part2_end)
		if('chrM' %in% junction$chr) junction <- junction[-which(junction$chr=='chrM'), ]
		if('MT' %in% junction$chr) junction <- junction[-which(junction$chr=='MT'), ]
		
		junRange <- GRanges(seqnames=junction$chr, ranges=IRanges(start=junction$part1_end, 
                end=junction$part2_sta), strand=junction$strand, id=junction$id,
                cov=junction$cov, part1_len=junction$part1_len, part2_len=junction$part2_len, 
                part1_sta=junction$part1_sta, part1_end=junction$part1_end, 
                part2_sta=junction$part2_sta, part2_end=junction$part2_end)
	
}
