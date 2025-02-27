

## Dereplicate read mapping data into unique alignements, eg. read-transcripts ("transfrags")

derep.alignements <- function(bam.filt) {
  #bam.filt
  dt <- bam.filt
  setDT(dt)
  
  # Sort the table by qname, seqnames, strand, start, and end
  dt <- dt[order(aln_ID, seqnames, strand, start, end)]
  # Create a column that captures the unique combination of exons for each qname
  dt[, exon_combination := toString(paste0(start, "-", end, '::', strand)), by = .(aln_ID)]
  # Create a unique TR_ID for each unique combination of exons (across qname, seqnames, and strand)
  dt[, TR_ID := .GRP, by = .(exon_combination, seqnames, strand)]
  
  return(dt)
}

bam.TR <- derep.alignements(bam.filt)
bam.TR[,TR_start := min(start), by=.(TR_ID)]
bam.TR[,TR_end   := max(end),   by=.(TR_ID)]
aln.TR <- unique(bam.TR[,.(aln_ID, TR_ID)])

# Count the number of unique qname for each TR_ID and sample combination
# This ensures that if a transcript with multiple exons appears more than once for a single qname, it's still counted as one
TR.counts     <- bam.TR[, .(count=uniqueN(qname)), by = .(TR_ID, sample)]
TR.sum.counts <- TR.counts[,.(sum_count=sum(count)), by=.(TR_ID)][order(sum_count, decreasing=T)]

# spreaded count table
TR.counts.sp  <- dcast(TR.counts, TR_ID ~ sample, value.var = 'count', fill = 0 )

## unique transcripts
TR.uni <- unique(bam.TR[,.(seqnames, start, end, strand, TR_ID, exon_combination, prime3, prime5)])
# Add exon_rank column
TR.uni[, exon_rank := 1:.N, by = .(TR_ID)]
# Add exon_ID column
TR.uni[, exon_ID   := .GRP, by = .(seqnames, strand, start, end)]

TR.uni <- TR.uni[order(TR_ID, seqnames, strand, start, end)]

## read counts 
read.counts   <- bam.TR[, .(count=uniqueN(qname)), by = .( sample)]


## summaries
TR.uni[,.N,by=.(seqnames, strand)]

TR.uni[,.N,by=.(seqnames, strand)] %>% spread(strand, N)


