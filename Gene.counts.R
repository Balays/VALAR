
##
TR.gene.ov.all    <- fread(paste0(outdir, '/TR.gene.ov.all.tsv'), na.strings = '')
TR.gene.ov.counts <- fread(paste0(outdir, '/TR.gene.ov.counts.tsv'), na.strings = '')
norm_cov_summary <- fread(paste0(outdir, '/norm.cov.summary.tsv'), na.strings = '')



## gather (melt)
TR.gene.ov.counts.gt <- melt(TR.gene.ov.counts, id.vars=1:12, variable.name = "sample", value.name = "TR.count") # measure.vars = metafilt$sample,


#### Include adapter-based filtering
if(adapters == 'both') {
  TR.gene.ov.counts.gt <- TR.gene.ov.counts.gt[correct_tss == T & correct_tes == T, .(TR.count = sum(TR.count)), 
                                               by=.(transcript_id, seqnames, strand, gene, gene_cluster, Kinetic_class, transcript_start, transcript_end, prime5, prime3, sample)]
  
} else if(adapters == 'either') {
  TR.gene.ov.counts.gt <- TR.gene.ov.counts.gt[correct_tss == T | correct_tes == T, .(TR.count = sum(TR.count)), 
                                               by=.(transcript_id, seqnames, strand, gene, gene_cluster, Kinetic_class, transcript_start, transcript_end, prime5, prime3, sample)]
  
} else if(adapters == 'any') {
  TR.gene.ov.counts.gt <- TR.gene.ov.counts.gt[, .(TR.count = sum(TR.count)), 
                                               by=.(transcript_id, seqnames, strand, gene, gene_cluster, Kinetic_class, transcript_start, transcript_end, prime5, prime3, sample)]
  
} else if(adapters == 'prime5') {
  TR.gene.ov.counts.gt <- TR.gene.ov.counts.gt[correct_tss == T, .(TR.count = sum(TR.count)), 
                                               by=.(transcript_id, seqnames, strand, gene, gene_cluster, Kinetic_class, transcript_start, transcript_end, prime5, prime3, sample)]
  
} else if(adapters == 'prime3') {
  TR.gene.ov.counts.gt <- TR.gene.ov.counts.gt[correct_tes == T, .(TR.count = sum(TR.count)), 
                                               by=.(transcript_id, seqnames, strand, gene, gene_cluster, Kinetic_class, transcript_start, transcript_end, prime5, prime3, sample)]
  
}


#### Important !! Merge with average coverage (and metadata) for count normalization !
TR.gene.ov.counts.gt <- merge(TR.gene.ov.counts.gt, norm_cov_summary, by=c('seqnames', 'sample'), all=T)






#### Frequency of genes on transcripts (not counts)
gene.TR.freq <- TR.gene.ov.all[,.N,by=.(seqnames, strand, gene, gene_cluster, Kinetic_class)][order(N, decreasing = T)]


#### Gene counts (abs) - summarise TR counts to get gene and gene cluster counts
gene.sample_count <- TR.gene.ov.counts.gt[, .(gene_count = sum(TR.count)),  by=.(seqnames, strand, gene, gene_cluster, Kinetic_class, sample)]

#gene.sample_count <- unique(TR.gene.ov.counts.gt[,.(seqnames, strand, gene, gene_cluster, Kinetic_class, sample, gene_count)])

gene.sample_count.sp <- dcast(gene.sample_count, 
                              #sample, gene_count,
                              seqnames + strand + gene_cluster + gene + Kinetic_class ~ sample, value.var = 'gene_count', #fun.aggregate = sum,
                              fill=0)


#### Gene counts (normalized)
gene.sample_norm.count <- TR.gene.ov.counts.gt[,  .(gene_count = sum(TR.count)), by=.(seqnames, strand, gene, gene_cluster, Kinetic_class, sample, average_coverage)]
gene.sample_norm.count[, gene_count.norm := gene_count / average_coverage]

gene.sample_norm.count.sp <- dcast(gene.sample_norm.count, 
                                   #sample, gene_count,
                                   seqnames + strand + gene_cluster + gene + Kinetic_class ~ sample, value.var = 'gene_count.norm', #fun.aggregate = sum,
                                   fill=0)






#fwrite(gene.sample_norm.count.sp, paste0(outdir, '/gene.sample_norm.count.sp.tsv'), sep = '\t')
#fwrite(gene.sample_count.sp, paste0(outdir, '/gene.sample.count.sp.tsv'), sep = '\t')
#fwrite(TR.gene.ov.counts.gt,      paste0(outdir, '/TR.gene.ov.counts.gt.tsv'), sep = '\t')
