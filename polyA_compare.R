

files <- list.files('D:/data/PRV_3cell/rebasecall/mapped_v6_whost/LoRTIA/PK-15')

file.dt <- data.table(file = files, group = stri_extract_first_regex(files, '[0-9]*h_.'))


file.dt[,.N,group]





#### PolyA tails: LoRTIA vs Nanopolish

bam.sc.correct[,.N,sample]





bam.sc.correct.filt <- bam.sc.correct[grepl('out_sorted', sample)]
bam.sc.correct.filt[, group  := gsub('_dRNA.*', '_dRNA', sample)]
bam.sc.correct.filt[, sample := group]
bam.sc.correct.filt[, group  := NULL]


polya.nanopolish <- fread('dRNA_polyA_tail_lengths.tsv.gz')
polya.nanopolish[, sample := gsub('_unmapped', '', sample)]
setnames(polya.nanopolish, 'pt', 'Nanopolish_polya_length')


polya.dt <- merge(bam.sc.correct.filt, polya.nanopolish,  by=c('sample', 'qname')) #, all=T)

polya.dt[,LoRTIA_polya_score := fifelse(strand == '+', score.r3, score.l3)]

polya.scores <- melt(polya.dt, 
                     id.vars       = c('qname', 'seqnames', 'strand', 'sample', 'tag', 'polya_sc_len'),
                     measure.vars  = c('LoRTIA_polya_length', 'Nanopolish_polya_length'), 
                     variable.name = 'score_type', value.name = 'polya_score')

#polya.scores[325328,]

ggplot(polya.dt, 
       aes(LoRTIA_polya_score, polya_score)) + 
  geom_boxplot(aes(fill = score_type))


polya.lengths <- melt(polya.dt, 
                     id.vars       = c('qname', 'seqnames', 'strand', 'sample', 'tag', 'LoRTIA_polya_score'),
                     measure.vars  = c('polya_sc_len', 'Nanopolish_polya_length'), 
                     variable.name = 'score_type', value.name = 'polya_score')


ggplot(polya.lengths, 
       aes(as.factor(LoRTIA_polya_score), polya_score)) + 
  geom_boxplot(aes(fill = score_type)) + 
  facet_nested(rows=vars(tag), scales = 'free', independent = T)



ggplot(polya.dt, 
       aes(LoRTIA_polya_score, Nanopolish_polya_length)) + 
  geom_point(aes(color = strand)) + 
  facet_nested(rows=vars(sample),  cols = vars(strand))

