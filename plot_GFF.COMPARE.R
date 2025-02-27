




### summarise read counts of reference transcripts

## in each sample
TR.ref.count.sample <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(count)),
                                     by=.(seqnames, strand, cmp_ref, class_code, sample, rep, hpi, Time, cell_line, group)]

## in sample groups
TR.ref.count.group <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(count)),
                                     by=.(seqnames, strand, cmp_ref, class_code, hpi, Time, cell_line, group)]

ggp <- ggplot( 
  TR.ref.count.sample[class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')],
  aes(x=cell_line, y=sum_count, fill=cell_line)) +
  geom_point(position = position_jitter(width = 0.2), size = 0.1) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_wrap(~class_code + as.factor(Time), scales='free', ncol=6)

ggsave(paste0(outdir, '/', 'gffcompare.TR.ref.count.boxplot.jpg'), ggp, width = 24, height = 16)




### summarise NORMALISED read counts of reference transcripts

## in each sample
TR.ref.norm.count.sample <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(norm_count)),
                                     by=.(seqnames, strand, cmp_ref, class_code, sample, rep, hpi, Time, cell_line, group)]

## in sample groups
TR.ref.norm.count.group <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(norm_count)),
                                     by=.(seqnames, strand, cmp_ref, class_code, hpi, Time, cell_line, group)]


ggp <- ggplot( 
  TR.ref.norm.count.sample[class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')],
  aes(x=cell_line, y=sum_count, fill=cell_line)) +
  geom_point(position = position_jitter(width = 0.2), size = 0.1) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_wrap(~class_code + as.factor(Time), scales='free', ncol=6)

ggsave(paste0(outdir, '/', 'gffcompare.TR.ref.norm.count.boxplot.jpg'), ggp, width = 24, height = 16)





### summarise read counts of transcript categories

## in each sample
TR.class.count.sample <- TR.gff.compare.merged.TR.counts.gt[,
                                                            .(sum_count = sum(count)),
                                                            by=.(seqnames, strand, class_code, sample, rep, hpi, Time, cell_line, group)]

## in each sample - filtered classes only
ggp <- ggplot( 
  TR.class.count.sample[class_code %in% c('~', '=', 'c', 'k', 'o', 'x', 'j')],
  aes(x=rep, y=sum_count, fill=class_code)) +
  geom_col(position = 'dodge') +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_nested(cols=vars(as.factor(Time), cell_line),
               rows=vars(class_code), scales = 'free')

ggp
ggsave(paste0(outdir, '/', 'gffcompare.TR.class.filt.sum.jpg'), ggp, width = 21, height = 16)


## in all of the sample groups
TR.class.count.sum <- TR.gff.compare.merged.TR.counts.gt[,
                                                         .(sum_count = sum(count)),
                                                         by=.(seqnames, strand, class_code, hpi, Time, cell_line, group)]


ggp <- ggplot( 
  TR.class.count.sum,
  aes(x=class_code, y=sum_count, fill=class_code)) +
  geom_col(position = 'stack') +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_grid(cols = vars(cell_line))

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.cell_lines.jpg'), ggp, width = 18, height = 15)


ggp <- ggplot( 
  TR.class.count.sum,
  aes(x=as.factor(Time), y=sum_count, fill=class_code)) +
  geom_col(position = 'stack', color = 'black') +
  #scale_fill_cosmic() +
  #coord_flip() +
  theme_bw() + 
  facet_nested_wrap(~as.factor(Time) + cell_line, nrow=1, scales= 'free')

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.cell_lines.stack.jpg'), ggp, width = 18, height = 15)

ggp <- ggplot( 
  TR.class.count.sample,
  aes(x=as.factor(sample), y=sum_count, fill=class_code)) +
  geom_col(position = 'stack', color = 'black') +
  #scale_fill_cosmic() +
  #coord_flip() +
  theme_bw() + 
  facet_nested_wrap(~as.factor(Time) + cell_line, nrow=1, scales= 'free')

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.samples.stack.jpg'), ggp, width = 24, height = 15)


### summarise NORMALISED read counts of transcript categories

## in each sample
TR.class.norm.count.sample <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(norm_count)),
                                     by=.(seqnames, strand, class_code, sample, rep, hpi, Time, cell_line, group)]

## in each sample - filtered classes only
ggp <- ggplot( 
  TR.class.norm.count.sample[class_code %in% c('~', '=', 'c', 'k', 'o', 'x', 'j')],
  aes(x=rep, y=sum_count, fill=class_code)) +
  geom_col(position = 'dodge') +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_nested(cols=vars(as.factor(Time), cell_line),
               rows=vars(class_code), scales = 'free')

#ggp
ggsave(paste0(outdir, '/', 'gffcompare.TR.class.filt.sum.norm.jpg'), ggp, width = 21, height = 16)



ggp <- ggplot( 
  TR.class.norm.count.sample,
  aes(x=as.factor(sample), y=sum_count, fill=class_code)) +
  geom_col(position = 'stack', color = 'black') +
  #scale_fill_cosmic() +
  #coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90)) +
  facet_nested_wrap(~as.factor(Time) + cell_line, nrow=1, scales= 'free')

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.samples.stack.norm.jpg'), ggp, width = 24, height = 15)


## in all of the sample groups
TR.class.norm.count.sum <- 
  TR.gff.compare.merged.TR.counts.gt[, .(sum_count = sum(norm_count)),
                                     by=.(seqnames, strand, class_code, hpi, Time, cell_line, group)]


ggp <- ggplot( 
  TR.class.norm.count.sum,
  aes(x=class_code, y=sum_count, fill=class_code)) +
  geom_col(position = 'stack') +
  #scale_fill_cosmic() +
  coord_flip() +
  theme_bw() + 
  facet_grid(cols = vars(cell_line))

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.cell_lines.norm.jpg'), ggp, width = 18, height = 15)


ggp <- ggplot( 
  TR.class.norm.count.sum,
  aes(x=cell_line, y=sum_count, fill=class_code)) +
  geom_col(position = 'stack', color = 'black') +
  #scale_fill_cosmic() +
  #coord_flip() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90)) +
  facet_nested_wrap(~as.factor(Time) + cell_line, nrow=1, scales= 'free')

ggsave(paste0(outdir, '/', 'gffcompare.TR.class.sum.cell_lines.stack.norm.jpg'), ggp, width = 18, height = 12)



#####


## counts of Transcript classes in each time point


##### Ditsrubution of distances (from ref transcript) in each sample / group, considering read counts



## Class frequency of TR-reads (without read count)
ggp <- ggplot(best.gff.compare.ref.TR.class.freq, 
              aes(class_code, N, fill=class_code)) +
  geom_col(position = 'dodge') +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = 'none') #+ facet_grid(cols=vars(class_code))

ggsave(paste0(outdir, '/', 'gffcompare.class.frequency.jpg'), ggp, width = 9, height = 12)





#### DISTANCE PLOTS #### 

## Absoulte distance of best hits according to class for each transcript
ggb <- ggplot(best.merged.result_gff.compare,
              aes(x=cmp_ref, y=abs_distance, color=class_code, fill=class_code
              )) + 
  geom_boxplot() +
  theme_bw() +
  coord_flip() +
  facet_wrap(~class_code, nrow=1)


ggsave(paste0(outdir, '/gffcompare.abs.distance.TR.class.boxplot.jpg'), ggb, height = 32, width = 24, limitsize = F)  


## prime5 distance of best hits according to class and strand 
ggp <- ggplot(TR.gff.compare.merged.TR[ seqnames   %in% seqs.to.keep       ## filter low baundnance chromosomes
                                        # & !class_code %in% c('x', 'u', 's')]  ## filter unclassifed and antisense categories
                                        # &   class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')
                                        ,],
              aes(class_code, distance_prime5) ) +
  #geom_boxplot() +
  geom_jitter(aes(color=class_code), size=0.1) +
  ylim(c(NA, 5000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(cols=vars(strand))

ggsave(paste0(outdir, '/', 'gffcompare.prime5.distance.class.strand.jitterplot.jpg'), ggp, width = 18, height = 15)


## prime3 distance of best hits according to class
ggp <- ggplot(TR.gff.compare.merged.TR[ seqnames   %in% seqs.to.keep       ## filter low baundnance chromosomes
                                        # & !class_code %in% c('x', 'u', 's')]  ## filter unclassifed and antisense categories
                                        # &   class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')
                                        ,],
              aes(class_code, distance_prime3) ) +
  #geom_boxplot() +
  geom_jitter(aes(color=class_code), size=0.1) +
  ylim(c(NA, 5000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  facet_grid(cols=vars(strand))

ggsave(paste0(outdir, '/', 'gffcompare.prime3.distance.class.strand.jitterplot.jpg'), ggp, width = 18, height = 15)


## prime5 and prime3 distance of best hits according to class
ggp5 <- ggplot(TR.gff.compare.merged.TR[ seqnames   %in% seqs.to.keep       ## filter low baundnance chromosomes
                                         # & !class_code %in% c('x', 'u', 's')]  ## filter unclassifed and antisense categories
                                         # &   class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')
                                         ,],
               aes(class_code, distance_prime5) ) +
  #geom_boxplot() +
  geom_jitter(aes(color=class_code), size=0.1) +
  ylim(c(NA, 5000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90), legend.position = 'none') # + facet_grid(cols=vars(strand))

ggp3 <- ggplot(TR.gff.compare.merged.TR[ seqnames   %in% seqs.to.keep       ## filter low baundnance chromosomes
                                         # & !class_code %in% c('x', 'u', 's')]  ## filter unclassifed and antisense categories
                                         # &   class_code %in% c('~', '=', 'c', 'k', 'm', 'n', 'j')
                                         ,],
               aes(class_code, distance_prime3) ) +
  #geom_boxplot() +
  geom_jitter(aes(color=class_code), size=0.1) +
  ylim(c(NA, 5000)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90), legend.position = 'none') # + facet_grid(cols=vars(strand))

ggp <- cowplot::plot_grid(ggp5, ggp3, ncol = 2)
ggsave(paste0(outdir, '/', 'gffcompare.prime5.prime3.distance.class.jitterplot.jpg'), ggp, width = 18, height = 15)





