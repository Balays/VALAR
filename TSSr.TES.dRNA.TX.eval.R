

##
outdir <- config$outdir

### dRNA TransFrags Results
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv.gz"), na.strings = '')


### TSSr Results from dRNA reads (combined clusters)
TSSr_clusters_uni <- fread(paste0(outdir, '/TSSr.TES.dRNA.uni.Cluster.Results.txt'))
TSSr_clusters_uni[, gene := paste0('TC_', cluster)]

ggplot(TSSr_clusters_uni) +
  #geom_histogram(bins=50) +
  geom_col(aes(x=support, y=score)) # + facet_grid(rows = vars(support))


ggplot(TSSr_clusters_uni) +
  geom_histogram(aes(TC.width), bins=50)




## Merge TSS clusters with dRNA reads
#### ####

cols_to_group <- setdiff(colnames(TR.gff.compare.merged.TR.counts.gt), c(metacols, "count", 'rep', 'total', 'contig_size', "average_coverage", "norm_count"))
TR.gff.compare.uni <- unique(TR.gff.compare.merged.TR.counts.gt[,..cols_to_group])
dup(TR.gff.compare.uni$transcript_id)

TR.gff.compare.uni[,prime5 := fifelse(strand == '+', start, end)]
TR.gff.compare.uni[,prime3 := fifelse(strand == '+', end, start)]

TR.gff.compare.uni[,TR.prime5.win.start := prime5]
TR.gff.compare.uni[,TR.prime5.win.end   := prime5]
TR.gff.compare.uni[,TR.prime3.win.start := prime3]
TR.gff.compare.uni[,TR.prime3.win.end   := prime3]

## window for overlaps
# +/- 10 bp
ov.win <- 10
# approx 1% increase in transfrag count per increase in window

##
DTx <- TR.gff.compare.uni[,.(seqnames, strand, start = TR.prime3.win.start, end = TR.prime3.win.end, transcript_id)]
DTy <- TSSr_clusters_uni[, .(seqnames, strand, start = TC.start - ov.win, end = TC.end + ov.win, width = TC.width, gene, support, score, consensus_peak)]


dRNA.TES.TF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1,
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
dRNA.TES.TF.OV.dup <- dRNA.TES.TF.OV[transcript_id %in% dup(unique(dRNA.TES.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dRNA.TES.TF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
dRNA.TES.TF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
dRNA.TES.TF.OV.dup <- dRNA.TES.TF.OV.dup[min_abs_dis == abs(dis)]

dRNA.TES.TF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
dRNA.TES.TF.OV.dup <- dRNA.TES.TF.OV.dup[order(score, decreasing = T)]
dRNA.TES.TF.OV.dup <- dRNA.TES.TF.OV.dup[!duplicated(transcript_id)]

dRNA.TES.TF.OV <- dRNA.TES.TF.OV[!transcript_id %in% dup(unique(dRNA.TES.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dRNA.TES.TF.OV <- rbind(dRNA.TES.TF.OV, dRNA.TES.TF.OV.dup)



length(unique(dRNA.TES.TF.OV$transcript_id))
length(unique(dRNA.TES.TF.OV$transcript_id)) / length(unique(TR.gff.compare.uni$transcript_id))


# 90% percent of dRNA TransFrags were supported by TSSr Clusters from the dRNA reads (+/- 10 bp)


TR.gff.compare.uni <- merge(dRNA.TES.TF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, dRNA.TES_TC_ID=gene, support, score, consensus_peak, transcript_id)],
                            by.x=c('seqnames', 'strand', 'transcript_id'),
                            TR.gff.compare.uni, by.y=c('seqnames', 'strand', 'transcript_id'), all.y=T)



## Calculate significance
## -->> FROM dRNA DATA
data <- unique( TSSr_clusters_uni[, .(seqnames, strand, dRNA.TES_TC_ID=gene, score, support)] )



# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'dRNA significance' based on the calculated thresholds
data[,dRNA.TES_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                           fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                     score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                           fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                           'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + dRNA.TES_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(dRNA.TES_TC_significance, support)]
support_data

## add dRNA_significanxe to transcript table
TR.gff.compare.uni <- merge(TR.gff.compare.uni, data[,.(dRNA.TES_TC_ID, dRNA.TES_TC_significance)], by='dRNA.TES_TC_ID', all.x=T)

dRNA.TR.support.freq <- TR.gff.compare.uni[, .(ratio=.N/nrow(TR.gff.compare.uni)), by = .(dRNA.TES_TC_significance) ][order(dRNA.TES_TC_significance)]


dRNA.fr.dt <- TR.gff.compare.uni[,.N,by=.(dRNA.TES_TC_significance, dRNA.TES_TC_ID)]


## ennyi TransFrag-okat NEM támaszott alá a dRNA
TR.gff.compare.uni[ is.na(dRNA.TES_TC_ID), .N]

## ennyi TransFrag-okat támaszott alá a dRNA
TR.gff.compare.uni[!is.na(dRNA.TES_TC_ID), .N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a dRNA
TR.gff.compare.uni[!is.na(dRNA.TES_TC_ID) & class_code == '=', .N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a dRNA PONTOSAN
TR.gff.compare.uni[!is.na(dRNA.TES_TC_ID) & consensus_peak == prime3 & class_code == '=', .N]


## ennyi REF TX-et, támaszott alá a dRNA
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(dRNA.TES_TC_ID), cmp_ref]))

## ennyi REF TX-et, támaszott alá a dRNA PONTOSAN
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(dRNA.TES_TC_ID) & consensus_peak == prime3, cmp_ref]))


#### ####
##


## Merge TSS clusters with reference annotation
#### ####
TR.merged.data <- TR.Ref.data

## window for overlaps
# +/- 10 bp
ov.win <- 10
# approx 1% increase in transfrag count per increase in window

DTx <- TR.merged.data[,.(seqnames, strand, start = prime3.TR,     end = prime3.TR,   transcript_id, exon_number)]
DTy <- TSSr_clusters_uni[, .(seqnames, strand, start = TC.start - ov.win, end = TC.end + ov.win, width = TC.width, gene, support, score, consensus_peak)]

dRNA.TES.REF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
dRNA.TES.REF.OV.dup <- dRNA.TES.REF.OV[transcript_id %in% dup(unique(dRNA.TES.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dRNA.TES.REF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
dRNA.TES.REF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
dRNA.TES.REF.OV.dup <- dRNA.TES.REF.OV.dup[min_abs_dis == abs(dis)]

dRNA.TES.REF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
dRNA.TES.REF.OV.dup <- dRNA.TES.REF.OV.dup[order(score, decreasing = T)]
dRNA.TES.REF.OV.dup <- dRNA.TES.REF.OV.dup[!duplicated(transcript_id)]

dRNA.TES.REF.OV <- dRNA.TES.REF.OV[!transcript_id %in% dup(unique(dRNA.TES.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dRNA.TES.REF.OV <- rbind(dRNA.TES.REF.OV, dRNA.TES.REF.OV.dup)

length(unique(dRNA.TES.REF.OV$transcript_id))
length(unique(dRNA.TES.REF.OV$transcript_id)) / length(unique(TR.merged.data$transcript_id))

# 90% percent of Reference Transcripts were supported by TSSr Clusters from the dRNA reads (+/- 10 bp)


TR.merged.data <- merge(TR.merged.data,
                        by.x=c('seqnames', 'strand', 'transcript_id', 'exon_number'),
                        dRNA.TES.REF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, dRNA.TES_TC_ID=gene, support, score, consensus_peak, transcript_id, exon_number)],
                        by.y=c('seqnames', 'strand', 'transcript_id', 'exon_number'), all.x=T)



### Calculate significance
## -->> FROM dRNA DATA
data <- unique( TSSr_clusters_uni[, .(seqnames, strand, dRNA.TES_TC_ID=gene, score, support)] )


# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'dRNA significance' based on the calculated thresholds
data[,dRNA.TES_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                                      fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                                score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                                              fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                                                      'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + dRNA.TES_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(dRNA.TES_TC_significance, support)]
support_data


## add dRNA_significanxe to transcript table
TR.merged.data <- merge(TR.merged.data, data[,.(dRNA.TES_TC_ID, dRNA.TES_TC_significance)], by='dRNA.TES_TC_ID', all.x=T)


dRNA.fr.dt <- TR.merged.data[,.N,by=.(dRNA.TES_TC_significance, dRNA.TES_TC_ID)]


dRNA.ref.support.freq <- TR.merged.data[, .(ratio=.N/nrow(TR.merged.data)), by = .(dRNA.TES_TC_significance) ][order(dRNA.TES_TC_significance)]


TR.merged.data[, .N, by = .(dRNA.TES_TC_significance) ][order(dRNA.TES_TC_significance)]


ggplot(TR.merged.data, aes(dRNA.TES_TC_significance)) +
  geom_bar(aes(fill=dRNA.TES_TC_significance), color='black') + 
  theme_bw() +
  labs(title = 'dRNA support of Reference Transcripts',
       y = 'Number of Refernce Transcripts')




####


#### Add dRNA significance (from all TR's read count) TSSr Table
dRNA.sig            <- unique(TR.gff.compare.uni[,.(seqnames, strand, dRNA.TES_TC_ID, dRNA.TES_TC_significance)])
TSSr_clusters_uni    <- merge(TSSr_clusters_uni, dRNA.sig, by.x=c('seqnames', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'dRNA.TES_TC_ID'), all.x=T)


#### Finalize tables
dRNA.support.freq       <- rbind(dRNA.TR.support.freq[,source := 'Kakuk_et_al'], dRNA.ref.support.freq[,source := 'Torma_et_al'])
TR.gff.compare.uni.dRNA <- TR.gff.compare.uni
TR.Ref.data.dRNA        <- TR.merged.data


#### Write out
fwrite(TR.Ref.data.dRNA,        paste0(outdir, '/TR.Ref.data.TES.dRNA.tsv'),        sep = '\t')
fwrite(TR.gff.compare.uni.dRNA, paste0(outdir, '/TR.gff.compare.uni.TES.dRNA.tsv'), sep = '\t')
fwrite(dRNA.support.freq,       paste0(outdir, '/TES.dRNA.support.freq.tsv'),       sep = '\t')
fwrite(TSSr_clusters_uni,        paste0(outdir, '/TES.dRNA.TSSr_clusters_uni.tsv'),  sep = '\t')


#### Read back
TR.Ref.data.dRNA        <- fread(paste0(outdir, '/TR.Ref.data.TES.dRNA.tsv'),        na.strings = '')
TR.gff.compare.uni.dRNA <- fread(paste0(outdir, '/TR.gff.compare.uni.TES.dRNA.tsv'), na.strings = '')
dRNA.support.freq       <- fread(paste0(outdir, '/TES.dRNA.support.freq.tsv'),       na.strings = '')
TSSr_clusters_uni        <- fread(paste0(outdir, '/TES.dRNA.TSSr_clusters_uni.tsv'),  na.strings = '')



#### Ennyi elég

####


