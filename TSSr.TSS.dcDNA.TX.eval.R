

##
outdir <- config$outdir

### dcDNA TransFrags Results
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv.gz"), na.strings = '')


### TSSr Results from dcDNA reads (combined clusters)
TSSr_clusters_uni <- fread(paste0(outdir, '/TSSr.TSS.dcDNA.uni.Cluster.Results.txt'))
TSSr_clusters_uni[, gene := paste0('TC_', cluster)]

ggplot(TSSr_clusters_uni) +
  #geom_histogram(bins=50) +
  geom_col(aes(x=support, y=score)) # + facet_grid(rows = vars(support))


ggplot(TSSr_clusters_uni) +
  geom_histogram(aes(TC.width), bins=50)




## Merge TSS clusters with dcDNA reads
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
DTx <- TR.gff.compare.uni[,.(seqnames, strand, start = TR.prime5.win.start, end = TR.prime5.win.end, transcript_id)]
DTy <- TSSr_clusters_uni[, .(seqnames, strand, start = TC.start - ov.win, end = TC.end + ov.win, width = TC.width, gene, support, score, consensus_peak)]


dcDNA.TSS.TF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1,
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
dcDNA.TSS.TF.OV.dup <- dcDNA.TSS.TF.OV[transcript_id %in% dup(unique(dcDNA.TSS.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dcDNA.TSS.TF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
dcDNA.TSS.TF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
dcDNA.TSS.TF.OV.dup <- dcDNA.TSS.TF.OV.dup[min_abs_dis == abs(dis)]

dcDNA.TSS.TF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
dcDNA.TSS.TF.OV.dup <- dcDNA.TSS.TF.OV.dup[order(score, decreasing = T)]
dcDNA.TSS.TF.OV.dup <- dcDNA.TSS.TF.OV.dup[!duplicated(transcript_id)]

dcDNA.TSS.TF.OV <- dcDNA.TSS.TF.OV[!transcript_id %in% dup(unique(dcDNA.TSS.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dcDNA.TSS.TF.OV <- rbind(dcDNA.TSS.TF.OV, dcDNA.TSS.TF.OV.dup)



length(unique(dcDNA.TSS.TF.OV$transcript_id))
length(unique(dcDNA.TSS.TF.OV$transcript_id)) / length(unique(TR.gff.compare.uni$transcript_id))


# 90% percent of dcDNA TransFrags were supported by TSSr Clusters from the dcDNA reads (+/- 10 bp)


TR.gff.compare.uni <- merge(dcDNA.TSS.TF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, dcDNA.TSS_TC_ID=gene, support, score, consensus_peak, transcript_id)],
                            by.x=c('seqnames', 'strand', 'transcript_id'),
                            TR.gff.compare.uni, by.y=c('seqnames', 'strand', 'transcript_id'), all.y=T)



## Calculate significance
## -->> FROM dcDNA DATA
data <- unique( TSSr_clusters_uni[, .(seqnames, strand, dcDNA.TSS_TC_ID=gene, score, support)] )



# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'dcDNA significance' based on the calculated thresholds
data[,dcDNA.TSS_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                           fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                     score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                           fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                           'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + dcDNA.TSS_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(dcDNA.TSS_TC_significance, support)]
support_data

## add dcDNA_significanxe to transcript table
TR.gff.compare.uni <- merge(TR.gff.compare.uni, data[,.(dcDNA.TSS_TC_ID, dcDNA.TSS_TC_significance)], by='dcDNA.TSS_TC_ID', all.x=T)

dcDNA.TR.support.freq <- TR.gff.compare.uni[, .(ratio=.N/nrow(TR.gff.compare.uni)), by = .(dcDNA.TSS_TC_significance) ][order(dcDNA.TSS_TC_significance)]


dcDNA.fr.dt <- TR.gff.compare.uni[,.N,by=.(dcDNA.TSS_TC_significance, dcDNA.TSS_TC_ID)]


## ennyi TransFrag-okat NEM támaszott alá a dcDNA
TR.gff.compare.uni[is.na(dcDNA.TSS_TC_ID),.N]

## ennyi TransFrag-okat támaszott alá a dcDNA
TR.gff.compare.uni[!is.na(dcDNA.TSS_TC_ID),.N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a dcDNA
TR.gff.compare.uni[!is.na(dcDNA.TSS_TC_ID) & class_code == '=', .N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a dcDNA PONTOSAN
TR.gff.compare.uni[!is.na(dcDNA.TSS_TC_ID) & consensus_peak == prime5 & class_code == '=', .N]


## ennyi REF TX-et, támaszott alá a dcDNA
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(dcDNA.TSS_TC_ID), cmp_ref]))

## ennyi REF TX-et, támaszott alá a dcDNA PONTOSAN
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(dcDNA.TSS_TC_ID) & consensus_peak == prime5, cmp_ref]))


#### ####
##


## Merge TSS clusters with reference annotation
#### ####
TR.merged.data <- TR.Ref.data

## window for overlaps
# +/- 10 bp
ov.win <- 10
# approx 1% increase in transfrag count per increase in window

DTx <- TR.merged.data[,.(seqnames, strand, start = prime5.TR,     end = prime5.TR,   transcript_id, exon_number)]
DTy <- TSSr_clusters_uni[, .(seqnames, strand, start = TC.start - ov.win, end = TC.end + ov.win, width = TC.width, gene, support, score, consensus_peak)]

dcDNA.TSS.REF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
dcDNA.TSS.REF.OV.dup <- dcDNA.TSS.REF.OV[transcript_id %in% dup(unique(dcDNA.TSS.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dcDNA.TSS.REF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
dcDNA.TSS.REF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
dcDNA.TSS.REF.OV.dup <- dcDNA.TSS.REF.OV.dup[min_abs_dis == abs(dis)]

dcDNA.TSS.REF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
dcDNA.TSS.REF.OV.dup <- dcDNA.TSS.REF.OV.dup[order(score, decreasing = T)]
dcDNA.TSS.REF.OV.dup <- dcDNA.TSS.REF.OV.dup[!duplicated(transcript_id)]

dcDNA.TSS.REF.OV <- dcDNA.TSS.REF.OV[!transcript_id %in% dup(unique(dcDNA.TSS.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
dcDNA.TSS.REF.OV <- rbind(dcDNA.TSS.REF.OV, dcDNA.TSS.REF.OV.dup)

length(unique(dcDNA.TSS.REF.OV$transcript_id))
length(unique(dcDNA.TSS.REF.OV$transcript_id)) / length(unique(TR.merged.data$transcript_id))

# 90% percent of Reference Transcripts were supported by TSSr Clusters from the dcDNA reads (+/- 10 bp)


TR.merged.data <- merge(TR.merged.data,
                        by.x=c('seqnames', 'strand', 'transcript_id', 'exon_number'),
                        dcDNA.TSS.REF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, dcDNA.TSS_TC_ID=gene, support, score, consensus_peak, transcript_id, exon_number)],
                        by.y=c('seqnames', 'strand', 'transcript_id', 'exon_number'), all.x=T)



### Calculate significance
## -->> FROM dcDNA DATA
data <- unique( TSSr_clusters_uni[, .(seqnames, strand, dcDNA.TSS_TC_ID=gene, score, support)] )


# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'dcDNA significance' based on the calculated thresholds
data[,dcDNA.TSS_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                                      fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                                score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                                              fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                                                      'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + dcDNA.TSS_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(dcDNA.TSS_TC_significance, support)]
support_data


## add dcDNA_significanxe to transcript table
TR.merged.data <- merge(TR.merged.data, data[,.(dcDNA.TSS_TC_ID, dcDNA.TSS_TC_significance)], by='dcDNA.TSS_TC_ID', all.x=T)


dcDNA.fr.dt <- TR.merged.data[,.N,by=.(dcDNA.TSS_TC_significance, dcDNA.TSS_TC_ID)]


dcDNA.ref.support.freq <- TR.merged.data[, .(ratio=.N/nrow(TR.merged.data)), by = .(dcDNA.TSS_TC_significance) ][order(dcDNA.TSS_TC_significance)]


TR.merged.data[, .N, by = .(dcDNA.TSS_TC_significance) ][order(dcDNA.TSS_TC_significance)]


ggplot(TR.merged.data, aes(dcDNA.TSS_TC_significance)) +
  geom_bar(aes(fill=dcDNA.TSS_TC_significance), color='black') + 
  theme_bw() +
  labs(title = 'dcDNA support of Reference Transcripts',
       y = 'Number of Refernce Transcripts')




####


#### Add dcDNA significance (from all TR's read count) TSSr Table
dcDNA.sig            <- unique(TR.gff.compare.uni[,.(seqnames, strand, dcDNA.TSS_TC_ID, dcDNA.TSS_TC_significance)])
TSSr_clusters_uni   <- merge(TSSr_clusters_uni, dcDNA.sig, by.x=c('seqnames', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'dcDNA.TSS_TC_ID'), all.x=T)

#### Finalize tables
dcDNA.support.freq       <- rbind(dcDNA.TR.support.freq[,source := 'Kakuk_et_al'], dcDNA.ref.support.freq[,source := 'Torma_et_al'])
TR.gff.compare.uni.dcDNA <- TR.gff.compare.uni
TR.Ref.data.dcDNA        <- TR.merged.data


#### Write out
fwrite(TR.Ref.data.dcDNA,        paste0(outdir, '/TR.Ref.data.TSS.dcDNA.tsv'),        sep = '\t')
fwrite(TR.gff.compare.uni.dcDNA, paste0(outdir, '/TR.gff.compare.uni.TSS.dcDNA.tsv'), sep = '\t')
fwrite(dcDNA.support.freq,       paste0(outdir, '/TSS.dcDNA.support.freq.tsv'),       sep = '\t')
fwrite(TSSr_clusters_uni,        paste0(outdir, '/TSS.dcDNA.TSSr_clusters_uni.tsv'),  sep = '\t')

#### Read back
TR.Ref.data.dcDNA        <- fread(paste0(outdir, '/TR.Ref.data.TSS.dcDNA.tsv'),        na.strings = '')
TR.gff.compare.uni.dcDNA <- fread(paste0(outdir, '/TR.gff.compare.uni.TSS.dcDNA.tsv'), na.strings = '')
dcDNA.support.freq       <- fread(paste0(outdir, '/TSS.dcDNA.support.freq.tsv'),       na.strings = '')
TSSr_clusters_uni        <- fread(paste0(outdir, '/TSS.dcDNA.TSSr_clusters_uni.tsv'),  na.strings = '')



#### Ennyi elég

####


