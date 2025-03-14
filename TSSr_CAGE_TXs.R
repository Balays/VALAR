

##
outdir <- config$outdir

### dcDNA TransFrags Results
TR.gff.compare.merged.TR.counts.gt <- fread(paste0(config$outdir, "/TR.gff.compare.merged.TR.counts.gt.tsv"), na.strings = '')


### TSSr Results from dcDNA reads (combined clusters)
TSSr_clusters_uni <- fread(paste0(CAGE_config$outdir, '/TSSr.CAGE.uni.Cluster.Results.txt'))
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


CAGE.TSS.TF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1,
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
CAGE.TSS.TF.OV.dup <- CAGE.TSS.TF.OV[transcript_id %in% dup(unique(CAGE.TSS.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
CAGE.TSS.TF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
CAGE.TSS.TF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
CAGE.TSS.TF.OV.dup <- CAGE.TSS.TF.OV.dup[min_abs_dis == abs(dis)]

CAGE.TSS.TF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
CAGE.TSS.TF.OV.dup <- CAGE.TSS.TF.OV.dup[order(score, decreasing = T)]
CAGE.TSS.TF.OV.dup <- CAGE.TSS.TF.OV.dup[!duplicated(transcript_id)]

CAGE.TSS.TF.OV <- CAGE.TSS.TF.OV[!transcript_id %in% dup(unique(CAGE.TSS.TF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
CAGE.TSS.TF.OV <- rbind(CAGE.TSS.TF.OV, CAGE.TSS.TF.OV.dup)


length(unique(CAGE.TSS.TF.OV$transcript_id))
length(unique(CAGE.TSS.TF.OV$transcript_id)) / length(unique(TR.gff.compare.uni$transcript_id))


# 90% percent of dcDNA TransFrags were supported by TSSr Clusters from the dcDNA reads (+/- 10 bp)


TR.gff.compare.uni <- merge(CAGE.TSS.TF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, CAGE_TC_ID=gene, support, score, consensus_peak, transcript_id)],
                            by.x=c('seqnames', 'strand', 'transcript_id'),
                            TR.gff.compare.uni, by.y=c('seqnames', 'strand', 'transcript_id'), all.y=T)



## Calculate significance
calc.CAGE.sig <- 'CAGE' ## OR: 'dcDNA'

if ( calc.CAGE.sig == 'CAGE' ) {

  ## -->> FROM CAGE DATA
  data <- unique( TSSr_clusters_uni[, .(seqnames, strand, CAGE_TC_ID=gene, score, support)] )

} else if ( calc.CAGE.sig == 'dcDNA' ) {  ### dont use this !

  ## -->> FROM dcDNA DATA
  data <- unique( TR.gff.compare.uni[, .(seqnames, strand, CAGE_TC_ID, score, support)] )

}


# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'CAGE significance' based on the calculated thresholds
data[,CAGE_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                           fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                     score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                           fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                           'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + CAGE_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(CAGE_TC_significance, support)]
support_data

## add CAGE_significanxe to transcript table
TR.gff.compare.uni <- merge(TR.gff.compare.uni, data[,.(CAGE_TC_ID, CAGE_TC_significance)], by='CAGE_TC_ID', all.x=T)

CAGE.TR.support.freq <- TR.gff.compare.uni[, .(ratio=.N/nrow(TR.gff.compare.uni)), by = .(CAGE_TC_significance) ][order(CAGE_TC_significance)]


cage.fr.dt <- TR.gff.compare.uni[,.N,by=.(CAGE_TC_significance, CAGE_TC_ID)]


## ennyi TransFrag-okat NEM támaszott alá a CAGE
TR.gff.compare.uni[is.na(CAGE_TC_ID),.N]

## ennyi TransFrag-okat támaszott alá a CAGE
TR.gff.compare.uni[!is.na(CAGE_TC_ID),.N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a CAGE
TR.gff.compare.uni[!is.na(CAGE_TC_ID) & class_code == '=', .N]

## ennyi olyan TransFrag-ot, ami REF TX-el megegyező ("="), támaszott alá a CAGE PONTOSAN
TR.gff.compare.uni[!is.na(CAGE_TC_ID) & consensus_peak == prime5 & class_code == '=', .N]


## ennyi REF TX-et, támaszott alá a CAGE
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(CAGE_TC_ID), cmp_ref]))

## ennyi REF TX-et, támaszott alá a CAGE PONTOSAN
length(unique(TR.gff.compare.uni[class_code == '=' & !is.na(CAGE_TC_ID) & consensus_peak == prime5, cmp_ref]))


#### ####
##


## Merge TSS clusters with reference annotation
#### ####
TR.merged.data <- TR.Ref.data

## window for overlaps
# +/- 10 bp
ov.win <- 9
# approx 1% increase in transfrag count per increase in window

DTx <- TR.merged.data[,.(seqnames, strand, start = prime5.TR,     end = prime5.TR,   transcript_id, exon_number)]
DTy <- TSSr_clusters_uni[, .(seqnames, strand, start = TC.start - ov.win, end = TC.end + ov.win, width = TC.width, gene, support, score, consensus_peak)]

CAGE.TSS.REF.OV <-
  foverlaps2(DTx=DTx,
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )

## Select the closest cluster for transcripts that could be assigned to more than one cluster
CAGE.TSS.REF.OV.dup <- CAGE.TSS.REF.OV[transcript_id %in% dup(unique(CAGE.TSS.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
CAGE.TSS.REF.OV.dup[, dis := fifelse(strand == '+', i.start - consensus_peak, i.end - consensus_peak)]
CAGE.TSS.REF.OV.dup[, min_abs_dis := min(abs(dis)), by=transcript_id ]
CAGE.TSS.REF.OV.dup <- CAGE.TSS.REF.OV.dup[min_abs_dis == abs(dis)]

CAGE.TSS.REF.OV.dup[, dis := NULL][, min_abs_dis := NULL]
CAGE.TSS.REF.OV.dup <- CAGE.TSS.REF.OV.dup[order(score, decreasing = T)]
CAGE.TSS.REF.OV.dup <- CAGE.TSS.REF.OV.dup[!duplicated(transcript_id)]


CAGE.TSS.REF.OV <- CAGE.TSS.REF.OV[!transcript_id %in% dup(unique(CAGE.TSS.REF.OV[, .(transcript_id, gene )])[,transcript_id]), ]
CAGE.TSS.REF.OV <- rbind(CAGE.TSS.REF.OV, CAGE.TSS.REF.OV.dup)


length(unique(CAGE.TSS.REF.OV$transcript_id))
length(unique(CAGE.TSS.REF.OV$transcript_id)) / length(unique(TR.merged.data$transcript_id))

# 90% percent of Reference Transcripts were supported by TSSr Clusters from the dcDNA reads (+/- 10 bp)





#### ITT JAROK :

## To-DO:

# have to check wether the CAGE also better with new clustering parameters !
# Assign Cluster 3-primes




TR.merged.data <- merge(TR.merged.data,
                        by.x=c('seqnames', 'strand', 'transcript_id', 'exon_number'),
                        CAGE.TSS.REF.OV[,.(seqnames, strand, TC.start=start, TC.end=end, CAGE_TC_ID=gene, support, score, consensus_peak, transcript_id, exon_number)],
                        by.y=c('seqnames', 'strand', 'transcript_id', 'exon_number'), all.x=T)



### Calculate significance
calc.CAGE.sig <- 'CAGE' ## OR: 'dcDNA'

if ( calc.CAGE.sig == 'CAGE' ) {
  
  ## -->> FROM CAGE DATA
  data <- unique( TSSr_clusters_uni[, .(seqnames, strand, CAGE_TC_ID=gene, score, support)] )
  
} else if ( calc.CAGE.sig == 'dcDNA' ) {   ### dont use this !
  
  ## -->> FROM dcDNA DATA
  data <- unique( TR.merged.data[, .(seqnames, strand, CAGE_TC_ID, score, support)] )
  
}


# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
#support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) #+ 10

# Categorize 'CAGE significance' based on the calculated thresholds
data[,CAGE_TC_significance := fifelse(support <= support_percentiles[1] | score <= score_percentiles[1], "*",
                                      fifelse(support >  support_percentiles[1] & support < support_percentiles[2] | 
                                                score >  score_percentiles[1]   & score <  score_percentiles[2], "**",
                                              fifelse(support >= support_percentiles[2] | score >= score_percentiles[2], "***",
                                                      'NA'))) ]


ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) +
  facet_wrap(~support + CAGE_TC_significance, nrow=3, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(CAGE_TC_significance, support)]
support_data


## add CAGE_significanxe to transcript table
TR.merged.data <- merge(TR.merged.data, data[,.(CAGE_TC_ID, CAGE_TC_significance)], by='CAGE_TC_ID', all.x=T)


cage.fr.dt <- TR.merged.data[,.N,by=.(CAGE_TC_significance, CAGE_TC_ID)]


CAGE.ref.support.freq <- TR.merged.data[, .(ratio=.N/nrow(TR.merged.data)), by = .(CAGE_TC_significance) ][order(CAGE_TC_significance)]


TR.merged.data[, .N, by = .(CAGE_TC_significance) ][order(CAGE_TC_significance)]


ggplot(TR.merged.data, aes(CAGE_TC_significance)) +
  geom_bar(aes(fill=CAGE_TC_significance), color='black') + 
  theme_ipsum() +
  labs(title = 'CAGE support of Reference Transcripts',
       y = 'Number of Refernce Transcripts')




####


#### Add CAGE significance (from all TR's read count) to CAGEfighter table 
cage.sig            <- unique(TR.gff.compare.uni[,.(seqnames, strand, CAGE_TC_ID, CAGE_TC_significance)])
TSSr_clusters_uni   <- merge(TSSr_clusters_uni, cage.sig, by.x=c('seqnames', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'CAGE_TC_ID'), all.x=T)

#### Finalize tables
CAGE.support.freq       <- rbind(CAGE.TR.support.freq[,source := 'Kakuk_et_al'], CAGE.ref.support.freq[,source := 'Torma_et_al'])
TR.gff.compare.uni.CAGE <- TR.gff.compare.uni
TR.Ref.data.CAGE        <- TR.merged.data


#### Write out
fwrite(TR.Ref.data.CAGE,        paste0(outdir, '/TR.Ref.data.CAGE.tsv'),        sep = '\t')
fwrite(TR.gff.compare.uni.CAGE, paste0(outdir, '/TR.gff.compare.uni.CAGE.tsv'), sep = '\t')
fwrite(CAGE.support.freq,       paste0(outdir, '/CAGE.support.freq.tsv'),       sep = '\t')
fwrite(TSSr_clusters_uni,       paste0(outdir, '/CAGE.TSSr_clusters_uni.tsv'),  sep = '\t')

#### Read back
TR.Ref.data.CAGE        <- fread(paste0(outdir, '/TR.Ref.data.CAGE.tsv'),        na.strings = '')
TR.gff.compare.uni.CAGE <- fread(paste0(outdir, '/TR.gff.compare.uni.CAGE.tsv'), na.strings = '')
CAGE.support.freq       <- fread(paste0(outdir, '/CAGE.support.freq.tsv'),       na.strings = '')
TSSr_clusters_uni       <- fread(paste0(outdir, '/CAGE.TSSr_clusters_uni.tsv'),  na.strings = '')



#### Ennyi elég

####


