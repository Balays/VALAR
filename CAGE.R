
##
### Import clusters
CAGE_clusters <- fread('all.tagClusters.txt')
colnames(CAGE_clusters)[2] <- 'seqnames'

CAGE_clusters[,orientation := fifelse(strand == '+', 1, 0)]
CAGE_clusters[,gene := paste0('CAGE_cluster_', .GRP), by=.(cluster)]
CAGE_clusters[,strand := factor(strand, levels = c('+', '-', '*'))]

CAGE_clusters[,prime5 := fifelse(strand == '+', start, end)][,prime3 := fifelse(strand == '+', end, start)]
CAGE_clusters[,TSS.win.start := start]
CAGE_clusters[,TSS.win.end   := end]
CAGE_clusters[,TSS.win.size  := abs(TSS.win.end - TSS.win.start) + 1]

### Import counts
## TSSr Results
CAGE_counts <- rbind(
  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 1)),
  data.table(readxl::read_xlsx('ALL.samples.TSS.orient.xlsx', 2))
)

colnames(CAGE_counts)[1] <- 'seqnames'
## melt
CAGE_counts <- melt(CAGE_counts, 1:3, value.name = 'count', variable.name = 'sample')
## Exclude unkown samples and sum
CAGE_counts <- CAGE_counts[!grepl('unknown', sample), ]
CAGE_counts <- CAGE_counts[!grepl('sum', sample), ]
## cast
CAGE_counts <- dcast.data.table(CAGE_counts, ... ~ sample, value.var = 'count')
##
CAGE_counts[,start := pos][, end := pos]

### Merge
CAGE_m <- foverlaps2(CAGE_counts, CAGE_clusters, minoverlap=1)

CAGE_m[,start := pos]
CAGE_m[,end   := pos]

CAGE_m[,i.start := NULL]
CAGE_m[,i.end   := NULL]
CAGE_m[,width_x := NULL]
CAGE_m[,width_y := NULL]
CAGE_m[,overlap_size:= NULL]
 
CAGE_m <- melt(CAGE_m, 1:19, value.name = 'count', variable.name = 'sample')

cagefr.clust <- CAGE_m
cagefr.clust[, score   := as.integer(tags)]
cagefr.clust[, support := sum(!is.0(count)), by=.(sample, gene)]

cagefr.clust[,width := TSS.win.size]
cagefr.clust[,TSS.win.size := NULL]

CAGE_m <- dcast.data.table(cagefr.clust, ... ~ sample, value.var = 'count', fill=0)

cagefr.clust[,sample := NULL] [,count := NULL][, pos := NULL]
cagefr.clust[,start := TSS.win.start]
cagefr.clust[,end   := TSS.win.end]

cagefr.clust <- unique(cagefr.clust)


##
library(rtracklayer)


## fastq read counts
pattern <- '*.fastq'
fastq_dir   <- 'D:/data/PRV_3cell/CAGE/fastq'
fastq_files <- list.files(fastq_dir, pattern, full.names = T)
fastq_seqlengths <- purrr::map(fastq_files, fastq.seqlengths)
names(fastq_seqlengths) <- gsub(pattern, '', gsub(paste0(fastq_dir, '/'), '', fastq_files))
fastq_readcounts <- lapply(fastq_seqlengths, length)
fastq_readcounts <- data.frame(sample = names(fastq_readcounts), read_count = t(as.data.frame(fastq_readcounts)), row.names = NULL)
fastq_readcounts

CAGEfightR <- fread('D:/data/PRV_3cell/CAGE/Gmail/CAGEPRV1.tsv')[,-1]
CAGEfightR[,orientation := fifelse(strand == '+', 1, 0)]
CAGEfightR[,gene := paste0('CAGE_cluster_', .GRP), by=.(thick.names)]
CAGEfightR[,strand := factor(strand, levels = c('+', '-', '*'))]

cagefr.clust <- CAGEfightR
cagefr.clust[,prime5 := fifelse(strand == '+', start, end)][,prime3 := fifelse(strand == '+', end, start)]
cagefr.clust[,TSS.win.start := start]
cagefr.clust[,TSS.win.end   := end]



#TR.CAGE <- data.table(as.data.frame(readGFF('../CAGE/Gmail/LT934125.1-2 (1).gff3')))
TR.CAGE <- fread('D:/data/PRV_3cell/CAGE/Gmail/LT934125.1-2 (1).gff3')
TR.CAGE <- TR.CAGE[,-12]
colnames(TR.CAGE) <- c('seqnames', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'category', 'gene', 'CAGE')
TR.CAGE[,CAGE := as.logical(CAGE)]
TR.CAGE[type == 'mRNA', .N, by =.(CAGE)]


ggplot(cagefr.clust) +
  #geom_histogram(bins=50) +
  geom_col(aes(x=support, y=score)) # + facet_grid(rows = vars(support))


ggplot(cagefr.clust) +
  geom_histogram(aes(width), bins=50)
  


##### Merge TSS clusters with dcDNA reads

cols_to_group <- setdiff(colnames(TR.gff.compare.merged.TR.counts.gt), c(metacols, "count", 'rep', 'total', 'contig_size', "average_coverage", "norm_count"))
TR.gff.compare.uni <- unique(TR.gff.compare.merged.TR.counts.gt[,..cols_to_group])
dup(TR.gff.compare.uni$transcript_id)

TR.gff.compare.uni[,prime5 := fifelse(strand == '+', start, end)]
TR.gff.compare.uni[,prime3 := fifelse(strand == '+', end, start)]

TR.gff.compare.uni[,TR.prime5.win.start := prime5]
TR.gff.compare.uni[,TR.prime5.win.end   := prime5]
TR.gff.compare.uni[,TR.prime3.win.start := prime3]
TR.gff.compare.uni[,TR.prime3.win.end   := prime3]



##
DTx <- TR.gff.compare.uni[,.(seqnames, strand, start = TR.prime5.win.start, end = TR.prime5.win.end, transcript_id)]
DTy <- cagefr.clust[,.(seqnames, strand, start = TSS.win.start, end = TSS.win.end, width, gene, support, score, thick.start, thick.end)]

CAGE.TR.OV <-
  foverlaps2(DTx=DTx, 
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )


length(unique(CAGE.TR.OV$transcript_id))

TR.gff.compare.uni <- merge(CAGE.TR.OV[,.(seqnames, strand, CAGE.cluster.start=start, CAGE.cluster.end=end, CAGE_ID=gene, support, score, thick.start, thick.end, transcript_id)], 
                            by.x=c('seqnames', 'strand', 'transcript_id'),
                            TR.gff.compare.uni, by.y=c('seqnames', 'strand', 'transcript_id'), all=T)



## Calculate significance 
calc.CAGE.sig <- 'CAGE' ## OR: 'dcDNA'

if ( calc.CAGE.sig == 'CAGE' ) {
  
  ## -->> FROM CAGE DATA
  data <- unique( cagefr.clust[, .(seqnames, strand, gene, score, support)] )
} else if ( calc.CAGE.sig == 'dcDNA' ) {
  
  ## -->> FROM dcDNA DATA
  data <- unique( TR.gff.compare.uni[, .(seqnames, strand, transcript_id, score, support)] )
  
}


# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
support_percentiles <- c(3,5) #quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE) + 10

# Categorize 'CAGE significance' based on the calculated thresholds
data <- data %>%
  mutate(CAGE_significance = case_when(
    support <= support_percentiles[1] & score <= score_percentiles[1] ~ "*",
    (support > support_percentiles[1] & support <= support_percentiles[2]) | 
      (score > score_percentiles[1] & score <= score_percentiles[2]) ~ "**",
    support > support_percentiles[2] | score > score_percentiles[2] ~ "***",
    TRUE ~ NA_character_
  ))

data[,CAGE_significance := fifelse(support < 3, '*', CAGE_significance)]

ggplot(data) +
  theme_bw() +
  geom_histogram(aes(score), bins=50) + 
  facet_wrap(~support + CAGE_significance, nrow=2, scales='free')


support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(CAGE_significance, support)]
support_data

## add CAGE_significanxe to transcript table
TR.gff.compare.uni <- merge(TR.gff.compare.uni, data[,.(gene, CAGE_significance)], by.x='CAGE_ID', by.y='gene', all=T)

CAGE.TR.support.freq <- TR.gff.compare.uni[, .(ratio=.N/nrow(TR.gff.compare.uni)), by = .(CAGE_significance) ][order(CAGE_significance)]


cage.fr.dt <- TR.gff.compare.uni[,.N,by=.(CAGE_significance, CAGE_ID)]


## ezeket TX-féleséget NEM támaszott alá a CAGE
TR.gff.compare.uni[is.na(CAGE_ID),.N]

## ezeket TX-féleséget támaszott alá a CAGE
TR.gff.compare.uni[!is.na(CAGE_ID),.N]

## ezeket az REF TX-eket ("=") támaszotta alá a CAGE cluster
TR.gff.compare.uni[!is.na(CAGE_ID) & class_code == '=', .N]

## ezeket az REF TX-eket ("=") támaszotta alá a CAGE cluster SZŰRÉSSEL (score > 3.5, előfordulás > 1)
#TR.gff.compare.uni[!is.na(CAGE_ID) & score > 3.5 & support > 1 & class_code == '=', .N]

## ezeket az REF TX-eket ("=") támaszotta alá a CAGE pontosan
TR.gff.compare.uni[!is.na(CAGE_ID) & thick.start == prime5 & class_code == '=', .N]

## ezeket az REF TX-eket ("=") támaszotta alá a CAGE pontosan SZŰRÉSSEL (score > 3.5, előfordulás > 1)
#TR.gff.compare.uni[!is.na(CAGE_ID) & thick.start == prime5 & score > 4 & support > 1 & class_code == '=', .N]

## ezeket az REF TX-eket ("=") NEM támaszotta alá a CAGE
TR.gff.compare.uni[is.na(CAGE_ID) & class_code == '=', .N]


## ezeket az REF-hez hasonló TX-eket ("~") támaszotta alá a CAGE
TR.gff.compare.uni[!is.na(CAGE_ID) & class_code == '~', .N]
## ezeket az REF-hez hasonló TX-eket ("~") támaszotta alá a CAGE
#TR.gff.compare.uni[!is.na(CAGE_ID) & score > 3.6 & support > 1 & class_code == '~', .N]


length(unique(TR.gff.compare.uni[class_code == '=', cmp_ref]))

fwrite(TR.gff.compare.uni, paste0(outdir, '/TR.gff.compare.uni.tsv'), sep = '\t')



##### Merge with reference annotation
TR.merged.data <- TR.Ref.data

## Merge
DTx <- TR.merged.data[,.(seqnames, strand, start = prime5.TR,     end = prime5.TR,   transcript_id, exon_number)]
DTy <- cagefr.clust  [,.(seqnames, strand, start = TSS.win.start, end = TSS.win.end, width, gene, support, score, thick.start, thick.end)]

CAGE.REF.OV <-
  foverlaps2(DTx=DTx, 
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )


length(unique(CAGE.REF.OV$transcript_id))

TR.merged.data <- merge(TR.merged.data, 
                        by.x=c('seqnames', 'strand', 'transcript_id', 'exon_number'),
                        CAGE.REF.OV[,.(seqnames, strand, CAGE.cluster.start=start, CAGE.cluster.end=end, CAGE_ID=gene, support, score, thick.start, thick.end, transcript_id, exon_number)], 
                        by.y=c('seqnames', 'strand', 'transcript_id', 'exon_number'), all.x=T)



### Calculate significance

data <- unique( TR.merged.data[, .(seqnames, strand, start.TR, end.TR, transcript_id, score, support)] )

# Calculate the 50th and 75th percentiles for 'support' and 'score'
support_percentiles <- quantile(data$support, probs = c(0.5, 0.75), na.rm = TRUE)
score_percentiles   <- quantile(data$score,   probs = c(0.5, 0.75), na.rm = TRUE)

# Categorize 'CAGE significance' based on the calculated thresholds
data <- data %>%
  mutate(CAGE_significance = case_when(
    support <= support_percentiles[1] & score <= score_percentiles[1] ~ "*",
    (support > support_percentiles[1] & support <= support_percentiles[2]) | 
      (score > score_percentiles[1] & score <= score_percentiles[2]) ~ "**",
    support > support_percentiles[2] | score > score_percentiles[2] ~ "***",
    TRUE ~ NA_character_
  ))

support_data <- data[,.(min=min(score), max=max(score), mean=mean(score)), by=.(CAGE_significance, support)]
support_data

# Categorize 'CAGE significance' based on score only
#data[,CAGE_significance := fifelse(score <= score_percentiles[1], '*', 
#                                   fifelse(score > score_percentiles[1] & score <= score_percentiles[2], '**', 
#                                           fifelse()))]

TR.merged.data <- merge(TR.merged.data, data, by=c('seqnames', 'strand', 'start.TR', 'end.TR', 'transcript_id', 'score', 'support'))
#TR.merged.data[,CAGE := fifelse(is.na(CAGE_ID), F, T)]


CAGE.ref.support.freq <- TR.merged.data[, .(ratio=.N/nrow(TR.merged.data)), by = .(CAGE_significance) ][order(CAGE_significance)]


#### Add CAGE significance (from all TR's read count) to CAGEfighter table
cage.sig       <- unique(TR.gff.compare.uni[,.(seqnames, strand, CAGE_ID, CAGE_significance)])
cagefr.clust   <- merge(cagefr.clust, cage.sig, by.x=c('seqnames', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'CAGE_ID'), all.x=T)

#### Finalize tables
CAGE.support.freq       <- rbind(CAGE.TR.support.freq[,source := 'Kakuk_et_al'], CAGE.ref.support.freq[,source := 'Torma_et_al'])
TR.gff.compare.uni.CAGE <- TR.gff.compare.uni
TR.Ref.data.CAGE        <- TR.merged.data


#### Write out
fwrite(TR.Ref.data.CAGE,        paste0(outdir, '/TR.Ref.data.CAGE.tsv'),        sep = '\t')
fwrite(TR.gff.compare.uni.CAGE, paste0(outdir, '/TR.gff.compare.uni.CAGE.tsv'), sep = '\t')    
fwrite(CAGE.support.freq,       paste0(outdir, '/CAGE.support.freq.tsv'),       sep = '\t')
fwrite(cagefr.clust,            paste0(outdir, '/CAGE.freq.clust.tsv'),         sep = '\t')

#### Read back
TR.Ref.data.CAGE        <- fread(paste0(outdir, '/TR.Ref.data.CAGE.tsv'), na.strings = '')
TR.gff.compare.uni.CAGE <- fread(paste0(outdir, '/TR.gff.compare.uni.CAGE.tsv'), na.strings = '')
CAGE.support.freq       <- fread(paste0(outdir, '/CAGE.support.freq.tsv'), na.strings = '')
cagefr.clust            <- fread(paste0(outdir, '/CAGE.freq.clust.tsv'), na.strings = '')



#### Ennyi elég

####


