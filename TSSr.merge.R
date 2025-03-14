


## Merge dcDNA TSSr Clustering results and CAGE TSSr clustering results with TransFrags


#
CAGE_TSSr_clusters  <- fread(paste0(CAGE_config$outdir, '/TSSr.CAGE.all.Cluster.Results.txt'))
#
dcDNA_TSSr_clusters <- fread(paste0(config$outdir,      '/TSSr.dcDNA.all.Cluster.Results.txt'))


outdir <- config$outdir

#### dcDNA TSSr
TR.Ref.data.dcDNA        <- fread(paste0(outdir, '/TR.Ref.data.dcDNA.tsv'),        na.strings = '')
TR.gff.compare.uni.dcDNA <- fread(paste0(outdir, '/TR.gff.compare.uni.dcDNA.tsv'), na.strings = '')
dcDNA.support.freq       <- fread(paste0(outdir, '/dcDNA.support.freq.tsv'),       na.strings = '')
dcDNA_TSSr_clusters_uni  <- fread(paste0(outdir, '/dcDNA.TSSr_clusters_uni.tsv'),  na.strings = '')

#### CAGE TSSr
TR.Ref.data.CAGE        <- fread(paste0(outdir, '/TR.Ref.data.CAGE.tsv'),        na.strings = '')
TR.gff.compare.uni.CAGE <- fread(paste0(outdir, '/TR.gff.compare.uni.CAGE.tsv'), na.strings = '')
CAGE.support.freq       <- fread(paste0(outdir, '/CAGE.support.freq.tsv'),       na.strings = '')
CAGE_TSSr_clusters_uni  <- fread(paste0(outdir, '/CAGE.TSSr_clusters_uni.tsv'),  na.strings = '')




#TR.Ref.data.TSSr.comb <- copy(TR.Ref.data.CAGE)
cnames <- c('TC.start', 'TC.end', 'support', 'score', 'consensus_peak')
setnames(TR.Ref.data.CAGE,  cnames, paste0('CAGE_',  cnames))
setnames(TR.Ref.data.dcDNA, cnames, paste0('dcDNA_', cnames))

TR.Ref.data.TSSr.comb <- merge(TR.Ref.data.CAGE, TR.Ref.data.dcDNA,
                               by=colnames(TR.Ref.data.CAGE)[2:17])

ggVennDiagram::ggVennDiagram(list(
  REF_TR = TR.Ref.data$transcript_id,
  CAGE_TC_supported   = TR.Ref.data.CAGE[!is.na(CAGE_TC_ID), transcript_id],
  dcDNA_TC_supported  = TR.Ref.data.dcDNA[!is.na(dcDNA_TC_ID), transcript_id]
  ))



### TransFrags
#TR.Ref.data.TSSr.comb <- copy(TR.Ref.data.CAGE)
cnames <- c('TC.start', 'TC.end', 'support', 'score', 'consensus_peak')
setnames(TR.gff.compare.uni.CAGE,  cnames, paste0('CAGE_',  cnames))
setnames(TR.gff.compare.uni.dcDNA, cnames, paste0('dcDNA_', cnames))

TransFrags.TSSr.comb <- merge(TR.gff.compare.uni.CAGE, TR.gff.compare.uni.dcDNA,
                               by=colnames(TR.gff.compare.uni.CAGE)[c(2:4,10:27)])


ggVennDiagram::ggVennDiagram(list(
  TransFrag = TR.gff.compare.uni$transcript_id,
  CAGE_TC_supported   = TR.gff.compare.uni.CAGE[!is.na(CAGE_TC_ID),   transcript_id],
  dcDNA_TC_supported  = TR.gff.compare.uni.dcDNA[!is.na(dcDNA_TC_ID), transcript_id]
))

### Non-spliced Transfrags ONly !!!
TR_EX           <- dcast(EX.uni, seqnames + strand + TR_ID ~ exon_rank, value.var = 'exon_ID')
non.spliced.TFs <- TR_EX[is.na(`2`), TR_ID]

TransFrags.TSSr.comb <- TransFrags.TSSr.comb[transcript_id %in% non.spliced.TFs, ]


### Validate TESs
# make start and end as prime3 for joining
TransFrags.TSSr.comb[,start := prime3]
TransFrags.TSSr.comb[,end   := prime3]

# "Valid" prime3 sites form reference transcripts
#valid.prime3[,TES := start + 10]
#fwrite(valid.prime3, 'valid.prime3.tsv', sep = '\t')

# had to correct, as some TESs were very close ! -->> cehck !!!
valid.prime3 <- fread('valid.prime3.tsv')

prime3.TF.ov <- foverlaps2(TransFrags.TSSr.comb, valid.prime3, by.x=c('seqnames', 'strand', 'start', 'end'), by.y=c('seqnames', 'strand', 'start', 'end'), minoverlap = 1)
prime3.TF.ov <- unique(prime3.TF.ov[,.(seqnames, strand, transcript_id, valid_tes, TES)])

# Merge with TransFrags
TransFrags.TSSr.comb <- merge(TransFrags.TSSr.comb, prime3.TF.ov, by=c('seqnames', 'strand', 'transcript_id'), all.x=T)
TransFrags.TSSr.comb[is.na(valid_tes), valid_tes := F]

# correct start and end
TransFrags.TSSr.comb[, start  := fifelse(strand == '+', prime5, prime3)]
TransFrags.TSSr.comb[, end    := fifelse(strand == '+', prime3, prime5)]

prime3.valid.corr.freq <- TransFrags.TSSr.comb[,.N, by=.(valid_tes)]

ggvf <- ggplot(prime3.valid.corr.freq) +
  geom_col(aes(x=valid_tes, y=N, fill=valid_tes), color='black') +
  theme_bw() +
  ggtitle('TransFrags 3-prime end result, according to ref mRNA TES')

ggsave(file.path(outdir, 'TF_TES_validation.jpg'), ggvf, height = 12, width = 9)




### Annotate new Transcripts from TSS-validated TransFrags

## TransFrags that are not equal to refernce
ggVennDiagram::ggVennDiagram(list(
  TransFrag = TransFrags.TSSr.comb$transcript_id,
  CAGE_TC_supported   = TransFrags.TSSr.comb[!is.na(CAGE_TC_ID)   & class_code != '=',  transcript_id],
  dcDNA_TC_supported  = TransFrags.TSSr.comb[!is.na(dcDNA_TC_ID)  & class_code != '=',  transcript_id]
))

## TransFrags that are have at least '**' 
ggVennDiagram::ggVennDiagram(list(
  TransFrag = TransFrags.TSSr.comb$transcript_id,
  CAGE_TC_supported   = TransFrags.TSSr.comb[!is.na(CAGE_TC_ID)   & class_code != '=' & CAGE_TC_significance  != '*',  transcript_id],
  dcDNA_TC_supported  = TransFrags.TSSr.comb[!is.na(dcDNA_TC_ID)  & class_code != '=' & dcDNA_TC_significance != '*',  transcript_id]
))

## TransFrags that share 3-prime with previous annotation -->> significant drop
ggVennDiagram::ggVennDiagram(list(
  TransFrag = TransFrags.TSSr.comb$transcript_id,
  CAGE_TC_supported   = TransFrags.TSSr.comb[!is.na(CAGE_TC_ID)   & class_code != '=' & CAGE_TC_significance  != '*' & valid_tes == T,  transcript_id],
  dcDNA_TC_supported  = TransFrags.TSSr.comb[!is.na(dcDNA_TC_ID)  & class_code != '=' & dcDNA_TC_significance != '*' & valid_tes == T,  transcript_id]
))


TransFrags.TSSr.filt <-  TransFrags.TSSr.comb[  (!is.na(CAGE_TC_ID)  & CAGE_TC_significance   %in% c('**', '***'))
                                              | (!is.na(dcDNA_TC_ID) & dcDNA_TC_significance  %in% c('**', '***'))
                                              #& class_code != '='  & valid_tes == T
                                              ,  ]

TransFrags.TSSr.filt[,.N, by=.(CAGE_TC_significance, dcDNA_TC_significance)]

TransFrags.TSSr.filt <-  TransFrags.TSSr.filt[class_code != '='  & valid_tes == T, ]

TransFrags.TSSr.filt[,.N, by=.(CAGE_TC_significance, dcDNA_TC_significance, valid_tes)]



### Merge TransFrags with TSSr group-wise cluster results (for the peaks)
TransFrags.TSSr.filt.counts.gt <- merge(TR.gff.compare.merged.TR.counts.gt[,.(seqnames, strand, transcript_id, sample, rep, group, cell_line, hpi, Time, count, norm_count)], 
                                        TransFrags.TSSr.filt,
                                        by= c('seqnames', 'strand', 'transcript_id'),
                                        all.y = T)

TransFrags.TSSr.filt.counts.gt[, cell_line := gsub('-', '_', cell_line)]

## Merge with dcDNA clusters
TransFrags.TSSr.filt.counts.gt  <- merge(TransFrags.TSSr.filt.counts.gt, 
                                         dcDNA_TSSr_clusters[,.(cluster, seqnames = chr, strand, cell_line, hpi, 
                                                                cluster.start = start, cluster.end = end, dcDNA_TC_ID = paste0('TC_', cluster),
                                                                dominant_tss, tags, tags.dominant_tss,  q_0.1, q_0.9, interquantile_width, shape.score, cluster_width)],
                                         #by.y = c('seqnames', 'strand', 'cluster', 'group'),
                                         by = c('seqnames', 'strand', 'dcDNA_TC_ID', 'cell_line', 'hpi'),
                                         all.x = T )

# Caluclate distance of TransFrags to peak of TC                      
TransFrags.TSSr.filt.counts.gt[, TSS_dcDNA_distance := dominant_tss - prime5, by= transcript_id]

# Reanme columns
cnames <- c('cluster.start', 'cluster.end', 'dominant_tss', 'tags', 'tags.dominant_tss',  'q_0.1', 'q_0.9', 'interquantile_width', 'shape.score', 'cluster_width')
setnames(TransFrags.TSSr.filt.counts.gt,  cnames, paste0('dcDNA_', cnames))


## Merge with dcDNA clusters
TransFrags.TSSr.filt.counts.gt  <- merge(TransFrags.TSSr.filt.counts.gt, 
                                         CAGE_TSSr_clusters[,.(cluster, seqnames = chr, strand, cell_line, hpi, 
                                                                cluster.start = start, cluster.end = end, CAGE_TC_ID = paste0('TC_', cluster),
                                                                dominant_tss, tags, tags.dominant_tss,  q_0.1, q_0.9, interquantile_width, shape.score, cluster_width)],
                                         #by.y = c('seqnames', 'strand', 'cluster', 'group'),
                                         by = c('seqnames', 'strand', 'CAGE_TC_ID', 'cell_line', 'hpi'),
                                         all.x = T )

# Caluclate distance of TransFrags to peak of TC                      
TransFrags.TSSr.filt.counts.gt[, TSS_CAGE_distance := dominant_tss - prime5, by= transcript_id]

# Reanme columns
setnames(TransFrags.TSSr.filt.counts.gt,  cnames, paste0('CAGE_', cnames))



# Check distances
TransFrags.TSSr.filt.counts.gt[, .N, TSS_dcDNA_distance]
ggplot(TransFrags.TSSr.filt.counts.gt, aes(TSS_dcDNA_distance) ) + geom_histogram(bins=100) + geom_vline(xintercept = c(-25, 25))

TransFrags.TSSr.filt.counts.gt[, .N, TSS_CAGE_distance]
ggplot(TransFrags.TSSr.filt.counts.gt, aes(TSS_CAGE_distance) ) + geom_histogram(bins=100) + geom_vline(xintercept = c(-25, 25))



## TransFrags that share 5-primes with CAGE or dcDNA TSS cluster peak in the respective sample group with a +/- 10 bp wobble
ggVennDiagram::ggVennDiagram(list(
  TransFrag           = unique(TransFrags.TSSr.filt.counts.gt$transcript_id),
  CAGE_TC_supported   = unique(TransFrags.TSSr.filt.counts.gt[abs(TSS_CAGE_distance) <= 10,   transcript_id]),
  dcDNA_TC_supported  = unique(TransFrags.TSSr.filt.counts.gt[abs(TSS_dcDNA_distance) <= 10,  transcript_id])
))


TransFrags.TSSr.filt.counts.gt <- TransFrags.TSSr.filt.counts.gt[abs(TSS_CAGE_distance) <= 10 | abs(TSS_dcDNA_distance) <= 10, ]
TransFrags.TSSr.filt.counts.gt[is.na(TSS_CAGE_distance) & is.na(TSS_dcDNA_distance)]

# The TSS here will thus be the peak of the TC (either in the CAGE ir dcDNA - which is closer)

TransFrags.TSSr.filt.counts.gt[, TSS := fifelse(is.na(TSS_dcDNA_distance), CAGE_dominant_tss, 
                                          fifelse(is.na(TSS_CAGE_distance),  dcDNA_dominant_tss,
                                            fifelse(abs(TSS_dcDNA_distance) < abs(TSS_CAGE_distance), dcDNA_dominant_tss, CAGE_dominant_tss
                                                    ))), by=transcript_id]

TransFrags.TSSr.filt.counts.gt[,.N, by=.(TSS, TSS_dcDNA_distance, TSS_CAGE_distance)]



TransFrags.TSSr.filt.counts.gt <- unique(TransFrags.TSSr.filt.counts.gt[,.(seqnames, strand, start, end, TSS, TES,
                                                                transcript_id, cmp_ref, class_code, 
                                                                group, hpi, Time, cell_line, group, sample, count)])

TransFrags.TSSr.filt.group.counts <- TransFrags.TSSr.filt.counts.gt[,.(sum_count = sum(count)), by=.(seqnames, strand, TSS, TES, group, cell_line, hpi, Time)]

## Filter for read_count
TransFrags.TSSr.filt.group.counts <- TransFrags.TSSr.filt.group.counts[ sum_count >= 5, ]


NOV.TXs.uni <- TransFrags.TSSr.filt.group.counts[,.(sum_count = sum(sum_count)), by = .(seqnames, strand, TSS, TES) ]
















