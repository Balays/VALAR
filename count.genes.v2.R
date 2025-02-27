
##
## Requires:
#     TR.counts.sp
#     gene.clusters.all
#     TR.EX
#     viral.mrna
#     viral.ref




#### Import and format data ####

### Import gene clusters
setDT(gene.clusters.all)
## has to be unique !

## CTO-L, as this is the only one which share the TSS with another gene (UL21 in this case)
if(exists('genes_to_filter')) {
  gene.clusters.all <- gene.clusters.all[ ! gene %in% c(genes_to_filter), ]
} else (message('No genes were filtered out for counting!'))

##
gene_clusters     <- unique(gene.clusters.all[,gene_cluster])
genes             <- unique(gene.clusters.all[,gene])
gene.clusters.uni <- unique(gene.clusters.all[,.(seqnames, strand, gene, gene_cluster, Kinetic_class)])


message('The following genes will be counted: \n', paste(genes, collapse = '\n'))
##

## IMPORTANT !! cluster and individual gene TESs

message("The following genes TES will be counted seprately from the cluster's TES, \n
         as they have differnet canonic TES than the cluster TES: ")
print(gene.clusters.all[cluster_TES != TES.canonic, .(gene, gene_cluster, cluster_TES, TES.canonic)])

gene.clusters.all[, gene_cluster := fifelse(cluster_TES != TES.canonic, gene, gene_cluster)]


## start has to be larger than end
stopifnot(nrow(gene.clusters.all[TSS.win.start > TSS.win.end,]) < 1)
stopifnot(nrow(gene.clusters.all[TES.win.start > TES.win.end,]) < 1)


### Make .gff out of gene cluster table
TES.df <- gene.clusters.all[,.(seqnames, start=TES.win.start, end=TES.win.end, source='Torma_etal_2020', type='TES_window', phase='.', strand, Name=gene, ID=gene, gene_cluster)]

TSS.df <- gene.clusters.all[,.(seqnames, start=TSS.win.start, end=TSS.win.end, source='Torma_etal_2020', type='TSS_window', phase='.', strand, Name=gene, ID=gene, gene_cluster)]

gene.clusters.gff <- rbind(TES.df, TSS.df)

#TSS.df[,start := ifelse(strand == '+', TSS, TES)]
#TSS.df[,end   := ifelse(strand == '+', TES, TSS)]


### exclude spliced genes -> just for transcript-level counting !
#gene.clusters <- gene.clusters.all[!gene %in% spliced.genes, ]
gene.clusters <- gene.clusters.all

## spliced Transcripts
check.spliced.TRs <- check.spliced.TRs
if(check.spliced.TRs) { 
  spliced.TRs <- unique(viral.mrna[exon_number > 1, .(transcript_id, gene_id)])  
  message('The following spliced transcripts will be counted to their parent gene from the annotation: \n')
  print(spliced.TRs) 
  spliced.TRs <- unique(viral.mrna[exon_number > 1, .(transcript_id)])  
}



### Checking: gene column has to be unique here!
stopifnot(nrow(gene.clusters) == luniq(gene.clusters[,gene]) )

### Checking: are there overlaps between the TSS windows of the genes
gene.cluster.TSS.ov <-
  foverlaps2(DTx=gene.clusters[,.(seqnames, strand, gene, start=TSS.win.start, end=TSS.win.end)], 
             DTy=gene.clusters[,.(seqnames, strand, gene, start=TSS.win.start, end=TSS.win.end)],
             by=c('seqnames', 'strand', 'start', 'end'), type=c('any'), minoverlap = 1 )
## the TES windows overlap of cousre because the poly clusters
## but the TSSs should not overlap
gene.cluster.TSS.ov[gene != i.gene]
stopifnot(nrow(gene.cluster.TSS.ov[gene != i.gene]) == 0 )



### Import the Read-Transcripts 
TR.reads <- unique(TR.EX[ # type == 'transcript'
        , .(seqnames, strand, transcript_id, 
            start  = start.transcript, end = end.transcript, transcript_start = start.transcript, transcript_end = end.transcript)])
if(luniq(TR.reads$transcript_id) != nrow(TR.reads)) {message('there are duplicate transcript ids! maybe these are exons?')}

TR.reads[,prime3 := ifelse(strand == '+', end, start)][, prime5 := ifelse(strand == '+', start, end)]

TR.reads[,TR.prime5.win.start := prime5]
TR.reads[,TR.prime5.win.end   := prime5]
TR.reads[,TR.prime3.win.start := prime3]
TR.reads[,TR.prime3.win.end   := prime3]

#### #####
##



## Start the Workflow !
#### Find their overlaps ####

### 1. Based on the TSS
DTx <- TR.reads[,.(seqnames, strand, start = TR.prime5.win.start, end = TR.prime5.win.end, transcript_id)]
DTy <- gene.clusters[,.(seqnames, strand, start = TSS.win.start, end = TSS.win.end, gene)]

TR.gene.ov <-
 foverlaps2(DTx=DTx, 
            DTy=DTy,
            by=c('seqnames', 'strand', 'start', 'end'),
            #by.x=c('seqnames', 'strand', '', ''),
            #by.y=c('seqnames', 'strand', '',	''),
            type=c('within'), minoverlap = 1
           )
stopifnot(nrow(TR.gene.ov[i.start == i.end,]) == nrow(TR.gene.ov[,]))

TR.gene.TSS.ov <- TR.gene.ov[,.(seqnames, strand, 
                                gene.TSS.ov = gene, gene.TSS.win.start=start, gene.TSS.win.end=end, 
                                transcript_id, TR.prime5=i.start)]



### 2. Based on the TES
DTx <- TR.reads[,.(seqnames, strand, start = TR.prime3.win.start, end = TR.prime3.win.end, transcript_id)]
DTy <- gene.clusters[,.(seqnames, strand, start = TES.win.start, end = TES.win.end, gene)]

TR.gene.ov <-
  foverlaps2(DTx=DTx, 
             DTy=DTy,
             by=c('seqnames', 'strand', 'start', 'end'),
             #by.x=c('seqnames', 'strand', '', ''),
             #by.y=c('seqnames', 'strand', '',	''),
             type=c('within'), minoverlap = 1
  )

stopifnot(nrow(TR.gene.ov[i.start == i.end,]) == nrow(TR.gene.ov[,]))

TR.gene.TES.ov <- TR.gene.ov[,.(seqnames, strand, 
                                gene.TES.ov = gene, gene.TES.win.start=start, gene.TES.win.end=end, 
                                transcript_id, TR.prime3=i.start)]

### check
#TR.gene.TES.ov[transcript_id == '99990']
#TR.gene.TSS.ov[transcript_id == '99990']
#TR.reads[transcript_id == '99990']

### Merge the overlaps
TR.gene.ov <- merge(TR.gene.TSS.ov, TR.gene.TES.ov, by=c('seqnames', 'strand', 'transcript_id'))
  
### Find the overlap between the two
message('Number of transcripts that overlap with bothe the TSS and the TES of the gene: \n')
nrow(TR.gene.ov[ gene.TES.ov ==  gene.TSS.ov ])
message('Total number of transcripts: \n')
nrow(TR.gene.ov)

## Filter those, where both the TSS and TES were assigned to the same gene
TR.gene.ov <- TR.gene.ov[ gene.TES.ov ==  gene.TSS.ov ]
TR.gene.ov[,gene := gene.TSS.ov]
TR.gene.ov[,gene.TES.ov := NULL]
TR.gene.ov[,gene.TSS.ov := NULL]

TR.gene.ov <- unique(TR.gene.ov)

##
dups <- dup(TR.gene.ov[,transcript_id])
if(length(dups != 0)) {
  dups
  stop('these transcripts overlap with more than one gene, based on both TSS and TES.')
} else {
  message('OK! no transcripts that overlap with more than one gene, based on both TSS and TES.')
}

## 

## merge back gene cluster info
TR.gene.ov <- merge(gene.clusters[,.(seqnames, strand, gene, gene_cluster, Kinetic_class)],
                                TR.gene.ov, 
                                by.x=c('seqnames', 'strand', 'gene'),
                                by.y=c('seqnames', 'strand', 'gene'), all=F)


TR.gene.ov.all <- TR.gene.ov



#### Transcripts based on TES overlap only
tr.genes <- unique(TR.gene.ov.all[,transcript_id])

TR.gene.TES.ov <- TR.gene.TES.ov[!transcript_id %in% tr.genes]

### see whcih genes overlap with the transcripts, based on the TES
TR.gene.TES.ov <- merge(gene.clusters[,.(seqnames, strand, gene, gene_cluster, Kinetic_class)],
                        TR.gene.TES.ov, 
                        by.x=c('seqnames', 'strand', 'gene'),
                        by.y=c('seqnames', 'strand', 'gene.TES.ov'), all=F)

TR.gene.TES.ov <- merge(TR.reads[,.(seqnames, strand, transcript_id, prime5, prime3)],
                        TR.gene.TES.ov, 
                        by.x=c('seqnames', 'strand', 'transcript_id'),
                        by.y=c('seqnames', 'strand', 'transcript_id'), all=F)

setnames(TR.gene.TES.ov, old=c('prime3', 'prime5'), new=c('TR.prime3', 'TR.prime5'))

## associate them to clusters only
TR.gene.TES.ov <- 
  unique(TR.gene.TES.ov[, .(seqnames, strand, transcript_id, gene_cluster, gene.TES.win.start, gene.TES.win.end, TR.prime3, TR.prime5)])

## merge
if(length(intersect(TR.gene.TES.ov$transcript_id, TR.gene.ov.all$transcript_id)) == 0 ) {message('OK')} else {message('NOT OK')}
TR.gene.ov.all <- rbindlist(list(TR.gene.ov.all, TR.gene.TES.ov), fill=T)
stopifnot(length(dup(TR.gene.ov.all$transcript_id)) == 0) 



#### Transcripts based on TSS overlap only
#tr.gene.clusters <- unique(TR.gene.ov.all[,transcript_id])
#TR.gene.TSS.ov <- TR.gene.TSS.ov[!transcript_id %in% tr.gene.clusters]
tr.genes <- unique(TR.gene.ov.all[,transcript_id])
TR.gene.TSS.ov <- TR.gene.TSS.ov[!transcript_id %in% tr.genes]

TR.gene.TSS.ov <- merge(TR.reads[,.(seqnames, strand, transcript_id, prime5, prime3)],
                        TR.gene.TSS.ov, 
                        by.x=c('seqnames', 'strand', 'transcript_id', 'prime5'),
                        by.y=c('seqnames', 'strand', 'transcript_id', 'TR.prime5'), all=F)

setnames(TR.gene.TSS.ov, old=c('prime3', 'prime5', 'gene.TSS.ov'), new=c('TR.prime3', 'TR.prime5', 'gene'))

## add kinetic class info from gene clusters DT
TR.gene.TSS.ov <- merge(unique(gene.clusters[,.(seqnames, strand, gene, Kinetic_class)]),
                        TR.gene.TSS.ov, 
                        by.x=c('seqnames', 'strand', 'gene'),
                        by.y=c('seqnames', 'strand', 'gene'), all=F)

## merge
if(length(intersect(TR.gene.TSS.ov$transcript_id, TR.gene.ov.all$transcript_id)) == 0 ) {message('OK')} else {message('NOT OK')}
TR.gene.ov.all <- rbindlist(list(TR.gene.ov.all, TR.gene.TSS.ov), fill=T)
stopifnot(length(dup(TR.gene.ov.all$transcript_id)) == 0) 



#### SPLICED !!! -->> omit if desired or if theres no gff-compare result

### = transcripts of spliced genes from gff compare 
## Use the Itereative gff.compare 'best' results
if(check.spliced.TRs) {
  TR.polyc <- best.merged.result_gff.compare[  ## TR.gff.compare[type == 'transcript' &
    abs_junc_distance.tr < thresh.eq.junc &
    equal.5_prime.tr == T & equal.3_prime.tr == T & 
    strand == strand.ref &
    #class_code == '=' &
    cmp_ref %in% spliced.TRs, .(seqnames, start, end, strand, transcript_id, cmp_ref)
    ]
  
  TR.polyc[,prime3 := ifelse(strand == '+', end,   start)]
  TR.polyc[,prime5 := ifelse(strand == '+', start, end)]
  setnames(TR.polyc, old=c('start', 'end', 'cmp_ref'), new=c('transctipt_start', 'transcript_end', 'gene'))
  TR.polyc[,gene_cluster := gene] 
  TR.polyc <- merge(TR.polyc, unique(gene.clusters.all[,.(seqnames, strand, gene, Kinetic_class)]), by=c('seqnames', 'strand', 'gene'), all.x=T, allow.cartesian=TRUE)
  #### !!!!
  
  #### -->> USE THIS !!! 
  ### = transcripts of spliced genes from gff compare
  ## Use the Itereative gff.compare 'best' results
  TR.splice <- best.merged.result_gff.compare[  ## TR.gff.compare[type == 'transcript' &
    abs_junc_distance.tr < thresh.eq.junc &
      equal.5_prime.tr == T & equal.3_prime.tr == T & 
      strand == strand.ref &
      #class_code == '=' &
      cmp_ref %in% spliced.TRs, .(seqnames, start, end, strand, transcript_id, cmp_ref)
  ]
  
  TR.splice[,prime3 := ifelse(strand == '+', end,   start)]
  TR.splice[,prime5 := ifelse(strand == '+', start, end)]
  
  ## merge ref transcripts with their gene and ORF
  TR.splice <- merge(TR.splice, TR.genes, 
                     by.x=c('seqnames', 'strand', 'cmp_ref'), by.y=c('seqnames', 'strand', 'transcript_id'), 
                     all.x=T, allow.cartesian=TRUE)
  
  #setnames(TR.splice, old=c('start', 'end', 'cmp_ref'), new=c('transctipt_start', 'transcript_end', 'gene'))
  #TR.splice[,gene_cluster := gene] 
  
  TR.splice <- merge(gene.clusters.uni, 
                     TR.splice,
                     by.x=c('seqnames', 'strand', 'gene'), by.y=c('seqnames', 'strand', 'gene_id'), 
                     all.y=T, allow.cartesian=TRUE)
  ###
  setnames(TR.splice, old=c('start', 'end'), new=c('transctipt_start', 'transcript_end'))  ####
  
  cnames <- colnames(TR.polyc)
  TR.polyc   <- TR.splice[,..cnames]
  
  TR.gene.ov.all <- TR.gene.ov.all[!transcript_id %in% TR.polyc$transcript_id, ]
  TR.gene.ov.all <- rbindlist(list(TR.gene.ov.all, TR.polyc), fill=T)
  
  stopifnot(length(dup(TR.gene.ov.all$transcript_id)) == 0) 
} else {
  message('spliced transcripts were not treated differently.')
}

#### OK!




#### Merge back to include TRs that have no such overlaps
TR.gene.ov.all <- merge(TR.gene.ov.all[, .(seqnames, strand, transcript_id, gene, gene_cluster, Kinetic_class)],
                        TR.reads[,.(seqnames, strand, transcript_id, transcript_start, transcript_end, prime5, prime3)],
                        by=c('seqnames', 'strand', 'transcript_id'), all=T)
stopifnot(length(dup(TR.gene.ov.all$transcript_id)) == 0) 




#### Merge w TR counts (including apater info)
TR.gene.ov.counts <- merge(TR.gene.ov.all, 
                           TR.counts.sp,
                           by.x=c('transcript_id'), by.y='TR_ID')




#### Write results
export.gff3(gene.clusters.gff, paste0(virus, '.genes.TSS.TES.gff3'))
fwrite(TR.gene.ov.all,         paste0(outdir, '/TR.gene.ov.all.tsv'), sep = '\t')
fwrite(TR.gene.ov.counts,      paste0(outdir, '/TR.gene.ov.counts.tsv'), sep = '\t')




#### !!! Consider adapters!
adapters <- 'any'
source('Gene.counts.R')
fwrite(gene.sample_count.sp, paste0(res.dir, '/gene.sample.count.sp.tsv'), sep = '\t')
gene.sample_count.sp <- fread(paste0(res.dir, '/gene.sample.count.sp.tsv'), na.strings = '')


