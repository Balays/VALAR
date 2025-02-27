

#### Export reads in  .gff format as an input for StringTie
require(rtracklayer)

TR.gff  <- merge(EX.uni, TR.sum.counts, by='TR_ID' )


TR.gff[,TR_ID:=as.character(TR_ID)]
TR.gff[,exon_rank:=as.character(exon_rank)]
colnames(TR.gff)[c(1,9)] <- c('transcript_id', 'exon_number')

TR.exon <- TR.gff[,type:='exon']
TR.exon <- TR.exon[,.(seqnames, strand, type, start, end, transcript_id, exon_number)]

TR.mrna <- TR.gff[,transcript_start:=min(start), by=.(transcript_id)][,transcript_end:=max(end), by=.(transcript_id)][,type:='transcript']
TR.mrna <- TR.mrna[,.(seqnames, strand, type, transcript_start, transcript_end, transcript_id)]
colnames(TR.mrna)[4:5] <- c('start', 'end')
TR.mrna <- unique(TR.mrna)

TR.gff <- rbind(TR.exon, TR.mrna, fill=T)


### merge READ exons and transcripts
## transcripts
DT.tr <- TR.gff[type == 'transcript']
setnames(DT.tr, old=c('start', 'end'), new=c('start.transcript', 'end.transcript'))

DT.tr  <- DT.tr[, .(seqnames, strand, start.transcript, end.transcript, transcript_id)]


## exons
DT.ex <- TR.gff[type == 'exon', .(seqnames, strand, start, end, transcript_id)]
## Add prime5 and prime3 columns
DT.ex[, prime5 := ifelse(strand == '+', start, end)]
DT.ex[, prime3 := ifelse(strand == '+', end, start)]
## order to calculate exon number
DT.ex <- DT.ex[order(seqnames, strand, transcript_id, 
                     fifelse(strand == '+', prime5, -prime5))]
## Assign exon numbers
DT.ex <- DT.ex[, exon_number := 1:.N, by = .(transcript_id)]
## exon count per transcript
DT.ex[, exon_count := max(exon_number), by = .(transcript_id)]


## merge
TR.EX  <- merge(DT.ex, DT.tr, 
                by=c('seqnames', 'strand', 'transcript_id'),
                suffixes=c('', '.transcript'))


##
fwrite(TR.EX, paste0(outdir, '/TR.EX.tsv'))
##
rtracklayer::export.gff2(TR.gff, TR.reads.gfffile)


## only for multiple cell lines

if (multiCellLines) {
  
  trs.PC12 <- unique(TR.counts[grepl('PC-12', sample) & count > 1,TR_ID])
  trs.C6   <- unique(TR.counts[grepl('C6', sample) & count > 1,TR_ID])
  trs.PK15 <- unique(TR.counts[grepl('PK-15', sample) & count > 1,TR_ID])
  
  
  rtracklayer::export.gff2(TR.gff, TR.reads.gfffile )
  
  
  TR.gff.PC12 <- TR.gff[transcript_id %in% trs.PC12,]
  rtracklayer::export.gff2(TR.gff.PC12, paste0(outdir, '/PC12.TR.reads.gff2'))
  
  TR.gff.C6  <- TR.gff[transcript_id %in% trs.C6,]
  rtracklayer::export.gff2(TR.gff.C6,  paste0(outdir, '/C6.TR.reads.gff2'))
  
  TR.gff.PK15 <- TR.gff[transcript_id %in% trs.PK15,]
  rtracklayer::export.gff2(TR.gff.PK15,  paste0(outdir, '/PK15.TR.reads.gff2'))
  
}